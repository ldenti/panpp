#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include "fmd_simple.hpp"
#include "usage.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread)

std::string interval2str(const FMDPosition &i) {
  return std::to_string(i.forward_start) + "," +
         std::to_string(i.reverse_start) + "," +
         std::to_string(i.end_offset + 1);
}

void assemble(std::vector<std::pair<int, int>> &specifics) {
  std::vector<std::pair<int, int>> assembled_specifics;
  int i = specifics.size() - 1;
  while (i >= 0) {
    int j;
    for (j = i - 1; j >= 0; --j) {
      if (specifics[j + 1].second < specifics[j].first) {
        // non-overlapping
        assembled_specifics.push_back(
            std::make_pair(specifics[i].first, specifics[j + 1].second));
        i = j;
        break;
      }
    }
    if (j < 0) {
      assembled_specifics.push_back(
          std::make_pair(specifics[i].first, specifics[0].second));
      i = j;
    }
  }
  specifics = assembled_specifics;
}

std::vector<std::pair<int, int>> ping_pong_search(const FMD *index, char *seq,
                                                  int l) {
  std::vector<std::pair<int, int>> result;
  int begin = l - 1;
  uint8_t c;
  FMDPosition sai;
  while (begin >= 0) {
    c = seq[begin];
    sai = index->getCharPosition(c);
    spdlog::debug("Starting from {} ({}): {}", int2char[seq[begin]], begin,
                  interval2str(sai));
    // Backward search. Find a mismatching sequence. Stop at first mismatch.
    int bmatches = 0;

    while (!sai.isEmpty() && begin > 0) {
      --begin;
      c = seq[begin];
      ++bmatches;
      sai = index->extend(sai, c, true);
      spdlog::debug("Backward extending with {} ({}): {}", int2char[seq[begin]],
                    begin, interval2str(sai));
    }
    // prefix is a match
    if (begin == 0 && !sai.isEmpty())
      break;

    spdlog::debug("Mismatch {} ({}). Matches: {}", int2char[seq[begin]], begin,
                  std::to_string(bmatches));

    // Forward search
    int end = begin;
    int fmatches = 0;
    c = seq[end];
    sai = index->getCharPosition(c);
    spdlog::debug("Starting from {} ({}): {}", int2char[seq[end]], end,
                  interval2str(sai));
    while (!sai.isEmpty()) {
      ++end;
      c = seq[end];
      ++fmatches;
      sai = index->extend(sai, c, false);
      spdlog::debug("Forward extending with {} ({}): {}", int2char[c], end,
                    interval2str(sai));
    }
    spdlog::debug("Mismatch {} ({}). Matches: {}", int2char[c], end, fmatches);

    // add solution
    // TODO: to this better
    spdlog::debug("Adding [{}, {}]", begin, end);
    result.push_back(std::make_pair(begin, end));

    // prepare for next round
    if (begin == 0)
      break;
    begin = end - 1;
    // TODO: add relaxed and overlap (?)
  }

  return result;
}

int main_pingpong(int argc, char **argv) {
  bool verbose = false;
  int c;
  int flank = 0;
  int threads = 1;
  bool assemble_flag = true;
  bool fx_out = false;
  while ((c = getopt(argc, argv, "f:@:xavh")) >= 0) {
    switch (c) {
    case 'f':
      flank = atoi(optarg);
      continue;
    case 'x':
      fx_out = true;
      continue;
    case '@':
      threads = atoi(optarg);
      continue;
    case 'a':
      assemble_flag = false;
      continue;
    case 'v':
      spdlog::set_level(spdlog::level::debug);
      verbose = true;
      continue;
    case 'h':
      std::cerr << SEARCH_USAGE_MESSAGE;
      return 0;
    default:
      std::cerr << SEARCH_USAGE_MESSAGE;
      return 1;
    }
  }

  char *index_prefix = argv[optind++];
  char *query_path = argv[optind++];

  spdlog::info("Loading index..");
  const FMD *fmd = new FMD(index_prefix, false);
  FMDPosition range;
  if (verbose) {
    for (int i = 1; i < 5; ++i) {
      range = fmd->getCharPosition(i);
      spdlog::debug("{}: {},{} ({})", int2char[i], range.forward_start,
                    range.reverse_start, range.end_offset + 1);
    }
  }
  gzFile fp = gzopen(query_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  std::vector<std::pair<int, int>> result;
  while ((l = kseq_read(seq)) >= 0) {
    spdlog::info("Checking {}..", seq->name.s);
    seq_char2nt6(l, seq->seq.s);
    result = ping_pong_search(fmd, seq->seq.s, l);
    for (auto &r : result) {
      r.first = std::max(0, r.first - flank);
      r.second = std::min(r.second + flank, l - 1);
    }
    if (assemble_flag)
      assemble(result);

    bool is_first = true;
    for (auto &r : result) {
      if (fx_out) {
        // TODO
        // std::cout << (is_first ? seq->name.s : "*") << "\t" << r.first <<
        // "\t"
        //           << r.second - r.first + 1 << "\t" << 0 << std::endl;
      } else {
        std::cout << (is_first ? seq->name.s : "*") << "\t" << r.first << "\t"
                  << r.second - r.first + 1 << "\t" << 0 << std::endl;
        is_first = false;
      }
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  spdlog::info("Done.");

  return 0;
}

int main_exact(int argc, char **argv) {
  char *index_prefix = argv[1];
  char *query_path = argv[2];

  const CSA::RLCSA *rlcsa = new CSA::RLCSA(index_prefix, false);
  CSA::pair_type range;
  gzFile fp = gzopen(query_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    seq_char2nt6(l, seq->seq.s);
    range = rlcsa->count(seq->seq.s);
    std::cout << seq->name.s << ": "
              << range.first + rlcsa->getNumberOfSequences() << ","
              << range.second + rlcsa->getNumberOfSequences() << std::endl;
  }
  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}

int main_search(int argc, char **argv) {
  bool verbose = false;
  int c;
  while ((c = getopt(argc, argv, "vh")) >= 0) {
    switch (c) {
    case 'v':
      spdlog::set_level(spdlog::level::debug);
      verbose = true;
      continue;
    case 'h':
      std::cerr << "ERROR";
      return 0;
    default:
      std::cerr << "ERROR";
      return 1;
    }
  }

  char *index_prefix = argv[optind++];
  char *query = argv[optind++];
  seq_char2nt6(strlen(query), query);

  const CSA::RLCSA *rlcsa = new CSA::RLCSA(index_prefix, false);

  if (verbose) {
    unsigned char *bwt = rlcsa->readBWT();
    uint size = rlcsa->getSize() + rlcsa->getNumberOfSequences();
    for (uint i = 0; i < size; ++i)
      if (bwt[i] == 0)
        std::cerr << '|';
      else
        std::cerr << int2char[bwt[i]];
    std::cerr << std::endl;
  }

  rlcsa->printInfo();

  auto x = rlcsa->getBWTRange();
  std::cerr << x.first << "," << x.second << std::endl;

  auto r = rlcsa->count(query);
  std::cerr << r.first + rlcsa->getNumberOfSequences() << ","
            << r.second + rlcsa->getNumberOfSequences() << std::endl;

  delete rlcsa;

  return 0;
}
