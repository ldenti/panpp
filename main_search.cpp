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

struct fxentry_t {
  char *name;
  char *seq;
  char *qual;
  int l;
  int size = 32768;

  fxentry_t() {
    name = (char *)malloc(size);
    seq = (char *)malloc(size);
    qual = (char *)malloc(size);
    l = 0;
  }

  void set(char *_name, char *_seq, char *_qual, int _l) {
    l = _l;
    if (l > size) {
      size = l + 1;
      seq = (char *)realloc(seq, size * sizeof(char));
      qual = (char *)realloc(qual, size * sizeof(char));
    }
    strcpy(name, _name);
    strcpy(seq, _seq);
    if (_qual == NULL)
      qual = NULL;
    else
      strcpy(qual, _qual);
  }
};

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

void ping_pong_search(const FMD *index, fxentry_t &fxe,
                      std::vector<std::pair<int, int>> &result) {
  char *seq = fxe.seq;
  int l = fxe.l;
  seq_char2nt6(l, seq);

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
}

int load_batch(kseq_t *seq, std::vector<fxentry_t> &input) {
  int l = 0;
  uint b = 0;
  while (b < input.size() && (l = kseq_read(seq)) >= 0) {
    if (seq->qual.l == 0)
      input[b].set(seq->name.s, seq->seq.s, NULL, l);
    else
      input[b].set(seq->name.s, seq->seq.s, seq->qual.s, l);
    ++b;
  }
  return b;
}

void dump(fxentry_t &fxe, std::vector<std::pair<int, int>> output,
          bool fx_out) {
  bool is_first = true;
  for (const auto &o : output) {
    if (fx_out) {
      std::string seq(fxe.seq, o.first, o.second - o.first + 1);
      std::transform(seq.cbegin(), seq.cend(), seq.begin(),
                     [](unsigned char c) { return int2char[c]; });
      if (fxe.qual == NULL) {
        std::cout << ">" << fxe.name << ":" << o.first << "-" << o.second
                  << "\n"
                  << seq << std::endl;
      } else {
        std::string qual(fxe.qual, o.first, o.second - o.first + 1);
        std::cout << "@" << fxe.name << ":" << o.first << "-" << o.second
                  << "\n"
                  << seq << "\n"
                  << "+"
                  << "\n"
                  << qual << std::endl;
      }
    } else {
      std::cout << (is_first ? fxe.name : "*") << "\t" << o.first << "\t"
                << o.second - o.first + 1 << std::endl;
      is_first = false;
    }
  }
}

int main_pingpong(int argc, char **argv) {
  (void)(argc); // suppress unused parameter warning

  int c;
  int flank = 0;
  int bsize = 10000;
  int threads = 1;
  bool assemble_flag = true;
  bool fx_out = false;
  // bool verbose = false;
  while ((c = getopt(argc, argv, "f:@:b:xavh")) >= 0) {
    switch (c) {
    case 'f':
      flank = atoi(optarg);
      continue;
    case 'b':
      bsize = atoi(optarg);
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
      // verbose = true;
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
  const FMD *index = new FMD(index_prefix, false);

  spdlog::info("Allocating space..");
  FMDPosition range;

  std::vector<fxentry_t> input(bsize);
  std::vector<std::vector<std::pair<int, int>>> output(bsize);

  gzFile fp = gzopen(query_path, "r");
  kseq_t *seq = kseq_init(fp);

  int curr_bsize = 0;
  int b = 0; // iterator for the batch
  int n = 1; // current batch
  while (true) {
    spdlog::info("Loading batch {}..\r", n);
    curr_bsize = load_batch(seq, input);
    if (curr_bsize <= 0)
      break;
    spdlog::info("Analyzing batch {}..\r", n);
#pragma omp parallel for num_threads(threads)
    for (b = 0; b < curr_bsize; ++b) {
      ping_pong_search(index, input[b], output[b]);
      for (auto &r : output[b]) {
        r.first = std::max(0, r.first - flank);
        r.second = std::min(r.second + flank, input[b].l - 1);
      }
      if (assemble_flag)
        assemble(output[b]);
    }
    spdlog::info("Dumping batch {}..\r", n);
    for (b = 0; b < curr_bsize; ++b)
      dump(input[b], output[b], fx_out);

    for (auto &o : output)
      o.clear();
    ++n;
  }
  spdlog::info("Dumped {} batches.", n - 1);

  kseq_destroy(seq);
  gzclose(fp);

  return 0;
}

int main_exact(int argc, char **argv) {
  (void)(argc); // suppress unused parameter warning

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
  (void)(argc); // suppress unused parameter warning

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

// TODO: Move this to test subcommand
// if (verbose) {
//   for (int i = 1; i < 5; ++i) {
//     range = fmd->getCharPosition(i);
//     spdlog::debug("{}: {},{} ({})", int2char[i], range.forward_start,
//                   range.reverse_start, range.end_offset + 1);
//   }
// }
