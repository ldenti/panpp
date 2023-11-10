#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include "usage.hpp"

KSEQ_INIT(gzFile, gzread)

int main_exact(int argc, char **argv) {
  char *index_prefix = argv[1];
  char *query_path = argv[2];

  const CSA::RLCSA *rlcsa = new CSA::RLCSA(index_prefix, false);
  CSA::pair_type range;
  gzFile fp = gzopen(query_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    range = rlcsa->count(seq->seq.s);
    std::cerr << seq->name.s << ": "
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
      std::cerr << SEARCH_USAGE_MESSAGE;
      return 0;
    default:
      std::cerr << SEARCH_USAGE_MESSAGE;
      return 1;
    }
  }

  char *index_prefix = argv[optind++];
  char *query = argv[optind++];

  const CSA::RLCSA *rlcsa = new CSA::RLCSA(index_prefix, false);

  if (verbose) {
    unsigned char *bwt = rlcsa->readBWT();
    uint size = rlcsa->getSize() + rlcsa->getNumberOfSequences();
    for (uint i = 0; i < size; ++i)
      if (bwt[i] == 0)
        std::cerr << '|';
      else
        std::cerr << bwt[i];
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
