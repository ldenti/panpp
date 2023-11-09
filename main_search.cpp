#include <cstdint>
#include <getopt.h>
#include <iostream>

#include "rlcsa.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

int main_search(int argc, char **argv) {
  bool verbose = false;
  // TODO: improve CLI
  int c;
  while ((c = getopt(argc, argv, "vh")) >= 0) {
    switch (c) {
    // case '@':
    //   threads = atoi(optarg);
    //   continue;
    case 'v':
      spdlog::set_level(spdlog::level::debug);
      verbose = true;
      continue;
    case 'h':
      // cerr << INDEX_USAGE_MESSAGE;
      return 0;
    default:
      // cerr << INDEX_USAGE_MESSAGE;
      return 1;
    }
  }
  char *index_prefix = argv[optind++];
  char *query = argv[optind++];

  const CSA::RLCSA *rlcsa = new CSA::RLCSA(index_prefix, false);

  if (verbose) {
    unsigned char *bwt = rlcsa->readBWT();
    uint size = rlcsa->getSize() + rlcsa->getNumberOfSequences();
    for (int i = 0; i < size; ++i)
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
