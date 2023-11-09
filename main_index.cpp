#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "rlcsa_builder.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

// #include "argument_parser.hpp"

KSEQ_INIT(gzFile, gzread)

// using namespace CSA;

// const int MAX_THREADS = 64;

// static unsigned char seq_nt6_table[128] = {
//     0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
//     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
//     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
//     5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
//     5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
//     5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

// void seq_char2nt6(int l, unsigned char *s) {
//   for (int i = 0; i < l; ++i)
//     s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
// }

uint64_t get_size_fa(const char *fa_path) {
  gzFile fp = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fp);
  uint64_t tot_l = 0;
  int l;
  while ((l = kseq_read(seq)) >= 0)
    tot_l += l + 1;
  kseq_destroy(seq);
  gzclose(fp);
  return tot_l;
}

void concatenate_fa(const char *fa_path, uint64_t size, unsigned char *data) {
  gzFile fp = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  uint64_t curr_l = 0;
  while ((l = kseq_read(seq)) >= 0) {
    // seq_char2nt6(l, (unsigned char *)seq->seq.s);
    memmove(data + curr_l, seq->seq.s, l);
    curr_l += l + 1;
    spdlog::debug("Moved {} chars (out of {})", curr_l, size);
  }
  kseq_destroy(seq);
  gzclose(fp);
}

int main_index(int argc, char **argv) {
  uint sample_rate = 0; // (1 << 31);
  int block_size = 32;
  int threads = 1;
  std::string index_prefix = "RLCSA";

  // TODO: improve CLI
  int c;
  while ((c = getopt(argc, argv, "i:@:vh")) >= 0) {
    switch (c) {
    case 'i':
      index_prefix = optarg;
      continue;
    case '@':
      threads = atoi(optarg);
      continue;
    case 'v':
      spdlog::set_level(spdlog::level::debug);
      continue;
    case 'h':
      // cerr << INDEX_USAGE_MESSAGE;
      return 0;
    default:
      // cerr << INDEX_USAGE_MESSAGE;
      return 1;
    }
  }
  if (argc - optind < 1) {
    // cerr << INDEX_USAGE_MESSAGE;
    return 1;
  }
  std::vector<std::string> inputs;
  while (optind < argc)
    inputs.push_back(argv[optind++]);

  spdlog::info("Initializing..");
  CSA::RLCSABuilder builder(block_size, sample_rate, 0, threads, NULL);

  uint64_t max_size = ((uint64_t)1 << 32) - 1;
  unsigned char *data =
      (unsigned char *)calloc(max_size, sizeof(unsigned char));
  for (const std::string &fa_path : inputs) {
    spdlog::info("Reading {}..", fa_path);
    uint64_t size = get_size_fa(fa_path.c_str());
    concatenate_fa(fa_path.c_str(), size, data);
    // Debug
    // for (uint i = 0; i < size; ++i)
    //   if (data[i] == 0)
    //     std::cerr << '|';
    //   else
    //     std::cerr << data[i];
    // std::cerr << std::endl;

    spdlog::info("Indexing {}..", fa_path);
    CSA::RLCSA *index =
        new CSA::RLCSA(data, size, block_size, sample_rate, threads, false);
    spdlog::info("Storing {}..", fa_path);
    if (index != 0 && index->isOk())
      index->writeTo(fa_path); // TODO: store to tmp dir

    spdlog::info("Merging {}..", fa_path);
    builder.insertFromFile(fa_path, data);

    delete index;
  }
  free(data);

  spdlog::info("Storing full index to {}..", index_prefix);
  CSA::RLCSA *rlcsa = builder.getRLCSA();
  if (rlcsa != 0 && rlcsa->isOk()) {
    std::cout << std::endl;
    rlcsa->printInfo();
    rlcsa->reportSize(true);
    rlcsa->writeTo(index_prefix);
  }
  delete rlcsa;

  return 0;
}
