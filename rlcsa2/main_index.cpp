#include <assert.h>
#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#include "kseq.h"
#include "rlcsa.h"
#include "rlcsa_builder.h"
#include "spdlog/spdlog.h"

#include "usage.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread)

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

void concatenate_fa(const char *fa_path, uint64_t size, unsigned char *data,
                    unsigned char *data_r) {
  gzFile fp = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  uint64_t curr_l = 0;
  int i;
  while ((l = kseq_read(seq)) >= 0) {
    seq_char2nt6(l, seq->seq.s);

    memmove(data + curr_l, seq->seq.s, l);

    // reverse
    for (i = 0; i < (l >> 1); ++i) {
      int tmp = seq->seq.s[l - 1 - i];
      tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
      seq->seq.s[l - 1 - i] = (seq->seq.s[i] >= 1 && seq->seq.s[i] <= 4)
                                  ? 5 - seq->seq.s[i]
                                  : seq->seq.s[i];
      seq->seq.s[i] = tmp;
    }
    if (l & 1)
      seq->seq.s[i] = (seq->seq.s[i] >= 1 && seq->seq.s[i] <= 4)
                          ? 5 - seq->seq.s[i]
                          : seq->seq.s[i];
    memmove(data_r + curr_l, seq->seq.s, l);

    curr_l += l;

    data[curr_l] = '\0';
    data_r[curr_l] = '\0';

    ++curr_l;
    spdlog::debug("Moved {} chars (out of {})", curr_l, size);
  }
  kseq_destroy(seq);
  gzclose(fp);
}

int main_index(int argc, char **argv) {
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
      std::cerr << INDEX_USAGE_MESSAGE;
      return 0;
    default:
      std::cerr << INDEX_USAGE_MESSAGE;
      return 1;
    }
  }
  if (argc - optind < 1) {
    std::cerr << INDEX_USAGE_MESSAGE;
    return 1;
  }
  std::vector<std::string> inputs;
  while (optind < argc)
    inputs.push_back(argv[optind++]);

  spdlog::info("Initializing..");
  CSA::RLCSABuilder builder(block_size, 0, threads, NULL);

  uint64_t max_size = ((uint64_t)1 << 32) - 1;
  unsigned char *data =
      (unsigned char *)calloc(max_size, sizeof(unsigned char));
  unsigned char *data_r =
      (unsigned char *)calloc(max_size, sizeof(unsigned char));
  for (const std::string &fa_path : inputs) {
    spdlog::info("Reading {}..", fa_path);
    uint64_t size = get_size_fa(fa_path.c_str());
    assert(size < max_size);
    concatenate_fa(fa_path.c_str(), size, data, data_r);
    // Debug
    // for (uint i = 0; i < size; ++i)
    //   if (data[i] == 0)
    //     std::cerr << '|';
    //   else
    //     std::cerr << data[i];
    // std::cerr << std::endl;

    spdlog::info("Indexing {}..", fa_path);
    CSA::RLCSA *index = new CSA::RLCSA(data, size, block_size, threads);
    spdlog::info("Storing {}..", fa_path);
    if (index != 0 && index->isOk())
      index->writeTo(fa_path); // TODO: store to tmp dir
    delete index;

    spdlog::info("Merging {}..", fa_path);
    builder.insertFromFile(fa_path, data);

    spdlog::info("Indexing {} (R)..", fa_path);
    index = new CSA::RLCSA(data_r, size, block_size, threads);
    spdlog::info("Storing {} (R)..", fa_path);
    char *rev = (char *)malloc(fa_path.size() + 2);
    strcpy(rev, fa_path.c_str());
    strcat(rev, ".R");
    if (index != 0 && index->isOk())
      index->writeTo(rev); // TODO: store to tmp dir
    delete index;
    spdlog::info("Merging {} (R)..", fa_path);
    builder.insertFromFile(rev, data_r);
  }
  free(data);

  spdlog::info("Storing full index to {}..", index_prefix);
  CSA::RLCSA *rlcsa = builder.getRLCSA();
  if (rlcsa != 0 && rlcsa->isOk()) {
    std::cout << std::endl;
    // rlcsa->printInfo();
    // rlcsa->reportSize(true);
    rlcsa->writeTo(index_prefix);
  }
  delete rlcsa;

  spdlog::info("Done.");

  return 0;
}
