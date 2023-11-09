#include <cstdint>
#include <iostream>
#include <zlib.h>

#include "kseq.h"
#include "misc/utils.h"
#include "rlcsa.h"
#include "rlcsa_builder.h"

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

uint get_size_fa(char *fa_path) {
  gzFile fp = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fp);
  uint tot_l = 0;
  int l;
  while ((l = kseq_read(seq)) >= 0)
    tot_l += l + 1;
  kseq_destroy(seq);
  gzclose(fp);
  return tot_l;
}

void concatenate_fa(char *fa_path, uint size, unsigned char *data) {
  gzFile fp = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  uint curr_l = 0;
  while ((l = kseq_read(seq)) >= 0) {
    // seq_char2nt6(l, (unsigned char *)seq->seq.s);
    memmove(data + curr_l, seq->seq.s, l);
    curr_l += l + 1;
    // TODO: better progression
    std::cerr << "Moved " << curr_l << " chars (out of " << size << ")"
              << std::endl;
  }
  kseq_destroy(seq);
  gzclose(fp);
}

// uint index_fa(char *fa_path, CSA::usint threads) {
//   gzFile fp = gzopen(fa_path, "r");
//   kseq_t *seq = kseq_init(fp);

//

//   std::cerr << "Concatenating.." << std::endl;
//   int l;
//   uint tot_l = 0;
//   while ((l = kseq_read(seq)) >= 0)
//     tot_l += l + 1;
//   gzclose(fp);
//   fp = gzopen(fa_path, "r");
//   seq = kseq_init(fp);
//   // tot_l *= 2;
//   std::cerr << "Allocating " << tot_l << " chars..";
//   unsigned char *data = (unsigned char *)calloc(tot_l, sizeof(unsigned
//   char)); std::cerr << "Done." << std::endl; uint curr_l = 0; while ((l =
//   kseq_read(seq)) >= 0) {
//     // seq_char2nt6(l, (unsigned char *)seq->seq.s);
//     std::cerr << "Moving " << l << " to " << curr_l << ".." << std::endl;
//     memmove(data + curr_l, seq->seq.s, l);
//     curr_l += l + 1;
//     // TODO: better progression
//     std::cerr << "Moved " << curr_l << " chars (out of " << tot_l << ")"
//               << std::endl;
//   }
//   kseq_destroy(seq);
//   gzclose(fp);

//

//   return tot_l;
// }

int main(int argc, char **argv) {
  char *fa_path = argv[1];
  int threads = 4;
  uint sample_rate = 0; // (1 << 31);
  int block_size = 32;

  // TODO: improve CLI

  // Phase 1: indexing each fa independently
  // TODO: cycle over files
  int size = get_size_fa(fa_path);
  unsigned char *data = (unsigned char *)calloc(size, sizeof(unsigned char));
  concatenate_fa(fa_path, size, data);
  // // Debug
  // for (uint i = 0; i < size; ++i)
  //   if (data[i] == 0)
  //     std::cerr << '|';
  //   else
  //     std::cerr << data[i];
  // std::cerr << std::endl;

  // std::cerr << "Indexing.." << std::endl;
  CSA::RLCSA *index =
      new CSA::RLCSA(data, size, block_size, sample_rate, threads, false);
  std::cerr << "Storing.." << std::endl;
  if (index != 0 && index->isOk()) {
    index->writeTo(fa_path); // TODO: store to tmp dir
    // total_size += size;
  }
  // delete index;

  CSA::RLCSABuilder builder(block_size, sample_rate, 0, threads, NULL);

  // Phase 2: merging
  // TODO: cycle over files

  builder.insertFromFile(fa_path, data);
  builder.insertFromFile(fa_path, data);

  CSA::RLCSA *rlcsa = builder.getRLCSA();
  // if (rlcsa != 0 && rlcsa->isOk()) {
  //   rlcsa->printInfo();
  //   rlcsa->reportSize(true);
  //   rlcsa->writeTo("XXX");
  //   // size = index->getSize() / (double)CSA::MEGABYTE;
  // }
  // delete rlcsa;

  // increment
  // std::vector<CSA::usint> end_markers;
  // CSA::usint* ranks = this->getRanks(sequence, length, end_markers);
  // if(delete_sequence) { delete[] sequence; }

  // double mark = readTimer();
  // parallelSort(ranks, ranks + length);
  // #pragma omp parallel for schedule(static)
  // for(usint i = 0; i < length; i++) { ranks[i] += i + 1; }
  // this->sort_time += readTimer() - mark;

  // this->mergeRLCSA(increment, ranks, length);

  // TODO
  // double megabytes = size / (double)CSA::MEGABYTE;
  // std::cout << "Indexed " << megabytes << " megabytes in " << total_time
  //             << " seconds (" << (megabytes / total_time) << " MB/s)."
  //             << std::endl;
  //   std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  //   std::cout << std::endl;

  // // TODO
  // // double stop = readTimer();
  // // std::cout << megabytes << " megabytes indexed in " << (stop - start)
  // //           << " seconds (" << (megabytes / (stop - start)) << " MB/s)."
  // //           << std::endl;
  // // if (do_merge) {
  // //   std::cout << "Build time:    " << build_time + builder.getBuildTime()
  // //             << " seconds" << std::endl;
  // //   std::cout << "Search time:   " << builder.getSearchTime() << "
  // seconds"
  // //             << std::endl;
  // //   std::cout << "Sort time:     " << builder.getSortTime() << " seconds"
  // //             << std::endl;
  // //   std::cout << "Merge time:    " << builder.getMergeTime() << " seconds"
  // //             << std::endl;
  // // }
  // // std::cout << "Memory usage:  " << memoryUsage() << " kB" << std::endl;
  // // std::cout << std::endl;

  // const CSA::RLCSA *rlcsa = new CSA::RLCSA("XXX", false);

  // unsigned char *bwt = index->readBWT();
  // for (int i = 0; i < size; ++i)
  //   if (bwt[i] == 0)
  //     std::cerr << '|';
  //   else
  //     std::cerr << bwt[i];
  // std::cerr << std::endl;

  // bwt = rlcsa->readBWT();
  // for (int i = 0; i < size * 2; ++i)
  //   if (bwt[i] == 0)
  //     std::cerr << '|';
  //   else
  //     std::cerr << bwt[i];
  // std::cerr << std::endl;

  std::cerr << rlcsa->getSize() << " vs " << index->getSize() << std::endl;
  std::cerr << rlcsa->getNumberOfSequences() << " vs "
            << index->getNumberOfSequences() << std::endl;

  auto x = rlcsa->getBWTRange();
  std::cerr << x.first << "," << x.second << std::endl;
  x = index->getBWTRange();
  std::cerr << x.first << "," << x.second << std::endl;

  index->printInfo();
  rlcsa->printInfo();

  delete index;
  delete rlcsa;
  delete data;

  return 0;
}
