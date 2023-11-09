#include <cstdint>
#include <iostream>
// #include <zlib.h>

// #include "kseq.h"
#include "misc/utils.h"
#include "rlcsa.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

// #include "argument_parser.hpp"

// KSEQ_INIT(gzFile, gzread)

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

int main_search(int argc, char **argv) {
  // parse_arguments(argc, argv);

  // if (opt::verbose)
  //   spdlog::set_level(spdlog::level::debug);

  spdlog::info("Initializing..");

  char *index_prefix = argv[1];
  char *query = argv[2];
  // int threads = 4;

  const CSA::RLCSA *rlcsa = new CSA::RLCSA(index_prefix, false);

  // unsigned char *bwt = rlcsa->readBWT();
  // for (int i = 0; i < size; ++i)
  //   if (bwt[i] == 0)
  //     std::cerr << '|';
  //   else
  //     std::cerr << bwt[i];
  // std::cerr << std::endl;

  std::cerr << rlcsa->getSize() << std::endl;
  std::cerr << rlcsa->getNumberOfSequences() << std::endl;

  auto x = rlcsa->getBWTRange();
  std::cerr << x.first << "," << x.second << std::endl;

  rlcsa->printInfo();

  auto r = rlcsa->count(query);
  std::cerr << r.first << "," << r.second << std::endl;

  delete rlcsa;

  return 0;
}
