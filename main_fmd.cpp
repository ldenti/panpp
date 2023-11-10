#include <assert.h>
#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#include "kseq.h"
#include "mrope.h"
#include "rld0.h"
#include "rle.h"
#include "spdlog/spdlog.h"

#include "utils.hpp"

KSEQ_INIT(gzFile, gzread)

int main_fmdindex(int argc, char *argv[]) {
  spdlog::info("Initializing..");
  // add CLI to include reverse
  char *fa_path = argv[1];

  // hardcoded parameters
  uint64_t m = (uint64_t)(.97 * 10 * 1024 * 1024 * 1024) +
               1; // batch size for multi-string indexing
  int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES,
      so = MR_SO_RCLO;
  int thr_min =
      100; // switch to single thread when < 100 strings remain in a batch

  // the index
  mrope_t *mr = mr_init(max_nodes, block_len, so);
  mr_thr_min(mr, thr_min);

  // Parsing the input sample
  spdlog::info("Iterating over sequences..");
  gzFile fp = gzopen(fa_path, "rb");
  kseq_t *ks = kseq_init(fp);
  kstring_t buf = {0, 0, 0}; // buffer, will contain the concatenation
  int l;
  uint8_t *s;
  int i;
  while ((l = kseq_read(ks)) >= 0) {
    spdlog::info("{}", ks->name.s);
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

    // Reverse the sequence
    for (i = 0; i < (l >> 1); ++i) {
      int tmp = s[l - 1 - i];
      s[l - 1 - i] = s[i];
      s[i] = tmp;
    }

    // Add forward to buffer
    kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);

    // // Add reverse to buffer
    // for (i = 0; i < (l >> 1); ++i) {
    //   int tmp = s[l - 1 - i];
    //   tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
    //   s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    //   s[i] = tmp;
    // }
    // if (l & 1)
    //   s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    // kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);

    if (buf.l >= m) {
      mr_insert_multi(mr, buf.l, (uint8_t *)buf.s, 1);
      buf.l = 0;
    }
  }

  if (buf.l) { // last batch
    mr_insert_multi(mr, buf.l, (uint8_t *)buf.s, 1);
  }

  free(buf.s);
  kseq_destroy(ks);
  gzclose(fp);

  // FMD format
  mritr_t itr;
  const uint8_t *block;
  rld_t *e = 0;
  rlditr_t di;
  e = rld_init(6, 3);
  rld_itr_init(e, &di, 0);
  mr_itr_first(mr, &itr, 1);
  while ((block = mr_itr_next_block(&itr)) != 0) {
    const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
    while (q < end) {
      int c = 0;
      int64_t l;
      rle_dec1(q, c, l);
      rld_enc(e, &di, l, c);
    }
  }
  rld_enc_finish(e, &di);
  char *index_path = (char *)malloc(strlen(fa_path) + 4);
  strcpy(index_path, fa_path);
  strcat(index_path, ".fmd");
  rld_dump(e, index_path);

  mr_destroy(mr);

  return 0;
}

rldintv_t backward_search(rld_t *index, const uint8_t *P, int p2) {
  rldintv_t sai; // rldintv_t is the struct used to store a SA interval.
  fm6_set_intv(index, P[p2], sai);
  // std::cerr << p2 << " " << int2char[P[p2]] << std::endl;
  while (sai.x[2] != 0 && p2 > 0) {
    --p2;
    // std::cerr << p2 << " " << int2char[P[p2]] << std::endl;
    rldintv_t osai[6];
    rld_extend(index, &sai, osai, 1); // 1: backward, 0: forward
    sai = osai[P[p2]];
  }
  return sai;
}

int main_fmdexact(int argc, char *argv[]) {
  char *index_path = argv[1];
  char *query_path = argv[2];

  // spdlog::info("Restoring index..");
  rld_t *index = rld_restore(index_path);

  // auto range = rlcsa->getBWTRange();
  // std::cerr << range.first << "," << range.second << std::endl;

  gzFile fp = gzopen(query_path, "r");
  kseq_t *ks = kseq_init(fp);
  int l;
  int i;
  uint8_t *s;
  rldintv_t range;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;
    // change encoding
    for (i = 0; i < l; ++i)
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;

    range = backward_search(index, s, l - 1);
    std::cout << ks->name.s << ": " << range.x[0] << ","
              << range.x[0] + range.x[2] - 1 << std::endl;
  }
  kseq_destroy(ks);
  gzclose(fp);
  return 0;
}
