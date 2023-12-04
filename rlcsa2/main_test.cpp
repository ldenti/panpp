#include <assert.h>
#include <cstdint>
#include <getopt.h>
#include <iostream>
#include <string>
#include <zlib.h>

#include "kseq.h"

#include "fmd.h"
#include "rlcsa.h"
#include "rlcsa_builder.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

/**
 * Return the "reverse" complement of a single character.
 */
#define fm6_comp(a) ((a) >= 1 && (a) <= 4 ? 5 - (a) : (a))

#define fm6_set_intv(e, c, ik)                                                 \
  ((ik).x[0] = (e)->cnt[(int)(c)],                                             \
   (ik).x[2] = (e)->cnt[(int)(c) + 1] - (e)->cnt[(int)(c)],                    \
   (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
    5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

static const std::string int2char = "$ACGTN";

static inline int kputsn(const char *p, int l, kstring_t *s) {
  if (s->l + l + 1 >= s->m) {
    char *tmp;
    s->m = s->l + l + 2;
    kroundup32(s->m);
    if ((tmp = (char *)realloc(s->s, s->m)))
      s->s = tmp;
    else
      return EOF;
  }
  memcpy(s->s + s->l, p, l);
  s->l += l;
  s->s[s->l] = 0;
  return l;
}

static void seq_char2nt6(int l, char *s) {
  for (int i = 0; i < l; ++i)
    s[i] =
        seq_nt6_table[(uint8_t)s[i]]; // s[i] < 128 ? seq_nt6_table[s[i]] : 5;
}

/**
 * Return true if a character is a valid DNA base, and false otherwise.
 */
static inline bool isBase(uint64_t input) { return input >= 1 && input <= 5; }

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

void concatenate_fa(const char *fa_path, unsigned char *data,
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
  }
  kseq_destroy(seq);
  gzclose(fp);
}

int main(int argc, char **argv) {
  (void)(argc); // suppress unused parameter warning

  char *fa_path = argv[1];
  char *query_path = argv[2];

  int block_size = 64;
  int threads = 1;

  uint64_t max_size = ((uint64_t)1 << 32) - 1;
  unsigned char *data =
      (unsigned char *)calloc(max_size, sizeof(unsigned char));
  unsigned char *data_r =
      (unsigned char *)calloc(max_size, sizeof(unsigned char));
  cerr << "Reading " << fa_path << ".." << endl;

  uint64_t size = get_size_fa(fa_path);
  assert(size < max_size);
  concatenate_fa(fa_path, data, data_r);
  // Debug
  for (uint i = 0; i < size; ++i)
    cerr << int2char[data[i]];
  cerr << endl;
  for (uint i = 0; i < size; ++i)
    cerr << int2char[data_r[i]];
  cerr << endl;

  cerr << "Indexing " << fa_path << ".." << endl;
  CSA::RLCSA *index = new CSA::RLCSA(data, size, block_size, threads);

  if (index != 0 && index->isOk()) {
    cerr << "Storing " << fa_path << ".." << endl;
    index->writeTo(fa_path); // TODO: store to tmp dir
  } else {
    cerr << "Index is corrupted. Halting." << endl;
    return 1;
  }
  delete index;

  cerr << "Indexing " << fa_path << " (R).." << endl;
  index = new CSA::RLCSA(data_r, size, block_size, threads);
  char *rev = (char *)malloc(strlen(fa_path) + 2);
  strcpy(rev, fa_path);
  strcat(rev, ".R");

  if (index != 0 && index->isOk()) {
    cerr << "Storing " << fa_path << " (R).." << endl;
    index->writeTo(rev); // TODO: store to tmp dir
  } else {
    cerr << "Index is corrupted. Halting." << endl;
    return 1;
  }
  delete index;

  cerr << "Loading " << fa_path << ".." << endl;
  const CSA::RLCSA *rlcsa = new CSA::RLCSA(fa_path);

  unsigned char *bwt = rlcsa->readBWT();
  size = rlcsa->getSize() + rlcsa->getNumberOfSequences();
  for (uint i = 0; i < size; ++i)
    if (bwt[i] == 0)
      std::cerr << '|';
    else
      std::cerr << int2char[bwt[i]];
  std::cerr << std::endl;

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

  CSA::RLCSABuilder builder(block_size, 0, threads, NULL);
  cerr << "Merging " << fa_path << ".." << endl;
  builder.insertFromFile(fa_path, data);
  cerr << "Merging " << rev << ".." << endl;
  builder.insertFromFile(rev, data_r);
  cerr << "Storing full index to "
       << "AAA"
       << ".." << endl;
  rlcsa = builder.getRLCSA();
  if (rlcsa != 0 && rlcsa->isOk())
    rlcsa->writeTo("AAA");

  fp = gzopen(query_path, "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    seq_char2nt6(l, seq->seq.s);
    range = rlcsa->count(seq->seq.s);
    std::cout << seq->name.s << ": "
              << range.first + rlcsa->getNumberOfSequences() << ","
              << range.second + rlcsa->getNumberOfSequences() << std::endl;
  }
  kseq_destroy(seq);
  gzclose(fp);

  std::cout << std::endl;
  const FMD *fmd = new FMD("AAA");
  bwt = fmd->readBWT();
  size = fmd->getSize() + fmd->getNumberOfSequences();
  for (uint i = 0; i < size; ++i)
    if (bwt[i] == 0)
      std::cerr << '|';
    else
      std::cerr << int2char[bwt[i]];
  std::cerr << std::endl;
  FMDPosition fmd_range;

  for (int i = 1; i < 5; ++i) {
    fmd_range = fmd->getCharPosition(i);
    std::cout << i << ": " << fmd_range.forward_start << ","
              << fmd_range.reverse_start << "," << fmd_range.getLength()
              << std::endl;
  }
  std::cout << std::endl;
  fp = gzopen(query_path, "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    seq_char2nt6(l, seq->seq.s);
    fmd_range = fmd->fmdCount(seq->seq.s);
    std::cout << seq->name.s << ": "
              << fmd_range.forward_start + fmd->getNumberOfSequences() << ","
              << fmd_range.reverse_start + fmd->getNumberOfSequences() << ","
              << fmd_range.getLength() << std::endl;
  }
  kseq_destroy(seq);
  gzclose(fp);

  cerr << "Done." << endl;

  return 0;
}
