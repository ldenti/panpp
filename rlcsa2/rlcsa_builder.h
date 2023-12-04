#ifndef RLCSA_BUILDER_H
#define RLCSA_BUILDER_H

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <omp.h>

#include "misc/utils.h"
#include "rlcsa.h"

namespace CSA {

class RLCSABuilder {
public:
  // We can optionally specify a starting RLCSA, so we can start out with a
  // big index without having to insertFromFile it (which requires the
  // original file to still be around).
  RLCSABuilder(usint _block_size, usint _buffer_size, usint _threads = 1,
               RLCSA *_index = NULL);
  ~RLCSABuilder();

  // Use this if you have already built an index for the file.
  void insertFromFile(const std::string &base_name, uchar *sequence);

  // User must free the index. Builder no longer contains it.
  RLCSA *getRLCSA();

  bool isOk();

private:
  RLCSA *index;

  usint block_size;
  usint buffer_size;

  usint threads;

  uchar *buffer;
  usint chars;

  bool ok;

  void flush();
  void reset();

  void addRLCSA(RLCSA *increment, uchar *sequence, usint length);
  void setRLCSA(RLCSA *new_index);
  void mergeRLCSA(RLCSA *increment, usint *ranks, usint length);

  usint *getRanks(uchar *sequence, usint length,
                  std::vector<usint> &end_markers);

  // These are not allowed.
  RLCSABuilder();
  RLCSABuilder(const RLCSABuilder &);
  RLCSABuilder &operator=(const RLCSABuilder &);
};

} // namespace CSA

#endif // RLCSA_BUILDER_H
