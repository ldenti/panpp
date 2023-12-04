#ifndef RLCSA_H
#define RLCSA_H

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <omp.h>
#include <vector>

#include "alphabet.h"
#include "bits/deltavector.h"
#include "misc/utils.h"
#include "suffixarray.h"

#include "sdsl/rle_vector.hpp"

namespace CSA {

class RLCSA {
  friend class RLCSABuilder;

public:
  //--------------------------------------------------------------------------
  //  CONSTRUCTION
  //--------------------------------------------------------------------------

  static const usint ENDPOINT_BLOCK_SIZE = 16;

  explicit RLCSA(const std::string &base_name);

  /*
    Build RLCSA for multiple sequences, treating each \0 as an end marker.
    There must be nonzero characters between the \0s, and the last
    character must also be \0.
    FIXME Crashes if bytes >= 4 GB.
  */
  RLCSA(uchar *data, usint bytes, usint block_size, usint threads);

  // Destroys contents of index and increment.
  RLCSA(RLCSA &index, RLCSA &increment, usint *positions, usint block_size,
        usint threads = 1);
  ~RLCSA();

  void writeTo(const std::string &base_name) const;

  inline bool isOk() const { return this->ok; }

  //--------------------------------------------------------------------------
  //  QUERIES
  //--------------------------------------------------------------------------

  // These queries use SA ranges.

  // Returns the closed range containing the matches.
  pair_type count(const std::string &pattern) const;

  // Used when merging CSAs.
  void reportPositions(uchar *data, usint length, usint *positions) const;

  // Returns the BWT of the collection including end of sequence markers.
  uchar *readBWT() const;
  uchar *readBWT(pair_type range) const;

  //   //--------------------------------------------------------------------------
  //   //  SUPPORT FOR EXTERNAL MODULES: POSITIONS
  //   //--------------------------------------------------------------------------

  //   inline usint LF(usint sa_index, usint c) const {
  //     if (c >= CHARS) {
  //       return this->data_size + this->number_of_sequences;
  //     }
  //     if (this->array[c] == 0) {
  //       if (c < this->alphabet->getFirstChar()) {
  //         return this->number_of_sequences - 1;
  //       }
  //       return this->alphabet->cumulative(c) + this->number_of_sequences - 1;
  //     }
  //     this->convertToBWTIndex(sa_index);

  //     PsiVector::Iterator iter(*(this->array[c]));
  //     return this->LF(sa_index, c, iter);
  //   }

  //   //   inline void convertToSAIndex(usint &bwt_index) const {
  //   //     bwt_index -= this->number_of_sequences;
  //   //   }
  //   inline void convertToBWTIndex(usint &sa_index) const {
  //     sa_index += this->number_of_sequences;
  //   }

  //   //--------------------------------------------------------------------------
  //   //  SUPPORT FOR EXTERNAL MODULES: RANGES
  //   //--------------------------------------------------------------------------

  pair_type getSARange() const;
  // pair_type getBWTRange() const;
  pair_type getCharRange(usint c) const;

  void convertToBWTRange(pair_type &sa_range) const;
  void convertToSARange(pair_type &bwt_range) const;
  // void convertToSARange(std::vector<pair_type> &bwt_ranges) const;

  // This is an unsafe function that does not check its parameters.
  pair_type LF(pair_type bwt_range, usint c) const;

  //   //--------------------------------------------------------------------------
  //   //  REPORTING
  //   //--------------------------------------------------------------------------

  inline usint getSize() const { return this->data_size; }
  inline usint getTextSize() const { return this->end_points->getSize(); }
  inline usint getNumberOfSequences() const {
    return this->number_of_sequences;
  }
  //   //   inline usint getBlockSize() const {
  //   //     return
  //   this->array[this->alphabet->getFirstChar()]->getBlockSize();
  //   //   }

  //   //   // Returns the size of the data structure.
  //   //   usint reportSize(bool print = false) const;

  //   //   void printInfo() const;

  //--------------------------------------------------------------------------
  //  INTERNAL VARIABLES
  //--------------------------------------------------------------------------

protected:
  bool ok;
  usint data_size;

  sdsl::rle_vector<64> *array[CHARS];
  sdsl::rank_support_rle<1, 64> *ranks_supp[CHARS];

  Alphabet *alphabet;

  usint number_of_sequences;
  DeltaVector *end_points;

  //   //--------------------------------------------------------------------------
  //   //  INTERNAL VERSIONS OF BASIC OPERATIONS
  //   //--------------------------------------------------------------------------

  inline usint LF(usint bwt_index, usint c) const {
    return this->alphabet->cumulative(c) + this->number_of_sequences +
           this->ranks_supp[c]->rank(bwt_index + 1) - 1;
  }

  //   //--------------------------------------------------------------------------
  //   //  INTERNAL STUFF
  //   //--------------------------------------------------------------------------

  void mergeEndPoints(RLCSA &index, RLCSA &increment);

  void buildRLCSA(uchar *data, usint *ranks, usint bytes, usint block_size,
                  usint threads);

  // Removes structures not necessary for merging.
  void strip();

  // These are not allowed.
  RLCSA();
  RLCSA(const RLCSA &);
  RLCSA &operator=(const RLCSA &);
};

} // namespace CSA

#endif // RLCSA_H
