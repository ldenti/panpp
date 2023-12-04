#include "rlcsa.h"

namespace CSA {

/*
This function merges two vectors using marked positions.
The original vectors are deleted.
*/

sdsl::rle_vector<64> *mergeVectors(sdsl::rle_vector<64> *first,
                                   sdsl::rle_vector<64> *second,
                                   usint *positions, usint n, usint size,
                                   usint block_size) {
  if ((first == 0 && second == 0) || positions == 0) {
    return 0;
  }

  // std::cerr << n << " " << size << std::endl;

  sdsl::rle_vector_builder<64> encoder(size);
  sdsl::rle_vector<64> ff = *first;
  sdsl::rle_vector<64> ss = *second;

  // for (int i = 0; i < ff.size(); ++i)
  //   std::cerr << ff[i];
  // std::cerr << std::endl;

  // for (int i = 0; i < ss.size(); ++i)
  //   std::cerr << ss[i];
  // std::cerr << std::endl;

  uint64_t curr_pos_on_f = 0;
  uint64_t curr_pos_on_s = 0;
  uint64_t curr_run_value = 0;
  uint64_t curr_run_len = 0;
  uint64_t curr_position = 0;
  if (positions[0] == 0) {
    // First from second
    curr_run_value = ss[0];
    ++curr_pos_on_s;
    ++curr_position;
  } else {
    // First from first
    curr_run_value = ff[0];
    ++curr_pos_on_f;
  }
  curr_run_len = 1;

  int curr_value = 0;
  uint64_t curr_run_start = 0;
  for (uint64_t p = 1; p < size; ++p) {
    if (p == positions[curr_position]) {
      // from second
      curr_value = ss[curr_pos_on_s];
      ++curr_position;
      ++curr_pos_on_s;
    } else {
      curr_value = ff[curr_pos_on_f];
      ++curr_pos_on_f;
    }
    if (curr_run_value != curr_value) {
      if (curr_run_value == 0) {
        curr_run_start = p;
        curr_run_value = curr_value;
        curr_run_len = 1;
      } else {
        // store run
        encoder.set(curr_run_start, curr_run_len);
        curr_run_value = curr_value;
        curr_run_len = 1;
      }
    } else {
      ++curr_run_len;
    }
  }
  if (curr_run_value == 1) {
    // store run
    encoder.set(curr_run_start, curr_run_len);
  }

  sdsl::rle_vector<64> *x = new sdsl::rle_vector<64>(encoder);
  // sdsl::rle_vector<64> xx = *x;
  // for (int i = 0; i < xx.size(); ++i)
  //   std::cerr << xx[i];
  // std::cerr << std::endl;

  return x;
}

RLCSA::RLCSA(const std::string &base_name)
    : ok(false), alphabet(0), end_points(0) {
  for (usint c = 0; c < CHARS; ++c)
    this->array[c] = 0;

  std::string in_name = base_name + ".ab";
  std::ifstream infile(in_name.c_str(), std::ios_base::binary);
  // if (!infile) {
  //   std::cerr << "RLCSA: Error opening Psi array file!" << std::endl;
  //   return;
  // }

  usint distribution[CHARS];
  infile.read((char *)distribution, CHARS * sizeof(usint));
  this->alphabet = new Alphabet(distribution);
  this->data_size = this->alphabet->getDataSize();
  infile.close();

  for (usint c = 0; c < CHARS; c++) {
    if (this->alphabet->hasChar(c)) {
      in_name = base_name + ".array" + std::to_string(c);
      infile.open(in_name.c_str(), std::ios_base::binary);
      this->array[c] = new sdsl::rle_vector<64>();
      this->array[c]->load(infile);
      infile.close();

      this->ranks_supp[c] = new sdsl::rank_support_rle<1, 64>(this->array[c]);
      // FIXME: build rank also when building RLCSA, not
      // only when loading from file
    }
  }

  in_name = base_name + ".ep";
  infile.open(in_name.c_str(), std::ios_base::binary);
  this->end_points = new DeltaVector(infile);
  this->number_of_sequences = this->end_points->getNumberOfItems();
  infile.close();

  this->ok = true;
}

RLCSA::RLCSA(uchar *data, usint bytes, usint block_size, usint threads)
    : ok(false), alphabet(0), end_points(0) {
  for (usint c = 0; c < CHARS; ++c)
    this->array[c] = 0;

  this->buildRLCSA(data, 0, bytes, block_size, threads);
}

RLCSA::RLCSA(RLCSA &index, RLCSA &increment, usint *positions, usint block_size,
             usint threads)
    : ok(false), alphabet(0), end_points(0) {
  for (usint c = 0; c < CHARS; c++) {
    this->array[c] = 0;
  }

  if (!index.isOk() || !increment.isOk()) {
    return; // Fail silently. Actual error has already been reported.
  }
  if (positions == 0) {
    std::cerr << "RLCSA: Positions for insertions not available!" << std::endl;
    return;
  }

  index.strip();
  increment.strip();

  // Build character tables etc.
  usint distribution[CHARS];
  for (usint c = 0; c < CHARS; c++) {
    distribution[c] =
        index.alphabet->countOf(c) + increment.alphabet->countOf(c);
  }
  this->alphabet = new Alphabet(distribution);
  this->data_size = this->alphabet->getDataSize();
  this->number_of_sequences =
      index.number_of_sequences + increment.number_of_sequences;

  // Merge end points, SA samples, and Psi.
  usint psi_size = this->data_size + this->number_of_sequences;
  bool should_be_ok = true;

  omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic, 1)
  for (int c = -1; c < (int)CHARS; c++) {
    if (c == -1) {
      this->mergeEndPoints(index, increment);
    } else if (this->alphabet->hasChar(c) != 0) {
      std::cerr << "Merging " << c << ": "
                << index.data_size + index.number_of_sequences << "+"
                << increment.data_size + increment.number_of_sequences << " = "
                << psi_size << std::endl;
      this->array[c] =
          mergeVectors(index.array[c], increment.array[c], positions,
                       increment.data_size + increment.number_of_sequences,
                       psi_size, block_size);
      index.array[c] = 0;
      increment.array[c] = 0;

      if (this->array[c] == 0) {
        std::cerr << "RLCSA: Merge failed for vectors " << c << "!"
                  << std::endl;
        should_be_ok = false;
      }
      this->ranks_supp[c] = new sdsl::rank_support_rle<1, 64>(this->array[c]);
    }
  }

  this->ok = should_be_ok;
}

RLCSA::~RLCSA() {
  for (usint c = 0; c < CHARS; ++c) {
    delete this->array[c];
    this->array[c] = 0;
  }
  delete this->alphabet;
  this->alphabet = 0;
  delete this->end_points;
  this->end_points = 0;
}

//--------------------------------------------------------------------------

void RLCSA::writeTo(const std::string &base_name) const {
  std::string out_name = base_name + ".ab";
  std::ofstream ofile(out_name.c_str(), std::ios_base::binary);
  //   if (!ofile) {
  //     std::cerr << "RLCSA: Error creating Psi array file!" << std::endl;
  //     return;
  //   }
  this->alphabet->writeTo(ofile);
  ofile.close();

  out_name = base_name + ".ep";
  ofile.open(out_name.c_str(), std::ios_base::binary);
  //   if (!ofile) {
  //     std::cerr << "RLCSA: Error creating Psi array file!" << std::endl;
  //     return;
  //   }
  this->end_points->writeTo(ofile);
  ofile.close();

  for (usint c = 0; c < CHARS; ++c) {
    if (this->array[c] != 0) {
      out_name = base_name + ".array" + std::to_string(c);
      ofile.open(out_name.c_str(), std::ios_base::binary);
      this->array[c]->serialize(ofile);
      ofile.close();
    }
  }
}

//--------------------------------------------------------------------------

pair_type RLCSA::count(const std::string &pattern) const {
  if (pattern.length() == 0) {
    return this->getSARange();
  }

  // for (int b = 0; b < 6; ++b)
  //   std::cerr << b << ": " << this->alphabet->cumulative(b) << std::endl;

  std::string::const_reverse_iterator iter = pattern.rbegin();
  pair_type index_range = this->getCharRange((uchar)*iter);

  if (isEmpty(index_range)) {
    return index_range;
  }

  for (++iter; iter != pattern.rend(); ++iter) {
    index_range = this->LF(index_range, (uchar)*iter);

    if (isEmpty(index_range)) {
      return EMPTY_PAIR;
    }
  }

  this->convertToSARange(index_range);

  return index_range;
}

//--------------------------------------------------------------------------

void RLCSA::reportPositions(uchar *data, usint length, usint *positions) const {
  if (data == 0 || length == 0 || positions == 0) {
    return;
  }

  usint current = this->number_of_sequences - 1;
  positions[length] = current; // "immediately after current"
  for (sint i = (sint)(length - 1); i >= 0; i--) {
    usint c = (usint)data[i];
    if (this->array[c] != 0) {
      current = this->LF(current, c);
    } else {
      if (c < this->alphabet->getFirstChar()) // No previous characters either.
      {
        current = this->number_of_sequences - 1;
      } else {
        current = this->alphabet->cumulative(c) - 1 + this->number_of_sequences;
      }
    }
    positions[i] = current; // "immediately after current"
  }
}

//--------------------------------------------------------------------------

pair_type RLCSA::getSARange() const {
  return pair_type(0, this->data_size - 1);
}

pair_type RLCSA::getCharRange(usint c) const {
  if (c >= CHARS) {
    return EMPTY_PAIR;
  }
  pair_type index_range = this->alphabet->getRange(c);
  this->convertToBWTRange(index_range);
  return index_range;
}

void RLCSA::convertToSARange(pair_type &bwt_range) const {
  bwt_range.first -= this->number_of_sequences;
  bwt_range.second -= this->number_of_sequences;
}

void RLCSA::convertToBWTRange(pair_type &sa_range) const {
  sa_range.first += this->number_of_sequences;
  sa_range.second += this->number_of_sequences;
}

pair_type RLCSA::LF(pair_type range, usint c) const {
  if (c >= CHARS || this->array[c] == 0) {
    return EMPTY_PAIR;
  }

  usint start = this->alphabet->cumulative(c) + this->number_of_sequences - 1;
  range.first = start + ranks_supp[c]->rank(range.first + 1 - 1) + 1;
  range.second = start + ranks_supp[c]->rank(range.second + 1);

  return range;
}

//--------------------------------------------------------------------------

uchar *RLCSA::readBWT() const {
  return this->readBWT(
      pair_type(0, this->data_size + this->number_of_sequences - 1));
}

uchar *RLCSA::readBWT(pair_type range) const {
  if (isEmpty(range) ||
      range.second >= this->data_size + this->number_of_sequences) {
    return 0;
  }

  usint n = length(range);

  uchar *bwt = new uchar[n];
  memset(bwt, 0, n);

  for (int p = 0; p < n; ++p) {
    for (usint c = 0; c < CHARS; ++c) {
      if (this->array[c] != 0) {
        sdsl::rle_vector<64> v = *this->array[c];
        if (v[p] == 1)
          bwt[p] = c;
      }
    }
  }

  return bwt;
}

// --------------------------------------------------------------------------

void RLCSA::mergeEndPoints(RLCSA &index, RLCSA &increment) {
  DeltaEncoder *endings = new DeltaEncoder(RLCSA::ENDPOINT_BLOCK_SIZE);

  DeltaVector::Iterator index_iter(*(index.end_points));
  DeltaVector::Iterator increment_iter(*(increment.end_points));

  endings->setBit(index_iter.select(0));
  for (usint i = 1; i < index.number_of_sequences; i++) {
    endings->setBit(index_iter.selectNext());
  }
  usint sum = index.end_points->getSize();
  delete index.end_points;
  index.end_points = 0;

  endings->setBit(sum + increment_iter.select(0));
  for (usint i = 1; i < increment.number_of_sequences; i++) {
    endings->setBit(sum + increment_iter.selectNext());
  }
  sum += increment.end_points->getSize();
  delete increment.end_points;
  increment.end_points = 0;

  this->end_points = new DeltaVector(*endings, sum);
  delete endings;
}

//--------------------------------------------------------------------------

void RLCSA::strip() {
  for (usint c = 0; c < CHARS; c++) {
    if (this->array[c] != 0) {
      // TODO: delete and set to 0
      // delete this->array[c];
      // delete this->ranks[c];
    }
  }
  this->end_points->strip();
}

//--------------------------------------------------------------------------

void RLCSA::buildRLCSA(uchar *data, usint *ranks, usint bytes, usint block_size,
                       usint threads) {

  threads = std::max(threads, (usint)1);
  omp_set_num_threads(threads);

  // Determine the number of sequences and mark their end points.
  DeltaEncoder endings(RLCSA::ENDPOINT_BLOCK_SIZE);

  this->number_of_sequences = 0;
  usint marker = 0;
  usint padding = 0, chars_encountered = 0;

  for (usint i = 0; i < bytes; i++) {
    if (data[i] == 0) {
      if (i == marker) {
        break;
      } // Empty sequence.
      this->number_of_sequences++;
      marker = i + 1;
      usint pos = chars_encountered + padding - 1;
      endings.setBit(pos);
      padding = ((pos + 1) / 1) * 1 - chars_encountered;
    } else {
      ++chars_encountered;
    }
  }

  if (this->number_of_sequences == 0 || marker != bytes) {
    std::cerr << "RLCSA: Collection must consist of 0-terminated nonempty "
                 "sequences !"
              << std::endl;
    return;
  }
  this->end_points = new DeltaVector(endings, chars_encountered + padding);

  // Build character tables etc.
  usint distribution[CHARS];
  for (usint c = 0; c < CHARS; c++) {
    distribution[c] = 0;
  }
  for (usint i = 0; i < bytes; i++) {
    distribution[(usint)data[i]]++;
  }
  distribution[0] = 0; // \0 is an end marker
  this->alphabet = new Alphabet(distribution);
  this->data_size = this->alphabet->getDataSize();

  // Build suffix array.
  short_pair *sa = 0;
  if (ranks == 0) {
    sa = simpleSuffixSort(data, bytes, this->number_of_sequences, threads);
  } else {
    sa = simpleSuffixSort(ranks, bytes, threads);
  }

// Build Psi.
#pragma omp parallel for schedule(static)
  for (usint i = 0; i < bytes; i++) {
    sa[i].first = sa[(sa[i].first + 1) % bytes].second;
  }

// Build RLCSA.
#pragma omp parallel for schedule(dynamic, 1)
  for (usint c = 0; c < CHARS; c++) {
    if (!(this->alphabet->hasChar(c))) {
      this->array[c] = 0;
      continue;
    }

    short_pair *curr =
        sa + this->alphabet->cumulative(c) + this->number_of_sequences;
    short_pair *limit = curr + this->alphabet->countOf(c);

    sdsl::rle_vector_builder<64> encoder(this->data_size +
                                         this->number_of_sequences);
    pair_type run((*curr).first, 1);
    ++curr;

    for (; curr < limit; ++curr) {
      if ((*curr).first == run.first + run.second) {
        run.second++;
      } else {
        encoder.set(run.first, run.second);
        run = pair_type((*curr).first, 1);
      }
    }
    encoder.set(run.first, run.second);
    // encoder.flush();

    this->array[c] = new sdsl::rle_vector<64>(encoder);
  }
  delete[] sa;

  this->ok = true;
}

//--------------------------------------------------------------------------

} // namespace CSA
