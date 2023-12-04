#include "rlcsa_builder.h"

namespace CSA {

RLCSABuilder::RLCSABuilder(usint _block_size, usint _buffer_size,
                           usint _threads, RLCSA *_index)
    : block_size(_block_size), buffer_size(_buffer_size), threads(_threads),
      buffer(0) {
  this->reset();

  if (_index != NULL) {
    this->setRLCSA(_index);
  }
}

RLCSABuilder::~RLCSABuilder() {
  delete this->index;
  delete[] this->buffer;
}

//--------------------------------------------------------------------------

void RLCSABuilder::insertFromFile(const std::string &base_name, uchar *data) {
  if (!this->ok) {
    return;
  }

  omp_set_num_threads(this->threads);

  this->flush();

  RLCSA *increment = new RLCSA(base_name);
  usint data_size = increment->getSize() + increment->getNumberOfSequences();

  this->addRLCSA(increment, data, data_size);
}

//--------------------------------------------------------------------------

RLCSA *RLCSABuilder::getRLCSA() {
  if (this->chars > 0) {
    this->flush();
  }

  RLCSA *temp = this->index;
  this->reset();

  return temp;
}

// --------------------------------------------------------------------------

void RLCSABuilder::flush() {
  if (this->buffer == 0 || this->chars == 0) {
    return;
  }

  RLCSA *temp =
      new RLCSA(this->buffer, this->chars, this->block_size, this->threads);
  this->addRLCSA(temp, this->buffer, this->chars);

  this->chars = 0;
  this->buffer = new uchar[this->buffer_size];
}

void RLCSABuilder::reset() {
  this->index = 0;

  if (this->buffer_size != 0) {
    delete[] this->buffer;
    this->buffer = new uchar[this->buffer_size];
  }
  this->chars = 0;

  this->ok = true;
}

//--------------------------------------------------------------------------

void RLCSABuilder::addRLCSA(RLCSA *increment, uchar *sequence, usint length) {
  if (this->index == 0) {
    this->setRLCSA(increment);
    return;
  }

  std::vector<usint> end_markers;
  usint *ranks = this->getRanks(sequence, length, end_markers);

  parallelSort(ranks, ranks + length);
#pragma omp parallel for schedule(static)
  for (usint i = 0; i < length; i++) {
    ranks[i] += i + 1;
  }
  this->mergeRLCSA(increment, ranks, length);
}

void RLCSABuilder::setRLCSA(RLCSA *new_index) {
  this->index = new_index;
  this->ok &= this->index->isOk();
}

void RLCSABuilder::mergeRLCSA(RLCSA *increment, usint *ranks, usint length) {
  RLCSA *merged = new RLCSA(*(this->index), *increment, ranks, this->block_size,
                            this->threads);
  delete[] ranks;
  delete this->index;
  delete increment;
  this->index = merged;

  this->ok &= this->index->isOk();
}

// --------------------------------------------------------------------------

usint *RLCSABuilder::getRanks(uchar *sequence, usint length,
                              std::vector<usint> &end_markers) {

  usint sequences = 0;
  for (usint i = 0; i < length; i++) {
    if (sequence[i] == 0) {
      end_markers.push_back(i);
      sequences++;
    }
  }

  usint *ranks = new usint[length];

  usint chunk = std::max((usint)1, sequences / (8 * this->threads));
#pragma omp parallel for schedule(dynamic, chunk)
  for (usint i = 0; i < sequences; i++) {
    usint begin = (i > 0 ? end_markers[i - 1] + 1 : 0);
    this->index->reportPositions(sequence + begin, end_markers[i] - begin,
                                 ranks + begin);
  }

  this->index->strip();

  return ranks;
}

//--------------------------------------------------------------------------

} // namespace CSA
