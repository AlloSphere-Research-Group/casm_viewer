#ifndef PROCESSORS_HPP
#define PROCESSORS_HPP

#include <vector>

#include "al/io/al_File.hpp"
#include "tinc/DataScript.hpp"
#include "tinc/DeferredComputation.hpp"

using namespace al;
using namespace tinc;

class DiffGenerator : public DataScript {
public:
  // Parameters

  std::string id = "DiffGenerator";

  std::string datasetPath;
  std::string condition;
};

#endif // PROCESSORS_HPP
