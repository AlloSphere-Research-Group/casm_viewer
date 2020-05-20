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

  std::map<std::string, std::string> configuration() {

    std::map<std::string, std::string> config;

    config["dataset_path"] = datasetPath;

    std::string prim_path;
    if (File::exists(datasetPath + "/prim_labels.json")) {
      prim_path = datasetPath + "/prim_labels.json";
    } else if (File::exists(datasetPath + "/prim.json")) {
      prim_path = datasetPath + "/prim.json";
    } else if (File::exists(datasetPath + "../prim_labels.json")) {
      prim_path = datasetPath + "../prim_labels.json";
    } else if (File::exists(datasetPath + "../prim.json")) {
      prim_path = datasetPath + "../prim.json";
    }
    config["prim_path"] = prim_path;
    config["condition"] = condition;

    config["__output_path"] = outputDirectory();

    return config;
  }
};

#endif // PROCESSORS_HPP
