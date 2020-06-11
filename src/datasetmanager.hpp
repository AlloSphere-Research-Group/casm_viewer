#ifndef INCLUDE_DATASETMANAGER_HPP
#define INCLUDE_DATASETMANAGER_HPP

#include <map>
#include <string>
#include <vector>

#include "tinc/ComputationChain.hpp"
#include "tinc/ImageDiskBuffer.hpp"
#include "tinc/ParameterSpaceDimension.hpp"
#include "tinc/VASPReader.hpp"

#include "al/ui/al_Parameter.hpp"

#include "slice.hpp"

#include "nlohmann/json.hpp"

using json = nlohmann::json;

using namespace tinc;

class DatasetManager {
public:
  bool mRunProcessors{false};

  VASPReader reader;
  std::mutex mDataLock;

  bool mRunComputation{true};

  // Dataset metadata
  std::vector<std::string>
      mAvailableAtomsJson; // Atoms in result.json <comp(XX)> field
  std::vector<std::string>
      mVacancyAtoms; // Vacancy atoms from prim.json or prim_labels.json
  std::vector<std::string>
      mShowAtomElements; // Labels from prim.json or prim_labels.json

  // Dataset metadata
  std::map<std::string, std::pair<float, float>>
      mDataRanges; // Ranges of the data across all the dataset

  std::string mGlobalRoot;
  std::string mLoadedDataset;

  std::string mTitle;
  json mCurrentBasis;

  std::string metaText;

  std::map<std::string, size_t> mCurrentLoadedIndeces;

  // Template and diffs for time history computation
  std::map<std::string, std::vector<float>> mTemplatePositions;
  std::vector<float> mEmptyTemplate;
  json mDiffs;

  std::vector<std::pair<Vec3f, Vec3f>>
      mHistory; // From ->to (first, second) of pair

  // TINC Computation
  CacheManager cacheManager;

  ComputationChain sampleComputationChain{ComputationChain::PROCESS_ASYNC,
                                          "SampleComputation"};

  ComputationChain atomPositionChain{"AtomPositionComputation"};

  DataScript labelProcessor{"AtomLabelProcessor"};
  DataScript graphGenerator{"GraphGenerator"};
  DataScript decompressTrajectory{"DecompressTrajectory"};

  // TINC Buffers.
  //  BufferManager<std::map<std::string, std::vector<float>>>
  //  positionBuffers{8};

  struct s1 {
    int8_t basis_index;
    int8_t occupancy_dof;
    //    float x, y, z;
  };
  BufferManager<std::vector<s1>> occupationData{8};
  std::vector<float> templateData;

  // ----------------

  // These are internal parameters to propagate data and triggering
  // computation, not for direct user control
  ParameterString mRootPath{"rootPath"};
  ParameterString mCurrentDataset{"currentDataset"}; // sub directory
  ParameterString currentGraphName{"currentGraphName", "internal", ""};
  ParameterString currentPoscarName{"currentPoscarName", "internal", ""};
  ParameterBool processing{"processing", "internal", false};

  std::map<std::string, bool>
      mParameterIsVariable; // Parameter space has changes

  std::string
      mConditionsParameter; // Parameter that maps to conditions.X directories
  std::vector<std::string>
      mParameterForSubDir; // Parameter value determines the subdirectory
  std::map<std::string, ParameterSpaceDimension *> mParameterSpaces;

  // Plot axes
  ParameterMenu mPlotYAxis{"PlotYAxis"};
  ParameterMenu mPlotXAxis{"PlotXAxis"};

  ParameterMenu mAtomOfInterest{"AtomOfInterest", "", 0};

  // Functions --------------
  DatasetManager();

  void initializeComputation();

  void setPythonBinary(std::string pythonBinaryPath);

  void setPythonScriptPath(std::string pythonScriptPath);

  std::string buildRootPath();
  std::string fullConditionPath();

  void initRoot();

  void analyzeDataset();

  // Run processing scripts across all the dataset prior to starting application
  // TODO we should run this asynchronously
  void preProcessDataset();

  bool valid();

  std::string currentDataset();

  std::string getSubDir();

  // Processing functions for atom positions
  bool loadDiff(int timeIndex);
  void processTemplatePositions();
  //  void loadFromPOSCAR();
  void computeNewSample();

  void updateText() {
    // Meta data texts
    metaText = "Global Root: " + mGlobalRoot + "\n";

    metaText += "Root: " + mRootPath.get() + "\n";
    metaText += "Dataset: " + mCurrentDataset.get() + "\n";
    metaText += "Subdir: " + getSubDir() + "\n";
    metaText += "Condition Param: " + mConditionsParameter + " condition: " +
                std::to_string(
                    mParameterSpaces[mConditionsParameter]->getCurrentIndex()) +
                "\n";

    metaText += " ----- Parameters -----\n";
    for (auto param : mParameterSpaces) {
      if (param.second->size() > 2) {
        metaText += param.first + " : " + param.second->getCurrentId() + "\n";
      } else if (param.second->size() == 1) {
        // For spaces with a single value, the id will be ./ so show the value
        metaText += param.first + " : " +
                    std::to_string(param.second->getCurrentValue()) + "\n";
      }
    }
    metaText += " ----- Data -----\n";
    for (auto compData : getCurrentCompositions()) {
      metaText +=
          compData.first + " = " + std::to_string(compData.second) + "\n";
    }
    metaText += "Current POSCAR :";
    metaText += labelProcessor.outputFile();
  }

  std::vector<std::string> getDataNames();

  typedef std::map<std::string, std::vector<std::string>> SpeciesLabelMap;

  SpeciesLabelMap getAvailableSpecies();

  std::vector<std::pair<std::string, float>> getCurrentCompositions() {
    std::vector<std::pair<std::string, float>> compositions;
    std::string str = readJsonResultsFile(mCurrentDataset.get(), getSubDir());
    if (str.size() > 0) {
      // Read the results file. This tells us which atoms are available
      auto resultsJson = json::parse(str);

      std::string comp_prefix = "<comp_n(";
      for (json::iterator it = resultsJson.begin(); it != resultsJson.end();
           ++it) {
        std::string key = it.key();
        //                // Look for available comp_n atom names
        if (key.compare(0, comp_prefix.size(), comp_prefix) == 0) {
          if (mConditionsParameter != "") {
            auto index =
                mParameterSpaces[mConditionsParameter]->getCurrentIndex();
            if (it.value().size() > index) {
              compositions.push_back({key, it.value()[index].get<float>()});
            } else {
              std::cerr
                  << "Warning: conditions parameter index out of range in "
                     "json results"
                  << std::endl;
            }
          } else {
          }
        }
      }
    }
    return compositions;
  }

  // Processing functions ----------------------------------------------

  al::Vec3f findSlicingDir(std::vector<al::Vec3f> &elem_pos_as_vec3,
                           float sampleRadius) {
    al::Vec3f layerDir = find_slicing_dir(
        elem_pos_as_vec3, elem_pos_as_vec3.size() / 5, sampleRadius, 8);
    //        std::cout << "Layer Direction " << layerDir.get()<< std::endl;
    return layerDir;
  }

  void readParameterSpace();

protected:
  // Look in the inner directory, then go up to the data root to try to find
  std::string findJsonFile(std::string datasetId, std::string subDir,
                           std::string fileName);

  inline std::string readJsonFile(std::string datasetId, std::string subDir,
                                  std::vector<std::string> fileNames) {
    std::string str;

    std::string foundFileName;
    for (auto name : fileNames) {
      if ((foundFileName = findJsonFile(datasetId, subDir, name)).size() > 0) {
        break;
      }
    }

    if (foundFileName.size() == 0) {
      return "";
    }

    std::ifstream f(foundFileName);

    if (f.fail()) {
      return "";
    }

    f.seekg(0, std::ios::end);
    str.reserve(f.tellg());
    f.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(f)),
               std::istreambuf_iterator<char>());

    return str;
  }

  std::string readJsonResultsFile(std::string datasetId, std::string subDir) {
    return readJsonFile(datasetId, subDir, {"results.json"});
  }

  std::string readJsonPrimLabelsFile(std::string datasetId,
                                     std::string subDir) {
    return readJsonFile(datasetId, subDir, {"prim_labels.json", "prim.json"});
  }

private:
};

#endif // INCLUDE_DATASETMANAGER_HPP
