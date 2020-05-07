#ifndef INCLUDE_DATASETMANAGER_HPP
#define INCLUDE_DATASETMANAGER_HPP

#include <map>
#include <string>
#include <vector>

#include "tinc/VASPReader.hpp"

#include "al/ui/al_Parameter.hpp"

#include "tinc/ParameterSpace.hpp"

#include "processors.hpp"
#include "slice.hpp"

#include "nlohmann/json.hpp"

using json = nlohmann::json;

class BoundingBox_ {
public:
  float minx, miny, minz;
  float maxx, maxy, maxz;
  inline void reset() {
    minx = FLT_MAX;
    maxx = -FLT_MAX;
    miny = FLT_MAX;
    maxy = -FLT_MAX;
    minz = FLT_MAX;
    maxz = -FLT_MAX;
  }
  inline void addPoint(al::Vec3f &pos) {
    if (pos.x > maxx) {
      maxx = pos.x;
    } else if (pos.x < minx) {
      minx = pos.x;
    }
    if (pos.y > maxy) {
      maxy = pos.y;
    } else if (pos.y < miny) {
      miny = pos.y;
    }
    if (pos.z > maxz) {
      maxz = pos.z;
    } else if (pos.z < minz) {
      minz = pos.z;
    }
  }
};

class DatasetManager {
public:
  bool mRunProcessors{false};

  VASPReader reader;
  std::mutex mDataLock;

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

  std::map<std::string, size_t> mCurrentLoadedIndeces;

  // Template and diffs for time history computation
  std::map<std::string, std::vector<float>> mTemplatePositions;
  std::vector<float> mEmptyTemplate;
  json mDiffs;

  std::vector<std::pair<Vec3f, Vec3f>>
      mHistory; // From ->to (first, second) of pair

  // External Processors
  CacheManager cacheManager;
  AtomLabelProcessor labelProcessor;
  GraphGenerator graphGenerator;
  TemplateGenerator templateGen;
  DiffGenerator diffGen;

  BufferManager<std::map<std::string, std::vector<float>>> positionBuffers{8};

  bool graphProcessing{false};

  // These are internal parameters to propagate data and triggering computation,
  // not for direct user control
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
  std::map<std::string, ParameterSpace *> mParameterSpaces;

  // Functions --------------
  DatasetManager();

  void setPythonBinary(std::string pythonBinaryPath);

  void setPythonScriptPath(std::string pythonScriptPath);

  std::string buildRootPath();

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
  void loadFromPOSCAR();
  void getAtomPositions();

  std::vector<std::string> getDataNames();

  void generateGraph(std::string xData, std::string yData,
                     std::string datasetId, bool multi);

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
