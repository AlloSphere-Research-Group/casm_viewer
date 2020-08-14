#ifndef INCLUDE_DATASETMANAGER_HPP
#define INCLUDE_DATASETMANAGER_HPP

#include <map>
#include <string>
#include <vector>

#include "tinc/ComputationChain.hpp"
#include "tinc/CppProcessor.hpp"
#include "tinc/ImageDiskBuffer.hpp"
#include "tinc/ParameterSpace.hpp"
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

  // Template and diffs for time history computation
  //  std::map<std::string, std::vector<float>> mTemplatePositions;
  //  std::vector<float> mEmptyTemplate;
  //  json mDiffs;

  std::vector<std::pair<Vec3f, Vec3f>>
      mHistory; // From ->to (first, second) of pair

  // TINC Computation

  ComputationChain sampleComputationChain{"SampleComputation"};

  ComputationChain atomPositionChain{"AtomPositionComputation"};

  ScriptProcessor trajectoryProcessor{"TrajectoryPreprocessor"};
  ScriptProcessor labelProcessor{"AtomLabelProcessor"};
  CppProcessor kmcOccupation{"KMCOccupationExtractor"};
  ScriptProcessor graphGenerator{"GraphGenerator"};

  // TINC Buffers.
  //  BufferManager<std::map<std::string, std::vector<float>>>
  //  positionBuffers{8};

  struct occupation_t {
    int8_t basis_index;
    int8_t occupancy_dof;
    //    float x, y, z;
  };
  BufferManager<std::vector<occupation_t>> occupationData{8};

  struct position_t {
    float x, y, z;
  };
  std::vector<position_t> templateData;

  typedef struct {
    std::vector<uint16_t> shell_sites;
    std::vector<uint8_t> occ_ref;
  } shell_site_t;

  std::vector<uint8_t> trajectoryData;
  std::map<std::string, shell_site_t> shellSiteMap;
  std::vector<int8_t> shellSiteTypes; // matches for current time step
  size_t numTimeSteps, numAtoms;

  // ----------------

  // These are internal parameters to propagate data and triggering
  // computation, not for direct user control
  //  ParameterString mRootPath{"rootPath"};
  ParameterString mCurrentDataset{"currentDataset"}; // sub directory
  ParameterString currentGraphName{"currentGraphName", "internal", ""};
  ParameterString currentPoscarName{"currentPoscarName", "internal", ""};
  ParameterBool processing{"processing", "internal", false};

  ParameterSpace mParameterSpace;

  // Plot axes
  ParameterMenu mPlotYAxis{"PlotYAxis"};
  ParameterMenu mPlotXAxis{"PlotXAxis"};

  ParameterMenu mAtomOfInterest{"AtomOfInterest", "", 0};
  ParameterChoice mShellSiteTypes{"ShellSiteTypes"};

  // Functions --------------
  DatasetManager();

  void initializeComputation();

  void setPythonBinary(std::string pythonBinaryPath);

  void setPythonScriptPath(std::string pythonScriptPath);

  std::string getGlobalRootPath();
  std::string fullConditionPath();

  void initDataset();
  void analyzeDataset();

  // Run processing scripts across all the dataset prior to starting application
  // TODO we should run this asynchronously
  void preProcessDataset();

  bool valid();

  std::string currentDataset();

  std::string getSubDir();

  // Processing functions for atom positions
  //  bool loadDiff(int timeIndex);
  //  void processTemplatePositions();
  //  void loadFromPOSCAR();
  void computeNewSample();

  void updateText();

  std::vector<std::string> getDataNames();

  void loadTrajectory();

  void loadShellSiteData();

  typedef std::map<std::string, std::vector<std::string>> SpeciesLabelMap;

  SpeciesLabelMap getAvailableSpecies();

  std::vector<std::pair<std::string, float>> getCurrentCompositions();
  std::vector<int8_t> getShellSiteTypes(size_t atomIndex);

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
