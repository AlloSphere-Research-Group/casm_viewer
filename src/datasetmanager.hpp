#ifndef INCLUDE_DATASETMANAGER_HPP
#define INCLUDE_DATASETMANAGER_HPP

#include <map>
#include <string>
#include <vector>

#include "tinc/DataPoolJson.hpp"
#include "tinc/DiskBufferImage.hpp"
#include "tinc/ParameterSpace.hpp"
#include "tinc/ProcessorCpp.hpp"
#include "tinc/ProcessorGraph.hpp"
#include "tinc/VASPReader.hpp"

#include "al/math/al_Vec.hpp"
#include "al/ui/al_Parameter.hpp"

#include "nlohmann/json.hpp"

using json = nlohmann::json;

using namespace tinc;

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
  json mCurrentBasis;

  std::string metaText;
  std::string lastError;

  std::vector<std::pair<al::Vec3f, al::Vec3f>>
      mHistory; // From ->to (first, second) of pair

  // TINC Computation

  ProcessorGraph sampleProcessorGraph{"SampleComputation"};

  ProcessorScript trajectoryProcessor{"TrajectoryPreprocessor"};

  ProcessorScript labelProcessor{"AtomLabelProcessor"};
  ProcessorScript graphGenerator{"GraphGenerator"};

  DataPoolJson dataPool{"resultsData", mParameterSpace, "slices"};
  DataPoolJson trajectoriesPool{"trajectories", mParameterSpace, "slices"};
  DataPoolJson neighborhoodPool{"neighborhood", mParameterSpace, "slices"};

  // -----------------

  struct occupation_t {
    int8_t basis_index;
    int8_t occupancy_dof;
  };
  BufferManager<std::vector<occupation_t>> occupationData{8};

  static const size_t maxPercolationTypes = 10;
  BufferManager<std::vector<uint16_t>> percolationData[maxPercolationTypes];

  struct position_t {
    float x, y, z;
  };
  std::vector<position_t> templateData;

  typedef struct {
    std::vector<uint16_t> shell_sites;
    std::vector<uint8_t> occ_ref;
    std::vector<uint8_t> flag;
  } shell_site_t;

  std::vector<uint8_t> trajectoryData;
  std::map<std::string, shell_site_t> shellSiteMap;
  std::vector<int8_t> shellSiteTypes; // matches for current time step
  size_t numTimeSteps, numAtoms;

  // ----------------

  // These are internal parameters to propagate data and triggering
  // computation, not for direct user control
  //  ParameterString mRootPath{"rootPath"};
  al::ParameterString mCurrentDataset{"currentDataset"}; // sub directory
  al::ParameterString currentGraphName{"currentGraphName", "internal", ""};
  al::ParameterString currentPoscarName{"currentPoscarName", "internal", ""};
  al::ParameterBool processing{"processing", "internal", false};

  ParameterSpace mParameterSpace{"casmParameters"};

  // Plot axes
  al::ParameterMenu mPlotYAxis{"PlotYAxis"};
  al::ParameterMenu mPlotXAxis{"PlotXAxis"};

  //
  al::ParameterMenu mAtomOfInterest{"AtomOfInterest", "", 0};
  al::ParameterChoice mShellSiteTypes{"ShellSiteTypes"};
  al::ParameterChoice mPercolationTypes{"PercolationTypes"};

  // Functions --------------
  DatasetManager();

  void initializeComputation();

  void setPythonBinary(std::string pythonBinaryPath);

  void setPythonScriptPath(std::string pythonScriptPath);

  std::string getGlobalRootPath();
  std::string fullConditionPath();

  bool initDataset();
  bool analyzeDataset();

  // Run processing scripts across all the dataset prior to starting application
  // TODO we should run this asynchronously
  void preProcessDataset();

  bool valid();

  std::string currentDataset();

  std::string getSubDir();

  // Processing functions for atom positions

  void updateText();

  std::vector<std::string> getDataNames();

  void loadTrajectory();

  void loadShellSiteData();

  typedef std::map<std::string, std::vector<std::string>> SpeciesLabelMap;

  SpeciesLabelMap getAvailableSpecies();

  std::vector<std::pair<std::string, float>> getCurrentCompositions();
  std::vector<int8_t> getShellSiteTypes(size_t timeIndex);

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

  void readPercolation();

private:
};

#endif // INCLUDE_DATASETMANAGER_HPP
