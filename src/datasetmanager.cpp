
#ifdef AL_WINDOWS
#define NOMINMAX
#include <Windows.h>
#undef far
#endif

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <regex>

#include "datasetmanager.hpp"

#include "tinc/DiskBufferNetCDFData.hpp"

using namespace al;

void replaceAll(std::string &s, const std::string &search,
                const std::string &replace) {
  for (size_t pos = 0;; pos += replace.length()) {
    // Locate the substring to replace
    pos = s.find(search, pos);
    if (pos == std::string::npos)
      break;
    // Replace by erasing and inserting
    s.erase(pos, search.length());
    s.insert(pos, replace);
  }
}

DatasetManager::DatasetManager() {
  // These parameters are used for data propagation, but shouldn't be set by the
  // user
  currentGraphName.setHint("hide", 1.0);
  currentPoscarName.setHint("hide", 1.0);

  initializeComputation();
}

void DatasetManager::initializeComputation() {
  // Configure parameter space
  mParameterSpace.parameterNameMap = {{"temperature", "T"},
                                      {"chempotA", "param_chem_pot(a)"},
                                      {"chempotB", "param_chem_pot(b)"},
                                      {"dir", "dir"}};

  // Parameter space callback
  mParameterSpace.onValueChange = [&](ParameterSpaceDimension *changedDimension,
                                      ParameterSpace *ps) {
    std::string condition;
    std::string folder;
    for (auto dim : ps->getDimensions()) {
      if (dim->getSpaceRepresentationType() == ParameterSpaceDimension::INDEX) {
        if (dim->getName() != "dir") {
          condition = std::to_string(dim->getCurrentIndex());
          break;
        } else {
          folder = dim->getCurrentId();
        }
      }
    }
    if (condition.size() == 0) {
      condition = "0";
    }

    if (changedDimension->getName() == "dir") {
      int32_t previousValue =
          dynamic_cast<Parameter *>(changedDimension->getParameterMeta())
              ->getPrevious();
      if (changedDimension->getCurrentId() !=
          changedDimension->idAt(
              changedDimension->getIndexForValue(previousValue))) {
        // Only reload if id has changed
        //      ps->stopSweep();

        loadTrajectory();
        loadShellSiteData();
        readPercolation();
        shellSiteTypes =
            getShellSiteTypes(ps->getDimension("time")->getCurrentIndex());

        mShellSiteTypes.set(shellSiteTypes);
        //      ps->sweep(sampleProcessorGraph, {"time"});
      }
    }

    sampleProcessorGraph.configuration["condition"] = condition;
    //    labelProcessor.configuration["dir"] = folder;
    if (ps->getDimension("time")) {
      sampleProcessorGraph.configuration["time"] =
          (int64_t)ps->getDimension("time")->getCurrentIndex();
    }

    // TODO is this necessary? Shouldn't this happen automatically?
    if (mParameterSpace.getDimension("time") &&
        mParameterSpace.getDimension("time")->size() > 0) {
      occupationData.doneWriting(occupationData.getWritable());

      for (auto &percData : percolationData) {
        percData.doneWriting(percData.getWritable());
      }
      shellSiteTypes =
          getShellSiteTypes(ps->getDimension("time")->getCurrentIndex());

      mShellSiteTypes.set(shellSiteTypes);
    } else {
      if (mRunProcessors) {
        labelProcessor.setInputDirectory(
            mParameterSpace.getRootPath() +
            mParameterSpace.getCurrentRelativeRunPath());
        mParameterSpace.runProcess(sampleProcessorGraph, {}, {});
        //      sampleProcessorGraph.process();
        if (graphGenerator.getOutputFileNames().size() > 0) {
          currentGraphName.set(graphGenerator.getOutputFileNames()[0]);
        }
        currentPoscarName.set(labelProcessor.getOutputDirectory() +
                              labelProcessor.getOutputFileNames()[0]);
      }
    }

    updateText();
  };

  trajectoryProcessor.setInputFileNames({"trajectory.json.gz"});
  trajectoryProcessor.setOutputFileNames({"trajectory.nc"});

  bool verbose = false;
  trajectoryProcessor.setVerbose(verbose);
  graphGenerator.setVerbose(verbose);
  labelProcessor.setVerbose(verbose);

  graphGenerator.setOutputFileNames({"graph.png"});
  graphGenerator.ignoreFail = true; // It doesn't matter if graph fails.

  graphGenerator.registerDependency(&mPlotYAxis);
  graphGenerator.registerDependency(&mPlotXAxis);

  graphGenerator.prepareFunction = [&]() {
    std::string datasetId = mCurrentDataset.get();
    std::string xLabel = mPlotXAxis.getCurrent();
    std::string yLabel = mPlotYAxis.getCurrent();

    std::vector<double> datax, datay;

    std::string highlightValue = "0.0";
    auto dim = mParameterSpace.getDimension(xLabel);
    if (dim) {
      size_t index = 0;
      while (index < dim->size()) {
        datax.push_back(dim->at(index));
        index += dim->getSpaceStride();
      }
      highlightValue = std::to_string(
          mParameterSpace.getDimension(xLabel)->getCurrentValue());
    } else {
      std::cerr << "Results file not suitable for graphing. Dimension "
                << xLabel << " not suitable." << std::endl;
      return false;
    }
    datay.resize(datax.size());
    auto ysliceSize =
        dataPool.readDataSlice(yLabel, xLabel, datay.data(), datay.size());

    if (datax.size() != datay.size()) {
      std::cerr << "ERROR graph axis size mismatch. " << xLabel << ":"
                << datax.size() << " " << yLabel << ":" << datay.size()
                << std::endl;
      return false;
    }

    if (mParameterSpace.parameterNameMap.find(xLabel) !=
        mParameterSpace.parameterNameMap.end()) {
      xLabel = mParameterSpace.parameterNameMap[xLabel];
    }

    std::string fullDatasetPath = File::conformPathToOS(
        getGlobalRootPath() + File::conformPathToOS(mCurrentDataset.get()));

    std::ofstream ofsx(fullDatasetPath + "cached_output/inx.bin",
                       std::ofstream::out | std::ios::binary);
    ofsx.write((char *)datax.data(), datax.size() * sizeof(double));
    ofsx.flush();
    ofsx.close();
    std::ofstream ofsy(fullDatasetPath + "cached_output/iny.bin",
                       std::ofstream::out | std::ios::binary);
    ofsy.write((char *)datay.data(), datay.size() * sizeof(double));
    ofsy.flush();
    ofsy.close();
    std::string subDir = getSubDir();
    // Check data range for whole dataset
    std::string str = readJsonResultsFile(datasetId, getSubDir());
    std::pair<float, float> dataRange(FLT_MAX, FLT_MIN);
    if (str.size() > 0) {
      auto resultsJson = json::parse(str);
      auto data = resultsJson[yLabel];
      for (json::iterator v = data.begin(); v != data.end(); v++) {
        if (v->get<float>() < dataRange.first) {
          dataRange.first = v->get<float>();
        } else if (v->get<float>() > dataRange.second) {
          dataRange.second = v->get<float>();
        }
      }
    }

    std::string yLabelSanitized = ProcessorScript::sanitizeName(yLabel);

    graphGenerator.configuration["title"] = mTitle;
    graphGenerator.configuration["dataset"] = datasetId;
    graphGenerator.configuration["temp_interest"] = highlightValue;
    graphGenerator.configuration["subdir"] = subDir;
    graphGenerator.configuration["xLabel"] = xLabel;
    graphGenerator.configuration["yLabel"] = yLabelSanitized;

    graphGenerator.configuration["miny"] = std::to_string(dataRange.first);
    graphGenerator.configuration["maxy"] = std::to_string(dataRange.second);
    graphGenerator.configuration["inxFile"] =
        std::string(fullDatasetPath + "cached_output/inx.bin");
    graphGenerator.configuration["inyFile"] =
        std::string(fullDatasetPath + "cached_output/iny.bin");

    return true;
  };

  // Connect processors into chains
  sampleProcessorGraph << labelProcessor << graphGenerator;

  // Configure parameter callbacks
  mPlotYAxis.registerChangeCallback([this](float value) {
    if (mPlotYAxis.get() != value) {
      mPlotYAxis.setNoCalls(value); // Force internal value, discard old value
      if (graphGenerator.process()) {
        if (graphGenerator.getOutputFileNames().size() > 0) {
          currentGraphName.set(graphGenerator.getOutputFileNames()[0]);
        }
      }
    }
  });

  mPlotXAxis.registerChangeCallback([this](float value) {
    if (mPlotXAxis.get() != value) {
      mPlotXAxis.setNoCalls(value); // Force internal value, discard old value
      if (graphGenerator.process()) {
        if (graphGenerator.getOutputFileNames().size() > 0) {
          currentGraphName.set(graphGenerator.getOutputFileNames()[0]);
        }
      }
    }
  });

  currentPoscarName.registerChangeCallback([&](std::string newName) {
    int ncid, retval;

    if ((retval = nc_open(newName.c_str(), NC_NOWRITE | NC_SHARE, &ncid))) {
      std::cerr << "Failed to open positions file: " << newName << std::endl;
      return /*false*/;
    }
    int varid;
    if ((retval = nc_inq_varid(ncid, "occupancy", &varid))) {
      return /*false*/;
    }

    nc_type xtypep;
    char name[32];
    int ndimsp;
    int dimidsp[32];
    int *nattsp = nullptr;
    if ((retval = nc_inq_var(ncid, varid, name, &xtypep, &ndimsp, dimidsp,
                             nattsp))) {
      return /*false*/;
    }

    size_t lenp;
    if ((retval = nc_inq_dimlen(ncid, dimidsp[0], &lenp))) {
      return /*false*/;
    }

    auto newPositionData = occupationData.getWritable();
    newPositionData->resize(lenp);

    if ((retval = nc_get_var(ncid, varid, newPositionData->data()))) {
      return /*false*/;
    }
    occupationData.doneWriting(newPositionData);

    if ((retval = nc_close(ncid))) {
      return /*false*/;
    }
  });

  // Datapools
  trajectoriesPool.registerDataFile("trajectories.nc", "time");
  trajectoriesPool.getAllPaths = [&](std::vector<std::string> fixedDimensions =
                                         std::vector<std::string>()) {
    std::vector<std::string> paths;
    for (auto id : mParameterSpace.getDimension("dir")->getSpaceIds()) {
      auto path = al::File::conformPathToOS(mParameterSpace.getRootPath()) +
                  al::File::conformPathToOS(id) + "/";
      for (auto ps : mParameterSpace.getDimensions()) {
        if (ps->getSpaceRepresentationType() ==
            ParameterSpaceDimension::INDEX) {
          std::string condition = std::to_string(ps->getCurrentIndex());
          path += "condition." + condition;
        }
      }
      if (path.size() > 0) {
        paths.push_back(path);
      }
    }
    return paths;
  };

  neighborhoodPool.registerDataFile("../cleared_path0.nc", "");
  neighborhoodPool.registerDataFile("../free_range_path0.nc", "");
  neighborhoodPool.getAllPaths = [&](std::vector<std::string> fixedDimensions =
                                         std::vector<std::string>()) {
    std::vector<std::string> paths;
    for (auto id : mParameterSpace.getDimension("dir")->getSpaceIds()) {
      auto path = al::File::conformPathToOS(mParameterSpace.getRootPath()) +
                  al::File::conformPathToOS(id) + "/";
      //          for (auto ps : mParameterSpace.getDimensions()) {
      //              if (ps->getSpaceType() ==
      //              ParameterSpaceDimension::INDEX)
      //              {
      //                  std::string condition =
      //                  std::to_string(ps->getCurrentIndex());
      //                  path += "condition." + condition;
      //              }
      //          }
      if (path.size() > 0) {
        paths.push_back(path);
      }
    }
    return paths;
  };

  mPercolationTypes.registerChangeCallback([&](auto val) {
    if (val != mPercolationTypes.get()) {
      for (auto &percData : percolationData) {
        percData.doneWriting(percData.getWritable());
      }
    }
  });
}

void DatasetManager::setPythonBinary(std::string pythonBinaryPath) {
  //      std::unique_lock<std::mutex> lk(mProcessingLock);
  graphGenerator.setCommand(pythonBinaryPath);
  labelProcessor.setCommand(pythonBinaryPath);
  trajectoryProcessor.setCommand(pythonBinaryPath);
}

void DatasetManager::setPythonScriptPath(std::string pythonScriptPath) {
  //      std::unique_lock<std::mutex> lk(mProcessingLock);
  graphGenerator.setScriptName(File::conformDirectory(pythonScriptPath) +
                               "graphing/plot.py");
  labelProcessor.setScriptName(File::conformDirectory(pythonScriptPath) +
                               "reassign_occs/reassign_occs.py");
  trajectoryProcessor.setScriptName(
      File::conformDirectory(pythonScriptPath) +
      "reassign_occs/extract_nc_from_trajectory_gzip.py");
}

std::string DatasetManager::getGlobalRootPath() {
  return File::conformPathToOS(mGlobalRoot);
}

std::string DatasetManager::fullConditionPath() {

  for (auto ps : mParameterSpace.getDimensions()) {
    if (ps->getSpaceRepresentationType() == ParameterSpaceDimension::INDEX) {
      std::string condition = std::to_string(ps->getCurrentIndex());
      return File::conformPathToOS(getGlobalRootPath() + mCurrentDataset.get() +
                                   "/" + getSubDir() + "/conditions." +
                                   condition + "/");
    }
  }
  return std::string();
}

void DatasetManager::readParameterSpace() {
  mParameterSpace.clear();
  mParameterSpace.setRootPath(File::conformPathToOS(
      getGlobalRootPath() + File::conformPathToOS(mCurrentDataset.get())));
  mParameterSpace.readFromNetCDF("parameter_space.nc");
  mParameterSpace.enableCache("ps_cache");
}

bool DatasetManager::initDataset() {
  bool ret = true;
  std::unique_lock<std::mutex> lk(mDataLock);
  mParameterSpace.stopSweep();

  // C:\Users\Andres\source\repos\casm_viewer\python\reassign_occs\extract_nc_from_trajectory_gzip.py

  std::string fullDatasetPath = File::conformPathToOS(
      getGlobalRootPath() + File::conformPathToOS(mCurrentDataset.get()));
  mLoadedDataset = fullDatasetPath;

  trajectoryProcessor.setInputDirectory(fullDatasetPath);
  trajectoryProcessor.process();

  if (!al::File::isDirectory(fullDatasetPath + "cached_output/")) {
    if (!al::Dir::make(fullDatasetPath + "cached_output/")) {
      std::cerr << "Error creating cache diretory: "
                << fullDatasetPath + "cached_output/" << std::endl;
    }
  }
  // Read template
  std::string templatePath = fullDatasetPath + "cached_output/template.nc";
  if (File::exists(templatePath)) {
    int retval, ncid, varid;
    if ((retval =
             nc_open(templatePath.c_str(), NC_NOWRITE | NC_SHARE, &ncid))) {
      lastError = "Error cached_output/opening template.nc";
      return false;
    }
    if ((retval = nc_inq_varid(ncid, "atoms_var", &varid))) {
      lastError = "Error getting 'atoms_var' from template.nc";
      return false;
    }

    nc_type xtypep;
    char name[32];
    int ndimsp;
    int dimidsp[32];
    int *nattsp = nullptr;
    if ((retval = nc_inq_var(ncid, varid, name, &xtypep, &ndimsp, dimidsp,
                             nattsp))) {

      lastError = "Error getting 'atoms_var' metadata from template.nc";
      return false;
    }

    size_t lenp;
    if ((retval = nc_inq_dimlen(ncid, dimidsp[0], &lenp))) {
      lastError = "Error getting 'atoms_var' length from template.nc";
      return false;
    }

    templateData.resize(lenp);

    if ((retval = nc_get_var(ncid, varid, templateData.data()))) {
      lastError = "Error getting 'atoms_var' data from template.nc";
      return false;
    }

    if ((retval = nc_close(ncid))) {
      lastError = "Error closing template.nc";
      return false;
    }
  }

  // ===
  {
    // FIXME hack to avoid loop that occurs from triggering these...
    labelProcessor.enabled = false;
    graphGenerator.enabled = false; // Graphing not working with kmc yet.
  }

  readParameterSpace();
  loadShellSiteData();
  readPercolation();
  ret &= analyzeDataset();
  loadTrajectory(); // The time parameter space is loaded here.

  std::string conditionDim;
  for (auto dim : mParameterSpace.getDimensions()) {
    if (dim->getSpaceRepresentationType() == ParameterSpaceDimension::INDEX &&
        dim->size() > 0) {
      conditionDim = dim->getName();
      std::cout << "Using dimension '" << dim->getName()
                << "' as condition dimension" << std::endl;
    }
  }

  std::vector<std::string> idDims;
  for (auto dim : mParameterSpace.getDimensions()) {
    if (dim->getSpaceRepresentationType() == ParameterSpaceDimension::ID) {
      idDims.push_back(dim->getName());
      std::cout << "Using dimension '" << dim->getName()
                << "' as subdir dimension" << std::endl;
      break;
    }
  }

  std::string pathTemplate;

  if (mParameterSpace.getDimension("dir") &&
      mParameterSpace.getDimension("dir")->size() > 0) {
    if (idDims.at(0) != mParameterSpace.getDimension("dir")->getName()) {
      pathTemplate = "%%dir%%";
    }
  }
  std::vector<std::string> multiIdDimensions;

  for (auto dim : mParameterSpace.getDimensions()) {
    if (dim->getCurrentIds().size() > 1) {
      multiIdDimensions.push_back(dim->getName());
    }
  }

  if (multiIdDimensions.size() > 0) {
    pathTemplate = "%%";
    for (auto dimName : multiIdDimensions) {
      pathTemplate += dimName + ",";
    }
    pathTemplate.resize(pathTemplate.size() - 1);
    pathTemplate += "%%conditions.%%" + conditionDim + "%%";
  } else if (idDims.size() > 0) {
    pathTemplate +=
        "%%" + idDims.at(0) + "%%conditions.%%" + conditionDim + "%%";
  } else {
    pathTemplate = "/conditions.%%" + conditionDim + "%%";
  }

  mParameterSpace.setCurrentPathTemplate(pathTemplate);

  std::vector<std::string> parameterSpaceNames;
  for (auto dim : mParameterSpace.getDimensions()) {
    if (dim->size() > 1) {
      parameterSpaceNames.push_back(dim->getName());
    }
  }

  mPlotXAxis.setElements(parameterSpaceNames);
  mPlotYAxis.setNoCalls(0);

  auto dataNames = getDataNames();
  mPlotYAxis.setElements(dataNames);
  std::string defaultYAxis = "<comp_n(" + mAtomOfInterest.getCurrent() + ")>";
  ptrdiff_t pos = std::find(dataNames.begin(), dataNames.end(), defaultYAxis) -
                  dataNames.begin();
  if (pos < (int)dataNames.size()) {
    mPlotYAxis.setNoCalls((int)pos);
  }

  // Configuration for processors
  graphGenerator.setRunningDirectory(fullDatasetPath);
  graphGenerator.setOutputDirectory(fullDatasetPath + "/cached_output");

  labelProcessor.setRunningDirectory(fullDatasetPath);
  labelProcessor.setOutputDirectory(fullDatasetPath + "/cached_output");
  labelProcessor.setOutputFileNames({"positions.nc"});
  labelProcessor.setInputFileNames({"final_state.json"});

  std::string globalRoot = getGlobalRootPath();

  std::vector<std::string> possiblePrims = {
      fullDatasetPath + "/prim_labels.json",
      fullDatasetPath + "/prim.json",
      globalRoot + "prim_labels.json",
      globalRoot + "prim.json",
      fullDatasetPath + "/../prim_labels.json",
      fullDatasetPath + "/../prim.json"};
  std::string prim_path;
  for (auto possiblePrim : possiblePrims) {
    if (File::exists(possiblePrim)) {
      prim_path = possiblePrim;
      break;
    }
  }
  labelProcessor.configuration["prim_path"] = prim_path;
  labelProcessor.configuration["dataset_path"] = fullDatasetPath;

  labelProcessor.configuration["template_path"] =
      File::conformPathToOS(fullDatasetPath + "/cached_output/template.nc");

  labelProcessor.setOutputDirectory(getGlobalRootPath() +
                                    mCurrentDataset.get() + "/cached_output");

  // If there is a time space, don't use labeling script. Data is ready
  if (mParameterSpace.getDimension("time") &&
      mParameterSpace.getDimension("time")->size() > 0) {
    labelProcessor.enabled = false;
    graphGenerator.enabled = false; // Graphing not working with kmc yet.
  } else {
    labelProcessor.enabled = true;
    graphGenerator.enabled = true;
  }
  for (auto ps : mParameterSpace.getDimensions()) {
    if (ps->getSpaceRepresentationType() == ParameterSpaceDimension::INDEX) {
      std::string condition = std::to_string(ps->getCurrentIndex());

      dataPool.registerDataFile("../results.json", ps->getName());
    }
  }
  dataPool.setCacheDirectory(fullDatasetPath + "slices");
  trajectoriesPool.setCacheDirectory(fullDatasetPath + "slices");
  neighborhoodPool.setCacheDirectory(fullDatasetPath + "slices");
  //  std::cout << fullDatasetPath + "slices" << std::endl;

  lk.unlock();
  return ret;
}

bool DatasetManager::analyzeDataset() {
  if (mCurrentDataset.get() == "") {
    std::cerr << "ERROR: No datasets found." << std::endl;
    return false;
  }
  std::string datasetId = mCurrentDataset.get();
  // mCurrentDataset might be more current than the internal parameter value
  // as this might be called from the parameter change callback, when the
  // internal value has not yet been updated

  mDataRanges.clear();
  mAvailableAtomsJson.clear();
  mVacancyAtoms.clear();
  mShowAtomElements.clear();

  std::string str = readJsonResultsFile(datasetId, getSubDir());
  if (str.size() > 0) {
    // Read the results file
    auto resultsJson = json::parse(str);

    for (json::iterator it = resultsJson.begin(); it != resultsJson.end();
         ++it) {
      std::string key = it.key();
      if (mDataRanges.find(key) == mDataRanges.end()) {
        mDataRanges[key] = std::pair<float, float>(FLT_MAX, FLT_MIN);
      }
    }
  }

  // Read prim_labels file. This tells us the atoms, the vacancies and the
  // labeling
  std::string primFileContents = readJsonPrimLabelsFile(datasetId, "");

  if (primFileContents.size() > 0) {
    auto primLabelsJson = json::parse(primFileContents);
    mCurrentBasis = primLabelsJson["basis"];
    for (auto basis : primLabelsJson["basis"]) {
      //                std::cout << basis.dump() <<std::endl;
      bool isVacancy = false;
      for (auto label : basis["occupant_dof"]) {
        std::string labelString = label;
        if (labelString == "Va") {
          isVacancy = true;
          continue; // The next label should be the atom name
        }
        mAvailableAtomsJson.push_back(labelString);
        //                  if (mDataRanges.find(labelString) ==
        //                  mDataRanges.end()) {
        //                      mDataRanges[labelString] = std::pair<float,
        //                      float>(FLT_MAX, FLT_MIN);
        //                  }
        if (std::find(mShowAtomElements.begin(), mShowAtomElements.end(),
                      labelString) == mShowAtomElements.end()) {
          mShowAtomElements.push_back(labelString);
        }
        if (isVacancy && std::find(mVacancyAtoms.begin(), mVacancyAtoms.end(),
                                   labelString) == mVacancyAtoms.end()) {
          mVacancyAtoms.push_back(labelString);
        }
      }
    }
    mTitle = primLabelsJson["title"].get<std::string>();

    // Match atom of interest to atom that can fill vacancies
    std::smatch match;
    std::regex atomNameRegex("[A-Z][a-z]*");
    auto availableSpecies = getAvailableSpecies();

    // Populate mAtomOfInterest
    std::vector<std::string> species;
    for (auto atomName : availableSpecies) {
      if (atomName.first != "Va") {
        // Use only first result
        if (std::regex_search(atomName.first, match, atomNameRegex)) {
          std::string result = match.str();
          if (std::find(species.begin(), species.end(), result) ==
              species.end()) {
            species.push_back(result);
          }
        }
      }
    }
    mAtomOfInterest.setElements(species);
    // Now set atom of interest
    for (auto atomName : availableSpecies) {
      if (atomName.first == "Va") {
        for (auto atom : atomName.second) {
          if (std::regex_search(atom, match, atomNameRegex)) {
            std::string result = match.str();
            auto availableAtoms = mAtomOfInterest.getElements();
            if (std::find(availableAtoms.begin(), availableAtoms.end(),
                          result) != availableAtoms.end()) {
              mAtomOfInterest.setCurrent(result, true);
              // Use only first result
              break;
            }
          }
        }
      }
    }
  } else {

    lastError = "failed to find prim labels file";
    std::cerr << "ERROR failed to find prim labels file" << std::endl;
    return false;
  }
  return true;
}

void DatasetManager::preProcessDataset() {
  // TODO implement preprocess for datasets
  std::string datasetId = mCurrentDataset.get();

#ifdef AL_BUILD_MPI
  int world_size, nameSize, proc_id;
  char processor_name[1024];

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Get_processor_name(processor_name, &nameSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

  vector<unsigned int> indeces;
  int numProcessesSlice;
  if (proc_id == 0) {
    indeces.reserve(mParameterSpaces["temperature"]->size());
    for (unsigned int i = 0; i < mParameterSpaces["temperature"]->size(); i++) {
      indeces.push_back(i);
    }
    numProcessesSlice = indeces.size() / world_size;
  }

  MPI_Bcast(&numProcessesSlice, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  unsigned int *nodeIndeces;
  nodeIndeces =
      (unsigned int *)malloc(numProcessesSlice * sizeof(unsigned int));

  MPI_Scatter(indeces.data(), numProcessesSlice, MPI::UNSIGNED, nodeIndeces,
              numProcessesSlice, MPI::UNSIGNED, 0, MPI_COMM_WORLD);
  std::cout << processor_name << "--" << numProcessesSlice
            << " node indeces start:::  " << nodeIndeces[0] << std::endl;
//        if (mParameterSpaces["chempot2"]->size() > 0) {
//            for (auto chempot: mParameterSpaces["chempotA"]->values()) {
//                for (auto chempot2:
//                mParameterSpaces["chempotB"]->values())
//                {
//                    for (int i = 0; i < numProcessesSlice; i++) {
//                        int temp = nodeIndeces[i];
//                        labelProcessor.setParams(chempot.first,
//                        chempot2.first, std::to_string(temp),
//                        datasetId); labelProcessor.processAsync();
//                    }
//                }
//            }
//            labelProcessor.waitForAsyncDone();
//        } else {
//            for (auto chempot: mParameterSpaces["chempotA"]->values()) {
//                for (int i = 0; i < numProcessesSlice; i++) {
//                    int temp = nodeIndeces[i];
//                    labelProcessor.setParams(chempot.first, "",
//                    std::to_string(temp), datasetId);
//                    labelProcessor.processAsync();
//                }
//            }
//        }
#else
  std::cout << "preProcessDataset() not supported without MPI" << std::endl;

#endif
}

bool DatasetManager::valid() {
  // TODO have a better validity check
  bool valid = mParameterSpace.getDimension("chempotA") &&
               mParameterSpace.getDimension("temperature");
  return valid;
}

std::string DatasetManager::currentDataset() { return mLoadedDataset; }

std::string DatasetManager::getSubDir() {

  bool isMultiId = false;

  for (auto dim : mParameterSpace.getDimensions()) {
    if (dim->getSpaceStride() > 1) {
      isMultiId = true;
      break;
    }
  }
  if (isMultiId) {
    return mParameterSpace.getCommonId();

  } else {
    for (auto ps : mParameterSpace.getDimensions()) {
      if (ps->getSpaceRepresentationType() == ParameterSpaceDimension::ID) {
        return ps->getCurrentId();
      }
    }
  }
  return "";
}

void DatasetManager::updateText() {
  // Meta data texts
  metaText = "Global Root: " + mGlobalRoot + "\n";
  metaText += "Dataset: " + mCurrentDataset.get() + "\n";
  metaText += "Subdir: " + getSubDir() + "\n";
  metaText += " ----- Data -----\n";
  for (auto compData : getCurrentCompositions()) {
    metaText += compData.first + " = " + std::to_string(compData.second) + "\n";
  }
}

std::vector<std::string> DatasetManager::getDataNames() {
  std::vector<std::string> names;
  for (auto entry : mDataRanges) {
    names.push_back(entry.first);
  }
  return names;
}

void DatasetManager::loadTrajectory() {
  auto trajectoryFile = fullConditionPath() + "trajectory.nc";
  //  std::cout << "trying " << trajectoryFile;
  if (File::exists(trajectoryFile)) {

    int retval, ncid, varid;
    if ((retval =
             nc_open(trajectoryFile.c_str(), NC_NOWRITE | NC_SHARE, &ncid))) {
      return /*false*/;
    }
    if ((retval = nc_inq_varid(ncid, "occupation_dofs", &varid))) {
      return /*false*/;
    }

    nc_type xtypep;
    char name[32];
    int ndimsp;
    int dimidsp[32];
    int *nattsp = nullptr;
    if ((retval = nc_inq_var(ncid, varid, name, &xtypep, &ndimsp, dimidsp,
                             nattsp))) {
      return /*false*/;
    }

    if ((retval = nc_inq_dimlen(ncid, dimidsp[0], &numTimeSteps))) {
      return /*false*/;
    }

    if ((retval = nc_inq_dimlen(ncid, dimidsp[1], &numAtoms))) {
      return /*false*/;
    }

    if (numAtoms != templateData.size()) {
      std::cerr << "ERROR trajectory data atoms mismatch" << std::endl;
    }

    trajectoryData.resize(numTimeSteps * numAtoms);

    if ((retval = nc_get_var(ncid, varid, trajectoryData.data()))) {
      return /*false*/;
    }

    // Read time steps -------------------------
    if ((retval = nc_inq_varid(ncid, "steps", &varid))) {
      return /*false*/;
    }

    if ((retval = nc_inq_var(ncid, varid, name, &xtypep, &ndimsp, dimidsp,
                             nattsp))) {
      return /*false*/;
    }
    size_t lenp;
    if ((retval = nc_inq_dimlen(ncid, dimidsp[0], &lenp))) {
      return /*false*/;
    }
    std::vector<uint32_t> timeValues;
    timeValues.resize(lenp);
    if ((retval = nc_get_var_uint(ncid, varid, timeValues.data()))) {
      return /*false*/;
    }

    //     Ensure organized
    int offset = 0;
    for (size_t i = 1; i < timeValues.size(); i++) {
      if (timeValues[i] < timeValues[i - 1]) {
        timeValues[i] += offset;
      }
      if (timeValues[i] < timeValues[i - 1]) {
        timeValues[i] -= offset;
        offset = timeValues[i - 1];
        timeValues[i] += offset;
        if (timeValues[i] == timeValues[i - 1]) {
          offset += 1;
          timeValues[i] += 1;
        }
      }
    }

    assert(mParameterSpace.getDimension("time"));
    mParameterSpace.getDimension("time")->clear();
    mParameterSpace.getDimension("time")->setSpaceValues(timeValues.data(),
                                                         timeValues.size());
    mParameterSpace.getDimension("time")->setSpaceRepresentationType(
        ParameterSpaceDimension::VALUE);
    mParameterSpace.getDimension("time")->conformSpace();

    if (numTimeSteps != mParameterSpace.getDimension("time")->size()) {
      std::cout << "ERROR: Time dimension mismatch!" << std::endl;
    }

    if ((retval = nc_close(ncid))) {
      return /*false*/;
    }
  } else {
    std::cout << "No trajectory file found" << std::endl;
    if (mParameterSpace.getDimension("time")) {
      mParameterSpace.getDimension("time")->clear();
      mParameterSpace.getDimension("time")->conformSpace();
    }
    trajectoryData.clear();
    auto w = occupationData.getWritable();
    w->clear();
    occupationData.doneWriting(w);
  }
}

void DatasetManager::readPercolation() {
  auto datasetPath =
      getGlobalRootPath() + mCurrentDataset.get() + "/" + getSubDir();
  auto percoFiles = datasetPath + "perco_files.json";
  std::ifstream percoFileStream;
  percoFileStream.open(percoFiles, std::ofstream::in);
  std::vector<std::string> siteTypes;
  if (percoFileStream.good()) {
    auto j = json::parse(percoFileStream);
    if (j.is_object()) {
      size_t counter = 0;
      for (auto filename : j["files"]) {
        int ncid, retval;
        int sizeid;
        size_t size;
        if ((retval =
                 nc_open((datasetPath + filename.get<std::string>()).c_str(),
                         NC_NOWRITE | NC_SHARE, &ncid))) {
          std::cerr << "Error opening percolation file: "
                    << datasetPath + filename.get<std::string>() << std::endl;
          return;
        }

        if (counter >= maxPercolationTypes) {
          std::cerr << "ERROR: Only " << maxPercolationTypes
                    << " percolation buffers available. Discarding: "
                    << filename.get<std::string>() << std::endl;
          continue;
        }
        siteTypes.push_back(filename.get<std::string>());
        auto perco_data = percolationData[counter].getWritable();
        if ((retval = nc_inq_dimid(ncid, "occupation_size", &sizeid))) {
          return;
        }
        if ((retval = nc_inq_dimlen(ncid, sizeid, &size))) {
          return;
        }

        perco_data->resize(size);
        int shell_sites_id;
        if ((retval = nc_inq_varid(ncid, "occupation_dof", &shell_sites_id))) {
          return;
        }

        if ((retval = nc_get_var(ncid, shell_sites_id, perco_data->data()))) {
          return;
        }
        nc_close(ncid);
        percolationData[counter].doneWriting(perco_data);

        counter++;
      }
    }
  }

  mPercolationTypes.setElements(siteTypes);
  mPercolationTypes.setNoCalls(0);
}

void DatasetManager::loadShellSiteData() {
  auto datasetPath =
      getGlobalRootPath() + mCurrentDataset.get() + "/" + getSubDir();
  auto shellSiteFiles = datasetPath + "shell_site_files.json";
  std::ifstream shellSitesFileStream;
  shellSitesFileStream.open(shellSiteFiles, std::ofstream::in);
  shellSiteMap.clear();
  std::vector<std::string> siteTypes;
  if (shellSitesFileStream.good()) {
    auto j = json::parse(shellSitesFileStream);
    if (j.is_object()) {
      for (auto filename : j["files"]) {

        int ncid, retval;
        int stepsid, shellsizeid;
        size_t steps, shellsize;
        if ((retval =
                 nc_open((datasetPath + filename.get<std::string>()).c_str(),
                         NC_NOWRITE | NC_SHARE, &ncid))) {
          std::cerr << "Error opening file: " << filename << std::endl;
          continue;
        }
        siteTypes.push_back(filename.get<std::string>());
        shellSiteMap[filename.get<std::string>()] = {};
        std::vector<uint16_t> &shell_sites =
            shellSiteMap[filename.get<std::string>()].shell_sites;
        std::vector<uint8_t> &occ_ref =
            shellSiteMap[filename.get<std::string>()].occ_ref;
        std::vector<uint8_t> &flag =
            shellSiteMap[filename.get<std::string>()].flag;
        if ((retval = nc_inq_dimid(ncid, "steps", &stepsid))) {
          continue;
        }
        if ((retval = nc_inq_dimlen(ncid, stepsid, &steps))) {
          continue;
        }
        if ((retval = nc_inq_dimid(ncid, "shell_size", &shellsizeid))) {
          continue;
        }
        if ((retval = nc_inq_dimlen(ncid, shellsizeid, &shellsize))) {
          continue;
        }

        shell_sites.resize(shellsize * steps);
        int shell_sites_id;
        if ((retval = nc_inq_varid(ncid, "shell_sites", &shell_sites_id))) {
          continue;
        }

        if ((retval = nc_get_var(ncid, shell_sites_id, shell_sites.data()))) {
          continue;
        }

        flag.resize(steps);
        int flag_id;
        if ((retval = nc_inq_varid(ncid, "flag", &flag_id))) {
          continue;
        }
        if ((retval = nc_get_var(ncid, flag_id, flag.data()))) {
          continue;
        }

        occ_ref.resize(shellsize);
        int occ_ref_id;
        if ((retval = nc_inq_varid(ncid, "occ_ref", &occ_ref_id))) {
          continue;
        }
        if ((retval = nc_get_var(ncid, occ_ref_id, occ_ref.data()))) {
          continue;
        }
        if ((retval = nc_close(ncid))) {
          continue;
        }
      }
    }
  }
  mShellSiteTypes.setElements(siteTypes);
  mShellSiteTypes.setNoCalls(0);
}

DatasetManager::SpeciesLabelMap DatasetManager::getAvailableSpecies() {
  SpeciesLabelMap labelMap;
  for (std::string atom : mAvailableAtomsJson) {
    labelMap[atom] = std::vector<std::string>();
  }
  for (std::string label : mShowAtomElements) {
    for (auto mapEntry : labelMap) {
      if (label.compare(0, mapEntry.first.size(), mapEntry.first) == 0) {
        labelMap[mapEntry.first].push_back(label);
        continue;
      }
    }
  }
  labelMap["Va"] = mVacancyAtoms;
  return labelMap;
}

std::vector<std::pair<std::string, float>>
DatasetManager::getCurrentCompositions() {
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

        for (auto ps : mParameterSpace.getDimensions()) {
          if (ps->getSpaceRepresentationType() ==
              ParameterSpaceDimension::INDEX) {
            if (ps->size() > 0) {
              auto index = ps->getCurrentIndex();
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
            break;
          }
        }
      }
    }
  }
  return compositions;
}

std::vector<int8_t> DatasetManager::getShellSiteTypes(size_t timeIndex) {
  std::vector<int8_t> matches;
  int8_t count = 0;

  for (auto &neighborhoodGroups : shellSiteMap) {
    std::vector<int16_t> sites;
    if (!mParameterSpace.getDimension("time") ||
        mParameterSpace.getDimension("time")->size() <= timeIndex) {
      continue;
    }
    if (neighborhoodGroups.second.shell_sites.size() <=
        timeIndex * neighborhoodGroups.second.occ_ref.size()) {
      continue;
    }
    auto *p =
        &neighborhoodGroups.second
             .shell_sites[timeIndex * neighborhoodGroups.second.occ_ref.size()];
    for (size_t i = 0; i < neighborhoodGroups.second.occ_ref.size(); i++) {
      if (*p != UINT16_MAX) {
        sites.push_back(*p);
      }
      p++;
    }

    assert(sites.size() == 0 ||
           sites.size() == neighborhoodGroups.second.occ_ref.size());
    if (sites.size() > 0) {
      // Previous algorithm that found matches from occupancy. This algorithm
      // has issues
      // because it doesn't include transpositions of the site shapes.
      //      auto siteIt = sites.begin();
      //      bool isMatch = true;
      //      for (auto occ_ref : neighborhoodGroups.second.occ_ref) {
      //        size_t index = numAtoms * timeIndex + *siteIt;
      //        auto occ_dof = trajectoryData[index];
      //        if (occ_dof != occ_ref) {
      //          isMatch = false;
      //          break;
      //        }
      //        siteIt++;
      //      }

      if (neighborhoodGroups.second.flag.size() > timeIndex) {
        if (neighborhoodGroups.second.flag[timeIndex] == 1) {
          matches.push_back(count);
        }
      }
    }
    count++;
  }
  return matches;
}

std::string DatasetManager::findJsonFile(std::string datasetId,
                                         std::string subDir,
                                         std::string fileName) {

  std::string prefix = getGlobalRootPath();
  if (prefix.size() > 0) {
    prefix += AL_FILE_DELIMITER_STR;
  }

  std::vector<std::string> paths = {
      prefix + datasetId + "/" + subDir + "/", prefix + datasetId + "/",
      getGlobalRootPath() + "/", prefix + datasetId,
      prefix + datasetId + "/../"};

  for (auto possiblePath : paths) {
    //            std::cout << "Looking for " << fileName << " in " <<
    //            possiblePath << std::endl;
    if (File::exists(possiblePath + fileName)) {
      return possiblePath + fileName;
    }
  }

  return std::string();
}
