
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

#include "tinc/NetCDFDiskBuffer.hpp"

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
  mParameterSpaces["temperature"] = new ParameterSpaceDimension("temperature");
  mParameterSpaces["chempotA"] = new ParameterSpaceDimension("chempotA");
  mParameterSpaces["chempotB"] = new ParameterSpaceDimension("chempotB");

  mParameterSpaces["chempotA"]->addConnectedParameterSpace(
      mParameterSpaces["chempotB"]);
  mParameterSpaces["chempotB"]->addConnectedParameterSpace(
      mParameterSpaces["chempotA"]);

  mParameterSpaces["time"] = new ParameterSpaceDimension("time");

  mParameterSpaces["temperature"]->parameter().set(300);

  // These parameters are used for data propagation, but shouldn't be set by the
  // user
  currentGraphName.setHint("hide", 1.0);
  currentPoscarName.setHint("hide", 1.0);

  initializeComputation();
}

void DatasetManager::initializeComputation() {

  cacheManager.registerProcessor(labelProcessor);
  cacheManager.registerProcessor(graphGenerator);

  //  setenv("HDF5_DISABLE_VERSION_CHECK", "2", 1);
  //  graphGenerator.verbose();
  labelProcessor.verbose();

  for (auto pspace : mParameterSpaces) {
    labelProcessor << pspace.second->parameter();
    graphGenerator << pspace.second->parameter();
  }
  atomPositionChain << labelProcessor;
  sampleComputationChain << atomPositionChain << graphGenerator;

  // Graph generator
  graphGenerator.prepareFunction = [&]() {
    std::string datasetId = mCurrentDataset.get();
    std::string xLabel = mPlotXAxis.getCurrent();
    std::string yLabel = mPlotYAxis.getCurrent();

    auto jsonText = readJsonResultsFile(mCurrentDataset, getSubDir());
    if (jsonText.size() > 0) {
      std::vector<double> datax, datay;
      auto resultsJson = json::parse(jsonText);

      if (resultsJson.find(yLabel) != resultsJson.end()) {
        try {
          datay.reserve(resultsJson[yLabel].size());
          datay.insert(datay.begin(), resultsJson[yLabel].begin(),
                       resultsJson[yLabel].end());
        } catch (std::exception &) {
        }
      }
      std::string highlightValue = "0.0";

      if (mParameterSpaces.find(xLabel) != mParameterSpaces.end()) {
        try {
          datax.reserve(mParameterSpaces[xLabel]->size());
          for (auto value : mParameterSpaces[xLabel]->values()) {
            datax.push_back(value);
          }
          highlightValue =
              std::to_string(int(mParameterSpaces[xLabel]->getCurrentValue()));
        } catch (std::exception &) {
        }
      }

      std::string fullDatasetPath = File::conformPathToOS(
          buildRootPath() + File::conformPathToOS(mCurrentDataset.get()));

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
        for (json::iterator v = data.begin(); v != data.end(); ++v) {
          if (v->get<float>() < dataRange.first) {
            dataRange.first = v->get<float>();
          } else if (v->get<float>() > dataRange.second) {
            dataRange.second = v->get<float>();
          }
        }
      }

      std::string yLabelSanitized = yLabel;
      replaceAll(yLabelSanitized, "<", "_");
      replaceAll(yLabelSanitized, ">", "_");
      replaceAll(yLabelSanitized, "(", "_");
      replaceAll(yLabelSanitized, ")", "_");

      graphGenerator.configuration["title"] = mTitle;
      graphGenerator.configuration["dataset"] = datasetId;
      graphGenerator.configuration["temp_interest"] = highlightValue;
      graphGenerator.configuration["subdir"] = subDir;
      graphGenerator.configuration["xLabel"] = xLabel;
      graphGenerator.configuration["yLabel"] = yLabelSanitized;

      graphGenerator.configuration["miny"] = std::to_string(dataRange.first);
      graphGenerator.configuration["maxy"] = std::to_string(dataRange.second);
      graphGenerator.configuration["inxFile"] =
          std::string("cached_output/inx.bin");
      graphGenerator.configuration["inyFile"] =
          std::string("cached_output/iny.bin");

      graphGenerator.setOutputDirectory(graphGenerator.runningDirectory() +
                                        "/cached_output");
      graphGenerator.setOutputFileNames({datasetId + "_" + highlightValue +
                                         "_" +
                                         DataScript::sanitizeName(subDir) +
                                         "_" + yLabelSanitized + "_graph.png"});
    } else {
      return false;
    }
    return true;
  };

  graphGenerator.registerDoneCallback([this](bool runOk) {
    if (runOk) {
      currentGraphName.set(graphGenerator.outputFile());
      //      std::cout << "Generated graph: " << graphGenerator.outputFile()
      //                << std::endl;
      //              graphProcessing = false;
    } else {
      std::cerr << "Graph generation failed." << std::endl;
    }
  });

  labelProcessor.prepareFunction = [&]() {
    std::string condition = std::to_string(
        mParameterSpaces[mConditionsParameter]->getCurrentIndex());

    std::string folder = File::conformDirectory(getSubDir());
    std::string root_path = buildRootPath();
    std::string template_pos_path = File::conformPathToOS(
        root_path + mCurrentDataset.get() + "/cached_output/template_POSCAR");
    // Try to find "template_POSCAR"
    if (!File::exists(template_pos_path)) {
      std::cerr << "Failed to find template file. Call it 'template_POSCAR'."
                << std::endl;
      return false;
    }
    std::string prim_path;
    if (File::exists(root_path + mCurrentDataset.get() + "/prim_labels.json")) {
      prim_path = root_path + mCurrentDataset.get() + "/prim_labels.json";
    } else if (File::exists(root_path + mCurrentDataset.get() + "/prim.json")) {
      prim_path = root_path + mCurrentDataset.get() + "/prim.json";
    } else if (File::exists(root_path + "prim_labels.json")) {
      prim_path = root_path + "prim_labels.json";
    } else if (File::exists(root_path + "prim.json")) {
      prim_path = root_path + "prim.json";
    }

    labelProcessor.setRunningDirectory(mLoadedDataset);

    labelProcessor.setOutputFileNames(
        {DataScript::sanitizeName(template_pos_path + mCurrentDataset.get() +
                                  "_" + folder + "_" + condition)});

    std::string conditionSubdir =
        File::conformPathToOS(folder + "conditions." + condition + "/");

    labelProcessor.configuration["dataset_path"] =
        File::conformPathToOS(root_path + mCurrentDataset.get());
    labelProcessor.configuration["prim_path"] = prim_path;

    std::string final_state_path =
        conditionSubdir + "final_state.json"; // final_state.json file
    labelProcessor.configuration["final_state_path"] = File::conformPathToOS(
        root_path + mCurrentDataset.get() + "/" + final_state_path);

    labelProcessor.configuration["template_path"] = File::conformPathToOS(
        root_path + mCurrentDataset.get() + "/cached_output/template_POSCAR");
    labelProcessor.configuration["condition"] = condition;
    return true;
  };

  labelProcessor.registerDoneCallback([&](bool ok) {
    if (ok) {
      currentPoscarName.set(labelProcessor.outputFile());

      // --------------netcdf ----------------------------

      int ncid, retval;

      std::string filename = currentPoscarName.get() + ".nc";
      /* Open the file. NC_NOWRITE tells netCDF we want read-only access
       * to the file.*/
      if ((retval = nc_open(filename.c_str(), NC_NOWRITE | NC_SHARE, &ncid))) {

        return /*false*/;
      }
      int varid;
      if ((retval = nc_inq_varid(ncid, "atoms_var", &varid))) {
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

      /* Read the data. */
      if ((retval = nc_get_var(ncid, varid, newPositionData->data()))) {
        return /*false*/;
      }
      occupationData.doneWriting(newPositionData);

      /* Close the file, freeing all resources. */
      if ((retval = nc_close(ncid))) {
        return /*false*/;
      }

      // ---------------------------------------------
    }
  });

  // Parameter spaces trigger computation

  for (auto &parameterSpace : mParameterSpaces) {
    parameterSpace.second->parameter().registerChangeCallback([&](float value) {
      if (parameterSpace.second->getCurrentId() !=
          parameterSpace.second->idAt(parameterSpace.second->getIndexForValue(
              value))) { // Only reload if id has changed

        // To have the internal value already changed for the following
        // functions. This throws away the previous value, but we don't need it
        // now.
        parameterSpace.second->parameter().setNoCalls(value);
        std::cout << value << " : " << parameterSpace.second->getCurrentId()
                  << "..."
                  << parameterSpace.second->idAt(
                         parameterSpace.second->getIndexForValue(value));
        if (mParameterSpaces["time"]->size() > 0) {
          if (parameterSpace.first == mConditionsParameter) {
            loadTrajectory();
          }

          computeNewSample();
        }
        updateText();
      }
    });
  }

  mPlotYAxis.registerChangeCallback([this](float value) {
    if (mPlotYAxis.get() != value) {
      graphGenerator.process();
    }
  });

  mPlotXAxis.registerChangeCallback([this](float value) {
    if (mPlotXAxis.get() != value) {
      graphGenerator.process();
    }
  });
}

void DatasetManager::setPythonBinary(std::string pythonBinaryPath) {
  //      std::unique_lock<std::mutex> lk(mProcessingLock);
  graphGenerator.setCommand(pythonBinaryPath);
  labelProcessor.setCommand(pythonBinaryPath);
}

void DatasetManager::setPythonScriptPath(std::string pythonScriptPath) {
  //      std::unique_lock<std::mutex> lk(mProcessingLock);
  graphGenerator.setScriptName(pythonScriptPath + "/graphing/plot.py");
  labelProcessor.setScriptName(File::conformDirectory(
      pythonScriptPath) + "/reassign_occs/reassign_occs.py");
}

std::string DatasetManager::buildRootPath() {
  return File::conformPathToOS(File::conformPathToOS(mGlobalRoot) +
                               File::conformPathToOS(mRootPath));
}

std::string DatasetManager::fullConditionPath() {
  std::string condition =
      std::to_string(mParameterSpaces[mConditionsParameter]->getCurrentIndex());
  return File::conformPathToOS(buildRootPath() + mCurrentDataset.get() + "/" +
                               getSubDir() + "/conditions." + condition + "/");
}

void DatasetManager::readParameterSpace() {
  std::string paramSpaceFile =
      File::conformPathToOS(buildRootPath() +
                            File::conformPathToOS(mCurrentDataset.get())) +
      "cached_output/_parameter_space.json";
  // Read parameter space json file
  std::string str;
  std::ifstream f(paramSpaceFile);
  std::map<std::string, std::string> parameterNameMap = {
      {"T", "temperature"},
      {"param_chem_pot(a)", "chempotA"},
      {"param_chem_pot(b)", "chempotB"}};
  if (!f.fail()) {
    mParameterForSubDir.clear();
    for (auto &paramSpace : mParameterSpaces) {
      paramSpace.second->clear();
    }

    f.seekg(0, std::ios::end);
    str.reserve(f.tellg());
    f.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(f)),
               std::istreambuf_iterator<char>());
    json j = json::parse(str);
    for (json::iterator it = j["parameters"].begin();
         it != j["parameters"].end(); ++it) {
      //          std::cout << it.key() << " : " << it.value() << std::endl;
      std::string key = it.key();
      if (parameterNameMap.find(key) != parameterNameMap.end()) {
        std::string mappedKey = parameterNameMap[key];
        mParameterForSubDir.push_back(mappedKey);
        for (auto temps : it.value()) {
          mParameterSpaces[mappedKey]->push_back(float(temps["value"]),
                                                 temps["dir"]);
        }
        if (mappedKey == "chempotA" || mappedKey == "chempotB") {
          mParameterSpaces[mappedKey]->sort();
        }
        //        // set to current value to clamp if needed
        //        mParameterSpaces[mappedKey]->parameter().set(
        //            mParameterSpaces[mappedKey]->parameter().get());
      }
    }

    for (json::iterator it = j["conditions"].begin();
         it != j["conditions"].end(); ++it) {
      std::string key = it.key();
      if (parameterNameMap.find(key) != parameterNameMap.end()) {
        key = parameterNameMap[key];
        for (auto temps : it.value()) {
          mParameterSpaces[key]->push_back(float(temps));
        }
        if (key == "chempotA" || key == "chempotB") {
          mParameterSpaces[key]->sort();
        }
        // set to current value to clamp if needed
        mParameterSpaces[key]->parameter().setNoCalls(
            mParameterSpaces[key]->parameter().get());
        if (mConditionsParameter.size() > 0) {
          std::cerr << "ERROR conditions parameter already set. Overwriting."
                    << std::endl;
        }
        mConditionsParameter = key;
      }
    }
    if (mConditionsParameter == "" &&
        mParameterSpaces["temperature"]->size() > 0) {
      // This is a fallback to support datasets with only a single condition.
      // In this case no conditions parameter is found, and therefore this would
      // be considered an invalid dataset...
      mConditionsParameter = "temperature";
    }

    //    Verify conditions parameter space from parameter space file against
    //    results.json
    std::string resultsFile =
        File::conformPathToOS(buildRootPath() +
                              File::conformPathToOS(mCurrentDataset.get())) +
        "/" + File::conformPathToOS(getSubDir()) + "results.json";
    if (File::exists(resultsFile)) {
      std::ifstream results(resultsFile);
      results.seekg(0, std::ios::end);
      str.reserve(results.tellg());
      results.seekg(0, std::ios::beg);

      str.assign((std::istreambuf_iterator<char>(results)),
                 std::istreambuf_iterator<char>());
      json resultsJson = json::parse(str);

      std::string conditionInJson = mConditionsParameter;
      for (auto mapEntry : parameterNameMap) {
        if (mapEntry.second == mConditionsParameter) {
          conditionInJson = mapEntry.first;
          break;
        }
      }

      auto conditionValues = resultsJson[conditionInJson];
      if (conditionValues.is_array()) {
        mParameterSpaces[mConditionsParameter]->clear();
        for (auto value : conditionValues) {
          mParameterSpaces[mConditionsParameter]->push_back(value.get<float>());
        }
      }
    }

    for (json::iterator it = j["internal_states"].begin();
         it != j["internal_states"].end(); ++it) {
      std::string key = it.key();
      json j2 = it.value();
      int timeOffset = 0;
      if (key == "time") {
        mParameterSpaces["time"]->reserve(mParameterSpaces["time"]->size() +
                                          j2.size());
        uint64_t previous = std::numeric_limits<uint64_t>::min();
        for (json::iterator values = j2.begin(); values != j2.end(); values++) {
          // TODO we should get the actual time values here (available in the
          // results file.
          if (mParameterSpaces["time"]->size() > 0 && previous > *values) {
            timeOffset = mParameterSpaces["time"]->at(
                mParameterSpaces["time"]->size() - 1);
          }
          mParameterSpaces["time"]->push_back(int(*values) + timeOffset);
          previous = int(*values);
        }
      }
    }
    // If there is a time space, don't use script. Data is already avaialble
    if (mParameterSpaces["time"]->size() > 0) {
      labelProcessor.enabled = false;
    } else {
      labelProcessor.enabled = true;
    }
  } else {
    std::cerr << "ERROR failed to find parameter space file: " << paramSpaceFile
              << std::endl;
  }
}

void DatasetManager::initRoot() {
  std::unique_lock<std::mutex> lk(mDataLock);

  std::string fullDatasetPath = File::conformPathToOS(
      buildRootPath() + File::conformPathToOS(mCurrentDataset.get()));
  mLoadedDataset = fullDatasetPath;

  // Update processors configuration
  cacheManager.setOutputDirectory(fullDatasetPath + "/cached_output");
  cacheManager.setRunningDirectory(fullDatasetPath);

  readParameterSpace();
  analyzeDataset();

  loadTrajectory();

  std::vector<std::string> parameterSpaceNames;
  parameterSpaceNames.push_back("temperature");

  mPlotXAxis.setElements(parameterSpaceNames);
  mPlotYAxis.set(0);

  auto dataNames = getDataNames();
  mPlotYAxis.setElements(dataNames);
  std::string defaultYAxis = "<comp_n(" + mAtomOfInterest.getCurrent() + ")>";
  ptrdiff_t pos = std::find(dataNames.begin(), dataNames.end(), defaultYAxis) -
                  dataNames.begin();
  if (pos < (int)dataNames.size()) {
    mPlotYAxis.set((int)pos);
  }

  std::string templatePath =
      fullDatasetPath + "cached_output/template_POSCAR.nc";
  if (File::exists(templatePath)) {
    int retval, ncid, varid;
    if ((retval =
             nc_open(templatePath.c_str(), NC_NOWRITE | NC_SHARE, &ncid))) {
      return /*false*/;
    }
    if ((retval = nc_inq_varid(ncid, "atoms_var", &varid))) {
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

    templateData.resize(lenp);

    /* Read the data. */
    if ((retval = nc_get_var(ncid, varid, templateData.data()))) {
      return /*false*/;
    }

    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid))) {
      return /*false*/;
    }
  }

  // Configure diff generator
  if (mParameterSpaces["time"]->size() > 0) {
    graphGenerator.enabled = false; // Graphing not working with kmc yet.
  } else {
    graphGenerator.enabled = true;
  }

  lk.unlock();
  // Update data
  // Not necessary as this is triggered automatically through change in
  // mLoadedDataset
  //  computeNewSample();
}

void DatasetManager::analyzeDataset() {
  if (mCurrentDataset.get() == "") {
    std::cerr << "ERROR: No datasets found." << std::endl;
    return;
  }
  std::string datasetId = mCurrentDataset.get();
  // mCurrentDataset might be more current than the internal parameter value
  // as this might be called from the parameter change callback, when the
  // internal value has not yet been updated

  mDataRanges.clear();
  mAvailableAtomsJson.clear();

  for (auto &paramSpace : mParameterSpaces) {
    paramSpace.second->lock();
    bool isVariable = false;
    if (paramSpace.second->size() > 0) {
      auto previous = paramSpace.second->at(0);
      for (auto &value : paramSpace.second->values()) {
        if (value != previous) {
          isVariable = true;
          break;
        }
      }
    }
    mParameterIsVariable[paramSpace.first] = isVariable;
    paramSpace.second->parameter().setHint(
        "hide", paramSpace.second->size() > 0 ? 0.0 : 1.0);
    paramSpace.second->unlock();
  }

  mVacancyAtoms.clear();
  mShowAtomElements.clear();

  std::string str = readJsonResultsFile(datasetId, getSubDir());
  if (str.size() > 0) {
    // Read the results file. This tells us which atoms are available
    auto resultsJson = json::parse(str);

    std::string comp_prefix = "<comp_n(";
    // Fill available atoms list from the available atoms in the <comp_n(XX)>
    // tag
    for (json::iterator it = resultsJson.begin(); it != resultsJson.end();
         ++it) {
      std::string key = it.key();
      //                // Look for available comp_n atom names
      //                if (key.compare(0, comp_prefix.size(), comp_prefix) ==
      //                0) {
      //                    assert(key.find(")>") != string::npos);
      //                    availableAtomsJson.push_back(key.substr(comp_prefix.size(),
      //                    key.find(")>") - comp_prefix.size()));
      //                }
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
    std::cerr << "ERROR failed to find prim labels file" << std::endl;
  }
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
  //                for (auto chempot2: mParameterSpaces["chempotB"]->values())
  //                {
  //                    for (int i = 0; i < numProcessesSlice; i++) {
  //                        int temp = nodeIndeces[i];
  //                        labelProcessor.setParams(chempot.first,
  //                        chempot2.first, std::to_string(temp), datasetId);
  //                        labelProcessor.processAsync();
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
  //        return false; // force testing
  // TODO have a better validity check
  bool valid = mParameterSpaces["chempotA"]->size() > 0 &&
               mParameterSpaces["temperature"]->size() > 0;
  return valid;
}

std::string DatasetManager::currentDataset() { return mLoadedDataset; }

std::string DatasetManager::getSubDir() {
  if (mParameterForSubDir.size() > 0) {
    return mParameterSpaces[mParameterForSubDir[0]]->getCurrentId();
  } else {
    return std::string();
  }
}

void DatasetManager::processTemplatePositions() {
  //  std::unique_lock<std::mutex> lk(mDataLock);
  auto templatePoscarName = labelProcessor.outputFile();

  // Load POSCAR data
  if (reader.loadFile(templatePoscarName)) {
    mTemplatePositions = reader.getAllPositions();
  } else {
    std::cerr << "ERROR creating template at time 0 ----------- " << std::endl;
  }
  // Load empty template
  VASPReader emptyTemplateReader;
  if (emptyTemplateReader.loadFile(
          File::conformPathToOS(buildRootPath() + mCurrentDataset.get() +
                                "/cached_output/template_POSCAR"))) {
    mEmptyTemplate = emptyTemplateReader.getElementPositions("X");
  }
}

bool DatasetManager::loadDiff(int timeIndex) {
  bool diffLoaded = false;
  //         Then accumulate all diffs until current time
  int targetIndex = mParameterSpaces["time"]->getCurrentIndex();

  //      std::cout << timeIndex << " --> " << targetIndex << std::endl;
  if (timeIndex < targetIndex) {

    mHistory.clear();
    for (int i = timeIndex; i != targetIndex; i++) {
      std::cout << "applying diff " << i << std::endl;
      std::pair<Vec3f, Vec3f> historyPoint;
      auto this_diff = mDiffs[i];
      auto this_diff_indeces = this_diff[0];
      auto this_diff_labels = this_diff[1];

      for (int change = 0; change < this_diff_indeces.size(); change++) {
        int changeOffset = this_diff_indeces[change].get<int>() * 4;
        if (changeOffset < mEmptyTemplate.size()) {
          Vec3f pos(mEmptyTemplate[changeOffset],
                    mEmptyTemplate[changeOffset + 1],
                    mEmptyTemplate[changeOffset + 2]);

          int index = -1;
          std::string label;
          for (auto &positions : mTemplatePositions) {
            if (positions.first !=
                this_diff_labels[change]) { // Look only for changes
              for (int curIndex = 0; curIndex < positions.second.size();
                   curIndex += 4) {
                //                  float distance = (pos -
                //                  Vec3f(positions.second[curIndex],
                //                                 positions.second[curIndex +
                //                                 1],
                //                                                positions.second[curIndex
                //                                                + 2])).mag();
                if ((pos - Vec3f(positions.second[curIndex],
                                 positions.second[curIndex + 1],
                                 positions.second[curIndex + 2]))
                        .mag() < 0.001f) {
                  label = positions.first;
                  index = curIndex;

                  //                      std::cout << "match " <<
                  //                      this_diff_indeces[change] << " " <<
                  //                      pos.x << "," << pos.y << std::endl;
                  break;
                }
              }
            }
            if (index != -1) {
              break;
            }
          }
          if (label == "") {
            std::cout << "ERROR: Label not found for index " << timeIndex
                      << " label " << this_diff_labels[change] << std::endl;
            return false;
          }
          if (this_diff_labels[change] == "Va") { // Remove atom
            historyPoint.first = pos;
            //                std::cout << "remove index " << index <<
            //                std::endl;
            assert(index >= 0);
            mTemplatePositions[this_diff_labels[change]].insert(
                mTemplatePositions[this_diff_labels[change]].begin(),
                mTemplatePositions[label].begin() + index,
                mTemplatePositions[label].begin() + index + 4);
            mTemplatePositions[label].erase(
                mTemplatePositions[label].begin() + index,
                mTemplatePositions[label].begin() + index + 4);
          } else { // Add Atom
            historyPoint.second = pos;
            //                std::cout << "add index " << index << std::endl;
            mTemplatePositions[this_diff_labels[change]].push_back(pos.x);
            mTemplatePositions[this_diff_labels[change]].push_back(pos.y);
            mTemplatePositions[this_diff_labels[change]].push_back(pos.z);
            mTemplatePositions[this_diff_labels[change]].push_back(0.0);
            mTemplatePositions[label].erase(
                mTemplatePositions[label].begin() + index,
                mTemplatePositions[label].begin() + index + 4);
          }
        } else {
          std::cerr << " ERROR: diff index greater than template size."
                    << std::endl;
          break;
        }
      }
      mHistory.push_back(historyPoint);
    }
    diffLoaded = true;
  } else if (timeIndex > targetIndex) {
  }
  return diffLoaded;
}

// void DatasetManager::loadFromPOSCAR() {}

void DatasetManager::computeNewSample() {
  if (mConditionsParameter == "" || !mRunProcessors) {
    return;
  }
  std::unique_lock<std::mutex> lk(mDataLock);

  sampleComputationChain.process();

  if (mParameterSpaces["time"]->size() > 0) {
    occupationData.doneWriting(occupationData.getWritable());
  }

  updateText();
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
  if (File::exists(trajectoryFile)) {

    int retval, ncid, varid;
    if ((retval =
             nc_open(trajectoryFile.c_str(), NC_NOWRITE | NC_SHARE, &ncid))) {
      return /*false*/;
    }
    if ((retval = nc_inq_varid(ncid, "occupation", &varid))) {
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

    /* Read the data. */
    if ((retval = nc_get_var(ncid, varid, trajectoryData.data()))) {
      return /*false*/;
    }

    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid))) {
      return /*false*/;
    }
  }
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

std::string DatasetManager::findJsonFile(std::string datasetId,
                                         std::string subDir,
                                         std::string fileName) {
  std::vector<std::string> paths = {
      buildRootPath() + "/" + datasetId + "/" + subDir + "/",
      buildRootPath() + "/" + datasetId + "/", buildRootPath() + "/"};

  for (auto possiblePath : paths) {
    //            std::cout << "Looking for " << fileName << " in " <<
    //            possiblePath << std::endl;
    if (File::exists(possiblePath + fileName)) {
      return possiblePath + fileName;
    }
  }

  return std::string();
}
