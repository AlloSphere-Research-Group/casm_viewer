
#ifdef AL_WINDOWS
#define NOMINMAX
#include <Windows.h>
#undef far
#endif

#include <iostream>
#include <regex>
#include <algorithm>
#include <memory>
#include <mutex>
#include <limits>
#include <queue>
#include <fstream>

#include "datasetmanager.hpp"

using namespace  al;

DatasetManager::DatasetManager() {
  mParameterSpaces["temperature"] = new ParameterSpace("temperature");
  mParameterSpaces["chempotA"] = new ParameterSpace("chempotA");
  mParameterSpaces["chempotB"] = new ParameterSpace("chempotB");

  mParameterSpaces["chempotA"]->addConnectedParameterSpace(mParameterSpaces["chempotB"]);
  mParameterSpaces["chempotB"]->addConnectedParameterSpace(mParameterSpaces["chempotA"]);

  mParameterSpaces["time"] = new ParameterSpace("time");

  mParameterSpaces["temperature"]->parameter().set(300);
#ifdef NDEBUG
  labelProcessor.verbose(true);
#endif
  // These parameters are used for data propagation, but shouldn't be set by the user
  currentGraphName.setHint("hide", 1.0);
  currentPoscarName.setHint("hide", 1.0);
#ifdef NDEBUG
  //        graphGenerator.verbose(true);
  //        mParameterSpaceProcessor.verbose(true);
  //        templateGen.verbose(true);
#endif

  cacheManager.registerProcessor(labelProcessor);
  cacheManager.registerProcessor(graphGenerator);
  cacheManager.registerProcessor(templateGen);

  labelProcessor.maxAsyncProcesses(10);
}

void DatasetManager::setPythonBinary(std::string pythonBinaryPath) {
  //      std::unique_lock<std::mutex> lk(mProcessingLock);
  graphGenerator.setCommand(pythonBinaryPath);
  labelProcessor.setCommand(pythonBinaryPath);
  templateGen.setCommand(pythonBinaryPath);
  diffGen.setCommand(pythonBinaryPath);
}

void DatasetManager::setPythonScriptPath(std::string pythonScriptPath) {
  //      std::unique_lock<std::mutex> lk(mProcessingLock);
  labelProcessor.pythonScriptsPath = pythonScriptPath;
  graphGenerator.pythonScriptsPath = pythonScriptPath;
  graphGenerator.setScriptName(pythonScriptPath + "/graphing/plot.py");
  templateGen.setScriptName(pythonScriptPath + "/reassign_occs/template_creator.py");
  diffGen.setScriptName(pythonScriptPath + "/reassign_occs/analyze_dmc.py");
}

std::string DatasetManager::buildRootPath() {
  return  File::conformPathToOS(File::conformPathToOS(mGlobalRoot) + File::conformPathToOS(mRootPath));
}

void DatasetManager::initRoot() {
  std::unique_lock<std::mutex> lk(mDataLock);

  std::string fullDatasetPath = File::conformPathToOS(buildRootPath() + File::conformPathToOS(mCurrentDataset.get()));
  mLoadedDataset = fullDatasetPath;

  if (mRunProcessors) {

    cacheManager.setOutputDirectory(fullDatasetPath + "/cached_output");
    cacheManager.setRunningDirectory(fullDatasetPath);

    labelProcessor.root_path = buildRootPath();

    // FIMXE this could be moved as there is no need to recompute on every load
    templateGen.process(templateGen.configuration(), true);

    diffGen.datasetPath = fullDatasetPath;
    diffGen.setRunningDirectory(fullDatasetPath);
    diffGen.setOutputDirectory(fullDatasetPath);
    //          diffGen.processAsync(templateGen.configuration(), false, [&](){ std::cout << "Computed diff for " << diffGen.datasetPath << std:: endl;});

  }
  std::string paramSpaceFile = File::conformPathToOS(buildRootPath() + File::conformPathToOS(mCurrentDataset.get())) + "cached_output/_parameter_space.json";
  // Read parameter space json file
  std::string str;
  std::ifstream f(paramSpaceFile);
  std::map<std::string, std::string> parameterNameMap = {{"T", "temperature"},
                                                         {"param_chem_pot(a)", "chempotA"},
                                                         {"param_chem_pot(b)", "chempotB"}};
  if (!f.fail()) {
    for (auto &paramSpace: mParameterSpaces) {
      paramSpace.second->lock();
      paramSpace.second->clear();
    }

    f.seekg(0, std::ios::end);
    str.reserve(f.tellg());
    f.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(f)),
               std::istreambuf_iterator<char>());
    json j = json::parse(str);
    for (json::iterator it = j["parameters"].begin(); it != j["parameters"].end(); ++it) {
      //          std::cout << it.key() << " : " << it.value() << std::endl;
      std::string key = it.key();
      if (parameterNameMap.find(key) != parameterNameMap.end()) {
        std::string mappedKey = parameterNameMap[key];
        mParameterForSubDir.push_back(mappedKey);
        for (auto temps: it.value()) {
          mParameterSpaces[mappedKey]->push_back(float(temps["value"]), temps["dir"]);
        }
        if (mappedKey == "chempotA" || mappedKey == "chempotB") {
          mParameterSpaces[mappedKey]->sort();
        }
        // set to current value to clamp if needed
        mParameterSpaces[mappedKey]->parameter().set(mParameterSpaces[mappedKey]->parameter().get());
      }
    }

    for (json::iterator it = j["conditions"].begin(); it != j["conditions"].end(); ++it) {
      std::string key = it.key();
      if (parameterNameMap.find(key) != parameterNameMap.end()) {
        key = parameterNameMap[key];
        for (auto temps: it.value()) {
          mParameterSpaces[key]->push_back(float(temps));
        }
        if (key == "chempotA" || key == "chempotB") {
          mParameterSpaces[key]->sort();
        }
        // set to current value to clamp if needed
        mParameterSpaces[key]->parameter().set(mParameterSpaces[key]->parameter().get());
        if (mConditionsParameter.size() > 0) {
          std::cerr << "ERROR conditions parameter already set. Overwriting." << std::endl;
        }
        mConditionsParameter = key;
      }
    }


    for (json::iterator it = j["internal_states"].begin(); it != j["internal_states"].end(); ++it) {
      std::string key = it.key();
      json j2 = it.value();
      int timeOffset = 0;
      if (key == "time") {
        for (json::iterator values = j2.begin(); values != j2.end(); values++) {
          // TODO we should get the actual time values here (available in the results file.
          if (mParameterSpaces["time"]->size() > 0 && mParameterSpaces["time"]->at(mParameterSpaces["time"]->size() -1) > int(*values) ) {
            timeOffset = mParameterSpaces["time"]->at(mParameterSpaces["time"]->size() - 1);
          }
          mParameterSpaces["time"]->push_back(int(*values) + timeOffset);
        }
      }
    }
    // If there is a time space, use different script
    if (mParameterSpaces["time"]->size() > 0) {
      labelProcessor.setScriptName(File::conformDirectory(labelProcessor.pythonScriptsPath) + "reassign_occs/result_extractor.py");
    } else {
      labelProcessor.setScriptName(File::conformDirectory(labelProcessor.pythonScriptsPath) + "reassign_occs/reassign_occs.py");
    }
    for (auto &paramSpace: mParameterSpaces) {
      paramSpace.second->unlock();
    }

  } else {
    std::cerr << "ERROR failed to find parameter space file: " << paramSpaceFile << std::endl;
  }
  analyzeDataset();
  lk.unlock();
  getAtomPositions();
}

void DatasetManager::analyzeDataset() {
  if(mCurrentDataset.get() == "") {
    std::cerr << "ERROR: No datasets found." << std::endl;
    return;
  }
  std::string datasetId = mCurrentDataset.get();; // mCurrentDataset might be more current than the internal parameter value as this might be called from the parameter change callback, when the internal value has not yet been updated

  mDataRanges.clear();
  mAvailableAtomsJson.clear();

  for (auto &paramSpace: mParameterSpaces) {
    paramSpace.second->lock();
    bool isVariable = false;
    if (paramSpace.second->size() > 0) {
      auto previous = paramSpace.second->at(0);
      for (auto &value: paramSpace.second->values()) {
        if (value.second != previous) {
          isVariable = true;
          break;
        }
      }
    }
    mParameterIsVariable[paramSpace.first] = isVariable;
    paramSpace.second->parameter().setHint("hide", paramSpace.second->size() > 0 ? 0.0:1.0);
    paramSpace.second->unlock();
  }

  mVacancyAtoms.clear();
  mShowAtomElements.clear();

  std::string str = readJsonResultsFile(datasetId, getSubDir());
  if (str.size() > 0) {
    // Read the results file. This tells us which atoms are available
    auto resultsJson = json::parse(str);

    std::string comp_prefix = "<comp_n(";
    // Fill available atoms list from the available atoms in the <comp_n(XX)> tag
    for (json::iterator it = resultsJson.begin(); it != resultsJson.end(); ++it) {
      std::string key = it.key();
      //                // Look for available comp_n atom names
      //                if (key.compare(0, comp_prefix.size(), comp_prefix) == 0) {
      //                    assert(key.find(")>") != string::npos);
      //                    availableAtomsJson.push_back(key.substr(comp_prefix.size(), key.find(")>") - comp_prefix.size()));
      //                }
      if (mDataRanges.find(key) == mDataRanges.end()) {
        mDataRanges[key] = std::pair<float, float>(FLT_MAX, FLT_MIN);
      }
    }
  }

  // Read prim_labels file. This tells us the atoms, the vacancies and the labeling
  std::string primFileContents = readJsonPrimLabelsFile(datasetId, "");
  if (primFileContents.size() > 0 ) {
    auto primLabelsJson = json::parse(primFileContents);
    for (auto basis: primLabelsJson["basis"]) {
      //                std::cout << basis.dump() <<std::endl;
      bool isVacancy = false;
      for (auto label: basis["occupant_dof"]) {
        std::string labelString = label;
        if (labelString == "Va") {
          isVacancy = true;
          continue; // The next label should be the atom name
        }
        mAvailableAtomsJson.push_back(labelString);
        //                  if (mDataRanges.find(labelString) == mDataRanges.end()) {
        //                      mDataRanges[labelString] = std::pair<float, float>(FLT_MAX, FLT_MIN);
        //                  }
        if (std::find(mShowAtomElements.begin(), mShowAtomElements.end(), labelString) == mShowAtomElements.end()) {
          mShowAtomElements.push_back(labelString);
        }
        if (isVacancy && std::find(mVacancyAtoms.begin(), mVacancyAtoms.end(), labelString) == mVacancyAtoms.end()) {
          mVacancyAtoms.push_back(labelString);
        }
      }
    }
    mTitle = primLabelsJson["title"].get<std::string>();
  } else {
    std::cerr << "ERROR failed to find prim labels file" << std::endl;
  }
}

void DatasetManager::preProcessDataset() {
  // TODO preprocess for all datasets
  std::string datasetId = mCurrentDataset.get();

#ifdef AL_BUILD_MPI
  int world_size,nameSize, proc_id;
  char processor_name[1024];

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Get_processor_name(processor_name, &nameSize);
  MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);

  vector<unsigned int> indeces;
  int numProcessesSlice;
  if (proc_id == 0) {
    indeces.reserve(mParameterSpaces["temperature"]->size());
    for (unsigned int i = 0; i < mParameterSpaces["temperature"]->size(); i++) {
      indeces.push_back(i);
    }
    numProcessesSlice = indeces.size()/ world_size;
  }


  MPI_Bcast(&numProcessesSlice, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  unsigned int *nodeIndeces;
  nodeIndeces = (unsigned int *) malloc(numProcessesSlice * sizeof(unsigned int));

  MPI_Scatter(indeces.data(), numProcessesSlice, MPI::UNSIGNED,
              nodeIndeces, numProcessesSlice, MPI::UNSIGNED,
              0, MPI_COMM_WORLD);
  std::cout << processor_name << "--" << numProcessesSlice <<  " node indeces start:::  " << nodeIndeces[0] << std::endl;
  //        if (mParameterSpaces["chempot2"]->size() > 0) {
  //            for (auto chempot: mParameterSpaces["chempotA"]->values()) {
  //                for (auto chempot2: mParameterSpaces["chempotB"]->values()) {
  //                    for (int i = 0; i < numProcessesSlice; i++) {
  //                        int temp = nodeIndeces[i];
  //                        labelProcessor.setParams(chempot.first, chempot2.first, std::to_string(temp), datasetId);
  //                        labelProcessor.processAsync();
  //                    }
  //                }
  //            }
  //            labelProcessor.waitForAsyncDone();
  //        } else {
  //            for (auto chempot: mParameterSpaces["chempotA"]->values()) {
  //                for (int i = 0; i < numProcessesSlice; i++) {
  //                    int temp = nodeIndeces[i];
  //                    labelProcessor.setParams(chempot.first, "", std::to_string(temp), datasetId);
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
  bool valid = mParameterSpaces["chempotA"]->size() > 0 && mParameterSpaces["temperature"]->size() > 0;
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

  std::unique_lock<std::mutex> lk(mDataLock);
  auto templatePoscarName = File::conformDirectory(mLoadedDataset) + labelProcessor.outputFile();

  // Load POSCAR data
  if (reader.loadFile(templatePoscarName)) {
    mTemplatePositions = reader.getAllPositions();
  } else {
    std::cerr << "ERROR creating template at time 0 ----------- " << std::endl;
  }
  mHistory.clear();
  // Load empty template
  VASPReader emptyTemplateReader;
  if (emptyTemplateReader.loadFile(File::conformPathToOS(buildRootPath() + mCurrentDataset.get() + "/template_POSCAR"))) {
    mEmptyTemplate = emptyTemplateReader.getElementPositions("X");
  }

}

bool DatasetManager::loadDiff(int timeIndex) {
  bool diffLoaded = false;
  //         Then accumulate all diffs until current time
  int targetIndex = mParameterSpaces["time"]->getCurrentIndex();

  //      std::cout << timeIndex << " --> " << targetIndex << std::endl;
  if (timeIndex < targetIndex) {
    for (int i = timeIndex; i != targetIndex; i ++) {
      std::cout << "applying diff " << i <<std::endl;
      std::pair<Vec3f, Vec3f> historyPoint;
      auto this_diff = mDiffs[i];
      auto this_diff_indeces = this_diff[0];
      auto this_diff_labels = this_diff[1];

      for (int change = 0; change < this_diff_indeces.size(); change++) {
        int changeOffset = this_diff_indeces[change].get<int>() * 4;
        if (changeOffset < mEmptyTemplate.size()) {
          Vec3f pos(mEmptyTemplate[changeOffset], mEmptyTemplate[changeOffset + 1], mEmptyTemplate[changeOffset + 2]);

          int index = -1;
          std::string label;
          for (auto &positions: mTemplatePositions) {
            if (positions.first != this_diff_labels[change]) { // Look only for changes
              for (int curIndex = 0; curIndex < positions.second.size(); curIndex += 4 ) {
                if (pos == Vec3f(positions.second[curIndex], positions.second[curIndex + 1], positions.second[curIndex + 2])) {
                  label = positions.first;
                  index = curIndex;

                  //                      std::cout << "match " << this_diff_indeces[change] << " " << pos.x << "," << pos.y << std::endl;
                  break;
                }
              }
            }
            if (index != -1) {
              break;
            }
          }
          assert(label != "");
          if (this_diff_labels[change] == "Va") { // Remove atom
            historyPoint.first = pos;
            //                std::cout << "remove index " << index << std::endl;
            assert(index >= 0);
            mTemplatePositions[this_diff_labels[change]].insert(
                  mTemplatePositions[this_diff_labels[change]].begin(),
                mTemplatePositions[label].begin() + index,
                mTemplatePositions[label].begin() + index + 4);
            mTemplatePositions[label].erase(mTemplatePositions[label].begin() + index, mTemplatePositions[label].begin() + index + 4);
          } else { // Add Atom
            historyPoint.second = pos;
            //                std::cout << "add index " << index << std::endl;
            mTemplatePositions[this_diff_labels[change]].push_back(pos.x);
            mTemplatePositions[this_diff_labels[change]].push_back(pos.y);
            mTemplatePositions[this_diff_labels[change]].push_back(pos.z);
            mTemplatePositions[this_diff_labels[change]].push_back(0.0);
            mTemplatePositions[label].erase(mTemplatePositions[label].begin() + index, mTemplatePositions[label].begin() + index + 4);
          }
        } else {
          std::cerr << " ERROR: diff index greater than template size." << std::endl;
          break;
        }
      }
      mHistory.push_back(historyPoint);
    }
    diffLoaded = true;
  }
  return diffLoaded;
}


void DatasetManager::loadFromPOSCAR() {
  mHistory.clear();
  std::string id = getSubDir();
  std::string condition =  std::to_string(mParameterSpaces[mConditionsParameter]->getCurrentIndex());
  if (mRunProcessors && mConditionsParameter != "") {
    int timeIndex = 0;
    if (mParameterSpaces["time"]->size() > 0) {
      timeIndex = mParameterSpaces["time"]->getCurrentIndex();
    }
    labelProcessor.setRunningDirectory(mLoadedDataset);
    labelProcessor.setParams(id, condition, mCurrentDataset, std::to_string(timeIndex));

    labelProcessor.processAsync(labelProcessor.configuration(),false,[&] (bool ok) {
      if (ok) {
        currentPoscarName.set(labelProcessor.outputFile());

        // Load POSCAR data
        std::cout << " ********  Before Read:" << currentPoscarName.get() << std::endl;
        if (reader.loadFile(currentPoscarName)) {
          std::cout << "*********** Read:" << currentPoscarName.get() << std::endl;
        } else {
          std::cout << "Cannot Read:" << currentPoscarName.get() << std::endl;
        }

        auto positions = positionBuffers.getWritable();
        *positions = reader.getAllPositions();

        positionBuffers.doneWriting(positions);
//        mTemplatePositions = *positions;
//        allPositions = mTemplatePositions;
      }
    });

  }
}

void DatasetManager::getAtomPositions() {

  if (mConditionsParameter == "") { return; } // FIXME hack to avoid crashes
  std::unique_lock<std::mutex> lk(mDataLock);

  std::string id = getSubDir();
  std::string condition =  std::to_string(mParameterSpaces[mConditionsParameter]->getCurrentIndex());
  // First chek if this has time steps. If it does, check to see if we have cached diff
  std::string fullconditionPath =  File::conformPathToOS(mGlobalRoot + mRootPath.get() + mCurrentDataset.get() + "/" + id + "/conditions." + condition + "/");

  if (mParameterSpaces["time"]->size() > 0 && File::exists(fullconditionPath + "time_diffs.json")) {
    // Generate an internal template file for data at time 0
    // Load new template if more than time has changed.
    bool onlyTimeChanged = true;
    for (auto &paramSpace: mParameterSpaces) {
      if (paramSpace.first == "time") {
        continue;
      }
      if (mCurrentLoadedIndeces[paramSpace.first] != paramSpace.second->getCurrentIndex()) {
        onlyTimeChanged = false;
        break;
      }
    }
    int timeIndex = mParameterSpaces["time"]->getCurrentIndex();
    if (!onlyTimeChanged) {
      timeIndex = 0;
      // Load diffs for this parameter space sample
      std::ifstream f(fullconditionPath + "time_diffs.json");
      std::string str;;
      if (!f.fail()) {

        f.seekg(0, std::ios::end);
        str.reserve(f.tellg());
        f.seekg(0, std::ios::beg);

        str.assign((std::istreambuf_iterator<char>(f)),
                   std::istreambuf_iterator<char>());

        mDiffs = json::parse(str);

      } else {
        std::cerr << "ERROR loading diff file" << std::endl;
      }
      //          std::string timeIndexId = mParameterSpaces["time"]->idAt(0);
      labelProcessor.setRunningDirectory(mLoadedDataset);
      labelProcessor.setParams(id, condition, mCurrentDataset, "0");
      labelProcessor.processAsync(labelProcessor.configuration(),false,[&] (bool ok) {
        if (ok) {
           processTemplatePositions();
        }
      });
    } else {
      timeIndex = mCurrentLoadedIndeces["time"];
    }

    bool diffLoaded = loadDiff(timeIndex);
    if (diffLoaded) {

      auto positions = positionBuffers.getWritable();
      *positions = mTemplatePositions;

      positionBuffers.doneWriting(positions);}
    else {
      loadFromPOSCAR();

    }
  } else { // Generate POSCAR file from template if no time space and no cached diffs
    loadFromPOSCAR();

  }
}

std::vector<std::string> DatasetManager::getDataNames() {
  std::vector<std::string> names;
  for (auto entry: mDataRanges) {
    names.push_back(entry.first);
  }
  return names;
}

void DatasetManager::generateGraph(std::string xData, std::string yData, std::string datasetId, bool multi) {
  // This function is called asynchronously through graphComputationFunction()
  // that in turn is called by the deferred computation facilities
  // Get graph data
  if (mRunProcessors) {
    //            std::unique_lock<std::mutex> lk(mDataLock);

    auto jsonText = readJsonResultsFile(mCurrentDataset, getSubDir());
    if (jsonText.size() > 0) {
      std::vector<double> datax, datay;
      auto resultsJson = json::parse(jsonText);

      if (resultsJson.find(yData) != resultsJson.end()) {
        try {

          datay.reserve(resultsJson[yData].size());
          datay.insert(datay.begin(), resultsJson[yData].begin(), resultsJson[yData].end());
        } catch (std::exception &) {

        }
      }

      std::string highlightValue = "0.0";

      if(mParameterSpaces.find(xData) != mParameterSpaces.end() ) {
        try {
          datax.reserve(mParameterSpaces[xData]->size());
          for (auto value: mParameterSpaces[xData]->values()) {
            datax.push_back(value.second);
          }
          highlightValue = std::to_string(int(mParameterSpaces[xData]->getCurrentValue()));
        } catch (std::exception &) {

        }
      }

      std::string fullDatasetPath = File::conformPathToOS(buildRootPath() + File::conformPathToOS(mCurrentDataset.get()));

      std::ofstream ofsx (fullDatasetPath + "cached_output/inx.bin", std::ofstream::out | std::ios::binary);
      ofsx.write((char *) datax.data(), datax.size() * sizeof(double));
      ofsx.flush();
      ofsx.close();
      std::ofstream ofsy (fullDatasetPath + "cached_output/iny.bin", std::ofstream::out | std::ios::binary);
      ofsy.write((char *) datay.data(), datay.size() * sizeof(double));
      ofsy.flush();
      ofsy.close();

      //                labelProcessor.setRunningDirectory("cached_output");

      graphGenerator.title = mTitle;
      graphGenerator.dataset = datasetId;
      graphGenerator.temp_interest = highlightValue;
      {
        std::string subDir = getSubDir();
        graphGenerator.multi = false;
        graphGenerator.subdir = subDir;
        graphGenerator.xLabel = xData;
        graphGenerator.yLabel = yData;

        std::string str = readJsonResultsFile(datasetId, getSubDir());
        std::pair<float, float> dataRange(FLT_MAX, FLT_MIN);
        if (str.size() > 0) {
          auto resultsJson = json::parse(str);
          auto data = resultsJson[yData];
          for (json::iterator v =data.begin(); v != data.end(); ++v) {
            if (v->get<float>() < dataRange.first) {
              dataRange.first = v->get<float>();
            } else if (v->get<float>() > dataRange.second) {
              dataRange.second = v->get<float>();
            }
          }

        }

        graphGenerator.miny = std::to_string(dataRange.first);
        graphGenerator.maxy = std::to_string(dataRange.second);
        graphGenerator.inxFile = "cached_output/inx.bin";
        graphGenerator.inyFile = "cached_output/iny.bin";
      }
      graphProcessing = true;
      graphGenerator.processAsync(graphGenerator.configuration(), true, [this](bool runOk) {
        if (runOk && graphProcessing && mRunProcessors) {
          currentGraphName.set(graphGenerator.outputFile());
          graphProcessing = false;
        }
      });
    }
  }
}

DatasetManager::SpeciesLabelMap DatasetManager::getAvailableSpecies() {
  SpeciesLabelMap labelMap;
  for (std::string atom: mAvailableAtomsJson) {
    labelMap[atom] = std::vector<std::string>();
  }
  for (std::string label: mShowAtomElements) {
    for (auto mapEntry: labelMap) {
      if (label.compare(0, mapEntry.first.size(), mapEntry.first) == 0) {
        labelMap[mapEntry.first].push_back(label);
        continue;
      }
    }
  }
  labelMap["Va"] = mVacancyAtoms;
  return labelMap;
}

std::string DatasetManager::findJsonFile(std::string datasetId, std::string subDir, std::string fileName) {

  std::vector<std::string> paths = {
    buildRootPath() + "/" + datasetId + "/" + subDir + "/",
    buildRootPath() + "/" + datasetId + "/",
    buildRootPath() + "/"
  };

  for (auto possiblePath : paths) {
    //            std::cout << "Looking for " << fileName << " in " << possiblePath << std::endl;
    if (File::exists(possiblePath + fileName)) {
      return possiblePath + fileName;
    }
  }

  return std::string();
}

std::string DatasetManager::readJsonFile(std::string datasetId, std::string subDir, std::vector<std::string> fileNames) {
  std::string str;

  std::string foundFileName;
  for (auto name: fileNames) {
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
