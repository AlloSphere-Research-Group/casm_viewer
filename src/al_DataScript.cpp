
#include <array>
#include <atomic>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>  // setprecision
#include <iostream>
#include <mutex>
#include <thread>
#include <utility>  // For pair

#include "al_DataScript.hpp"

#if defined(AL_OSX) || defined(AL_LINUX) || defined(AL_EMSCRIPTEN)
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#elif defined(AL_WINDOWS)
#include <Windows.h>
#include <direct.h>  // for _chdir() and _getcwd()
#define chdir _chdir
#define getcwd _getcwd

#define popen _popen
#define pclose _pclose
#include <fileapi.h>
#endif

using namespace al;

constexpr auto DATASCRIPT_META_FORMAT_VERSION = 0;

PushDirectory::PushDirectory(std::string directory, bool verbose)
    : mVerbose(verbose) {
  mDirectoryLock.lock();
  getcwd(previousDirectory, 512);
  chdir(directory.c_str());
  if (mVerbose) {
    std::cout << "Pushing directory: " << directory << std::endl;
  }
}

PushDirectory::~PushDirectory() {
  chdir(previousDirectory);
  if (mVerbose) {
    std::cout << "Setting directory back to: " << previousDirectory
              << std::endl;
  }
  mDirectoryLock.unlock();
}

// --------------------------------------------------

std::string DataScript::scriptFile(bool fullPath) {
  std::string scriptName =
      mScriptName;  // Set to this by default. Will be overriden if a flag has
                    // been set for FLAG_SCRIPT
  std::string scriptPath = mRunningDirectory;
  for (auto &flag : mFlags) {
    if (flag.type == FLAG_SCRIPT && scriptName == "") {
      scriptName = flag.flagText;
    } else if (flag.type == FLAG_INPUT_DIR && scriptPath == "" && fullPath) {
      //                scriptPath = flag.flagText;
      //                if (scriptPath.size() > 0 && scriptPath.back() != '/') {
      //                    scriptPath += "/";
      //                }
    }
  }
  return scriptName;
}

std::string DataScript::inputFile(bool fullPath, int index) {
  std::string inputName;
  std::string inputPath;
  for (auto &flag : mFlags) {
    if (flag.type == FLAG_INPUT_NAME) {
      if (index == 0) {
        inputName = flag.flagText;
      }
      index--;
    } else if (flag.type == FLAG_INPUT_DIR && inputPath == "" && fullPath) {
      inputPath = flag.flagText;
      if (inputPath.size() > 0 && inputPath.back() != '/') {
        inputPath += "/";
      }
    }
  }
  return inputPath + inputName;
}

std::string DataScript::outputFile(bool fullPath, int index) {
  std::string outputName;
  if (mOutputFileNames.size() > index) {
    outputName = mOutputFileNames[index];
  }
  if (fullPath) {
    return outputDirectory() + outputName;
  } else {
    return outputName;
  }
}

void DataScript::configure() {
  bool hasOutputFlag = false;
  for (auto &flag : mFlags) {
    if (flag.type == FLAG_OUTPUT_NAME) {
      hasOutputFlag = true;
      break;
    }
  }
  if (!hasOutputFlag && mOutputFileNames.size() > 0) {
    appendFlag(mOutputFileNames[0], FLAG_OUTPUT_NAME);
  }
}

bool DataScript::process(bool forceRecompute) {
  bool ok = true;
  std::unique_lock<std::mutex> lk(mProcessingLock);
  configure();
  if (needsRecompute() || forceRecompute) {
    std::string command = makeCommandLine();
    ok = runCommand(command);
    if (ok) {
      PushDirectory p(mRunningDirectory, mVerbose);
      writeMeta();
    }
  } else {
    if (mVerbose) {
      std::cout << "No need to update cache according to " << metaFilename()
                << std::endl;
    }
  }
  if (mDoneCallback) {
    mDoneCallback(ok);
  }
  return ok;
}

bool DataScript::process(std::map<std::string, std::string> options,
                         bool forceRecompute, std::string outputDir,
                         std::string outputName, std::string inputDir,
                         std::string inputName) {
  std::unique_lock<std::mutex> lk(mProcessingLock);
  if (mScriptName == "" || mScriptCommand == "") {
    std::cout << "ERROR: process() missing script name or script command."
              << std::endl;
    return false;
  }
  // Using the flags structure to store the data
  clearFlags();
  appendFlag(mScriptName, FLAG_SCRIPT);
  if (outputDir != "") {
    setOutputDirectory(outputDir);
    appendFlag(outputDir, FLAG_OUTPUT_DIR);
  } else if (options.find("__output_dir") != options.end()) {
    appendFlag(options["__output_dir"], FLAG_OUTPUT_DIR);
  }
  if (outputName != "") {
    setOutputFileNames({outputName});
  } else if (options.find("__output_name") != options.end()) {
    setOutputFileNames({options["__output_name"]});
  }
  if (inputDir != "") {
    appendFlag(inputDir, FLAG_INPUT_DIR);
  } else if (options.find("__input_dir") != options.end()) {
    appendFlag(options["__input_dir"], FLAG_INPUT_DIR);
    inputDir = options["__input_dir"];
  }
  if (inputName != "") {
    appendFlag(inputName, FLAG_INPUT_NAME);
  } else if (options.find("__input_name") != options.end()) {
    appendFlag(options["__input_name"], FLAG_INPUT_NAME);
  }

  using json = nlohmann::json;
  json j;
  j["__output_dir"] = outputDirectory();
  j["__output_name"] = outputFile(false);
  j["__input_dir"] = inputDir;
  j["__input_name"] = inputFile(false);
  for (auto &option : options) {
    j[option.first] = option.second;
  }

  std::string jsonFilename = "_" + sanitizeName(mRunningDirectory) +
                             std::to_string(long(this)) + "_config.json";
  if (mVerbose) {
    std::cout << "Writing json: " << jsonFilename << std::endl;
  }
  {
    PushDirectory p(mRunningDirectory, mVerbose);
    std::ofstream of(jsonFilename, std::ofstream::out);
    if (of.good()) {
      of << j.dump(4);
      of.close();
      if (!of.good()) {
        std::cout << "Error writing json file." << std::endl;
        return false;
      }
    } else {
      std::cout << "Error writing json file." << std::endl;
      return false;
    }
  }

  bool ok = true;
  if (needsRecompute() || forceRecompute) {
    std::string command =
        mScriptCommand + " \"" + mScriptName + "\" \"" + jsonFilename + "\"";
    ok = runCommand(command);
    if (ok) {
      writeMeta();
    }
  } else {
    if (mVerbose) {
      std::cout << "No need to update cache according to " << metaFilename()
                << std::endl;
    }
  }
  if (mDoneCallback) {
    mDoneCallback(ok);
  }
  return ok;
}

std::string DataScript::sanitizeName(std::string output_name) {
  std::replace(output_name.begin(), output_name.end(), '/', '_');
  std::replace(output_name.begin(), output_name.end(), '.', '_');
  std::replace(output_name.begin(), output_name.end(), ':', '_');
  std::replace(output_name.begin(), output_name.end(), '\\', '_');
  return output_name;
}

bool DataScript::processAsync(bool noWait,
                              std::function<void(bool)> doneCallback) {
  std::lock_guard<std::mutex> lk(mProcessingLock);
  configure();
  if (needsRecompute()) {
    std::string command = makeCommandLine();
    while (mNumAsyncProcesses.fetch_add(1) > mMaxAsyncProcesses) {
      mNumAsyncProcesses--;
      if (noWait) {
        return false;  // Async process not started
      }
      std::unique_lock<std::mutex> lk2(mAsyncDoneTriggerLock);
      std::cout << "Async waiting 2 " << mNumAsyncProcesses << std::endl;
      mAsyncDoneTrigger.wait(lk2);
      std::cout << "Async done waiting 2 " << mNumAsyncProcesses << std::endl;
    }
    //            std::cout << "Async " << mNumAsyncProcesses << std::endl;
    mAsyncThreads.emplace_back(std::thread([this, command, doneCallback]() {
      bool ok = runCommand(command);

      if (doneCallback) {
        doneCallback(ok);
      }
      mNumAsyncProcesses--;
      mAsyncDoneTrigger.notify_all();
      //                std::cout << "Async runner done" << mNumAsyncProcesses
      //                << std::endl;
    }));

    PushDirectory p(mRunningDirectory, mVerbose);
    writeMeta();
  } else {
    if (doneCallback) {
      doneCallback(true);
    }
  }
  return true;
}

bool DataScript::processAsync(std::map<std::string, std::string> options,
                              bool noWait,
                              std::function<void(bool)> doneCallback) {
  std::unique_lock<std::mutex> lk(mProcessingLock);
  configure();
  if (needsRecompute()) {
    std::string jsonFilename =
        "_" + sanitizeName(mRunningDirectory) + "_config.json";
    std::string command =
        mScriptCommand + " " + mScriptName + " " + jsonFilename + "";
    lk.unlock();
    while (mNumAsyncProcesses.fetch_add(1) > mMaxAsyncProcesses) {
      mNumAsyncProcesses--;
      if (noWait) {
        return false;  // Async process not started
      }
      std::unique_lock<std::mutex> lk2(mAsyncDoneTriggerLock);
      std::cout << "Async waiting " << mNumAsyncProcesses << std::endl;
      mAsyncDoneTrigger.wait(lk2);
      std::cout << "Async done waiting " << mNumAsyncProcesses << std::endl;
    }
    //            std::cout << "Async " << mNumAsyncProcesses << std::endl;

    if (mVerbose) {
      std::cout << "Starting asyc thread" << std::endl;
    }
    mAsyncThreads.emplace_back(
        std::thread([this, options, command, doneCallback]() {
          bool ok = process(options, true);

          PushDirectory p(mRunningDirectory, mVerbose);
          writeMeta();

          if (doneCallback) {
            doneCallback(ok);
          }
          mNumAsyncProcesses--;
          mAsyncDoneTrigger.notify_all();
          //                std::cout << "Async runner done" <<
          //                mNumAsyncProcesses << std::endl;
        }));
  } else {
    doneCallback(true);
  }
  return true;
}

bool DataScript::runningAsync() {
  if (mNumAsyncProcesses > 0) {
    return true;
  } else {
    return false;
  }
}

bool DataScript::waitForAsyncDone() {
  bool ok = true;
  for (auto &t : mAsyncThreads) {
    t.join();
  }
  return ok;
}

void DataScript::appendFlag(std::string flagText, al::FlagType type) {
  if (type == FLAG_OUTPUT_DIR) {
    setOutputDirectory(flagText);
  } else if (type == FLAG_OUTPUT_NAME) {
    mOutputFileNames.push_back(flagText);
  }
  mFlags.push_back({flagText, type});
}

std::string DataScript::makeCommandLine() {
  std::string commandLine = mScriptCommand + " ";
  for (auto flag : mFlags) {
    commandLine += flag.flagText + " ";
  }
  return commandLine;
}

bool DataScript::runCommand(const std::string &command) {
  PushDirectory p(mRunningDirectory, mVerbose);

  if (mVerbose) {
    std::cout << "DataScript command: " << command << std::endl;
  }
  std::array<char, 128> buffer{0};
  std::string output;
  FILE *pipe = popen(command.c_str(), "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  while (!feof(pipe)) {
    if (fgets(buffer.data(), 128, pipe) != nullptr) {
      output += buffer.data();
      if (mVerbose) {
        std::cout << buffer.data() << std::endl;
      }
    }
  }

  int returnValue = 0;
  // TODO: Put back return value checking. Currently crashing
  //        int returnValue = -1;
  if (!ferror(pipe)) {
    returnValue = pclose(pipe);
  } else {
    returnValue = -1;
  }

  if (mVerbose) {
    std::cout << "Script result: " << returnValue << std::endl;
    std::cout << output << std::endl;
  }
  return returnValue == 0;
}

void DataScript::writeMeta() {
  std::ofstream metaFileStream;

  metaFileStream.open(metaFilename(), std::ofstream::out);

  metaFileStream << DATASCRIPT_META_FORMAT_VERSION << std::endl;
  for (auto &flag : mFlags) {
    FlagType type = flag.type;
    if (type == FLAG_INPUT_NAME) {
      metaFileStream << (int)flag.type << std::endl;
      metaFileStream << std::setprecision(12) << modified(inputFile().c_str())
                     << std::endl;
    } else if (type == FLAG_SCRIPT) {
      metaFileStream << (int)flag.type << std::endl;
      metaFileStream << std::setprecision(12) << modified(scriptFile().c_str())
                     << std::endl;
    }
  }

  metaFileStream.close();
  if (mVerbose) {
    std::cout << "Wrote cache in: " << metaFilename() << std::endl;
  }
}

al_sec DataScript::modified(const char *path) const {
  struct stat s;
  if (::stat(path, &s) == 0) {
    // const auto& t = s.st_mtim;
    // return t.tv_sec + t.tv_usec/1e9;
    return s.st_mtime;
  }
  return 0.;
}

bool DataScript::needsRecompute() {
  std::ifstream metaFileStream;
  metaFileStream.open(metaFilename(), std::ofstream::in);

  if (metaFileStream.fail()) {
    if (mVerbose) {
      std::cout << "Failed to open metadata: Recomputing. " << metaFilename()
                << std::endl;
    }
    return true;
  }

  std::string line;
  std::getline(metaFileStream, line);
  if (line != std::to_string(DATASCRIPT_META_FORMAT_VERSION)) {
    if (mVerbose) {
      std::cout << "Metadata format mismatch. Forcing recompute" << std::endl;
    }
    metaFileStream.close();
    return true;
  }
  while (std::getline(metaFileStream, line)) {
    FlagType type;
    int typeInt;
    std::istringstream(line) >> typeInt;
    type = (FlagType)typeInt;
    if (type == FLAG_INPUT_NAME) {
      al_sec time;
      std::getline(metaFileStream, line);
      std::istringstream(line) >> time;
      al_sec lastModified = modified(inputFile().c_str());
      if (time != lastModified) {
        metaFileStream.close();
        return true;
      }
    } else if (type == FLAG_SCRIPT) {
      al_sec time;
      std::getline(metaFileStream, line);
      std::istringstream(line) >> time;
      al_sec lastModified = modified(scriptFile().c_str());
      if (time != lastModified) {
        metaFileStream.close();
        return true;
      }
    }
  }
  metaFileStream.close();
  if (!File::exists(outputFile())) {
    return true;
  }
  return false;
}

std::string DataScript::metaFilename() {
  std::string outPath;
  std::string outName;
  for (auto &flag : mFlags) {
    if (flag.type == FLAG_OUTPUT_DIR) {
      outPath = flag.flagText;
    }
    if (flag.type == FLAG_OUTPUT_NAME) {
      outName = flag.flagText;
    }
  }
  if (outName == "") {
    outName = outputFile(false);
  }
  if (outPath == "") {
    outPath = outputDirectory();
  }
  std::string metafilename = File::conformPathToOS(outPath) + outName + ".meta";
  return metafilename;
}
