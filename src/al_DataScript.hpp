
#ifndef INCLUDE_AL_DATASCRIPT
#define INCLUDE_AL_DATASCRIPT

#include <string>
#include <vector>
#include <mutex>

#include "al/util/ui/al_Parameter.hpp"
#include "al/core/io/al_File.hpp"
#include "al/util/ui/al_ParameterServer.hpp"

#include "json.hpp"

namespace al {

enum FlagType {
    FLAG_GENERIC = 0,
    FLAG_SCRIPT,  // The script to be run
    FLAG_INPUT_DIR,
    FLAG_INPUT_NAME,  // Should be relative to input dir. If file has changed, this forces recompute
    FLAG_OUTPUT_DIR,
    FLAG_OUTPUT_NAME
};

struct Flag {
    std::string flagText;
    FlagType type {FLAG_GENERIC};
};

static std::mutex mDirectoryLock; // Protects all instances of PushDirectory

class PushDirectory {
public:
  PushDirectory(std::string directory, bool verbose = false);

  ~PushDirectory();
private:
  char previousDirectory[512];
  bool mVerbose;

};

class Processor {
public:

  /**
   * @brief Set the directory for output files
   * @param outputDirectory
   */
  void setOutputDirectory(std::string outputDirectory) {
      mOutputDirectory = File::conformPathToOS(outputDirectory);
      std::replace(mOutputDirectory.begin(), mOutputDirectory.end(), '\\', '/');
      if (!File::isDirectory(mOutputDirectory)) {
          if (!Dir::make(mOutputDirectory)) {
              std::cout <<  "Unable to create cache directory:" << mOutputDirectory << std::endl;
          }
      }
  }

  /**
   * @brief Set the current directory for process to run in.
   * @param directory
   */
  void setRunningDirectory(std::string directory) {
      mRunningDirectory = File::conformPathToOS(directory);
      std::replace(mRunningDirectory.begin(), mRunningDirectory.end(), '\\', '/');
      if (!File::exists(mRunningDirectory)) {
          if (!Dir::make(mRunningDirectory)) {
              std::cout << "Error creating directory: " << mRunningDirectory << std::endl;
          }
      }
  }


  /**
   * @brief Set the names of output file
   * @param list of output file names.
   */
  void setOutputFileNames(std::vector<std::string> outputFiles) {
    mOutputFileNames.clear();
    for (auto fileName : outputFiles) {

      auto name = File::conformPathToOS(fileName);
      std::replace(mOutputDirectory.begin(), mOutputDirectory.end(), '\\', '/');
      // FIXME this is not being used everywhere it should be....
      mOutputFileNames.push_back(name);
    }
  }

protected:

  std::string mRunningDirectory;
  std::string mOutputDirectory { "cached_output/"};
  std::vector<std::string> mOutputFileNames;
};

/**
 * @brief The DataScript class
 *
 * You can use the DataScript class in two ways. First, you can use the
 * appendFlag() to add sequential flags to pass to the
 * script on the command line. The simplest usage looks like:
 *
 * DataScript ds("output_folder");
 * ds.setCommand("python");
 * ds.setRunningDirectory("/home/sweet/home");
 * ds.appendFlag("myscript.py", FLAG_SCRIPT);
 * ds.appendFlag("option");
 * ds.appendFlag("3.14159");
 * ds.process();
 *
 * You can override the configure() function for better control of the
 * passed arguments. This is a quick and easy way to set things up
 * and will work well when there is a small known set of parameters. You can
 * define any additional functions to get and set configuration options.
 *
 * @code
 *
 * class MyScript : public DataScript {
 * public:
 *     float value;
 *     void configure() override {
 *         setCommand("python");
 *         setRunningDirectory("/home/sweet/home");
 *         clearFlags();
 *         appendFlag("myscript.py", FLAG_SCRIPT);
 *         appendFlag("output_dir", FLAG_OUTPUT_DIR);
 *         appendFlag("option");
 *         appendFlag("3.14159");
 *         appendFlag(std::to_string(value));
 *     }
 * };
 *
 * void process() {
 *     DataScript ds;
 *     ds.value = 0.5;
 *     ds.process();
 *     std::string
 * }
 *
 * @endcode
 *
 * This will execute the command:
 *
 * @codeline
 * python myscript.py output_dir option 3.14159 0.5
 *
 * In the /home/sweet/home directory.
 *
 * For more complex scenarios, use the setConfiguration() function. This
 * will create a json file with the options that can then be read by the script
 * This provides the greatest flexibility and extensibility.
 */
class DataScript : public Processor {
public:
    DataScript(std::string outputDirectory = "cached_output/") {
        setOutputDirectory(outputDirectory);
    }

    virtual ~DataScript() {
       waitForAsyncDone();
    }

    /**
     * @brief Set the script's main command (e.g. python)
     */
    void setCommand(std::string command) {
        mScriptCommand = command;
    }

    /**
     * @brief Set name of script to be run
     * @param scriptName
     *
     * This name will be used unless a flag is set with type FLAG_SCRIPT
     */
    void setScriptName(std::string scriptName) {
        std::replace(mScriptName.begin(), mScriptName.end(), '\\', '/');
        mScriptName = scriptName;
    }


    std::string scriptFile(bool fullPath = false);

    std::string inputFile(bool fullPath = true, int index = 0);

    std::string outputFile(bool fullPath = true, int index = 0);

    /**
    * @brief configure internal values prior to run
    *
    * Override this function to prepare command line flags.
    * This function is called internally by process()
    */
    virtual void configure();

    bool process(bool forceRecompute = false);

    /**
     * @brief execute script with configuration options provided
     * @param options
     * @param forceRecompute
     * @param outputDir
     * @param outputName
     * @param inputDir
     * @param inputName
     * @return false if there was any internal error
     *
     * Alternate method for processing. This discards any flags that have been
     * added with appendFlags and will generate a json file from the options
     * argument whose name is passed to the python script.
     *
     * You must call setScriptName() prior to any call of this function.
     */
    bool process(std::map<std::string, std::string> options,
                 bool forceRecompute = false,
                 std::string outputDir = "", std::string outputName = "",
                 std::string inputDir = "", std::string inputName = ""
                 );

    /**
     * @brief Cleans a name up so it can be written to disk
     * @param output_name the source name
     * @return the cleaned up name
     *
     * This function will remove any characters not allowed by the operating
     * system like :,/,\ etc. And will remove any characters like '.' that
     * can confuse the parsing of the name on read.
     */
    std::string sanitizeName(std::string output_name);

    bool processAsync(bool noWait = false, std::function<void(bool)> doneCallback = nullptr);

    bool processAsync(std::map<std::string, std::string> options, bool noWait = false, std::function<void(bool)> doneCallback = nullptr);

    bool runningAsync();

    bool waitForAsyncDone();

    void maxAsyncProcesses(int num) {
        mMaxAsyncProcesses = num;
    }

    void clearFlags() { mFlags.clear();}

    void appendFlag(std::string flagText, FlagType type = FLAG_GENERIC);

    void verbose(bool verbose = true) {
        mVerbose = true;
    }

    void registerDoneCallback() { throw "Not implemented yet."; }

    std::string outputDirectory() { return mOutputDirectory;}

    std::string runningDirectory() { return mRunningDirectory;}

    DataScript &registerParameter(Parameter &param) { mParameters.push_back(&param); return *this;}

    DataScript &operator << (Parameter& newParam){ return registerParameter(newParam); }

    DataScript &registerParameterServer(ParameterServer &server) { mParamServer= &server; return *this;}

    DataScript &operator << (ParameterServer& server){ return registerParameterServer(server); }

    DataScript &operator << (std::string flag){ appendFlag(flag); return *this; }

protected:
    // These need to be accessible by the subclass
    std::vector<Flag> mFlags;
    std::vector<Parameter *> mParameters;
    ParameterServer *mParamServer;

private:
    std::string mScriptCommand {"/usr/bin/python3"};
    std::string mScriptName;

    std::mutex mProcessingLock;
    bool mVerbose;
    int mMaxAsyncProcesses {4};
    std::atomic<int> mNumAsyncProcesses {0};
    std::vector<std::thread> mAsyncThreads;
    std::thread mAsyncDoneThread;
    std::condition_variable mAsyncDoneTrigger;
    std::mutex mAsyncDoneTriggerLock;

    std::function<void(bool ok)> mDoneCallback;

    std::string makeCommandLine();

    bool runCommand(const std::string &command);

    void writeMeta();

    al_sec modified(const char *path) const;

    bool needsRecompute();

    std::string metaFilename();

};

class CacheManager {
public:

  void registerProcessor(Processor &processor) {
    mProcessors.push_back(&processor);
  }

  void setRunningDirectory(std::string outputDirectory) {
    for (auto *processor: mProcessors) {
      processor->setRunningDirectory(outputDirectory);
    }

  }

  void setOutputDirectory(std::string outputDirectory) {
    for (auto *processor: mProcessors) {
      processor->setOutputDirectory(outputDirectory);
    }

  }

  void clearCache() {
      throw "Mot implemented yet.";
      //TODO implement clear cache
  }

private:
  std::vector<Processor *> mProcessors;


};

class ParallelProcessor : public DataScript
{
public:
    // TODO is this worth finishing?

    void processSpace(const std::vector<std::vector<std::string> > &allVecs, size_t vecIndex, std::vector<size_t> indeces = std::vector<size_t>())
    {
        if (vecIndex >= allVecs.size()) {
            return;
        }
        if (indeces.size() == 0) {
            indeces.resize(allVecs.size(), 0);
        }
        std::vector<std::string> currentValues(allVecs.size());
        auto indecesIt = indeces.begin();
        auto currentValuesIt = currentValues.begin();
        for (auto allVecsIt = allVecs.begin(); allVecsIt != allVecs.end(); allVecsIt++ ){
            *currentValuesIt = (*allVecsIt)[*indecesIt];
            indecesIt++;
            currentValuesIt++;
        }
        for (size_t i=0; i<allVecs[vecIndex].size(); i++) {
            indeces[vecIndex] = i;
            auto value = allVecs[vecIndex][i];
//            labelProcessor.setParams(chempot, std::to_string(i), datasetId, mConfig.pythonScriptPath);
//            labelProcessor.process();
            processSpace(allVecs, vecIndex+1, indeces);
        }
    }

    void start(std::vector<std::vector<std::string>> parameterSpace,
               std::vector<std::pair<unsigned, unsigned>> parameterRanges) {
//        mParallelProcess = std::make_shared<std::thread>(
//                    [this](parameterSpace,parameterRanges) {
//                for (auto parameterValues: parameterSpace) {
//                auto interator =
//    }

//    });
    }

    void waitForEnd() {

    }

    void stop() {

    }


private:

    std::shared_ptr<std::thread> mParallelProcess;
    std::atomic<bool> mRunning;
};


} // namespace al

#ifdef AL_WINDOWS
#undef popen
#undef pclose
#endif

#endif // INCLUDE_AL_DATASCRIPT
