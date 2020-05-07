#ifndef PROCESSORS_HPP
#define PROCESSORS_HPP

#include <vector>

#include "al/io/al_File.hpp"
#include "tinc/DataScript.hpp"
#include "tinc/DeferredComputation.hpp"

using namespace al;
using namespace tinc;

class AtomLabelProcessor : public DataScript {
public:
  std::string pythonScriptsPath;
  std::string root_path;
  std::string scriptName = "reassign_occs/reassign_occs.py";

  AtomLabelProcessor() {}

  void setParams(std::string id, std::string condition, std::string datasetId,
                 std::string timeStep = "") {
    mId = id;
    mCond = condition;
    mDatasetId = datasetId;
    mTimeStep = timeStep;

    std::string folder = File::conformDirectory(mId);
    std::string template_pos_path = File::conformPathToOS(
        root_path + mDatasetId + "/cached_output/template_POSCAR");

    // Try to find "template_POSCAR"
    if (!File::exists(template_pos_path)) {
      std::cout << "Failed to find template file. Call it 'template_POSCAR'."
                << std::endl;
    }
    setOutputFileNames({sanitizeName(template_pos_path + mDatasetId + "_" +
                                     folder + "_" + mCond + "_" + mTimeStep)});
  }

  //    void setFolderNamingTemplate(std::string nameTemplate) {
  //        mTemplate = nameTemplate;
  //    }

  std::map<std::string, std::string> configuration() {

    std::string folder = File::conformDirectory(mId);
    std::string subdir =
        File::conformPathToOS(folder + "conditions." + mCond + "/");
    std::map<std::string, std::string> config;

    //        config["root_path"] = File::conformDirectory(root_path);
    config["dataset_path"] = File::conformPathToOS(root_path + mDatasetId);

    std::string prim_path;
    if (File::exists(root_path + mDatasetId + "/prim_labels.json")) {
      prim_path = root_path + mDatasetId + "/prim_labels.json";
    } else if (File::exists(root_path + mDatasetId + "/prim.json")) {
      prim_path = root_path + mDatasetId + "/prim.json";
    } else if (File::exists(root_path + "prim_labels.json")) {
      prim_path = root_path + "prim_labels.json";
    } else if (File::exists(root_path + "prim.json")) {
      prim_path = root_path + "prim.json";
    }
    config["prim_path"] = prim_path;

    std::string final_state_path =
        subdir + "final_state.json"; // final_state.json file
    config["final_state_path"] =
        File::conformPathToOS(root_path + mDatasetId + "/" + final_state_path);

    config["template_path"] = File::conformPathToOS(
        root_path + mDatasetId + "/cached_output/template_POSCAR");
    config["condition"] = mCond;
    config["time_step"] = mTimeStep;

    config["__output_name"] = outputFile(false);
    config["__output_path"] = outputDirectory();

    return config;
  }

private:
  std::string mId;
  std::string mMu2;
  std::string mCond;
  std::string mDatasetId;
  std::string mTimeStep;

  //  std::string mTemplate {"mu_$"};
  //  std::string mTemplateMarker {"$"};
};

class ParameterSpaceProcessor : public DataScript {
public:
  virtual void configure() override {
    clearFlags();
    appendFlag(scriptFile(false), FLAG_SCRIPT);
    appendFlag(runningDirectory(), FLAG_OUTPUT_DIR);
    appendFlag("cached_output/_parameter_space.json", FLAG_OUTPUT_NAME);
  }
};

class GraphGenerator : public DataScript {
public:
  // Parameters
  std::string pythonScriptsPath;
  std::string title;
  std::string dataset;
  std::string temp_interest;
  std::string subdir;

  // Single plot
  std::string xLabel;
  std::string yLabel;
  std::string miny;
  std::string maxy;
  std::string inxFile;
  std::string inyFile;

  // Multiplot
  std::string casm_path;
  std::string subpath1;
  std::string subpath2;
  std::string path_template;

  bool multi{false};

  GraphGenerator() {}

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

  std::map<std::string, std::string> configuration() {
    mFlags.clear();
    std::map<std::string, std::string> config;
    std::string yLabelSanitized = yLabel;
    replaceAll(yLabelSanitized, "<", "_");
    replaceAll(yLabelSanitized, ">", "_");
    replaceAll(yLabelSanitized, "(", "_");
    replaceAll(yLabelSanitized, ")", "_");

    config["__output_name"] = dataset + "_" + temp_interest + "_" +
                              sanitizeName(subdir) + "_" + yLabelSanitized +
                              "_graph.png";
    appendFlag(config["__output_name"], FLAG_OUTPUT_NAME);

    config["__output_path"] = outputDirectory();

    if (multi) {
      //            appendFlag(pythonScriptsPath + "/graphing/plot_multi.py",
      //            FLAG_SCRIPT); appendFlag(casm_path); appendFlag(subpath1);
      //            appendFlag(subpath2);
      //            appendFlag(path_template);
      //            appendFlag(temp_interest);
      //            appendFlag(title);
      //            appendFlag(outputDirectory(), FLAG_OUTPUT_DIR);
      //            appendFlag(title + "_" + temp_interest + "_graph-multi.png",
      //            FLAG_OUTPUT_NAME);
    } else {
      config["temp_interest"] = temp_interest;
      config["xLabel"] = xLabel;
      config["yLabel"] = yLabel;
      config["miny"] = miny;
      config["maxy"] = maxy;
      config["inxFile"] = inxFile;
      config["inyFile"] = inyFile;
    }
    return config;
  }

private:
};

class TemplateGenerator : public DataScript {
public:
  // Parameters
  std::string outName = "template_POSCAR";

  std::map<std::string, std::string> configuration() {

    std::map<std::string, std::string> config;
    if (File::exists(runningDirectory() + "prim.json")) {
      config["prim_path"] = "prim.json";
    } else {
      config["prim_path"] = "../prim.json";
    }
    if (File::exists(runningDirectory() + "cached_output/transfmat")) {
      config["transfmat"] = "cached_output/transfmat";
    } else {
      config["transfmat"] = "../transfmat";
    }

    config["__output_name"] = sanitizeName(outName);
    config["__output_path"] = outputDirectory();

    return config;
  }
};

class DiffGenerator : public DataScript {
public:
  // Parameters

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
