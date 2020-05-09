
#ifdef AL_WINDOWS
#define NOMINMAX
#include <Windows.h>
#undef far
#undef near
#undef DELETE
#endif

#include <condition_variable>

#include "datadisplay.hpp"
#include "imgui.h"

void DataDisplay::init() {
#ifdef AL_BUILD_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mWorldRank);
#endif

  mPickableManager << graphPickable << parallelPickable << perspectivePickable;

  graphPickable.bb.set(Vec3f(-1.2f, -0.9f, -0.05f), Vec3f(1.2f, 0.9f, 0.0f));
  graphPickable.pose = Pose(Vec3d(-2.0, 0.8, -.75), Quatf());
  graphPickable.scale = 0.2f;

  parallelPickable.bb.setCenterDim(Vec3f(), Vec3f(4, 4, 0.1f));
  parallelPickable.pose = Pose(Vec3d(-0.68, 1.1, -1.5), Quatf());
  parallelPickable.scale = 0.4f;

  perspectivePickable.pose =
      Pose(Vec3d(0.55, 0.35, -0.7), Quatf().fromEuler(75.0 * M_DEG2RAD, 0, 0));
  perspectivePickable.scale = 0.02f;
  //    perspectivePickable.scaleVec = mPerspectiveScale.get();
  // perspectivePickable.addChild(rh);
  // rh.size = 25;
  // rh.dr = 2.2;

  perspectivePickable.addChild(slicePickable);
  perspectivePickable.testChildren = false;
  perspectivePickable.containChildren = true;

  //        mPerspectiveScale.registerChangeCallback([this](float s){
  //        perspectivePickable.scaleVec.set(s); });
  mPerspectiveRotY.registerChangeCallback([this](float d) {
    perspectivePickable.pose = Pose(perspectivePickable.pose.get().pos(),
                                    Quatf().fromEuler(-d * M_DEG2RAD, 0, 0));
  });

  // Initialization of graphics objects
  // #INSTANCED_RENDERING: init
  //  {
  //    addSphere(orthoMesh0);
  //    orthoMesh0.update();
  //  }

  EasyFBOSetting setting;
  setting.mUseMipmap = true;
  setting.filterMin = GL_LINEAR_MIPMAP_LINEAR;
  setting.filterMag = GL_LINEAR_MIPMAP_LINEAR;

  // Settings for parallel projection rendering fbo
  fbo_iso.init(2048, 2048, setting);

  mGraphTexture.mipmap(true);

  int num_verts_added;
  Mat4f transform;

  // Prepare axis mesh
  // x
  num_verts_added = addCube(axis);
  transform.setIdentity();
  transform *= Matrix4f::rotation(M_PI / 2, 2, 0); // rotate from z to x
  transform *= Matrix4f::translation(0, 0, 0.3);
  transform *= Matrix4f::scaling(0.02, 0.02, 1);
  axis.transform(transform, axis.vertices().size() - num_verts_added);
  for (int i = 0; i < num_verts_added; i += 1) {
    axis.color(1, 0, 0);
  }

  // y
  num_verts_added = addCube(axis);
  transform.setIdentity();
  transform *= Matrix4f::rotation(M_PI / 2, 2, 1); // rotate from z to y
  transform *= Matrix4f::translation(0, 0, 0.3);
  transform *= Matrix4f::scaling(0.02, 0.02, 1);
  axis.transform(transform, axis.vertices().size() - num_verts_added);
  for (int i = 0; i < num_verts_added; i += 1) {
    axis.color(0, 1, 0);
  }

  // z
  num_verts_added = addCube(axis);
  transform.setIdentity();
  transform *= Matrix4f::translation(0, 0, 0.3);
  transform *= Matrix4f::scaling(0.02, 0.02, 1);
  axis.transform(transform, axis.vertices().size() - num_verts_added);
  for (int i = 0; i < num_verts_added; i += 1) {
    axis.color(0, 0, 1);
  }
  axis.update();

  // Parameters setup

  for (auto &parameterSpace : mDatasetManager.mParameterSpaces) {
    parameterSpace.second->parameter().registerChangeCallback([&](float value) {
      if (parameterSpace.second->getCurrentId() !=
          parameterSpace.second->idAt(parameterSpace.second->getIndexForValue(
              value))) { // Only reload if id has changed

        parameterSpace.second->parameter().setNoCalls(
            value); // To have the internal value already changed for the
                    // following functions.
        std::cout << value << " : " << parameterSpace.second->getCurrentId()
                  << "..."
                  << parameterSpace.second->idAt(
                         parameterSpace.second->getIndexForValue(value));
        mDatasetManager.getAtomPositions();
        updateText();
      }
    });
  }

  // We need to process the callbacks for this within the graphics thread.
  mDatasetManager.mCurrentDataset.setSynchronousCallbacks();
  mDatasetManager.mCurrentDataset.registerChangeCallback(
      [&](std::string value) {
        //        if (value != mDatasetManager.mCurrentDataset.get()) {
        // First force application of value
        mDatasetManager.mCurrentDataset.setLocking(value);
        initRootDirectory();
        resetSlicing();
        //        }
      });

  mShowAtoms.registerChangeCallback([this](uint16_t value) {
    if (mShowAtoms.get() != value) {
      mDatasetManager.getAtomPositions();
    }
  });

  atomrender.mSlicingPlaneThickness.registerChangeCallback([this](float v) {
    auto m = slicePickable.bb.max;
    m.z = v;
    slicePickable.bb.set(slicePickable.bb.min, m);
  });

  slicePickable.pose.registerChangeCallback([this](Pose pose) {
    atomrender.mSlicingPlanePoint.setNoCalls(pose.pos());
    // mSlicingPlaneNormal.setNoCalls();
  });

  backgroundColor.setHint("showAlpha", 1.0);

  mGridType.setElements({"square", "triangle"});
  mShowGrid.registerChangeCallback([this](float value) {
    if (value == 0.0f) {
      mGridType.setHint("hide", 1.0);
      mGridSpacing.setHint("hide", 1.0);
      mGridXOffset.setHint("hide", 1.0);
      mGridYOffset.setHint("hide", 1.0);
    } else {
      mGridType.setHint("hide", 0.0);
      mGridSpacing.setHint("hide", 0.0);
      mGridXOffset.setHint("hide", 0.0);
      mGridYOffset.setHint("hide", 0.0);
    }
  });
  mShowGrid = false;

  // TODO we don't need to do a full data load here, just recompute the graph
  mPlotYAxis.registerChangeCallback([this](float value) {
    if (mPlotYAxis.get() != value) {
      mDatasetManager.getAtomPositions();
    }
  });
  // TODO we don't need to do a full data load here, just recompute the graph
  mPlotXAxis.registerChangeCallback([this](float value) {
    if (mPlotXAxis.get() != value) {
      mDatasetManager.getAtomPositions();
    }
  });

  mDatasetManager.currentGraphName.setSynchronousCallbacks();
  mDatasetManager.currentGraphName.registerChangeCallback(
      [this](std::string value) {
        std::string fullDatasetPath = File::conformPathToOS(
            mDatasetManager.buildRootPath() +
            File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));

        std::cout << "loading graph " << value << " at " << fullDatasetPath
                  << std::endl;
        // New image module puts origin on top right
        //        if (mGraphTextureLock.try_lock()) {
        //            if (mGraphFilePathToLoad.size() > 0) {
        Image img(value);
        if (!img.loaded()) {
          static std::string lastFailed = "";
          if (lastFailed != value) {
#ifdef AL_BUILD_MPI
            char name[MPI_MAX_PROCESSOR_NAME];
            int len;
            int ret = MPI_Get_processor_name(name, &len);
            std::cout << name << ": ";
#endif
            cout << "failed to load image " << value << endl;
            lastFailed = value;
          }
          mGraphTexture.resize(4, 4);
        } else {
          mGraphTexture.resize(img.width(), img.height());
          mGraphTexture.submit(img.pixels(), GL_RGBA, GL_UNSIGNED_BYTE);

          mGraphTexture.filter(Texture::LINEAR_MIPMAP_LINEAR);
        }
      });

  mLabelFont.alignLeft();

  bundle << mDatasetManager.mRootPath;
  bundle << mDatasetManager.mCurrentDataset;
  bundle << mVisible;
  bundle << mDatasetManager.mParameterSpaces["temperature"]->parameter();
  bundle << mDatasetManager.mParameterSpaces["chempotA"]->parameter()
         << mDatasetManager.mParameterSpaces["chempotB"]->parameter();
  bundle << mDatasetManager.mParameterSpaces["time"]->parameter();
  bundle << atomrender.mAtomMarkerSize;

  bundle << mShowAtoms;
  bundle << mPlotXAxis;
  bundle << mPlotYAxis;
  bundle << mShowGraph << mShowParallel << mShowPerspective;

  bundle << mPerspectiveRotY;
  bundle << atomrender.mLayerSeparation;
  bundle << atomrender.mSlicingPlaneNormal;
  bundle << atomrender.mSliceRotationPitch << atomrender.mSliceRotationRoll;
  bundle << atomrender.mSlicingPlanePoint << atomrender.mSlicingPlaneThickness
         << mLayerScaling;
  bundle << mShowGrid << mGridType << mGridSpacing << mGridXOffset
         << mGridYOffset;

  bundle << atomrender.mShowRadius;
  bundle << mCumulativeTrajectory << mIndividualTrajectory;
  bundle << mBillboarding;
  bundle << mSmallLabel;
  bundle << mDrawLabels;
}

void DataDisplay::initRootDirectory() {
  mDatasetManager.mLoadedDataset = "";

  mDatasetManager.initRoot();

  DatasetManager::SpeciesLabelMap labelMap =
      mDatasetManager.getAvailableSpecies();

  std::vector<std::string> availableAtoms;
  std::vector<std::string> atomLabels;
  for (auto species : labelMap) {
    availableAtoms.push_back(species.first);
    if (species.first != "Va") {
      // Vacancy atom names will be there already so wee need to skip to avoid
      // duplication
      atomLabels.insert(atomLabels.end(), species.second.begin(),
                        species.second.end());
    }
  }
  mShowAtoms.setNoCalls(0);
  // Only update menu and selection if it has changed
  if (atomLabels != mShowAtoms.getElements()) {
    // Fill the atoms that can be showed
    mShowAtoms.setElements(atomLabels);
  }
  // Only udpate available atoms menu and selections if it has changed
  if (availableAtoms != mAtomOfInterest.getElements()) {
    mAtomOfInterest.setElements(availableAtoms);

    // Match atom of interest to atom that can fill vacancies
    std::smatch match;
    std::regex atomNameRegex("[A-Z][a-z]*");
    for (auto atomName : labelMap["Va"]) {
      // Use only first result
      if (std::regex_search(atomName, match, atomNameRegex)) {
        string result = match.str();
        ptrdiff_t pos = std::distance(
            availableAtoms.begin(),
            find(availableAtoms.begin(), availableAtoms.end(), result));
        if (pos < (int)availableAtoms.size()) {
          mAtomOfInterest.set((int)pos);
        }
        mShowAtoms.setElementSelected(atomName);
      }
    }
  }
  auto dataNames = mDatasetManager.getDataNames();
  mPlotYAxis.setElements(dataNames);
  string defaultYAxis = "<comp_n(" + mAtomOfInterest.getCurrent() + ")>";
  ptrdiff_t pos = std::find(dataNames.begin(), dataNames.end(), defaultYAxis) -
                  dataNames.begin();
  if (pos < (int)dataNames.size()) {
    mPlotYAxis.set((int)pos);
  }

  std::vector<std::string> parameterSpaceNames;
  parameterSpaceNames.push_back("temperature");

  mPlotXAxis.setElements(parameterSpaceNames);
  mPlotYAxis.set(0);

  updateText();
}

void DataDisplay::updateText() {
  // Meta data texts
  metaText = "Global Root: " + mDatasetManager.mGlobalRoot + "\n";

  metaText += "Root: " + mDatasetManager.mRootPath.get() + "\n";
  metaText += "Dataset: " + mDatasetManager.mCurrentDataset.get() + "\n";
  metaText += "Subdir: " + mDatasetManager.getSubDir() + "\n";
  metaText +=
      "Condition Param: " + mDatasetManager.mConditionsParameter +
      " condition: " +
      std::to_string(mDatasetManager
                         .mParameterSpaces[mDatasetManager.mConditionsParameter]
                         ->getCurrentIndex()) +
      "\n";

  metaText += " ----- Parameters -----\n";
  for (auto param : mDatasetManager.mParameterSpaces) {
    if (param.second->size() > 2) {
      metaText += param.first + " : " + param.second->getCurrentId() + "\n";
    } else if (param.second->size() == 1) {
      // For spaces with a single value, the id will be ./ so show the value
      metaText += param.first + " : " +
                  std::to_string(param.second->getCurrentValue()) + "\n";
    }
  }
  metaText += " ----- Data -----\n";
  for (auto compData : mDatasetManager.getCurrentCompositions()) {
    metaText += compData.first + " = " + std::to_string(compData.second) + "\n";
  }
  metaText += "Current POSCAR :";
  metaText += mDatasetManager.labelProcessor.outputFile();

  // Parameter text
  auto subDir = mDatasetManager.getSubDir();
  auto temperatureId =
      mDatasetManager.mParameterSpaces["temperature"]->getCurrentId();
  std::string timeId;
  if (mDatasetManager.mParameterSpaces["time"]->size() > 0) {
    timeId = mDatasetManager.mParameterSpaces["time"]->getCurrentId();
  }
  if (subDir.size() > 0) {
    mParamText = subDir;
  }
  if (temperatureId.size() > 0) {
    mParamText += " temp:" + temperatureId + " ";
  }
  if (timeId.size() > 0) {
    mParamText += " time:" + timeId;
  }
  if (mParamText.size() == 0) {
    mParamText = "Dataset unavailable";
  }
}

void DataDisplay::prepare(Graphics &g, Matrix4f &transformMatrix) {
  if (mNeedsProcessing) {
    mDatasetManager.mCurrentDataset.processChange();
    mNeedsProcessing = false;
  }
  if (mDatasetManager.mCurrentDataset.hasChange()) {
    // TODO display some message here that we are computing
    mNeedsProcessing = true;
  }
  updateDisplayBuffers();

  // Load new graph data if needed.
  mDatasetManager.currentGraphName.processChange();

  // History mesh displays individual movements from their actual positions
  mHistoryMesh.primitive(Mesh::TRIANGLES);
  mHistoryMesh.reset();

  float historyWidth = 0.65f;
  size_t counter = mDatasetManager.mHistory.size() - 1;
  for (auto historyPoint = mDatasetManager.mHistory.begin();
       historyPoint != mDatasetManager.mHistory.end(); historyPoint++) {
    HSV hsvColor(0.5f * float(counter) / mDatasetManager.mHistory.size(), 1.0,
                 1.0);
    Color c;

    // Assumes the plane's normal is the z-axis
    Vec3f orthogonalVec = (historyPoint->second - historyPoint->first)
                              .cross({0, 0, 1})
                              .normalize(historyWidth);
    Vec3f orthogonalVec2 = (historyPoint->second - historyPoint->first)
                               .cross({1, 0, 0})
                               .normalize(historyWidth);
    if (orthogonalVec2.mag() < 0.0001f) {
      orthogonalVec2 = (historyPoint->second - historyPoint->first)
                           .cross({0, 1, 0})
                           .normalize(historyWidth);
    }
    assert(orthogonalVec2.mag() > 0.0001f);
    ImGui::ColorConvertHSVtoRGB(hsvColor.h, hsvColor.s, hsvColor.v, c.r, c.g,
                                c.b);
    c.a = 0.35f;
    unsigned int previousSize = mHistoryMesh.vertices().size();
    mHistoryMesh.color(c);
    mHistoryMesh.vertex(historyPoint->first - orthogonalVec);
    mHistoryMesh.color(c);
    mHistoryMesh.vertex(historyPoint->second - orthogonalVec * 0.1f);
    mHistoryMesh.color(c);
    mHistoryMesh.vertex(historyPoint->first + orthogonalVec);
    mHistoryMesh.color(c);
    mHistoryMesh.vertex(historyPoint->second + orthogonalVec * 0.1f);

    mHistoryMesh.index(previousSize);
    mHistoryMesh.index(previousSize + 1);
    mHistoryMesh.index(previousSize + 2);

    mHistoryMesh.index(previousSize + 2);
    mHistoryMesh.index(previousSize + 3);
    mHistoryMesh.index(previousSize + 1);
    counter--;
  }

  // Trajectory mesh displays the cumulative trajectory
  mTrajectoryMesh.primitive(Mesh::TRIANGLES);
  mTrajectoryMesh.reset();
  Vec3f previousPoint(0, 0, 0);
  float previousMag = 0.0;

  float trajectoryWidth = 0.35f;
  counter = mDatasetManager.mHistory.size() - 1;
  for (auto historyPoint = mDatasetManager.mHistory.begin();
       historyPoint != mDatasetManager.mHistory.end(); historyPoint++) {
    HSV hsvColor(0.5f * float(counter) / mDatasetManager.mHistory.size(), 1.0,
                 1.0);
    Color c;

    // Assumes the plane's normal is the z-axis
    Vec3f thisMovement = historyPoint->second - historyPoint->first;
    Vec3f orthogonalVec =
        thisMovement.cross({0, 0, 1}).normalize(trajectoryWidth);
    ImGui::ColorConvertHSVtoRGB(hsvColor.h, hsvColor.s, hsvColor.v, c.r, c.g,
                                c.b);
    c.a = 0.8f;
    unsigned int previousSize = mTrajectoryMesh.vertices().size();
    if (thisMovement.mag() >
        fabs(mDataBoundaries.max.x - mDataBoundaries.min.x) /
            2.0f) { // Atom is wrapping around
      //              c = Color(0.8f, 0.8f, 0.8f, 1.0f);
      thisMovement = -thisMovement;
      thisMovement.normalize(previousMag);
    } else {
      previousMag = thisMovement.mag();
    }
    mTrajectoryMesh.color(c);
    mTrajectoryMesh.vertex(previousPoint - orthogonalVec);
    mTrajectoryMesh.color(c);
    mTrajectoryMesh.vertex(previousPoint + thisMovement - orthogonalVec * 0.2f);
    mTrajectoryMesh.color(c);
    mTrajectoryMesh.vertex(previousPoint + orthogonalVec);
    mTrajectoryMesh.color(c);
    mTrajectoryMesh.vertex(previousPoint + thisMovement + orthogonalVec * 0.2f);

    mTrajectoryMesh.index(previousSize);
    mTrajectoryMesh.index(previousSize + 1);
    mTrajectoryMesh.index(previousSize + 2);

    mTrajectoryMesh.index(previousSize + 2);
    mTrajectoryMesh.index(previousSize + 3);
    mTrajectoryMesh.index(previousSize + 1);
    previousPoint = previousPoint + thisMovement;
    counter--;
  }

  prepareParallelProjection(g, transformMatrix);
  //        g.scale(1.0/(mDataBoundaries.maxy - mDataBoundaries.miny));
}

void DataDisplay::draw(Graphics &g) { // Load data after drawing frame to allow
                                      // showing "in process" state
  if (mVisible == 0.0f) {
    return;
  }

  std::unique_lock<std::mutex> lk(mDrawLock);
  // g.depthTesting(false);
  gl::blendAdd();

  g.lens().eyeSep(0.2);

  if (mShowPerspective.get() == 1.0f) {
    drawPerspective(g);
  }

  // now draw graph and iso view 3d slabs
  // set view to identity, effectively putting coord in eye space
  // g.pushViewMatrix(Matrix4f::identity());
  g.texture();
  gl::blendTrans();

  if (mShowParallel.get() == 1.0f) {
    drawParallelProjection(g);
  }

  gl::blending(false);
  if (mShowGraph.get() == 1.0f) {
    drawGraph(g);
  }
  // g.popViewMatrix();
}

void DataDisplay::setFont(string name, float size) {
  std::unique_lock<std::mutex> lk(mDrawLock);
  if (File::exists(name)) {
    if (!mLabelFont.load(name.c_str(), size, 1024)) {
      std::cout << "Failed to load font: " << name << std::endl;
      if (!mLabelFont.load("C:\\Windows\\Fonts\\arial.ttf", size, 1024)) {
      }

      if (!mLabelFont.load(
              "/usr/share/fonts/truetype/ttf-bitstream-vera/Vera.ttf", size,
              1024)) {
        throw;
      }
    } else {
      std::cout << "Loaded font " << name << std::endl;
    }
  } else {
    std::cout << "Could not sync font: " << name << std::endl;
  }
}

void DataDisplay::dumpImages(string dumpPrefix) {
  std::unique_lock<std::mutex> lk(mDrawLock);
  std::string dumpDirectory = File::conformPathToOS(
      mDatasetManager.buildRootPath() + mDatasetManager.mCurrentDataset.get() +
      "/graphics");
  if (!File::exists(dumpDirectory)) {
    if (!Dir::make(dumpDirectory)) {
      std::cerr << "Failed to create directory: " << dumpDirectory << std::endl;
      return;
    }
  }

  std::string fullDatasetPath = File::conformPathToOS(
      mDatasetManager.buildRootPath() +
      File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));
  std::string graphFilename =
      File::conformPathToOS(dumpDirectory + "/" + dumpPrefix + "_graph.png");
  std::cout << "dumping "
            << fullDatasetPath + mDatasetManager.currentGraphName.get()
            << " -> " << graphFilename << std::endl;
  if (!File::copy(File::conformPathToOS(fullDatasetPath +
                                        mDatasetManager.currentGraphName.get()),
                  graphFilename)) {
    std::cerr << "ERROR copying png file." << std::endl;
  }

  auto allPositions = mDatasetManager.positionBuffers.get(false);
  File allPositionsFile(File::conformPathToOS(dumpDirectory + "/" + dumpPrefix +
                                              "_positions.csv"),
                        "w", true);
  std::string header = "element,x,y,z\n";
  allPositionsFile.write(header);
  for (auto elemPositions : *allPositions) {
    auto positionIt = elemPositions.second.begin();
    while (positionIt != elemPositions.second.end()) {
      std::string line = elemPositions.first + ",";
      line += std::to_string(*positionIt++) + ",";
      line += std::to_string(*positionIt++) + ",";
      line += std::to_string(*positionIt++) + "\n";
      positionIt++; // Ignore label?
      allPositionsFile.write(line);
    }
  }
  allPositionsFile.close();

  auto between_planes = [&](Vec3f &point, Vec3f &plane_point,
                            Vec3f &plane_normal, float &second_plane_distance) {
    Vec3f difference = point - plane_point;
    float proj = plane_normal.dot(difference);
    if (proj >= 0 && proj <= second_plane_distance) {
      return true;
    } else {
      return false;
    }
  };

  auto plane_point = atomrender.mSlicingPlanePoint.get();
  auto plane_normal = atomrender.mSlicingPlaneNormal.get().normalized();
  auto second_plane_distance = atomrender.mSlicingPlaneThickness.get();
  auto slicePositions = *mDatasetManager.positionBuffers.get(false);
  File slicePositionsFile(File::conformPathToOS(dumpDirectory + "/" +
                                                dumpPrefix +
                                                "_slice_positions.csv"),
                          "w", true);
  std::string slicePositionsHeader = "element,x,y,z\n";
  slicePositionsFile.write(slicePositionsHeader);

  auto selectedElements = mShowAtoms.getSelectedElements();
  for (auto atom : selectedElements) {
    for (size_t i = 0; i < slicePositions[atom].size(); i += 4) {
      auto x = slicePositions[atom][i];
      auto y = slicePositions[atom][i + 1];
      auto z = slicePositions[atom][i + 2];
      Vec3f pos(x, y, z);
      if (between_planes(pos, plane_point, plane_normal,
                         second_plane_distance)) {
        std::string line = atom + ",";
        line += std::to_string(x) + ",";
        line += std::to_string(y) + ",";
        line += std::to_string(z) + "\n";
        slicePositionsFile.write(line);
      }
    }
  }
  slicePositionsFile.close();

  json metadata;
  metadata["dataset"]["path"] = mDatasetManager.currentDataset();
  metadata["dataset"]["subdir"] = mDatasetManager.getSubDir();
  metadata["dataset"]["rootpath"] = mDatasetManager.buildRootPath();
  std::string condition = std::to_string(
      mDatasetManager.mParameterSpaces[mDatasetManager.mConditionsParameter]
          ->getCurrentIndex());
  metadata["dataset"]["condition"] = condition;
  for (auto space : mDatasetManager.mParameterSpaces) {
    if (space.second && space.second->size() > 0) {
      metadata["parameters"][space.first] = space.second->getAllCurrentIds()[0];
    }
  }

  metadata["SlicingPlanePoint"] = {plane_point.x, plane_point.y, plane_point.z};
  metadata["SliceNormal"] = {plane_normal.x, plane_normal.y, plane_normal.z};

  metadata["SliceThickness"] = atomrender.mSlicingPlaneThickness.get();
  std::ofstream metadatafile(dumpDirectory + "/" + dumpPrefix +
                             "_metadata.json");
  metadatafile << std::setw(4) << metadata;
}

void DataDisplay::updateDisplayBuffers() {
  if (mDatasetManager.positionBuffers.newDataAvailable()) {
    std::map<string, int> elementCounts;

    auto allPositions = mDatasetManager.positionBuffers.get();

    mDatasetManager.mCurrentLoadedIndeces.clear();
    for (auto &paramSpace : mDatasetManager.mParameterSpaces) {
      mDatasetManager.mCurrentLoadedIndeces[paramSpace.first] =
          paramSpace.second->getCurrentIndex();
    }

    vector<vector<float> *> elemPositions;

    // Now that the layer direction has been computed from the atom of interest,
    // we need to generate aligned data from the visible atoms.
    elemPositions.clear();
    auto visibleAtoms = mShowAtoms.getSelectedElements();
    for (auto &elementData : *allPositions) {
      if (std::find(visibleAtoms.begin(), visibleAtoms.end(),
                    elementData.first) != visibleAtoms.end()) {
        elemPositions.push_back(&(elementData.second));
      }
    }
    // also prepare layer normal direction aligned data
    size_t totalSize = 0;
    for (auto *elems : elemPositions) {
      totalSize += elems->size();
    }

    mDataBoundaries.resetInv();
    mAligned4fData.resize(totalSize * 4);
    auto outit = mAligned4fData.begin();
    // now fill the
    for (auto *elems : elemPositions) {
      auto it = elems->begin();
      while (it != elems->end()) {
        assert(outit != mAligned4fData.end());
        float &x = *it++;
        float &y = *it++;
        float &z = *it++;
        float &w = *it++;
        Vec3f vec(x, y, z);
        *outit++ = x;
        *outit++ = y;
        *outit++ = z;
        *outit++ = w;
        mDataBoundaries.includePoint(vec);
      }
    }

    if (mAligned4fData.size() > 0) {
      auto &b = mDataBoundaries;
      perspectivePickable.bb.set(Vec3f(b.min.x, b.min.y, b.min.z),
                                 Vec3f(b.max.x, b.max.y, b.max.z));
      slicePickable.bb.set(Vec3f(b.min.x, b.min.y, b.min.z),
                           Vec3f(b.max.x, b.max.y, (b.max.z - b.min.z) * 0.5f));
      // rh.pose.pos().set(perspectivePickable.bb.cen);
      atomrender.setDataBoundaries(b);
    }

    // Set active atoms and colors
    atomPropertiesProj.clear();
    vector<AtomProperties> atomPropertiesPersp;

    // TODO these colors should be exposed as a preference
    vector<Color> colorList = {
        Color(0.0, 1.0, 1.0, 1.0), Color(1.0, 1.0, 0.0, 1.0),
        Color(1.0, 0.0, 1.0, 1.0), Color(0.0, 1.0, 1.0, 1.0),
        Color(1.0, 1.0, 0.0, 1.0), Color(1.0, 0.0, 1.0, 1.0)};

    vector<Color> colorList2 = {
        Color(0.7f, 0.0f, 0.7f, 0.4f), Color(0.7f, 0.0f, 0.7f, 1.0f),
        Color(0.0f, 0.7f, 0.0f, 0.4f), Color(0.0f, 0.7f, 0.0f, 1.0f),
        Color(0.0f, 0.0f, 0.7f, 0.4f), Color(0.0f, 0.0f, 0.7f, 1.0f),
        Color(0.7f, 0.0f, 0.7f, 0.4f), Color(0.7f, 0.0f, 0.7f, 1.0f),
        Color(0.0f, 0.7f, 0.0f, 0.4f), Color(0.0f, 0.7f, 0.0f, 1.0f),
        Color(0.0f, 0.0f, 0.7f, 0.4f), Color(0.0f, 0.0f, 0.7f, 1.0f),
    };
    auto colorListIt = colorList.begin();
    auto colorList2It = colorList2.begin();

    auto selectedElements = mShowAtoms.getSelectedElements();
    for (auto atom : mShowAtoms.getElements()) {
      if (std::find(selectedElements.begin(), selectedElements.end(), atom) !=
          selectedElements.end()) {
        if (elementData.find(atom) != elementData.end()) {
          //                    std::cout << "Color for: " << atom << ":" <<
          //                    elementData[atom].color.r << " " <<
          //                    elementData[atom].color.g << " "<<
          //                    elementData[atom].color.b << " "  <<std::endl ;
          // Atom was matched in elements.ini file, so use those colors

          atomPropertiesProj.push_back(AtomProperties{
              atom, elementData[atom].radius, elementData[atom].color});

          atomPropertiesPersp.push_back(AtomProperties{
              atom, elementData[atom].radius, elementData[atom].color});
        } else { // Use defaults
          atomPropertiesProj.push_back(
              AtomProperties{atom, 1.0f, *colorListIt});

          atomPropertiesPersp.push_back(
              AtomProperties{atom, 1.0f, *colorList2It});
        }
      }
      colorListIt++;
      colorList2It++;
      colorList2It++;
    }

    // Apply colors to aligned data -----------
    list<Color> colors;
    mAtomData.clear();
    for (auto elem : *allPositions) {
      elementCounts[elem.first] = elem.second.size() / 4;
    }
    for (auto atomProps : atomPropertiesProj) {
      colors.push_back(atomProps.color);
      mAtomData.push_back(
          {elementCounts[atomProps.name], atomProps.drawScale, atomProps.name});
    }
    float hue = 0.0f;

    if (colors.size() > 0) {
      hue = rgb2hsv(colors.front().rgb()).h;
    }
    if (mAtomData.size() > 0) {
      auto atomDataIt = mAtomData.begin();
      int atomCounter = 0;
      for (size_t i = 0; i < mAligned4fData.size() / 4; i++) {
        //            assert(atomCountsIt != mAtomCounts.end());
        if (atomDataIt != mAtomData.end() && --atomCounter <= 0) {
          atomCounter = atomDataIt->counts;
          atomDataIt++;
          hue = rgb2hsv(colors.front().rgb()).h;
          colors.pop_front();
        }
        mAligned4fData[4 * i + 3] = hue;
      }
    }

    // Load image data ----------------------------

    if (mRunComputation) {
      std::string xData = mPlotXAxis.getCurrent();
      std::string yData = mPlotYAxis.getCurrent();
      mDatasetManager.generateGraph(
          xData, yData, mDatasetManager.mCurrentDataset.get(), false);
    }
  }

  updateParameterText();
}

void DataDisplay::prepareParallelProjection(Graphics &g,
                                            Matrix4f &transformMatrix) {
  //        std::unique_lock<std::mutex> lk(mDataLock);
  g.pushFramebuffer(fbo_iso);
  g.pushViewport(fbo_iso.width(), fbo_iso.height());
  g.pushModelMatrix();
  g.pushViewMatrix();
  g.pushProjMatrix();

  // note that translation in z direction will mess up layer clipping
  g.modelMatrix(transformMatrix);
  // projection matrix from clipping params
  // will only work with z axis aligned data and single viewport mode

  //        const float cameraZ = mDataBoundaries.maxz;

  float ar = float(fbo_iso.width()) / fbo_iso.height();
  //        float rangeSizeZ = mDataBoundaries.maxz - mDataBoundaries.minz;
  float centerX = (mDataBoundaries.max.x + mDataBoundaries.min.x) / 2.0f;
  float centerY = (mDataBoundaries.max.y + mDataBoundaries.min.y) / 2.0f;
  float rangeSizeX = mDataBoundaries.max.x - mDataBoundaries.min.x;
  float rangeSizeY = mDataBoundaries.max.y - mDataBoundaries.min.y;
  float maxrange = std::max(rangeSizeX, rangeSizeY);
  float padding = maxrange * 0.05f;
  float left = centerX - 0.5f * maxrange * mLayerScaling - padding;
  float right = centerX + 0.5f * maxrange * mLayerScaling + padding;
  float top = centerY + 0.5f * maxrange * mLayerScaling + padding;
  float bottom = centerY - 0.5f * maxrange * mLayerScaling - padding;
  //        const float near = mDataBoundaries.minz + (mNearClip * rangeSizeZ);
  //        const float farClip = mDataBoundaries.minz + (mFarClip *
  //        rangeSizeZ); float minxy = std::min(mDataBoundaries.minx,
  //        mDataBoundaries.miny);
  g.projMatrix(Matrix4f::ortho(ar * left, ar * right, bottom, top,
                               mDataBoundaries.min.z - 100,
                               mDataBoundaries.max.z + 100));

  bool mAlignData = true;
  double scalingFactor = 1.0 / (mDataBoundaries.max.y - mDataBoundaries.min.y);

  // std::cout << near << "..." << farClip <<std::endl;

  // save for later
  GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);
  gl::scissorTest(true);
  int w = fbo_iso.width() / 2;
  int h = fbo_iso.height() / 2;
  gl::polygonFill();
  gl::depthMask(true); // for axis rendering

  if (mSingleProjection.get() == 1.0f) {
    g.pushViewport(0, 0, 2 * w, 2 * h);
    gl::scissorArea(0, 0, 2 * w, 2 * h);
    g.clear(backgroundColor);

    g.pushViewMatrix();
    //            g.viewMatrix(getLookAt({ 0.0f, 0.0f, cameraZ },
    //            { 0.0f, 0.0f, 0.0f },
    //            { 0.0f, 1.0f, 0.0f }));
    g.viewMatrix(
        getLookAt(atomrender.mSlicingPlanePoint,
                  atomrender.mSlicingPlanePoint.get() -
                      atomrender.mSlicingPlaneNormal.get().normalized(),
                  {0.0f, 1.0f, 0.0f}));

    gl::blending(false);
    gl::depthTesting(true);
    g.meshColor();
    g.draw(axis);

    if (mShowGrid == 1.0f) {
      mGridMesh.reset();
      addRect(mGridMesh, 0, -2.0, 0.1f / (maxrange * mLayerScaling), 4.0f);
      mGridMesh.update();
      if (backgroundColor.get().luminance() > 0.5f) {
        g.color(0.2f);
      } else {
        g.color(0.7f);
      }
      for (int i = 0; i < 40; i++) {
        if (mGridType.getCurrent() == "square") {
          g.pushMatrix();
          g.translate(mGridXOffset, mGridYOffset + mGridSpacing * i, 0);
          g.scale(maxrange * mLayerScaling);
          g.rotate(90, 0, 0, 1);
          g.draw(mGridMesh);
          g.popMatrix();
          g.pushMatrix();
          g.translate(mGridXOffset + mGridSpacing * i, mGridYOffset, 0);
          g.scale(maxrange * mLayerScaling);
          g.draw(mGridMesh);
          g.popMatrix();
        } else if (mGridType.getCurrent() == "triangle") {
          float spacing = 0.86602540378f * mGridSpacing; // sin(60)
          g.pushMatrix();
          g.translate(mGridXOffset, mGridYOffset + spacing * i, 0);
          g.scale(maxrange * mLayerScaling);
          g.rotate(90, 0, 0, 1);
          g.draw(mGridMesh);
          g.popMatrix();
          g.pushMatrix();
          g.translate(mGridXOffset, mGridYOffset, 0);
          g.rotate(30, 0, 0, 1);
          g.translate(spacing * i, 0, 0);
          g.scale(maxrange * mLayerScaling);
          g.draw(mGridMesh);
          g.popMatrix();
          g.pushMatrix();
          g.translate(mGridXOffset, mGridYOffset, 0);
          g.rotate(-30, 0, 0, 1);
          g.translate(spacing * i, 0, 0);
          g.scale(maxrange * mLayerScaling);
          g.draw(mGridMesh);
          g.popMatrix();
        }
      }
    }

    gl::blending(true);
    gl::blendTrans();
    gl::depthTesting(false);
    g.pushMatrix();

    if (atomPropertiesProj.size() > 0) {
      // ----------------------------------------
      int cumulativeCount = 0;
      for (auto &data : mAtomData) {
        if (mAlignData) {
          //            instancing_mesh0.attrib_data(
          //              mAligned4fData.size() * sizeof(float),
          //              mAligned4fData.data(),
          //              mAligned4fData.size()/4
          //            );
          int count = data.counts;
          assert((int)mAligned4fData.size() >= (cumulativeCount + count) * 4);
          atomrender.instancing_mesh0.attrib_data(
              count * 4 * sizeof(float),
              mAligned4fData.data() + (cumulativeCount * 4), count);
          cumulativeCount += count;
          //                        std::cout << "Drawing " << counts << " of "
          //                        << cumulativeCount << std::endl;
        } else {
          //            auto& elemPositions = r0->getElementPositions(p0.name);
          //            instancing_mesh0.attrib_data(
          //              elemPositions.size() * sizeof(float),
          //              elemPositions.data(),
          //              elemPositions.size()/4
          //            );
        }
        // now draw data with custom shader
        g.shader(atomrender.instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        g.shader().uniform("layerSeparation", 1.0);
        // A scaling value of 4.0 found empirically...
        if (atomrender.mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", data.radius *
                                                atomrender.mAtomMarkerSize *
                                                atomrender.mMarkerScale);
        } else {
          g.shader().uniform("markerScale", atomrender.mAtomMarkerSize *
                                                atomrender.mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);

        g.shader().uniform("plane_point", atomrender.mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           atomrender.mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance",
                           atomrender.mSlicingPlaneThickness);

        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);

        atomrender.instancing_mesh0.draw();
      }
    }

    drawHistory(g);

    g.popMatrix();
    g.popViewMatrix();
    g.popViewport();
  } else {
    // (from x, y up), (right bottom)
    {
      g.pushViewport(w, 0, w, h);
      gl::scissorArea(w, 0, w, h);
      g.clear(backgroundColor);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({4, 0, 0}, {0, 0, 0}, {0, 1, 0}));

      gl::blending(false);
      gl::depthTesting(true);
      g.meshColor();
      g.draw(axis);

      gl::blending(true);
      gl::blendAdd();
      gl::depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(atomrender.instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        g.shader().uniform("layerSeparation", 1.0);

        if (atomrender.mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", atomrender.mAtomMarkerSize *
                                                atomrender.mMarkerScale);
        } else {
          g.shader().uniform("markerScale", atomrender.mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);

        g.shader().uniform("plane_point", atomrender.mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           atomrender.mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance",
                           atomrender.mSlicingPlaneThickness);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        atomrender.instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }

    // (from y, -z up), (left top)
    {
      g.pushViewport(0, h, w, h);
      gl::scissorArea(0, h, w, h);
      g.clear(0, 0.2f, 0);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({0, 4, 0}, {0, 0, 0}, {0, 0, -1}));

      gl::blending(false);
      gl::depthTesting(true);
      g.meshColor();
      g.draw(axis);

      gl::blendAdd();
      gl::depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(atomrender.instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        if (atomrender.mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", atomrender.mAtomMarkerSize *
                                                atomrender.mMarkerScale);
        } else {
          g.shader().uniform("markerScale", atomrender.mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);

        g.shader().uniform("plane_point", atomrender.mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           atomrender.mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance",
                           atomrender.mSlicingPlaneThickness);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        atomrender.instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }

    // (from z, y up), (left bottom)
    {
      g.pushViewport(0, 0, w, h);
      gl::scissorArea(0, 0, w, h);
      g.clear(0, 0, 0.2f);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({0, 0, 4}, {0, 0, 0}, {0, 1, 0}));

      gl::blending(false);
      gl::depthTesting(true);
      g.meshColor();
      g.draw(axis);

      gl::blending(true);
      gl::blendAdd();
      gl::depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(atomrender.instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        if (atomrender.mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", atomrender.mAtomMarkerSize *
                                                atomrender.mMarkerScale);
        } else {
          g.shader().uniform("markerScale", atomrender.mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);
        g.shader().uniform("plane_point", atomrender.mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           atomrender.mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance",
                           atomrender.mSlicingPlaneThickness);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        atomrender.instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }

    // diagonal view, right-top
    {
      g.pushViewport(w, h, w, h);
      gl::scissorArea(w, h, w, h);
      g.clear(0.2f, 0.2f, 0.2f);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({3, 3, 3}, {0, 0, 0}, {0, 1, 0}));

      gl::blending(false);
      gl::depthTesting(true);
      g.meshColor();
      g.draw(axis);

      gl::blending(true);
      gl::blendAdd();
      gl::depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(atomrender.instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        if (atomrender.mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", atomrender.mAtomMarkerSize *
                                                atomrender.mMarkerScale);
        } else {
          g.shader().uniform("markerScale", atomrender.mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);
        g.shader().uniform("plane_point", atomrender.mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           atomrender.mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance",
                           atomrender.mSlicingPlaneThickness);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        atomrender.instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }
  }

  // put back scissoring
  gl::scissorTest(last_enable_scissor_test);
  g.popProjMatrix();
  g.popViewMatrix();
  g.popModelMatrix();
  g.popFramebuffer();
  g.popViewport();
}

void DataDisplay::drawPerspective(Graphics &g) {
  //  g.camera(Viewpoint::ORTHO_FOR_2D);

  //  float rpd =
  //      getCurrentWindowPixelDensity();  // reciprocal of pixel density
  //  auto v = viewport();                 // viewport in framebuffer unit
  //  mViewStack.setIdentity();
  //  mProjStack.set(Matrix4f::ortho(v.l * rpd, v.w * rpd,  // left, right
  //                                 v.b * rpd, v.h * rpd,  // bottom, top
  //                                 0, 1                   // near, far
  //                                 ));

  //        Vec3d normVector = reader.getNormalizingVector();
  //        Vec3d centeringVector = reader.getCenteringVector();
  if (mAligned4fData.size() == 0) {
    return; // No data has been loaded
  }

  g.lens().far(2000);
  perspectivePickable.pushMatrix(g);

  // move to center
  // g.translate(-(mDataBoundaries.maxx -
  // mDataBoundaries.minx)/2,-(mDataBoundaries.maxy -
  // mDataBoundaries.miny)/2,-(mDataBoundaries.maxz - mDataBoundaries.minz)/2);
  // g.translate(0,0,-7 -(mDataBoundaries.minz+mDataBoundaries.maxz)/2);
  // g.scale(1/(mDataBoundaries.maxy - mDataBoundaries.miny));
  // g.rotate(mPerspectiveRotY.get(), 0, 1.0, 0);

  gl::blendAdd();
  gl::depthTesting(false);

  g.meshColor();

  g.pushMatrix();
  g.scale(20);
  g.draw(axis);
  g.popMatrix();
  // g.scale(mPerspectiveScale);

  //        float near = mDataBoundaries.minz + (mNearClip *
  //        (mDataBoundaries.maxz- mDataBoundaries.minz)); float farClip =
  //        mDataBoundaries.minz + (mFarClip * (mDataBoundaries.maxz-
  //        mDataBoundaries.minz));

  atomrender.draw(g, perspectivePickable.scale, mAtomData, mAligned4fData);

  gl::depthTesting(true);
  gl::blending(true);

  g.pushMatrix();
  // Update box mesh that marks near and far clipping
  // assumption
  //   1. data is aligned to z axis
  //   2. data range is [0:1] for all x, y, z
  //   3. modelview matrix for perspective view should have been applied
  boxMesh.reset();
  boxMesh.primitive(Mesh::LINES);
  // const float z[2] = {mNearClip.get(), mFarClip.get()};
  // double scalingFactor = 1.0/(mDataBoundaries.maxz - mDataBoundaries.minz);

  // Plane equation Ax + By + Cz - D = 0
  //        float D =  -mSlicingPlanePoint.get().dot(mSlicingPlaneNormal);
  // Line Equation

  const float x[2] = {mDataBoundaries.max.x, mDataBoundaries.min.x};
  const float y[2] = {mDataBoundaries.max.y, mDataBoundaries.min.y};

  const float z[2] = {0.0, (atomrender.mSlicingPlaneThickness)};

  for (int k = 0; k <= 1; k++) {
    for (int j = 0; j <= 1; j++) {
      for (int i = 0; i <= 1; i++) {
        boxMesh.vertex(x[i], y[j], z[k]);
      }
    }
  }

  static const int I[] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 1, 3,
                          4, 6, 5, 7, 0, 4, 1, 5, 2, 6, 3, 7};
  boxMesh.index(I, sizeof(I) / sizeof(*I), 0);
  boxMesh.update();
  g.shader().uniform("eye_sep", perspectivePickable.scale * g.lens().eyeSep() *
                                    g.eye() / 2.0f);

  g.shader().uniform("eye_sep",
                     /*perspectivePickable.scale * */ g.lens().eyeSep() *
                         g.eye() / 2.0f);

  glLineWidth(5);
  gl::polygonLine();
  g.color(0.8f, 0.8f, 1.0f, 0.9f);
  g.scale(1.0, 1.0, 1.0f + atomrender.mLayerSeparation);
  g.translate(atomrender.mSlicingPlanePoint.get());
  g.rotate(atomrender.mSliceRotationPitch * 360.0f / (M_2PI), 1.0, 0.0, 0.0);
  g.rotate(atomrender.mSliceRotationRoll * 360.0f / (M_2PI), 0.0, -1.0, 0.0);
  g.draw(boxMesh);
  gl::polygonFill();
  glLineWidth(1);

  //  -----------

  g.popMatrix();

  drawHistory(g);

  g.popMatrix(); // pickable

  gl::depthTesting(false);
  gl::blendAdd();
  g.pushMatrix();
  // TODO billboard this text
  perspectivePickable.updateAABB();

  if (mDrawLabels == 1.0f) {
    Vec3f textPos = perspectivePickable.aabb.min;
    //        textPos.x = perspectivePickable.aabb.cen.x;
    textPos.z = perspectivePickable.aabb.cen.z;
    float dist = textPos.mag() * 0.025f;
    g.translate(textPos);
    if (mSmallLabel == 1.0f) {
      g.scale(dist * 0.5f);
    } else {
      //            g.scale(dist * 0.5);
    }
    g.translate(0.0, -1.5, 0); // Push the text down a bit

    // Draw label
    Mesh mesh;
    mLabelFont.write(
        mesh, (mDatasetManager.currentDataset() + " " + mParamText).c_str(),
        1.0);
    g.texture();
    mLabelFont.tex.bind();
    g.draw(mesh);
    mLabelFont.tex.unbind();
  }

  g.popMatrix();

  perspectivePickable.drawBB(g);
  //        perspectivePickable.drawChildren(g);
  // if(perspectivePickable.hover || rh.hover[0] || rh.hover[1] || rh.hover[2]
  // ){
  //     g.depthTesting(false);
  //     perspectivePickable.drawChildren(g);
  //     g.depthTesting(true);
  // }
}

void DataDisplay::drawParallelProjection(Graphics &g) {
  float pixel_w = float(fbo_iso.width());
  float pixel_h = float(fbo_iso.height());
  float w = 2;
  float h = w * pixel_h / pixel_w;
  Mesh m;
  addTexQuad(m, w, h);
  iso_scene().filter(Texture::LINEAR_MIPMAP_LINEAR);
  iso_scene().generateMipmap(); // XXX this works.. confusing api needs help..
  iso_scene().bind();
  parallelPickable.pushMatrix(g);
  // g.pushMatrix();
  // g.translate(mParallelPose.get().pos() );
  Vec3f forward = parallelPickable.pose.get().pos();
  forward.normalize();
  if (mBillboarding.get() == 1.0f) {
    Quatf rot = Quatf::getBillboardRotation(-forward, Vec3f{0.0, 1.0f, 0.0f});
    g.rotate(rot);
  }
  gl::depthTesting(true);
  g.draw(m);
  // Draw label
  if (mDrawLabels == 1.0f) {
    gl::depthTesting(false);
    gl::blendAdd();
    g.translate(-1.0, -2.0, 0.0);
    if (mSmallLabel == 1.0f) {
      g.scale(0.15f);
    } else {
      g.scale(0.3f);
    }

    // Draw label
    Mesh mesh;
    mLabelFont.write(
        mesh, (mDatasetManager.currentDataset() + " " + mParamText).c_str(),
        1.0);
    g.texture();
    mLabelFont.tex.bind();
    g.draw(mesh);
    mLabelFont.tex.unbind();
  }
  g.popMatrix();
  parallelPickable.drawBB(g);
}

void DataDisplay::drawGraph(Graphics &g) {
  graphPickable.pushMatrix(g);

  Vec3f forward = graphPickable.pose.get().pos();
  forward.normalize();

  if (mBillboarding.get() == 1.0f) {
    Quatf rot = Quatf::getBillboardRotation(-forward, Vec3f{0.0, 1.0f, 0.0f});
    g.rotate(rot);
  }

  gl::depthTesting(true);
  g.quad(mGraphTexture, -1.2f, 0.9f, 2.4f, -1.8f);
  // Draw label

  if (mDrawLabels == 1.0f) {
    gl::depthTesting(false);
    gl::blendAdd();
    g.translate(-1.0, -1.0, 0.0);
    if (mSmallLabel == 1.0f) {
      g.scale(0.1f);
    } else {
      g.scale(0.2f);
    }

    // Draw label
    Mesh mesh;
    mLabelFont.write(
        mesh, (mDatasetManager.currentDataset() + " " + mParamText).c_str(),
        1.0);
    g.texture();
    mLabelFont.tex.bind();
    g.draw(mesh);
    mLabelFont.tex.unbind();
  }
  graphPickable.popMatrix(g);

  graphPickable.drawBB(g);
}
