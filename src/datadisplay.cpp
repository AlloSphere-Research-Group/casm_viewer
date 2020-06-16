
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
      // New values must be available for computeNewSample, throws away previous
      // value
      mShowAtoms.setNoCalls(value);
      // TODO we should have a less heavy function to turn atoms on and off
      mDatasetManager.computeNewSample();
    }
  });

  atomrender.mSlicingPlaneThickness.registerChangeCallback([this](float v) {
    auto m = slicePickable.bb.max;
    m.z = v;
    slicePickable.bb.set(slicePickable.bb.min, m);
  });

  atomrender.init();

  slicePickable.pose.registerChangeCallback([this](Pose pose) {
    atomrender.mSlicingPlanePoint.setNoCalls(pose.pos());
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
  bundle << mDatasetManager.mPlotXAxis;
  bundle << mDatasetManager.mPlotYAxis;
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
      // Vacancy atom names will be there already so we need to skip to avoid
      // duplication
      atomLabels.insert(atomLabels.end(), species.second.begin(),
                        species.second.end());
    }
  }
  // Only update menu and selection if it has changed
  if (atomLabels != mShowAtoms.getElements()) {
    // Fill the atoms that can be shown
    mShowAtoms.setElements(atomLabels);
    mShowAtoms.setNoCalls(0);
  }
  auto speciesOfInterest = mDatasetManager.mAtomOfInterest.getCurrent();
  // Show atoms that begin with atom of interest e.g. Na1 if atom of interest Na
  for (auto elementLabel : mShowAtoms.getElements()) {
    if (elementLabel.substr(0, speciesOfInterest.size()) == speciesOfInterest) {
      mShowAtoms.setElementSelected(elementLabel);
    }
  }
}

void DataDisplay::prepareHistoryMesh() {
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
}

void DataDisplay::prepare(Graphics &g, Matrix4f &transformMatrix) {

  if (mNeedsProcessing) {
    mDatasetManager.mCurrentDataset.processChange();
    mDatasetManager.currentGraphName.processChange();
    mNeedsProcessing = false;
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
  if (mDatasetManager.mCurrentDataset.hasChange() ||
      mDatasetManager.currentGraphName.hasChange()) {
    // Replace parameter text with temporary text:
    mParamText = "Processing ...";
    mNeedsProcessing = true;
    // Schedule processing of changes for next pass to allow drawing one frame.
  }

  if (mDatasetManager.occupationData.newDataAvailable()) {
    updateDisplayBuffers();
  }
  prepareHistoryMesh();
  prepareParallelProjection(g, transformMatrix);
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

// void DataDisplay::dumpImages(string dumpPrefix) {
//  std::unique_lock<std::mutex> lk(mDrawLock);
//  std::string dumpDirectory = File::conformPathToOS(
//      mDatasetManager.buildRootPath() + mDatasetManager.mCurrentDataset.get()
//      +
//      "/graphics");
//  if (!File::exists(dumpDirectory)) {
//    if (!Dir::make(dumpDirectory)) {
//      std::cerr << "Failed to create directory: " << dumpDirectory <<
//      std::endl; return;
//    }
//  }

//  std::string fullDatasetPath = File::conformPathToOS(
//      mDatasetManager.buildRootPath() +
//      File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));
//  std::string graphFilename =
//      File::conformPathToOS(dumpDirectory + "/" + dumpPrefix + "_graph.png");
//  std::cout << "dumping "
//            << fullDatasetPath + mDatasetManager.currentGraphName.get()
//            << " -> " << graphFilename << std::endl;
//  if (!File::copy(File::conformPathToOS(fullDatasetPath +
//                                        mDatasetManager.currentGraphName.get()),
//                  graphFilename)) {
//    std::cerr << "ERROR copying png file." << std::endl;
//  }

//  auto allPositions = mDatasetManager.occupationData.get(false);
//  File allPositionsFile(File::conformPathToOS(dumpDirectory + "/" + dumpPrefix
//  +
//                                              "_positions.csv"),
//                        "w", true);
//  std::string header = "element,x,y,z\n";
//  allPositionsFile.write(header);
//  for (auto &elemPositions : *allPositions) {
//    std::string line =
//        mDatasetManager.mCurrentBasis[elemPositions.basis_index]
//                                     [elemPositions.occupancy_dof] +
//        ",";
//    line += std::to_string(elemPositions.x) + ",";
//    line += std::to_string(elemPositions.y) + ",";
//    line += std::to_string(elemPositions.z) + "\n";
//    allPositionsFile.write(line);
//  }
//  allPositionsFile.close();

//  auto between_planes = [&](Vec3f &point, Vec3f &plane_point,
//                            Vec3f &plane_normal, float &second_plane_distance)
//                            {
//    Vec3f difference = point - plane_point;
//    float proj = plane_normal.dot(difference);
//    if (proj >= 0 && proj <= second_plane_distance) {
//      return true;
//    } else {
//      return false;
//    }
//  };

//  auto plane_point = atomrender.mSlicingPlanePoint.get();
//  auto plane_normal = atomrender.mSlicingPlaneNormal.get().normalized();
//  auto second_plane_distance = atomrender.mSlicingPlaneThickness.get();
//  File slicePositionsFile(File::conformPathToOS(dumpDirectory + "/" +
//                                                dumpPrefix +
//                                                "_slice_positions.csv"),
//                          "w", true);
//  std::string slicePositionsHeader = "element,x,y,z\n";
//  slicePositionsFile.write(slicePositionsHeader);
//  {
//    auto allPositions = mDatasetManager.occupationData.get(false);
//    for (auto &elemPositions : *allPositions) {

//      Vec3f pos(elemPositions.x, elemPositions.y, elemPositions.z);
//      if (between_planes(pos, plane_point, plane_normal,
//                         second_plane_distance)) {
//        std::string line =
//            mDatasetManager.mCurrentBasis[elemPositions.basis_index]
//                                         [elemPositions.occupancy_dof] +
//            ",";
//        line += std::to_string(elemPositions.x) + ",";
//        line += std::to_string(elemPositions.y) + ",";
//        line += std::to_string(elemPositions.z) + "\n";
//        allPositionsFile.write(line);
//      }
//    }
//  }
//  slicePositionsFile.close();

//  json metadata;
//  metadata["dataset"]["path"] = mDatasetManager.currentDataset();
//  metadata["dataset"]["subdir"] = mDatasetManager.getSubDir();
//  metadata["dataset"]["rootpath"] = mDatasetManager.buildRootPath();
//  std::string condition = std::to_string(
//      mDatasetManager.mParameterSpaces[mDatasetManager.mConditionsParameter]
//          ->getCurrentIndex());
//  metadata["dataset"]["condition"] = condition;
//  for (auto space : mDatasetManager.mParameterSpaces) {
//    if (space.second && space.second->size() > 0) {
//      metadata["parameters"][space.first] =
//      space.second->getAllCurrentIds()[0];
//    }
//  }

//  metadata["SlicingPlanePoint"] = {plane_point.x, plane_point.y,
//  plane_point.z}; metadata["SliceNormal"] = {plane_normal.x, plane_normal.y,
//  plane_normal.z};

//  metadata["SliceThickness"] = atomrender.mSlicingPlaneThickness.get();
//  std::ofstream metadatafile(dumpDirectory + "/" + dumpPrefix +
//                             "_metadata.json");
//  metadatafile << std::setw(4) << metadata;
//}

void DataDisplay::updateDisplayBuffers() {
  //  std::map<string, int> elementCounts;

  auto allPositions = mDatasetManager.occupationData.get();

  mDatasetManager.mCurrentLoadedIndeces.clear();
  for (auto &paramSpace : mDatasetManager.mParameterSpaces) {
    mDatasetManager.mCurrentLoadedIndeces[paramSpace.first] =
        paramSpace.second->getCurrentIndex();
  }

  auto curVisibleAtoms = mShowAtoms.getSelectedElements();
  // TODO these colors should be exposed as a preference

  vector<Color> colorList = {
      Color(0.0, 1.0, 0.0, 0.0), Color(0.0, 0.0, 0.0, 0.0),
      Color(0.0, 1.0, 1.0, 1.0), Color(1.0, 1.0, 0.0, 1.0),
      Color(1.0, 0.0, 1.0, 1.0), Color(0.0, 1.0, 1.0, 1.0),
      Color(1.0, 1.0, 0.0, 1.0), Color(1.0, 0.0, 1.0, 1.0)};

  //  vector<Color> colorList2 = {
  //      Color(0.7f, 0.0f, 0.7f, 0.4f), Color(0.7f, 0.0f, 0.7f, 1.0f),
  //      Color(0.0f, 0.7f, 0.0f, 0.4f), Color(0.0f, 0.7f, 0.0f, 1.0f),
  //      Color(0.0f, 0.0f, 0.7f, 0.4f), Color(0.0f, 0.0f, 0.7f, 1.0f),
  //      Color(0.7f, 0.0f, 0.7f, 0.4f), Color(0.7f, 0.0f, 0.7f, 1.0f),
  //      Color(0.0f, 0.7f, 0.0f, 0.4f), Color(0.0f, 0.7f, 0.0f, 1.0f),
  //      Color(0.0f, 0.0f, 0.7f, 0.4f), Color(0.0f, 0.0f, 0.7f, 1.0f),
  //  };

  mDataBoundaries.resetInv();
  mAligned4fData.clear();

  mAtomData.clear();

  auto colorIt = colorList.begin();
  for (auto visibleAtom : curVisibleAtoms) {
    if (elementData.find(visibleAtom) != elementData.end()) {
      mAtomData[visibleAtom] =
          tinc::AtomData{0, visibleAtom, elementData[visibleAtom].radius,
                         elementData[visibleAtom].color};
    } else {

      mAtomData[visibleAtom] = tinc::AtomData{0, visibleAtom, 1.0, *colorIt};
      colorIt++;
    }
  }

  if (mDatasetManager.mParameterSpaces["time"]->size() > 0) {
    // This is a Kinetic MC dataset.
    auto currentIndex =
        mDatasetManager.mParameterSpaces["time"]->getAllCurrentIndeces()[0];
    auto templateDataIt = mDatasetManager.templateData.begin();
    uint8_t *occupationPtr =
        &(mDatasetManager.trajectoryData
              .data()[mDatasetManager.numAtoms * currentIndex]);

    if (mDatasetManager.templateData.size() != mDatasetManager.numAtoms) {
      std::cerr << "Error template size mismatch" << std::endl;
    }

    uint64_t counter = 0;
    uint64_t atomsInBasis = mDatasetManager.templateData.size() /
                            mDatasetManager.mCurrentBasis.size();
    size_t basis_index = 0;
    while (templateDataIt != mDatasetManager.templateData.end()) {

      std::string atomName =
          mDatasetManager
              .mCurrentBasis[basis_index]["occupant_dof"][*occupationPtr];
      if (std::find(curVisibleAtoms.begin(), curVisibleAtoms.end(), atomName) !=
          curVisibleAtoms.end()) {
        mAligned4fData.push_back(templateDataIt->x);
        mAligned4fData.push_back(templateDataIt->y);
        mAligned4fData.push_back(templateDataIt->z);
        mAtomData[atomName].counts++;

        auto hue = rgb2hsv(mAtomData[atomName].color.rgb()).h;
        mAligned4fData.push_back(hue);
        Vec3f vec(templateDataIt->x, templateDataIt->y, templateDataIt->z);
        mDataBoundaries.includePoint(vec);
      }
      counter++;
      if (counter == atomsInBasis) {
        counter = 0;
        basis_index++;
      }
      occupationPtr++;
      templateDataIt++;
    }

  } else {
    // Dataset is Grand Canonical MC
    auto templateDataIt = mDatasetManager.templateData.begin();
    for (auto &atom : *allPositions) {
      std::string atomName =
          mDatasetManager.mCurrentBasis[atom.basis_index]["occupant_dof"]
                                       [atom.occupancy_dof];
      if (std::find(curVisibleAtoms.begin(), curVisibleAtoms.end(), atomName) !=
          curVisibleAtoms.end()) {
        mAligned4fData.push_back(templateDataIt->x);
        mAligned4fData.push_back(templateDataIt->y);
        mAligned4fData.push_back(templateDataIt->z);
        mAtomData[atomName].counts++;
      }
      templateDataIt++;
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
  //  atomPropertiesProj.clear();
  //  vector<AtomProperties> atomPropertiesPersp;

  //  auto colorListIt = colorList.begin();
  //  auto colorList2It = colorList2.begin();

  //  auto selectedElements = mShowAtoms.getSelectedElements();
  //  for (auto atom : mShowAtoms.getElements()) {
  //    if (std::find(selectedElements.begin(), selectedElements.end(), atom) !=
  //        selectedElements.end()) {
  //      if (elementData.find(atom) != elementData.end()) {
  //        //                    std::cout << "Color for: " << atom << ":" <<
  //        //                    elementData[atom].color.r << " " <<
  //        //                    elementData[atom].color.g << " "<<
  //        //                    elementData[atom].color.b << " "  <<std::endl
  //        ;
  //        // Atom was matched in elements.ini file, so use those colors

  //        atomPropertiesProj.push_back(AtomProperties{
  //            atom, elementData[atom].radius, elementData[atom].color});

  //        atomPropertiesPersp.push_back(AtomProperties{
  //            atom, elementData[atom].radius, elementData[atom].color});
  //      } else { // Use defaults
  //        atomPropertiesProj.push_back(AtomProperties{atom, 1.0f,
  //        *colorListIt});

  //        atomPropertiesPersp.push_back(
  //            AtomProperties{atom, 1.0f, *colorList2It});
  //      }
  //    }
  //    colorListIt++;
  //    colorList2It++;
  //    colorList2It++;
  //  }

  // Apply colors to aligned data -----------
  //  std::vector<Color> colors;
  //  mAtomData.clear();
  //  for (auto elem : *allPositions) {
  //    elementCounts[elem.first] = elem.second.size() / 4;
  //  }
  //  for (auto atomProps : atomPropertiesProj) {
  //    colors.push_back(atomProps.color);
  //    mAtomData.push_back(
  //        {elementCounts[atomProps.name], atomProps.drawScale,
  //        atomProps.name});
  //  }
  //  float hue = 0.0f;

  //  if (colors.size() > 0) {
  //    hue = rgb2hsv(colors.front().rgb()).h;
  //  }
  //  if (mAtomData.size() > 0) {
  //    auto atomDataIt = mAtomData.begin();
  //    int atomCounter = 0;
  //    for (size_t i = 0; i < mAligned4fData.size() / 4; i++) {
  //      //            assert(atomCountsIt != mAtomCounts.end());
  //      if (atomDataIt != mAtomData.end() && --atomCounter <= 0) {
  //        atomCounter = atomDataIt->counts;
  //        atomDataIt++;
  //        hue = rgb2hsv(colors.front().rgb()).h;
  //        colors.pop_front();
  //      }
  //      mAligned4fData[4 * i + 3] = hue;
  //    }
  //  }
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

  g.pushViewport(0, 0, 2 * w, 2 * h);
  gl::scissorArea(0, 0, 2 * w, 2 * h);
  g.clear(backgroundColor);

  g.pushViewMatrix();
  //            g.viewMatrix(getLookAt({ 0.0f, 0.0f, cameraZ },
  //            { 0.0f, 0.0f, 0.0f },
  //            { 0.0f, 1.0f, 0.0f }));
  g.viewMatrix(getLookAt(atomrender.mSlicingPlanePoint,
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

  /*if (atomPropertiesProj.size() > 0) */ {
    // ----------------------------------------
    int cumulativeCount = 0;
    for (auto &data : mAtomData) {

      atomrender.instancing_mesh0.attrib_data(
          data.second.counts * sizeof(float),
          mAligned4fData.data() + (cumulativeCount * 4),
          data.second.counts / 4);
      cumulativeCount += data.second.counts;
      // now draw data with custom shader
      g.shader(atomrender.instancing_mesh0.shader);
      g.update(); // sends modelview and projection matrices
      g.shader().uniform("dataScale", scalingFactor);
      g.shader().uniform("layerSeparation", 1.0);
      // A scaling value of 4.0 found empirically...
      if (atomrender.mShowRadius == 1.0f) {
        g.shader().uniform("markerScale", atomrender.mAtomMarkerSize *
                                              atomrender.mMarkerScale);
        //        g.shader().uniform("markerScale", data.radius *
        //                                              atomrender.mAtomMarkerSize
        //                                              *
        //                                              atomrender.mMarkerScale);
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

    drawHistory(g);

    g.popMatrix();
    g.popViewMatrix();
    g.popViewport();
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
