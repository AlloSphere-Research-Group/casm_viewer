
#ifdef AL_WINDOWS
#define NOMINMAX
#include <Windows.h>
#undef far
#undef near
#undef DELETE
#endif

#include "datadisplay.hpp"

#include "imgui.h"

void DataDisplay::init() {

#ifdef AL_BUILD_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mWorldRank);
#endif

  mPickableManager << graphPickable << parallelPickable << perspectivePickable;

  graphPickable.bb.set(Vec3f(-1.2f,-0.9f,-0.05f), Vec3f(1.2f,0.9f,0.0f));
  graphPickable.pose = Pose(Vec3d(-2.0, 0.8, -.75), Quatf());
  graphPickable.scale = 0.2f;

  parallelPickable.bb.setCenterDim(Vec3f(), Vec3f(4,4,0.1f));
  parallelPickable.pose = Pose(Vec3d(-0.68, 1.1, -1.5), Quatf());
  parallelPickable.scale = 0.4f;

  perspectivePickable.pose = Pose(Vec3d(0.55, 0.35, -0.7), Quatf().fromEuler(75.0*M_DEG2RAD,0,0));
  perspectivePickable.scale = 0.02f;
  //    perspectivePickable.scaleVec = mPerspectiveScale.get();
  // perspectivePickable.addChild(rh);
  // rh.size = 25;
  // rh.dr = 2.2;

  perspectivePickable.addChild(slicePickable);
  perspectivePickable.testChildren = false;
  perspectivePickable.containChildren = true;

  //        mPerspectiveScale.registerChangeCallback([this](float s){ perspectivePickable.scaleVec.set(s); });
  mPerspectiveRotY.registerChangeCallback([this](float d){
    perspectivePickable.pose = Pose(perspectivePickable.pose.get().pos(), Quatf().fromEuler(-d*M_DEG2RAD,0,0));
  });

  mMarkerScale = 0.01f;

  // Initialization of graphics objects
  // #INSTANCED_RENDERING: init
  {
    addSphere(instancing_mesh0.mesh, 1, 12, 6);
    instancing_mesh0.mesh.update();
    instancing_mesh0.init(instancing_vert, instancing_frag,
                          1, // location
                          4, // num elements
                          GL_FLOAT); // type

    instancing_shader.compile(instancing_vert, instancing_frag);

    addSphere(orthoMesh0);
    orthoMesh0.update();
  }

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
  for (int i = 0 ; i < num_verts_added; i += 1) {
    axis.color(1, 0, 0);
  }

  // y
  num_verts_added = addCube(axis);
  transform.setIdentity();
  transform *= Matrix4f::rotation(M_PI / 2, 2, 1); // rotate from z to y
  transform *= Matrix4f::translation(0, 0, 0.3);
  transform *= Matrix4f::scaling(0.02, 0.02, 1);
  axis.transform(transform, axis.vertices().size() - num_verts_added);
  for (int i = 0 ; i < num_verts_added; i += 1) {
    axis.color(0, 1, 0);
  }

  // z
  num_verts_added = addCube(axis);
  transform.setIdentity();
  transform *= Matrix4f::translation(0, 0, 0.3);
  transform *= Matrix4f::scaling(0.02, 0.02, 1);
  axis.transform(transform, axis.vertices().size() - num_verts_added);
  for (int i = 0 ; i < num_verts_added; i += 1) {
    axis.color(0, 0, 1);
  }
  axis.update();

  // Parameters setup

  for (auto &parameterSpace: mDatasetManager.mParameterSpaces) {
    parameterSpace.second->parameter().registerChangeCallback([&](float value){

      if (parameterSpace.second->getCurrentId()
          != parameterSpace.second->idAt(parameterSpace.second->getIndexForValue(value))) { // Only reload if id has changed

        parameterSpace.second->parameter().setNoCalls(value); // To have the internal value already changed for the following functions.
        std::cout << value << " : " << parameterSpace.second->getCurrentId() << "..."
                  << parameterSpace.second->idAt(parameterSpace.second->getIndexForValue(value));
        mDatasetManager.getAtomPositions();
//        this->requestDataLoad();
      }
    });
  }

  mDatasetManager.mCurrentDataset.registerChangeCallback([&](std::string value){
    if (value != mDatasetManager.mCurrentDataset.get()) {
      this->requestInitDataset();
//      this->requestDataLoad();
    }
  });

  mShowAtoms.registerChangeCallback([this](uint16_t value) {
    if (mShowAtoms.get() != value) {
//      this->requestDataLoad();
    }
  });

  mSliceRotationPitch.registerChangeCallback([this](float value) {
    mSlicingPlaneNormal.setNoCalls(Vec3f(sin(mSliceRotationRoll), cos(mSliceRotationRoll) *sin(value), cos(value)).normalize());
  });

  mSliceRotationRoll.registerChangeCallback([this](float value) {
    mSlicingPlaneNormal.setNoCalls(Vec3f(sin(value), cos(value) * sin(mSliceRotationPitch), cos(mSliceRotationPitch)).normalize());
  });

  mSlicingPlaneNormal.registerChangeCallback([this](Vec3f value) {
    value = value.normalized();
    float pitch = std::atan(value.y/value.z);
    float roll = std::atan(value.x/value.z);
    mSliceRotationPitch.setNoCalls(pitch);
    mSliceRotationRoll.setNoCalls(roll);
  });

  slicePickable.pose.registerChangeCallback([this](Pose pose) {
    mSlicingPlanePoint.setNoCalls(pose.pos());
    // mSlicingPlaneNormal.setNoCalls();
  });

  mSlicingPlaneDistance.registerChangeCallback([this](float v){
    auto m = slicePickable.bb.max;
    m.z = v;
    slicePickable.bb.set(slicePickable.bb.min, m);
  });

  backgroundColor.setHint("showAlpha", 1.0);

  mGridType.setElements({"square", "triangle"});
  mShowGrid.registerChangeCallback([this](float value){
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
  mPlotYAxis.registerChangeCallback([this](float value){
    if (mPlotYAxis.get() != value) {
//      this->requestDataLoad();
    }
  });
  // TODO we don't need to do a full data load here, just recompute the graph
  mPlotXAxis.registerChangeCallback([this](float value){
    if (mPlotXAxis.get() != value) {
//      this->requestDataLoad();
    }
  });

  mDatasetManager.currentGraphName.registerChangeCallback([this](std::string value) {
    std::string fullDatasetPath = File::conformPathToOS(mDatasetManager.buildRootPath()
                                                        + File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));

    std::cout << "loading graph " << value << " at " << fullDatasetPath <<std::endl;
    loadGraphTexture(value);
  });

  mLabelFont.align(TEXT_ALIGN::LEFT);

  bundle << mDatasetManager.mRootPath;
  bundle << mDatasetManager.mCurrentDataset;
  bundle << mVisible;
  bundle << mDatasetManager.mParameterSpaces["temperature"]->parameter();
  bundle << mDatasetManager.mParameterSpaces["chempotA"]->parameter() << mDatasetManager.mParameterSpaces["chempotB"]->parameter();
  bundle << mDatasetManager.mParameterSpaces["time"]->parameter();
  bundle << mAtomMarkerSize;
  //        bundle << mAtomOfInterest;
  bundle << mShowAtoms;
  bundle << mPlotXAxis;
  bundle << mPlotYAxis;
  bundle << mShowGraph << mShowParallel /*<< mShowSurface */ << mShowPerspective /*<< mSingleProjection*/;

  bundle << mLayerSeparation;
  bundle << mPerspectiveRotY;
  bundle << mSlicingPlaneNormal;
  bundle << mSliceRotationPitch << mSliceRotationRoll;
  bundle << mSlicingPlanePoint << mSlicingPlaneDistance << mLayerScaling;
  bundle << mShowGrid << mGridType << mGridSpacing << mGridXOffset << mGridYOffset;

  bundle << mShowRadius;
  bundle << mCumulativeTrajectory << mIndividualTrajectory;
  bundle << mBillboarding;
  bundle << mSmallLabel;
  bundle << mDrawLabels;

}

void DataDisplay::initRootDirectory() {
  mDatasetManager.mLoadedDataset = "";

  mDatasetManager.initRoot();

  DatasetManager::SpeciesLabelMap labelMap = mDatasetManager.getAvailableSpecies();

  std::vector<std::string> availableAtoms;
  std::vector<std::string> atomLabels;
  for(auto species: labelMap) {
    availableAtoms.push_back(species.first);
    if (species.first != "Va") {
      // Vacancy atom names will be there already so wee need to skip to avoid duplication
      atomLabels.insert(atomLabels.end(), species.second.begin(), species.second.end());
    }
  }
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
    for(auto atomName: labelMap["Va"]) {
      // Use only first result
      if (std::regex_search(atomName, match, atomNameRegex)) {
        string result = match.str();
        ptrdiff_t pos = std::distance(availableAtoms.begin(), find(availableAtoms.begin(), availableAtoms.end(), result));
        if (pos < (int) availableAtoms.size()) {
          mAtomOfInterest.set((int) pos);
        }
        mShowAtoms.setElementSelected(atomName);
      }
    }
  }
  auto dataNames = mDatasetManager.getDataNames();
  mPlotYAxis.setElements(dataNames);
  string defaultYAxis= "<comp_n(" + mAtomOfInterest.getCurrent() + ")>";
  ptrdiff_t pos = std::find(dataNames.begin(), dataNames.end(), defaultYAxis) - dataNames.begin();
  if (pos < (int) dataNames.size()) {
    mPlotYAxis.set((int) pos);
  }

  std::vector<std::string> parameterSpaceNames;
  parameterSpaceNames.push_back("temperature");

  mPlotXAxis.setElements(parameterSpaceNames);
  mPlotYAxis.set(0);
  //        if (mConfig.preProcess) {
  ////            std::cout << "Starting preprocess " << std::endl;
  //            preProcessDataset();
  //        }
//  requestDataLoad();
}

void DataDisplay::prepare(Graphics &g, Matrix4f &transformMatrix) {
  if (mProcessing) {
    if (mRequestInit) {
      initRootDirectory();
      mRequestInit = false;
      resetSlicing();
    } /*else if (mRequestLoad)  {
      loadCurrentData();
      mRequestLoad = false;
    }*/
    mProcessing = false;
    updateParameterText();
  }
  updateDisplayBuffers();

  if (/*mRequestLoad ||*/ mRequestInit)  {
    mProcessing = true;
  } // Schedule processing for the start of next frame

  if (!mProcessing) { // If not currently doing computation prepare for drawing
    // History mesh displays individual movements from their actual positions
    mHistoryMesh.primitive(Mesh::TRIANGLES);
    mHistoryMesh.reset();

    float historyWidth = 0.65f;
    size_t counter = mDatasetManager.mHistory.size() - 1;
    for (auto historyPoint = mDatasetManager.mHistory.begin(); historyPoint != mDatasetManager.mHistory.end(); historyPoint++) {

      HSV hsvColor(0.5f * float(counter)/mDatasetManager.mHistory.size(), 1.0, 1.0);
      Color c;

      // Assumes the plane's normal is the z-axis
      Vec3f orthogonalVec = (historyPoint->second - historyPoint->first).cross({0, 0, 1}).normalize(historyWidth);
      Vec3f orthogonalVec2 = (historyPoint->second - historyPoint->first).cross({1, 0, 0}).normalize(historyWidth);
      if (orthogonalVec2.mag() < 0.0001f) {
        orthogonalVec2 = (historyPoint->second - historyPoint->first).cross({0, 1, 0}).normalize(historyWidth);
      }
      assert(orthogonalVec2.mag() > 0.0001f);
      ImGui::ColorConvertHSVtoRGB(hsvColor.h, hsvColor.s, hsvColor.v, c.r, c.g, c.b);
      c.a = 0.35f;
      unsigned int previousSize = mHistoryMesh.vertices().size();
      mHistoryMesh.color(c);
      mHistoryMesh.vertex(historyPoint->first - orthogonalVec);
      mHistoryMesh.color(c);
      mHistoryMesh.vertex(historyPoint->second - orthogonalVec* 0.1f);
      mHistoryMesh.color(c);
      mHistoryMesh.vertex(historyPoint->first + orthogonalVec);
      mHistoryMesh.color(c);
      mHistoryMesh.vertex(historyPoint->second + orthogonalVec* 0.1f);

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
    Vec3f previousPoint(0,0,0);
    float previousMag = 0.0;

    float trajectoryWidth = 0.35f;
    counter = mDatasetManager.mHistory.size() - 1;
    for (auto historyPoint = mDatasetManager.mHistory.begin(); historyPoint != mDatasetManager.mHistory.end(); historyPoint++) {

      HSV hsvColor(0.5f * float(counter)/mDatasetManager.mHistory.size(), 1.0, 1.0);
      Color c;

      // Assumes the plane's normal is the z-axis
      Vec3f thisMovement = historyPoint->second - historyPoint->first;
      Vec3f orthogonalVec = thisMovement.cross({0, 0, 1}).normalize(trajectoryWidth);
      ImGui::ColorConvertHSVtoRGB(hsvColor.h, hsvColor.s, hsvColor.v, c.r, c.g, c.b);
      c.a = 0.8f;
      unsigned int previousSize = mTrajectoryMesh.vertices().size();
      if (thisMovement.mag() > fabs(mDataBoundaries.maxx - mDataBoundaries.minx)/2.0f) { // Atom is wrapping around
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
}

void DataDisplay::draw(Graphics &g) {        // Load data after drawing frame to allow showing "in process" state
  if (mVisible == 0.0f) {
    return;
  }

  std::unique_lock<std::mutex> lk(mDrawLock);
  // g.depthTesting(false);
  g.blendModeAdd();

  g.lens().eyeSep(0.2);

  if (mShowPerspective.get() == 1.0f) {
    drawPerspective(g);
  }

  // now draw graph and iso view 3d slabs
  // set view to identity, effectively putting coord in eye space
  // g.pushViewMatrix(Matrix4f::identity());
  g.texture();
  g.blendModeTrans();

  if (mShowParallel.get() == 1.0f) {
    drawParallelProjection(g);
  }

  g.blending(false);
  if (mShowGraph.get() == 1.0f) {
    drawGraph(g);
  }
  // g.popViewMatrix();
}

void DataDisplay::setFont(string name, float size) {
  std::unique_lock<std::mutex> lk(mDrawLock);
  if (File::exists(name)) {
    if (!mLabelFont.load(name, size)) {
      std::cout << "Failed to load font: " << name << std::endl;
      if (!mLabelFont.load("C:\\Windows\\Fonts\\arial.ttf")) {
      }

      if (!mLabelFont.load("/usr/share/fonts/truetype/ttf-bitstream-vera/Vera.ttf")) {
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
  std::string dumpDirectory = File::conformPathToOS(mDatasetManager.buildRootPath() + mDatasetManager.mCurrentDataset.get() + "/graphics");
  if (!File::exists(dumpDirectory)) {
    if (!Dir::make(dumpDirectory)) {
      std::cerr << "Failed to create directory: " << dumpDirectory <<std::endl;
      return;
    }
  }

  std::string fullDatasetPath = File::conformPathToOS(mDatasetManager.buildRootPath()
                                                      + File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));
  std::string graphFilename = File::conformPathToOS(dumpDirectory + "/" + dumpPrefix + "_graph.png");
  std::cout << "dumping " << fullDatasetPath + mDatasetManager.currentGraphName.get() << " -> " << graphFilename << std::endl;
  if (!File::copy( File::conformPathToOS(fullDatasetPath + mDatasetManager.currentGraphName.get()) ,
                   graphFilename)) {
    std::cerr << "ERROR copying png file." <<std::endl;
  }

  auto allPositions = mDatasetManager.positionBuffers.get(false);
  File allPositionsFile(File::conformPathToOS(dumpDirectory + "/" + dumpPrefix + "_positions.csv"), "w", true);
  std::string header = "element,x,y,z\n";
  allPositionsFile.write(header);
  for (auto elemPositions: *allPositions) {
    auto positionIt = elemPositions.second.begin();
    while(positionIt != elemPositions.second.end()) {
      std::string line = elemPositions.first + ",";
      line += std::to_string(*positionIt++) + ",";
      line += std::to_string(*positionIt++) + ",";
      line += std::to_string(*positionIt++) + "\n";
      positionIt++; // Ignore label?
      allPositionsFile.write(line);
    }
  }
  allPositionsFile.close();

}
