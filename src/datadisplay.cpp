
#ifdef AL_WINDOWS
#define NOMINMAX
#include <Windows.h>
#undef far
#undef near
#undef DELETE
#endif

#include "datadisplay.hpp"

#include "imgui.h"

#include <condition_variable>

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

  mMarkerScale = 0.01f;

  // Initialization of graphics objects
  // #INSTANCED_RENDERING: init
  {
    addSphere(instancing_mesh0.mesh, 1, 12, 6);
    instancing_mesh0.mesh.update();
    instancing_mesh0.init(instancing_vert, instancing_frag,
                          1,         // location
                          4,         // num elements
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
        //        this->requestDataLoad();
      }
    });
  }

  mDatasetManager.mCurrentDataset.registerChangeCallback(
      [&](std::string value) {
        if (value != mDatasetManager.mCurrentDataset.get()) {
          this->requestInitDataset();
        }
      });

  mShowAtoms.registerChangeCallback([this](uint16_t value) {
    if (mShowAtoms.get() != value) {
      mDatasetManager.getAtomPositions();
    }
  });

  mSliceRotationPitch.registerChangeCallback([this](float value) {
    mSlicingPlaneNormal.setNoCalls(Vec3f(sin(mSliceRotationRoll),
                                         cos(mSliceRotationRoll) * sin(value),
                                         cos(value))
                                       .normalize());
  });

  mSliceRotationRoll.registerChangeCallback([this](float value) {
    mSlicingPlaneNormal.setNoCalls(Vec3f(sin(value),
                                         cos(value) * sin(mSliceRotationPitch),
                                         cos(mSliceRotationPitch))
                                       .normalize());
  });

  mSlicingPlaneNormal.registerChangeCallback([this](Vec3f value) {
    value = value.normalized();
    float pitch = std::atan(value.y / value.z);
    float roll = std::atan(value.x / value.z);
    mSliceRotationPitch.setNoCalls(pitch);
    mSliceRotationRoll.setNoCalls(roll);
  });

  slicePickable.pose.registerChangeCallback([this](Pose pose) {
    mSlicingPlanePoint.setNoCalls(pose.pos());
    // mSlicingPlaneNormal.setNoCalls();
  });

  mSlicingPlaneDistance.registerChangeCallback([this](float v) {
    auto m = slicePickable.bb.max;
    m.z = v;
    slicePickable.bb.set(slicePickable.bb.min, m);
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
      //      this->requestDataLoad();

      mDatasetManager.getAtomPositions();
    }
  });
  // TODO we don't need to do a full data load here, just recompute the graph
  mPlotXAxis.registerChangeCallback([this](float value) {
    if (mPlotXAxis.get() != value) {

      mDatasetManager.getAtomPositions();
      //      this->requestDataLoad();
    }
  });

  mDatasetManager.currentGraphName.registerChangeCallback(
      [this](std::string value) {
        std::string fullDatasetPath = File::conformPathToOS(
            mDatasetManager.buildRootPath() +
            File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));

        std::cout << "loading graph " << value << " at " << fullDatasetPath
                  << std::endl;
        loadGraphTexture(value);
      });

  mLabelFont.alignLeft();

  bundle << mDatasetManager.mRootPath;
  bundle << mDatasetManager.mCurrentDataset;
  bundle << mVisible;
  bundle << mDatasetManager.mParameterSpaces["temperature"]->parameter();
  bundle << mDatasetManager.mParameterSpaces["chempotA"]->parameter()
         << mDatasetManager.mParameterSpaces["chempotB"]->parameter();
  bundle << mDatasetManager.mParameterSpaces["time"]->parameter();
  bundle << mAtomMarkerSize;
  //        bundle << mAtomOfInterest;
  bundle << mShowAtoms;
  bundle << mPlotXAxis;
  bundle << mPlotYAxis;
  bundle << mShowGraph << mShowParallel /*<< mShowSurface */
         << mShowPerspective /*<< mSingleProjection*/;

  bundle << mLayerSeparation;
  bundle << mPerspectiveRotY;
  bundle << mSlicingPlaneNormal;
  bundle << mSliceRotationPitch << mSliceRotationRoll;
  bundle << mSlicingPlanePoint << mSlicingPlaneDistance << mLayerScaling;
  bundle << mShowGrid << mGridType << mGridSpacing << mGridXOffset
         << mGridYOffset;

  bundle << mShowRadius;
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
  }
  updateDisplayBuffers();

  if (/*mRequestLoad ||*/ mRequestInit) {
    mProcessing = true;
  } // Schedule processing for the start of next frame

  if (!mProcessing) { // If not currently doing computation prepare for drawing
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
          fabs(mDataBoundaries.maxx - mDataBoundaries.minx) /
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
      mTrajectoryMesh.vertex(previousPoint + thisMovement -
                             orthogonalVec * 0.2f);
      mTrajectoryMesh.color(c);
      mTrajectoryMesh.vertex(previousPoint + orthogonalVec);
      mTrajectoryMesh.color(c);
      mTrajectoryMesh.vertex(previousPoint + thisMovement +
                             orthogonalVec * 0.2f);

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

void DataDisplay::draw(Graphics &g) { // Load data after drawing frame to allow
                                      // showing "in process" state
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

    mDataBoundaries.reset();
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
        mDataBoundaries.addPoint(vec);
      }
    }

    if (mAligned4fData.size() > 0) {
      auto &b = mDataBoundaries;
      perspectivePickable.bb.set(Vec3f(b.minx, b.miny, b.minz),
                                 Vec3f(b.maxx, b.maxy, b.maxz));
      slicePickable.bb.set(Vec3f(b.minx, b.miny, b.minz),
                           Vec3f(b.maxx, b.maxy, (b.maxz - b.minz) * 0.5f));
      // rh.pose.pos().set(perspectivePickable.bb.cen);
      mSlicingPlanePoint.setHint("maxx", b.maxx);
      mSlicingPlanePoint.setHint("minx", b.minx - (b.maxx));
      mSlicingPlanePoint.setHint("maxy", b.maxy);
      mSlicingPlanePoint.setHint("miny", b.miny - (b.maxy));
      mSlicingPlanePoint.setHint("maxz", b.maxz);
      mSlicingPlanePoint.setHint("minz", b.minz - (b.maxz));
      mSlicingPlaneDistance.min(mDataBoundaries.minz);
      mSlicingPlaneDistance.max(mDataBoundaries.maxz);
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

          atomPropertiesPersp.push_back(
              AtomProperties{atom, elementData[atom].radius,
                             elementData[atom].color, Graphics::FILL});
          atomPropertiesPersp.push_back(
              AtomProperties{atom, elementData[atom].radius,
                             elementData[atom].color, Graphics::LINE});
        } else { // Use defaults
          atomPropertiesProj.push_back(
              AtomProperties{atom, 1.0f, *colorListIt});

          atomPropertiesPersp.push_back(
              AtomProperties{atom, 1.0f, *colorList2It, Graphics::FILL});
          atomPropertiesPersp.push_back(
              AtomProperties{atom, 1.0f, *colorList2It, Graphics::LINE});
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
      mAtomData.push_back({elementCounts[atomProps.name], atomProps.drawScale});
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
  float centerX = (mDataBoundaries.maxx + mDataBoundaries.minx) / 2.0f;
  float centerY = (mDataBoundaries.maxy + mDataBoundaries.miny) / 2.0f;
  float rangeSizeX = mDataBoundaries.maxx - mDataBoundaries.minx;
  float rangeSizeY = mDataBoundaries.maxy - mDataBoundaries.miny;
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
                               mDataBoundaries.minz - 100,
                               mDataBoundaries.maxz + 100));

  bool mAlignData = true;
  double scalingFactor = 1.0 / (mDataBoundaries.maxy - mDataBoundaries.miny);

  // std::cout << near << "..." << farClip <<std::endl;

  // save for later
  GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);
  g.scissorTest(true);
  int w = fbo_iso.width() / 2;
  int h = fbo_iso.height() / 2;
  g.polygonMode(Graphics::FILL);
  g.depthMask(true); // for axis rendering

  if (mSingleProjection.get() == 1.0f) {
    g.pushViewport(0, 0, 2 * w, 2 * h);
    g.scissor(0, 0, 2 * w, 2 * h);
    g.clear(backgroundColor);

    g.pushViewMatrix();
    //            g.viewMatrix(getLookAt({ 0.0f, 0.0f, cameraZ },
    //            { 0.0f, 0.0f, 0.0f },
    //            { 0.0f, 1.0f, 0.0f }));
    g.viewMatrix(getLookAt(mSlicingPlanePoint,
                           mSlicingPlanePoint.get() -
                               mSlicingPlaneNormal.get().normalized(),
                           {0.0f, 1.0f, 0.0f}));

    g.blending(false);
    g.depthTesting(true);
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

    g.blending(true);
    g.blendModeTrans();
    g.depthTesting(false);
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
          instancing_mesh0.attrib_data(
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
        g.shader(instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        g.shader().uniform("layerSeparation", 1.0);
        // A scaling value of 4.0 found empirically...
        if (mShowRadius == 1.0f) {
          g.shader().uniform("markerScale",
                             data.radius * mAtomMarkerSize * mMarkerScale);
        } else {
          g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);

        g.shader().uniform("plane_point", mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance", mSlicingPlaneDistance);

        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);

        instancing_mesh0.draw();
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
      g.scissor(w, 0, w, h);
      g.clear(backgroundColor);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({4, 0, 0}, {0, 0, 0}, {0, 1, 0}));

      g.blending(false);
      g.depthTesting(true);
      g.meshColor();
      g.draw(axis);

      g.blending(true);
      g.blendModeAdd();
      g.depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        g.shader().uniform("layerSeparation", 1.0);

        if (mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale);
        } else {
          g.shader().uniform("markerScale", mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);

        g.shader().uniform("plane_point", mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance", mSlicingPlaneDistance);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }

    // (from y, -z up), (left top)
    {
      g.pushViewport(0, h, w, h);
      g.scissor(0, h, w, h);
      g.clear(0, 0.2f, 0);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({0, 4, 0}, {0, 0, 0}, {0, 0, -1}));

      g.blending(false);
      g.depthTesting(true);
      g.meshColor();
      g.draw(axis);

      g.blending(true);
      g.blendModeAdd();
      g.depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        if (mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale);
        } else {
          g.shader().uniform("markerScale", mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);

        g.shader().uniform("plane_point", mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance", mSlicingPlaneDistance);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }

    // (from z, y up), (left bottom)
    {
      g.pushViewport(0, 0, w, h);
      g.scissor(0, 0, w, h);
      g.clear(0, 0, 0.2f);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({0, 0, 4}, {0, 0, 0}, {0, 1, 0}));

      g.blending(false);
      g.depthTesting(true);
      g.meshColor();
      g.draw(axis);

      g.blending(true);
      g.blendModeAdd();
      g.depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        if (mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale);
        } else {
          g.shader().uniform("markerScale", mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);
        g.shader().uniform("plane_point", mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance", mSlicingPlaneDistance);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }

    // diagonal view, right-top
    {
      g.pushViewport(w, h, w, h);
      g.scissor(w, h, w, h);
      g.clear(0.2f, 0.2f, 0.2f);

      g.pushViewMatrix();
      g.viewMatrix(getLookAt({3, 3, 3}, {0, 0, 0}, {0, 1, 0}));

      g.blending(false);
      g.depthTesting(true);
      g.meshColor();
      g.draw(axis);

      g.blending(true);
      g.blendModeAdd();
      g.depthTesting(false);

      if (atomPropertiesProj.size() > 0) {
        // now draw data with custom shader
        g.shader(instancing_mesh0.shader);
        g.update(); // sends modelview and projection matrices
        g.shader().uniform("dataScale", scalingFactor);
        if (mShowRadius == 1.0f) {
          g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale);
        } else {
          g.shader().uniform("markerScale", mMarkerScale);
        }
        g.shader().uniform("is_line", 0.0f);
        g.shader().uniform("is_omni", 0.0f);
        g.shader().uniform("plane_point", mSlicingPlanePoint.get());
        g.shader().uniform("plane_normal",
                           mSlicingPlaneNormal.get().normalized());
        g.shader().uniform("second_plane_distance", mSlicingPlaneDistance);
        //                    g.shader().uniform("far_clip", farClip);
        //                    g.shader().uniform("near_clip", near);
        g.shader().uniform("clipped_mult", 0.0);
        instancing_mesh0.draw();
      }

      g.popViewMatrix();
      g.popViewport();
    }
  }

  // put back scissoring
  g.scissorTest(last_enable_scissor_test);
  g.popProjMatrix();
  g.popViewMatrix();
  g.popModelMatrix();
  g.popFramebuffer();
  g.popViewport();
}

void DataDisplay::drawPerspective(Graphics &g) {
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

  g.blending(true);
  g.blendModeAdd();
  g.depthTesting(false);

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

  int cumulativeCount = 0;
  // now draw data with custom shaderg.shader(instancing_mesh0.shader);
  g.shader(instancing_mesh0.shader);
  g.shader().uniform("dataScale",
                     1.0f / ((mDataBoundaries.maxy - mDataBoundaries.miny) *
                             perspectivePickable.scale));
  g.shader().uniform("layerSeparation", mLayerSeparation);
  g.shader().uniform("is_omni", 1.0f);
  g.shader().uniform("eye_sep", perspectivePickable.scale * g.lens().eyeSep() *
                                    g.eye() / 2.0f);
  // g.shader().uniform("eye_sep", g.lens().eyeSep() * g.eye() / 2.0f);
  g.shader().uniform("foc_len", g.lens().focalLength());

  g.shader().uniform("plane_point", mSlicingPlanePoint.get());
  g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
  g.shader().uniform("second_plane_distance", mSlicingPlaneDistance);
  //        g.shader().uniform("near_clip", near);
  //        g.shader().uniform("far_clip", farClip);
  g.shader().uniform("clipped_mult", 0.45);
  g.update();

  for (auto data : mAtomData) {
    if (mShowRadius == 1.0f) {
      g.shader().uniform("markerScale", data.radius * mAtomMarkerSize *
                                            mMarkerScale /
                                            perspectivePickable.scale);
      //                std::cout << data.radius << std::endl;
    } else {
      g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale /
                                            perspectivePickable.scale);
    }
    int count = data.counts;
    assert((int)mAligned4fData.size() >= (cumulativeCount + count) * 4);
    instancing_mesh0.attrib_data(count * 4 * sizeof(float),
                                 mAligned4fData.data() + (cumulativeCount * 4),
                                 count);
    cumulativeCount += count;

    g.polygonMode(Graphics::FILL);
    g.shader().uniform("is_line", 0.0f);
    instancing_mesh0.draw();

    //            g.shader().uniform("is_line", 1.0f);
    //            g.polygonMode(Graphics::LINE);
    //            instancing_mesh0.draw();
    //            g.polygonMode(Graphics::FILL);
  }

  g.depthTesting(true);
  g.blending(true);

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

  const float x[2] = {mDataBoundaries.maxx, mDataBoundaries.minx};
  const float y[2] = {mDataBoundaries.maxy, mDataBoundaries.miny};

  const float z[2] = {0.0, (mSlicingPlaneDistance)};

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
  g.polygonMode(Graphics::LINE);
  g.color(0.8f, 0.8f, 1.0f, 0.9f);
  g.scale(1.0, 1.0, 1.0f + mLayerSeparation);
  g.translate(mSlicingPlanePoint.get());
  g.rotate(mSliceRotationPitch * 360.0f / (M_2PI), 1.0, 0.0, 0.0);
  g.rotate(mSliceRotationRoll * 360.0f / (M_2PI), 0.0, -1.0, 0.0);
  g.draw(boxMesh);
  g.polygonMode(Graphics::FILL);
  glLineWidth(1);

  //  -----------

  g.popMatrix();

  drawHistory(g);

  g.popMatrix(); // pickable

  g.depthTesting(false);
  g.blending(true);
  g.blendModeAdd();
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
  g.depthTesting(true);
  g.draw(m);
  // Draw label
  if (mDrawLabels == 1.0f) {
    g.depthTesting(false);
    g.blending(true);
    g.blendModeAdd();
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
  // g.pushMatrix();
  // g.translate(mGraphPose.get().pos() );
  Vec3f forward = graphPickable.pose.get().pos();
  forward.normalize();

  if (mBillboarding.get() == 1.0f) {
    Quatf rot = Quatf::getBillboardRotation(-forward, Vec3f{0.0, 1.0f, 0.0f});
    g.rotate(rot);
  }
  // New image module puts origin on top right
  if (mGraphTextureLock.try_lock()) {
    if (mGraphFilePathToLoad.size() > 0) {

      Image img(mGraphFilePathToLoad.c_str());
      if (img.loaded()) {
        static bool messagePrinted = false;
        static std::string lastFailed;
        if (lastFailed != mGraphFilePathToLoad) {
          messagePrinted = false;
          lastFailed = mGraphFilePathToLoad;
        }
        if (!messagePrinted) {
#ifdef AL_BUILD_MPI
          char name[MPI_MAX_PROCESSOR_NAME];
          int len;
          int ret = MPI_Get_processor_name(name, &len);
          std::cout << name << ": ";
#endif
          cout << "failed to load image " << mGraphFilePathToLoad << endl;
          messagePrinted = true;
        }
        mGraphTexture.resize(4, 4);
        // Will try again on next frame...
      } else {
        //                         cout << "loaded image size: " <<
        //                         imageData.width << ", " << imageData.height
        //                         << endl;

        mGraphTexture.resize(img.width(), img.height());
        mGraphTexture.submit(img.pixels(), GL_RGBA, GL_UNSIGNED_BYTE);

        mGraphTexture.filter(Texture::LINEAR_MIPMAP_LINEAR);
        mGraphFilePath = mGraphFilePathToLoad;
        mGraphFilePathToLoad = "";
      }
    }
    mGraphTextureLock.unlock();
  }
  g.depthTesting(true);
  g.quad(mGraphTexture, -1.2f, 0.9f, 2.4f, -1.8f);
  // Draw label

  if (mDrawLabels == 1.0f) {
    g.depthTesting(false);
    g.blending(true);
    g.blendModeAdd();
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
  g.popMatrix();
  graphPickable.drawBB(g);
}
