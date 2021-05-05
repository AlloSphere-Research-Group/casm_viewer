
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
#include "modalsynth.hpp"

void DataDisplay::init() {
#ifdef AL_BUILD_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mWorldRank);
#endif

  for (size_t i = 0; i < graphCount; i++) {
    std::string name = "graph" + std::to_string(i + 1);
    std::string filename = "currentGraph" + std::to_string(i + 1) + ".png";
    imageDiskBuffers.emplace_back(
        std::make_unique<DiskBufferImage>(name, filename, "cachedGraph"));
    graphPickables.emplace_back(std::make_unique<PickableBB>(name));

    std::string graphName = "graph" + std::to_string(i + 1) + "name";
    currentGraphNames.push_back(
        std::make_unique<ParameterString>(graphName, "internal", ""));

    mPickableManager << *graphPickables.back();
    graphPickables.back()->bb.set(Vec3f(-1.2f, -0.9f, -0.05f),
                                  Vec3f(1.2f, 0.9f, 0.0f));
    graphPickables.back()->pose =
        Pose(Vec3d(-2.0 + 0.5 * (i + 1), 0.2, -.75), Quatf());
    graphPickables.back()->scale = 0.2f;

    mShowGraphs.emplace_back(std::make_unique<ParameterBool>(
        "ShowGraph" + std::to_string(i), "", 1));
  }
  synth.gain(5.0);

  // We need to process the callbacks for this within the graphics thread.
  mDatasetManager.mCurrentDataset.setSynchronousCallbacks(false);
  mDatasetManager.mCurrentDataset.registerChangeCallback(
      [&](std::string value) {
        // First force application of value
        lastError.set("");
        //        mDatasetManager.mCurrentDataset.setNoCalls(value);
        if (!initDataset()) {
          lastError.set(mDatasetManager.lastError);
        }
      });

  mDatasetManager.mShellSiteTypes.registerChangeCallback([&](uint16_t value) {
    if (value & 1) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 1500;
      voice->setFrequencies(voice->smallHandBell, 0.0004);
      voice->globalAmp = 30.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 1) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 1200;
      voice->setFrequencies(voice->xylo, 0.07);
      voice->globalAmp = 3.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 2) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 800;
      voice->setFrequencies(voice->smallHandBell, 0.0004);
      voice->globalAmp = 30.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 3) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 650;
      voice->setFrequencies(voice->xylo, 0.004);
      voice->globalAmp = 2.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 4) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 420;
      voice->setFrequencies(voice->smallHandBell, 0.0004);
      voice->globalAmp = 40.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 5) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 300;
      voice->setFrequencies(voice->xylo, 0.004);
      voice->globalAmp = 7.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 6) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 280;
      voice->setFrequencies(voice->smallHandBell, 0.0004);
      voice->globalAmp = 40.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 7) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 210;
      voice->setFrequencies(voice->xylo, 0.004);
      voice->globalAmp = 7.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 8) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 1100;
      voice->setFrequencies(voice->aluminiumBar, 0.004);
      voice->globalAmp = 6.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 9) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 800;
      voice->setFrequencies(voice->tubularBell, 0.0004);
      voice->globalAmp = 40.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 10) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 570;
      voice->setFrequencies(voice->aluminiumBar, 0.004);
      voice->globalAmp = 6.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 11) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 400;
      voice->setFrequencies(voice->tubularBell, 0.0004);
      voice->globalAmp = 30.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 12) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 3180;
      voice->setFrequencies(voice->aluminiumBar, 0.004);
      voice->globalAmp = 6.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 13) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 259;
      voice->setFrequencies(voice->tubularBell, 0.0004);
      voice->globalAmp = 30.0;
      synth.triggerOn(voice);
    }
    if (value & 1 << 14) {
      auto voice = synth.getVoice<ModalVoice>();
      voice->fundamentalFreq = 223;
      voice->setFrequencies(voice->aluminiumBar, 0.004);
      voice->globalAmp = 2.0;
      synth.triggerOn(voice);
    }
  });

  mShowAtoms.registerChangeCallback([this](uint16_t value) {
    if (mShowAtoms.get() != value) {
      // New values must be available for computeNewSample, throws away previous
      // value
      mShowAtoms.setNoCalls(value);
      mDatasetManager.currentPoscarName.set(
          mDatasetManager.labelProcessor.getOutputDirectory() +
          mDatasetManager.labelProcessor.getOutputFileNames()[0]);
    }
  });

  atomrender.mSlicingPlaneThickness.registerChangeCallback([this](float v) {
    auto m = slicePickable.bb.max;
    m.z = perspectivePickable.bb.min.z + v;
    slicePickable.bb.set(slicePickable.bb.min, m);
  });

  atomrender.mSlicingPlaneCorner.registerChangeCallback([this](Vec3f v) {
    auto p = Pose(v - slicePickable.bb.min, slicePickable.pose.get().quat());
    slicePickable.pose.setNoCalls(p);
  });

  atomrender.mSlicingPlaneNormal.setHint("hide", 0.0);

  atomrender.init();

  mPerspectiveRotY.registerChangeCallback([this](float d) {
    perspectivePickable.pose = Pose(perspectivePickable.pose.get().pos(),
                                    Quatf().fromEuler(-d * M_DEG2RAD, 0, 0));
  });

  currentSelection.registerChangeCallback([this](int32_t value) {
    auto pos = mDatasetManager.templateData[value];
    selectedPosition.x = pos.x;
    selectedPosition.y = pos.y;
    selectedPosition.z = pos.z;
  });

  slicePickable.pose.registerChangeCallback([this](Pose pose) {
    auto p = pose.pos() + slicePickable.bb.min;
    atomrender.mSlicingPlaneCorner.setNoCalls(p);
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
            mDatasetManager.getGlobalRootPath() +
            File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));

        std::cout << "loading graph " << value << " at " << fullDatasetPath
                  << std::endl;
        // New image module puts origin on top right
        //        if (mGraphTextureLock.try_lock()) {
        //            if (mGraphFilePathToLoad.size() > 0) {
        imageDiskBuffer.loadData(value);
      });

  for (size_t i = 0; i < graphCount; i++) {
    currentGraphNames[i]->setSynchronousCallbacks();
    currentGraphNames[i]->registerChangeCallback([&](std::string value) {
      std::string fullDatasetPath = File::conformPathToOS(
          mDatasetManager.getGlobalRootPath() +
          File::conformPathToOS(mDatasetManager.mCurrentDataset.get()));

      std::cout << "loading graph [" << i << "] " << value << " at "
                << fullDatasetPath << std::endl;
      // New image module puts origin on top right
      //        if (mGraphTextureLock.try_lock()) {
      //            if (mGraphFilePathToLoad.size() > 0) {
      imageDiskBuffers[i]->loadData(value);
    });
  }

  mLabelFont.alignLeft();

  //  bundle << mDatasetManager.mRootPath;
  bundle << mDatasetManager.mCurrentDataset;
  bundle << mVisible;
  bundle << atomrender.mDataScale;
  bundle << atomrender.mAtomMarkerSize;

  bundle << mShowAtoms;
  // Colors

  for (int i = 0; i < 16; i++) {
    mColorList.emplace_back(
        std::make_shared<ParameterColor>("atomColor" + std::to_string(i)));
    mColorList.back()->set(colorList[i % colorList.size()]);
    mColorList.back()->setHint("hsv", 1.0);
    bundle << *mColorList.back();
  }

  for (int i = 0; i < DatasetManager::maxPercolationTypes; i++) {

    mPercoColorList.emplace_back(
        std::make_shared<ParameterColor>("percoColor" + std::to_string(i)));
    mPercoColorList.back()->set(colorList[(i + 4) % colorList.size()]);
    mPercoColorList.back()->setHint("hsv", 1.0);
  }

  bundle << mDatasetManager.mPlotXAxis;
  bundle << mDatasetManager.mPlotYAxis;
  bundle << mShowPerspective;

  bundle << mDisplaySlicing;
  bundle << mDisplayGraphs;

  bundle << mPerspectiveRotY;
  for (auto *p : atomrender.parameters()) {
    bundle << p;
  }
  bundle << mShowGrid << mGridType << mGridSpacing << mGridXOffset
         << mGridYOffset;
  bundle << mBillboarding;
  bundle << mSmallLabel;
  bundle << mDrawLabels;

  bundle << mColorTrigger;

  // Pickable

  mPickableManager << graphPickable << parallelPickable << perspectivePickable;

  graphPickable.bb.set(Vec3f(-1.2f, -0.9f, -0.05f), Vec3f(1.2f, 0.9f, 0.0f));
  graphPickable.pose = Pose(Vec3d(-2.0, 0.2, -.75), Quatf());
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

  initGL();
}

void DataDisplay::initGL() {

  addSphere(mMarker, 1, 6, 6);
  mMarker.update();

  // Box mesh for slicing
  boxMesh.reset();
  boxMesh.primitive(Mesh::LINES);

  const float x[2] = {0, 1};
  const float y[2] = {0, 1};
  const float z[2] = {0, 1};

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

  EasyFBOSetting setting;
  setting.mUseMipmap = true;
  setting.filterMin = GL_LINEAR_MIPMAP_LINEAR;
  setting.filterMag = GL_LINEAR_MIPMAP_LINEAR;

  // Settings for parallel projection rendering fbo
  fbo_iso.init(2048, 2048, setting);

  mGraphTexture.mipmap(true);
  for (auto &graph : mGraphTextures) {
    graph.mipmap(true);
  }

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
}

bool DataDisplay::initDataset() {
  bool ret = true;
  mDatasetManager.mLoadedDataset = "";

  ret &= mDatasetManager.initDataset();

  currentSelection.max(mDatasetManager.templateData.size() - 1);

  DatasetManager::SpeciesLabelMap labelMap =
      mDatasetManager.getAvailableSpecies();

  std::vector<std::string> availableAtoms;
  std::vector<std::string> atomLabels;
  for (auto &species : labelMap) {
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
  for (auto &elementLabel : mShowAtoms.getElements()) {
    if (elementLabel.substr(0, speciesOfInterest.size()) == speciesOfInterest) {
      mShowAtoms.setElementSelected(elementLabel);
    }
  }
  for (size_t i = 0; i < mColorList.size(); i++) {
    if (i < mShowAtoms.getElements().size()) {
      mColorList[i]->setHint("hide", 0.0); // Show in GUI
      mColorList[i]->displayName(mShowAtoms.getElements()[i]);
    } else {
      mColorList[i]->setHint("hide", 1.0); // hide in GUI
    }
  }

  mDatasetManager.mParameterSpace.runProcess(
      mDatasetManager.sampleProcessorGraph, {}, {});
  if (mDatasetManager.graphGenerator.getOutputFileNames().size() > 0) {
    mDatasetManager.currentGraphName.set(
        mDatasetManager.graphGenerator.getOutputDirectory() +
        mDatasetManager.graphGenerator.getOutputFileNames()[0]);
  }
  size_t count = 0;
  for (auto &perc : mDatasetManager.mPercolationTypes.getElements()) {

    mPercoColorList[count]->displayName(perc);
    count++;
  }

  resetSlicing();
  return ret;
}

void DataDisplay::prepare(Graphics &g) {

  if (mNeedsProcessing) {
    mDatasetManager.mCurrentDataset.processChange();
    mDatasetManager.currentGraphName.processChange();
    mNeedsProcessing = false;
    // Parameter text
    auto subDir = mDatasetManager.getSubDir();
    std::string temperatureId;
    if (mDatasetManager.mParameterSpace.getDimension("temperature")) {
      temperatureId =
          mDatasetManager.mParameterSpace.getDimension("temperature")
              ->getCurrentId();
    }
    std::string timeId;
    if (mDatasetManager.mParameterSpace.getDimension("time")) {
      timeId =
          mDatasetManager.mParameterSpace.getDimension("time")->getCurrentId();
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

  if (mDatasetManager.occupationData.newDataAvailable() ||
      mColorTrigger.get() == true) {
    mColorTrigger.set(false);
    updateDisplayBuffers();
  }
  updatePercolationBuffers();
  prepareParallelProjection(g);
}

void DataDisplay::draw(Graphics &g) { // Load data after drawing frame to allow
                                      // showing "in process" state
  if (mVisible == 0.0f) {
    return;
  }

  std::unique_lock<std::mutex> lk(mDrawLock);
  g.blendAdd();

  if (mShowPerspective.get() == 1.0f) {
    drawPerspective(g);
  }

  g.texture();
  g.blendTrans();

  if (mDisplaySlicing.get() == 1.0) {
    drawParallelProjection(g);
  }

  g.blending(false);

  if (mDisplayGraphs.get() == 1.0) {
    drawGraph(g);
  }
}

void DataDisplay::setFont(string name, float size) {
  std::unique_lock<std::mutex> lk(mDrawLock);
  if (File::exists(name)) {
    if (!mLabelFont.load(name.c_str(), size, 1024)) {
      std::cout << "Failed to load font: " << name << std::endl;
#ifdef AL_WINDOWS
      if (!mLabelFont.load("C:\\Windows\\Fonts\\arial.ttf", size, 1024)) {
        std::cout
            << "Attempt to load alternate font: C:\\Windows\\Fonts\\arial.ttf"
            << std::endl;
      }
#elif defined(AL_LINUX)
      if (!mLabelFont.load(
              "/usr/share/fonts/truetype/ttf-bitstream-vera/Vera.ttf", size,
              1024)) {
        std::cout << "Attempt to load alternate font: "
                     "/usr/share/fonts/truetype/ttf-bitstream-vera/Vera.ttf"
                  << std::endl;
      }
#elif defined(AL_OSX)
      if (!mLabelFont.load("/Library/Fonts/Arial Unicode.ttf", size, 1024)) {
        std::cout << "Attempt to load alternate font: /Library/Fonts/Arial "
                     "Unicode.ttf"
                  << std::endl;
      }
#endif
    } else {
      std::cout << "Loaded font " << name << std::endl;
    }
  } else {
    std::cout << "Could not sync font: " << name << std::endl;
  }
}

void DataDisplay::nextTemp() {
  mDatasetManager.mParameterSpace.getDimension("temperature")->stepIncrement();
}

void DataDisplay::previousTemp() {
  mDatasetManager.mParameterSpace.getDimension("temperature")->stepDecrease();
}

void DataDisplay::nextChempot() {
  mDatasetManager.mParameterSpace.getDimension("chempotA")->stepIncrement();
}

void DataDisplay::previousChempot() {
  mDatasetManager.mParameterSpace.getDimension("chempotA")->stepDecrease();
}

void DataDisplay::nextChempot2() {
  if (mDatasetManager.mParameterSpace.getDimension("chempotB")) {
    mDatasetManager.mParameterSpace.getDimension("chempotB")->stepIncrement();
  }
}

void DataDisplay::previousChempot2() {
  if (mDatasetManager.mParameterSpace.getDimension("chempotB")) {
    mDatasetManager.mParameterSpace.getDimension("chempotB")->stepDecrease();
  }
}

void DataDisplay::nextTime() {
  if (mDatasetManager.mParameterSpace.getDimension("time")) {
    mDatasetManager.mParameterSpace.getDimension("time")->stepIncrement();
  }
}

void DataDisplay::previousTime() {
  if (mDatasetManager.mParameterSpace.getDimension("time")) {
    mDatasetManager.mParameterSpace.getDimension("time")->stepDecrease();
  }
}

void DataDisplay::nextLayer() { atomrender.nextLayer(); }

void DataDisplay::previousLayer() { atomrender.previousLayer(); }

void DataDisplay::resetSlicing() {
  atomrender.mSlicingPlaneCorner.set(slicePickable.bb.min);

  atomrender.mSlicingPlaneThickness = atomrender.dataBoundary.dim.z;
  atomrender.mSlicingPlaneSize =
      std::max(atomrender.dataBoundary.dim.x, atomrender.dataBoundary.dim.y);
  atomrender.mSliceRotationRoll.set(0);
  atomrender.mSliceRotationPitch.set(0);
}

void DataDisplay::updateDisplayBuffers() {
  //  std::map<string, int> elementCounts;

  auto curVisibleAtoms = mShowAtoms.getSelectedElements();

  mTemplateDataBoundaries.resetInv();
  mAligned4fData.clear();

  mAtomData.clear();

  auto colorIt = mColorList.begin();
  for (auto atom : mShowAtoms.getElements()) {
    for (auto visibleAtom : curVisibleAtoms) {
      if (visibleAtom == atom) {
        mAtomData[visibleAtom] =
            tinc::AtomData{0, visibleAtom, 1.0, *colorIt->get()};
      }
    }
    colorIt++;
    if (colorIt == mColorList.end()) {
      colorIt = mColorList.begin();
    }
  }

  if (mDatasetManager.mParameterSpace.getDimension("time") &&
      mDatasetManager.mParameterSpace.getDimension("time")->size() > 0) {
    // This is a Kinetic MC dataset.
    mDatasetManager.occupationData.get(); // To mark as consumed.
    auto currentIndex =
        mDatasetManager.mParameterSpace.getDimension("time")->getCurrentIndex();
    auto templateDataIt = mDatasetManager.templateData.begin();
    uint8_t *occupationPtr =
        &(mDatasetManager.trajectoryData
              .data()[mDatasetManager.numAtoms * currentIndex]);

    uint8_t *prevOccupationPtr = nullptr;
    if (currentIndex > 0) {
      prevOccupationPtr =
          &(mDatasetManager.trajectoryData
                .data()[mDatasetManager.numAtoms * (currentIndex - 1)]);
    }

    if (mDatasetManager.templateData.size() != mDatasetManager.numAtoms) {
      std::cerr << "Error template size mismatch" << std::endl;
    }

    atomAdded.clear();
    atomRemoved.clear();
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
        mTemplateDataBoundaries.includePoint(vec);
      }
      if (prevOccupationPtr) {
        std::string prevAtomName =
            mDatasetManager
                .mCurrentBasis[basis_index]["occupant_dof"][*prevOccupationPtr];
        if (std::find(curVisibleAtoms.begin(), curVisibleAtoms.end(),
                      atomName) != curVisibleAtoms.end() ||
            std::find(curVisibleAtoms.begin(), curVisibleAtoms.end(),
                      prevAtomName) != curVisibleAtoms.end()) {
          if (prevAtomName != atomName) {
            if (atomName == "Va") {
              atomRemoved.push_back(*templateDataIt);

              previousSelection.set(std::distance(
                  mDatasetManager.templateData.begin(), templateDataIt));
            } else {
              if (atomAdded.find(atomName) == atomAdded.end()) {
                atomAdded[atomName] = std::vector<DatasetManager::position_t>();
              }
              atomAdded[atomName].push_back(*templateDataIt);

              currentSelection.set(std::distance(
                  mDatasetManager.templateData.begin(), templateDataIt));
            }
          }
        }
      }

      counter++;
      if (counter == atomsInBasis) {
        counter = 0;
        basis_index++;
      }
      if (prevOccupationPtr) {
        prevOccupationPtr++;
      }
      occupationPtr++;
      templateDataIt++;
    }
  } else {
    // Dataset is Grand Canonical MC

    auto allPositions = mDatasetManager.occupationData.get();
    auto templateDataIt = mDatasetManager.templateData.begin();
    for (auto &atom : *allPositions) {
      if (templateDataIt == mDatasetManager.templateData.end()) {
        std::cout << "Template and occupation size mismatch" << std::endl;
        break;
      }
      std::string atomName =
          mDatasetManager.mCurrentBasis[atom.basis_index]["occupant_dof"]
                                       [atom.occupancy_dof];
      if (std::find(curVisibleAtoms.begin(), curVisibleAtoms.end(), atomName) !=
          curVisibleAtoms.end()) {
        mAligned4fData.push_back(templateDataIt->x);
        mAligned4fData.push_back(templateDataIt->y);
        mAligned4fData.push_back(templateDataIt->z);

        auto hue = rgb2hsv(mAtomData[atomName].color.rgb()).h;
        mAligned4fData.push_back(hue);
        mAtomData[atomName].counts++;

        Vec3f vec(templateDataIt->x, templateDataIt->y, templateDataIt->z);
        mTemplateDataBoundaries.includePoint(vec);
      }
      templateDataIt++;
    }
  }

  if (mAligned4fData.size() > 0) {
    auto &b = mTemplateDataBoundaries;
    perspectivePickable.bb.set(Vec3f(b.min.x, b.min.y, b.min.z),
                               Vec3f(b.max.x, b.max.y, b.max.z));
    slicePickable.bb.set(Vec3f(b.min.x, b.min.y, b.min.z),
                         Vec3f(b.max.x, b.max.y, (b.max.z - b.min.z) * 0.25f));
    // rh.pose.pos().set(perspectivePickable.bb.cen);
    atomrender.setDataBoundaries(b);
    atomrender.setPositions(mAligned4fData.data(), mAligned4fData.size());
  }
}

void DataDisplay::updatePercolationBuffers() {

  for (size_t i = 0; i < mDatasetManager.maxPercolationTypes; i++) {
    if (mDatasetManager.percolationData[i].newDataAvailable()) {
      std::string name = std::to_string(i);
      mPercolationData[i][name].radius = 1.5;

      mPercolationData[i][name].color = mPercoColorList[i]->get();

      mPercolationData[i][name].counts = 0;
      mAlignedPercolation4fData[i].clear();

      auto allPositions = mDatasetManager.percolationData[i].get();
      auto templateDataIt = mDatasetManager.templateData.begin();
      for (auto &atom : *allPositions) {
        if (atom == 1) {
          mAlignedPercolation4fData[i].push_back(templateDataIt->x);
          mAlignedPercolation4fData[i].push_back(templateDataIt->y);
          mAlignedPercolation4fData[i].push_back(templateDataIt->z);

          auto hue = rgb2hsv(mPercolationData[i][name].color.rgb()).h;
          mAlignedPercolation4fData[i].push_back(hue);
          mPercolationData[i][name].counts++;

          Vec3f vec(templateDataIt->x, templateDataIt->y, templateDataIt->z);
        }

        templateDataIt++;
        if (templateDataIt == mDatasetManager.templateData.end()) {
          std::cout << "Template and occupation size mismatch" << std::endl;
          break;
        }
      }
    }
  }
}

void DataDisplay::drawTrajectories(Graphics &g) {
  g.pushMatrix();
  g.meshColor();
  g.polygonFill();
  g.blending(true);
  g.blendAdd();
  //  g.translate((mDataBoundaries.max.x - mDataBoundaries.min.x) / 2.0f,
  //              (mDataBoundaries.max.y - mDataBoundaries.min.y) / 2.0f,
  //              (mDataBoundaries.max.z - mDataBoundaries.min.z) / 2.0f);

  mHistoryRender.update(0);
  mTrajRender.update(0);

  mHistoryRender.onProcess(g);
  mTrajRender.onProcess(g);

  g.popMatrix();
}

void DataDisplay::prepareParallelProjection(Graphics &g) {
  //        std::unique_lock<std::mutex> lk(mDataLock);
  g.pushFramebuffer(fbo_iso);
  g.pushViewport(fbo_iso.width(), fbo_iso.height());
  g.pushProjMatrix();

  // projection matrix from clipping params
  // will only work with z axis aligned data and single viewport mode

  //        const float cameraZ = mDataBoundaries.maxz;

  float ar = float(fbo_iso.width()) / fbo_iso.height();
  float maxrange =
      std::max(atomrender.dataBoundary.dim.x, atomrender.dataBoundary.dim.y);
  float padding = maxrange * 0.05f;

  float scaling = atomrender.mSlicingPlaneSize;

  float left =
      /* atomrender.mDataScale.get().x * */ (
          atomrender.mSlicingPlaneCorner.get().x) -
      padding;
  float right =
      //      atomrender.mDataScale.get().x *
      (atomrender.mSlicingPlaneCorner.get().x + atomrender.mSlicingPlaneSize) +
      padding;
  float top =
      //      atomrender.mDataScale.get().y *
      (atomrender.mSlicingPlaneCorner.get().y + atomrender.mSlicingPlaneSize) +
      padding;
  float bottom =
      /*  atomrender.mDataScale.get().y **/ (
          atomrender.mSlicingPlaneCorner.get().y) -
      padding;

  g.projMatrix(Matrix4f::ortho(ar * left, ar * right, bottom, top,
                               mTemplateDataBoundaries.min.z - 100,
                               mTemplateDataBoundaries.max.z + 100));

  // save for later
  GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);
  g.scissorTest(true);
  int w = fbo_iso.width();
  int h = fbo_iso.height();
  g.polygonFill();

  {
    g.pushViewport(0, 0, w, h);
    g.scissorArea(0, 0, w, h);
    g.clear(backgroundColor);

    g.pushViewMatrix();
    g.viewMatrix(getLookAt(atomrender.mSlicingPlaneCorner,
                           atomrender.mSlicingPlaneCorner.get() -
                               atomrender.mSlicingPlaneNormal.get(),
                           {0.0f, 1.0f, 0.0f}));

    g.blending(false);
    g.depthTesting(true);

    if (mShowGrid == 1.0f) {
      mGridMesh.reset();
      addRect(mGridMesh, 0, -2.0, 0.1f / (scaling), 4.0f);
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
          g.scale(scaling);
          g.rotate(90, 0, 0, 1);
          g.draw(mGridMesh);
          g.popMatrix();
          g.pushMatrix();
          g.translate(mGridXOffset + mGridSpacing * i, mGridYOffset, 0);
          g.scale(scaling);
          g.draw(mGridMesh);
          g.popMatrix();
        } else if (mGridType.getCurrent() == "triangle") {
          float spacing = 0.86602540378f * mGridSpacing; // sin(60)
          g.pushMatrix();
          g.translate(mGridXOffset, mGridYOffset + spacing * i, 0);
          g.scale(scaling);
          g.rotate(90, 0, 0, 1);
          g.draw(mGridMesh);
          g.popMatrix();
          g.pushMatrix();
          g.translate(mGridXOffset, mGridYOffset, 0);
          g.rotate(30, 0, 0, 1);
          g.translate(spacing * i, 0, 0);
          g.scale(scaling);
          g.draw(mGridMesh);
          g.popMatrix();
          g.pushMatrix();
          g.translate(mGridXOffset, mGridYOffset, 0);
          g.rotate(-30, 0, 0, 1);
          g.translate(spacing * i, 0, 0);
          g.scale(scaling);
          g.draw(mGridMesh);
          g.popMatrix();
        }
      }
    }

    g.blending(true);
    g.blendTrans();
    g.depthTesting(false);
    {
      g.pushMatrix();
      // ----------------------------------------
      int cumulativeCount = 0;
      g.shader(atomrender.instancingMesh.shader);
      g.update(); // sends modelview and projection matrices
      g.shader().uniform("dataScale", atomrender.mDataScale.get());
      g.shader().uniform("is_omni", 0.0f);
      g.shader().uniform("eye_sep", g.lens().eyeSep() * g.eye() / 2.0f);
      //      g.shader().uniform("foc_len", g.lens().focalLength());

      g.shader().uniform("alpha", atomrender.mAlpha);
      g.shader().uniform("markerScale", atomrender.mAtomMarkerSize);

      g.shader().uniform("plane_point", atomrender.mSlicingPlaneCorner.get());
      g.shader().uniform("plane_normal",
                         atomrender.mSlicingPlaneNormal.get().normalized());

      // TODO add rotation vectors

      g.shader().uniform("slice_size", atomrender.mSlicingPlaneSize.get());
      g.shader().uniform("second_plane_distance",
                         atomrender.mSlicingPlaneThickness);
      g.shader().uniform("clipped_mult", 0.0);

      g.shader().uniform("is_line", 0.0f);

      for (auto &data : mAtomData) {

        int count = data.second.counts;
        // now draw data with custom shader
        atomrender.instancingMesh.attrib_data(
            count * 4 * sizeof(float),
            mAligned4fData.data() + (cumulativeCount * 4), count);
        cumulativeCount += data.second.counts;
        atomrender.instancingMesh.draw();
      }

      g.popMatrix();
    }
    g.popViewMatrix();
    g.popViewport();
  }

  g.scissorTest(last_enable_scissor_test); // put back scissoring
  g.popProjMatrix();
  g.popViewport();
  g.popFramebuffer();
}

void DataDisplay::drawPerspective(Graphics &g) {
  if (mAligned4fData.size() == 0) {
    return; // No data has been loaded
  }

  g.lens().far(2000);

  {
    perspectivePickable.pushMatrix(g);

    g.blendAdd();
    g.depthTesting(false);

    g.meshColor();

    if (mShowAxes) {
      g.pushMatrix();
      g.scale(atomrender.dataBoundary.dim.x * 0.5);
      g.draw(axis);
      g.popMatrix();
    }

    drawTrajectories(g);

    {
      g.pushMatrix();
      atomrender.applyTransformations(g);
      {
        g.pushMatrix();
        g.polygonLine();
        g.scale(atomrender.dataBoundary.dim);
        g.color(0.3f, 0.3f, 0.9f, 0.6f);
        g.draw(boxMesh);
        g.popMatrix();
      }
      atomrender.onProcess(g);

      g.popMatrix();
    }

    {
      g.pushMatrix();
      atomrender.applyTransformations(g);
      // draw percolation data
      g.blendTrans();
      g.depthTesting(true);
      //      g.meshColor();
      auto size = mDatasetManager.mPercolationTypes.getElements().size();
      auto enabledPerc = mDatasetManager.mPercolationTypes.get();
      auto previousScale = atomrender.mAtomMarkerSize.get();

      // Temporarily override scale to atom render
      atomrender.mAtomMarkerSize.setNoCalls(previousScale *
                                            mPercoMarkerScale.get());
      for (size_t i = 0; i < size; i++) {
        if (enabledPerc & ((uint64_t)1 << i)) {
          atomrender.draw(g, mPercolationData[i], mAlignedPercolation4fData[i]);
        }
      }
      // Put back previous scale to atom render
      atomrender.mAtomMarkerSize.setNoCalls(previousScale);

      g.popMatrix();
    }

    g.depthTesting(true);
    g.blending(true);

    if (mDisplaySlicing.get() == 1.0) {
      g.pushMatrix();

      //      atomrender.applyTransformations(g);

      g.translate(atomrender.mSlicingPlaneCorner.get());
      g.rotate(atomrender.mSlicingPlaneQuat);
      g.scale(atomrender.mSlicingPlaneSize, atomrender.mSlicingPlaneSize,
              atomrender.mSlicingPlaneThickness.get());

      //      g.scale(perspectivePickable.scale);

      g.polygonLine();
      g.color(0.8f, 0.8f, 1.0f, 0.9f);
      g.draw(boxMesh);
      g.polygonFill();
      g.popMatrix();
    }

    if (mShowMarkers) { // Draw change markers
      g.pushMatrix();

      g.depthTesting(true);
      g.blending(true);
      //  g.scale(perspectivePickable.scale);
      g.polygonLine();
      for (auto added : atomAdded) {
        for (auto addedPos : added.second) {
          g.pushMatrix();
          g.translate(addedPos.x, addedPos.y, addedPos.z);
          g.color(mMarkerColor);
          g.scale(mMarkerScale);
          g.draw(mMarker);
          g.popMatrix();
        }
      }
      for (auto removed : atomRemoved) {
        g.pushMatrix();
        g.translate(removed.x, removed.y, removed.z);
        g.color(0.4, 0.4, 0.4);
        g.scale(1.7);
        g.draw(mMarker);
        g.popMatrix();
      }

      g.popMatrix();
      { // Draw selected atom
        g.depthTesting(true);
        g.blending(true);
        //  g.scale(perspectivePickable.scale);
        g.polygonLine();
        g.pushMatrix();
        g.translate(selectedPosition);
        g.color(0.2f, 1.f, 0.2f);
        g.scale(2.8);
        g.draw(mMarker);
        g.popMatrix();
      }
    }

    //    perspectivePickable.drawBB(g);
    perspectivePickable.popMatrix(g); // pickable
  }

  g.polygonFill();
  g.depthTesting(false);
  g.blendAdd();
  { // Draw label
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
  }

  //  g.popMatrix(); // pop pickable matrix
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
    g.blendAdd();
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

void DataDisplay::drawGraphPickable(Graphics &g, PickableBB *graphPickable,
                                    Texture *graphTexture, bool billboarding,
                                    bool drawLabels, bool smallLabel) {
  graphPickable->pushMatrix(g);

  if (billboarding) {
    Vec3f forward = graphPickable->pose.get().pos();
    forward.normalize();
    Quatf rot = Quatf::getBillboardRotation(-forward, Vec3f{0.0, 1.0f, 0.0f});
    g.rotate(rot);
  }

  g.depthTesting(true);
  g.quad(*graphTexture, -1.2f, 0.9f, 2.4f, -1.8f);
  // Draw label

  if (drawLabels) {
    g.depthTesting(false);
    g.blendAdd();
    g.translate(-1.0, -1.0, 0.0);
    if (smallLabel) {
      g.scale(0.1f);
    } else {
      g.scale(0.2f);
    }

    // Draw label
    //    Mesh mesh;
    //    mLabelFont.write(
    //      mesh, (mDatasetManager.currentDataset() + " " + mParamText).c_str(),
    //      1.0);
    //    g.texture();
    //    mLabelFont.tex.bind();
    //    g.draw(mesh);
    //    mLabelFont.tex.unbind();
  }
  graphPickable->popMatrix(g);

  graphPickable->drawBB(g);
}

void DataDisplay::drawGraph(Graphics &g) {

  // Update graph textures if new data available
  if (imageDiskBuffer.newDataAvailable()) {
    auto imageBuffer = imageDiskBuffer.get();
    mGraphTexture.resize(imageBuffer->width(), imageBuffer->height());
    mGraphTexture.submit(imageBuffer->pixels(), GL_RGBA, GL_UNSIGNED_BYTE);
    mGraphTexture.filter(Texture::LINEAR_MIPMAP_LINEAR);
  }
  {
    auto *graphTextPtr = mGraphTextures;
    for (auto &imageBuffer : imageDiskBuffers) {
      if (imageBuffer->newDataAvailable()) {
        auto buf = imageBuffer->get();
        graphTextPtr->resize(buf->width(), buf->height());
        graphTextPtr->submit(buf->pixels(), GL_RGBA, GL_UNSIGNED_BYTE);
        graphTextPtr->filter(Texture::LINEAR_MIPMAP_LINEAR);
      }
      graphTextPtr++;
    }

    drawGraphPickable(g, &graphPickable, &mGraphTexture,
                      mBillboarding.get() == 1.0, mDrawLabels.get() == 1.0,
                      mSmallLabel.get());
  }
  {
    auto *graphTextPtr = mGraphTextures;
    auto showGraphIt = mShowGraphs.begin();
    graphTextPtr = mGraphTextures;
    for (auto &pickable : graphPickables) {
      if ((*showGraphIt)->get() == 1.0) {
        drawGraphPickable(g, pickable.get(), graphTextPtr,
                          mBillboarding.get() == 1.0, mDrawLabels.get() == 1.0,
                          mSmallLabel.get());
      }
      showGraphIt++;
      graphTextPtr++;
    }
  }
}
