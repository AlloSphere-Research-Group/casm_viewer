
#include <memory>

#ifdef AL_WINDOWS
//#define NOMINMAX
#undef far
#endif

#include "al/app/al_DistributedApp.hpp"
#include "al/graphics/al_Image.hpp"
#include "al/sphere/al_SphereUtils.hpp"
#include "al/ui/al_FileSelector.hpp"
#include "al/ui/al_HtmlInterfaceServer.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/ui/al_ParameterGUI.hpp"
#include "al/ui/al_Pickable.hpp"
#include "al/ui/al_PresetHandler.hpp"
#include "al/ui/al_SequenceRecorder.hpp"

#include "al_ext/openvr/al_OpenVRDomain.hpp"

#include "tinc/DeferredComputation.hpp"
#include "tinc/PeriodicTask.hpp"
#include "tinc/TincClient.hpp"
#include "tinc/TincServer.hpp"
#include "tinc/vis/GUI.hpp"

#undef AL_BUILD_MPI

#include "datadisplay.hpp"

#ifdef AL_WINDOWS
#undef rank
#undef CIEXYZ
#endif

using namespace al;
using namespace std;

struct State {
  float zoom = 0.0;
};

// ---------------------------------------------------

struct ObjectTransformHandler : WindowEventHandler {
  float scale = 1;
  float mSpeed = 0.1f;
  Quatf quat{Quatf::identity()};
  Vec3f spin;
  Vec3f home{0.0f, 0.0f, 4.0f};
  Vec3f pos{0.0f, 0.0f, 0.0f};
  Vec3f vel{0.0f, 0.0f, 0.0f};

  void setHome(Vec3f newHome) {
    // TODO home should include camera pose too
    home = newHome;
    //   base_scale = newScale;
  }

  void reset() {
    quat = Quatf::identity();
    pos.set(home);
  }

  void step() {
    quat = Quatf().fromEuler(spin) * quat;
    pos += vel;
  }

  Matrix4f mat() {
    Matrix4f rot_mat;
    quat.toMatrix(rot_mat.elems());
    return Matrix4f::translation(pos) * rot_mat;
  }
};

class ScriptRunner {
public:
  void init() {
    kmcPercoTools.setScriptName("..\\python\\tinc_kmc_paths_perco.py");
    asyncWrapper.setProcessor(&kmcPercoTools);
    kmcPercoTools.enableJsonConfig(false);
  }

  void start() {
    auto curPath = al::File::currentPath();
    std::cout << curPath << std::endl;
    kmcPercoTools.setCommand(
        "start cmd.exe /K %userprofile%\\anaconda3\\Scripts\\conda.exe run "
        "--cwd " +
        curPath + " -n base python");
    asyncWrapper.process();
  }

  void stop() {}

  ProcessorScript kmcPercoTools{"kmcPercoTools"};
  ProcessorAsyncWrapper asyncWrapper{"asyncWrapper"};
};

class MyApp : public DistributedAppWithState<State> {
public:
  // Parameters and triggers for interaction commands
  Trigger resetView{"resetView"};
  Trigger stepXpos{"stepXpos"};
  Trigger stepXneg{"stepXneg"};
  Trigger stepYpos{"stepYpos"};
  Trigger stepYneg{"stepYneg"};
  Trigger stepZpos{"stepZpos"};
  Trigger stepZneg{"stepZneg"};
  Trigger mJumpLayerPos{"jumpLayerPos"};
  Trigger mJumpLayerNeg{"jumpLayerNeg"};
  Trigger stepTempPos{"stepTempPos"};
  Trigger stepTempNeg{"stepTempNeg"};
  Trigger stepChempotPos{"stepChempotPos"};
  Trigger stepChempotNeg{"stepChempotNeg"};
  Trigger stepChempot2Pos{"stepChempot2Pos"};
  Trigger stepChempot2Neg{"stepChempot2Neg"};
  Trigger stepTimePos{"stepTimePos"};
  Trigger stepTimeNeg{"stepTimeNeg"};
  Parameter pitchAngleStep{"PitchAngleStep", "AngleControl", 15, 0.0, 60.0};
  Trigger stepPitchAnglePos{"stepPitchAnglePos", "AngleControl"};
  Trigger stepPitchAngleNeg{"stepPitchAngleNeg", "AngleControl"};
  Parameter rollAngleStep{"RollAngleStep", "AngleControl", 15, 0.0, 60.0};
  Trigger stepRollAnglePos{"stepRollAnglePos", "AngleControl"};
  Trigger stepRollAngleNeg{"stepRollAngleNeg", "AngleControl"};
  Trigger ResetSlicing{"ResetSlicing"};
  Trigger mSaveGraphics{"SaveScreenshot"};
  Trigger mAlignTemperatures{"AlignTemperatures"};
  Trigger mAlignChempots{"AlignChempots"};

  ParameterBool mDebugScripts{"DebugScripts"};
  ParameterBool mComputeCache{"ComputeCache"};

  ParameterBool mAutoAdvance{"autoAdvance"};
  Parameter mAutoAdvanceFreq{"autoAdvanceFreq", "", 5, 0.25, 10.0};

  // File selection
  ParameterString mDataset{"dataset"};
  ParameterMenu mRecentDatasets{"recentDatasets"};
  ParameterString mWindowTitle{"windowTitle"};

  // 3D nav parameters
  Parameter Z{"Z", "", 2.8f, -15, 15};
  Parameter X{"X", "", 0.0, -2, 2};

  // Settings
  ParameterColor backgroundColor{"background", "",
                                 Color(0.12f, 0.12f, 0.12f, 1.0f)};
  ParameterColor sliceBackground{"sliceBackground", "",
                                 Color(0.0f, 0.0f, 0.0f, 1.0f)};
  ParameterString font{"font"};
  Parameter fontSize{"fontSize", "", 32, 2, 128};
  ParameterString pythonScriptsPath{"pythonScriptsPath"};
  ParameterString pythonBinary{"pythonBinary"};

  ParameterBool recallPositions{"recallPositions", "", 1.0};

  std::unique_ptr<PresetHandler> presetHandler;

  std::unique_ptr<PresetHandler> positionPresets;

  std::unique_ptr<PresetSequencer> sequencer;
  std::unique_ptr<SequenceRecorder> recorder;
  std::unique_ptr<PresetServer> presetServer;

  BundleGUIManager vdvBundle;
  FileSelector *mDatasetSelector;
  std::string mPreviousBrowseDir; // Previously selected dir

  std::mutex mScreenshotMutex;
  std::string mScreenshotPrefix;

  PeriodicTask mParameterPlayback;

  // TINC computation chains
  ProcessorGraph initRootProcessorGraph{"PrepareDataset"};

  ProcessorScript parameterSpaceProcessor{"ParameterSpaceProcessor"};
  ProcessorScript shellSiteFileAnalyzer{"ShellSiteFileAnalyzer"};
  ProcessorScript transfmatExtractor{"TransfmatExtractor"};
  ProcessorScript templateGen{"TemplateGenerator"};

  //  ScriptRunner percoTools;

  TincServer tincServer;
  TincClient tincClient;
  TincProtocol *tinc;

  // --------
  std::vector<DataDisplay *> dataDisplays;

  ObjectTransformHandler object_transform;
  bool showGui{true};

  PersistentConfig config;

#ifdef AL_EXT_OPENVR
  std::shared_ptr<OpenVRDomain> openVRDomain;
  VAOMesh mControllerMesh;
#endif

  // ------------- Application callbacks

  void onInit() override {
    presetHandler = std::make_unique<PresetHandler>();
    positionPresets = std::make_unique<PresetHandler>();
    sequencer = std::make_unique<PresetSequencer>();
    recorder = std::make_unique<SequenceRecorder>();
    if (isPrimary()) {
      presetServer = std::make_unique<PresetServer>();
    }

    // Create and configure displays
    int numDisplays = 2;

    for (int i = 0; i < numDisplays; i++) {
      dataDisplays.emplace_back(new DataDisplay);
      dataDisplays.back()->mDatasetManager.mRunProcessors = rank == 0;
      dataDisplays.back()->init();
      if (dataRoot.size() > 0) {
        dataDisplays.back()->mDatasetManager.mGlobalRoot = dataRoot;
      }

      DataDisplay *newDisplay = dataDisplays.back();
      newDisplay->atomrender.registerUpdateCallback([newDisplay](bool ok) {
        auto b = newDisplay->atomrender.getDataBoundaries();
        newDisplay->perspectivePickable.bb.set(
            Vec3f(b.min.x, b.min.y, b.min.z), Vec3f(b.max.x, b.max.y, b.max.z));
        newDisplay->slicePickable.bb.set(
            Vec3f(b.min.x, b.min.y, b.min.z),
            Vec3f(b.max.x, b.max.y, (b.max.z - b.min.z) * 0.25f));
        newDisplay->resetSlicing();
      });
    }

    char *szBuff;
    szBuff = std::getenv("USERPROFILE");
    if (szBuff) {
      mPreviousBrowseDir = szBuff;
    } else {
      // mPreviousBrowseDir = "~/";
    }

    // Initialize default view.
    int counter = 0;
    for (auto *disp : dataDisplays) {
      disp->graphPickable.pose.setPos(Vec3d(-2.0 + (counter * 0.5), 0.8, -.75));
      disp->graphPickable.scale = 0.2f;
      disp->parallelPickable.pose.setPos(
          Vec3d(-0.68 + (counter * 0.5), 1.1, -1.5));
      disp->parallelPickable.scale = 0.4f;
      disp->perspectivePickable.pose.setPos(
          Vec3d(0.55 + (counter * 0.5), 0.35, -0.7));
      disp->perspectivePickable.scale = 0.02f;

      disp->mShowPerspective = true;
      disp->mBillboarding = false;
      disp->mSmallLabel = true;
      disp->mVisible = true;
    }

    for (size_t i = 1; i < dataDisplays.size(); i++) {
      dataDisplays[i]->mVisible = false;
    }

    // ----------
    loadConfiguration();
    initializeComputation();
    registerParameterCallbacks();

    setupParameterServer();
    readElementsIni();
    for (auto *display : dataDisplays) {
      *presetHandler << display->bundle;
      *presetHandler << display->graphPickable.bundle
                     << display->parallelPickable.bundle
                     << display->perspectivePickable.bundle;

      *positionPresets << display->graphPickable.bundle;
      *positionPresets << display->graphPickable.bundle
                       << display->parallelPickable.bundle
                       << display->perspectivePickable.bundle;
      for (auto &graph : display->graphPickables) {
        *positionPresets << graph->bundle;
      }
      vdvBundle << display->bundle;
      presetHandler->skipParameter(
          display->mDatasetManager.mCurrentDataset.getFullAddress());
    }

    // Add parameters that are not part of the bundle
    *presetHandler << Z << X; //

    if (isPrimary()) {
      if (presetServer) {
        *presetServer << *presetHandler;
      }
      *sequencer << *presetHandler;
      *recorder << *presetHandler;
      parameterServer().notifyAll();
      backgroundColor.set(
          backgroundColor.get()); // Set current value to update on renderers
      sliceBackground.set(
          sliceBackground.get()); // Set current value to update on renderers
    }

    fps(60);
    title("Casm Viewer");
    dimensions(1200, 800);

    backgroundColor.processChange();
    sliceBackground.processChange();
    pythonScriptsPath.processChange();
    pythonBinary.processChange();
  }

  virtual void onCreate() override {
    imguiInit();

    // disable nav control mouse drag to look
    navControl().useMouse(false);

    font.processChange(); // Loads font
    // use object control for model matrix
    defaultWindow().append(object_transform);
    object_transform.setHome(Vec3f(0, 0, 0));
    object_transform.reset();

    if (isPrimary()) {
      nav().pos(0, 1.5, 5);
      nav().setHome();
      navControl().useMouse(false);
    } else if (!hasCapability(CAP_OMNIRENDERING)) {
      // i.e. desktop replica
      nav().pos(0, 1.5, 5);
    }

    // For distributed running in allo infrastructure
    if (!isPrimary()) {
      std::string simulatorIp = "atari.1g";
      parameterServer().requestAllParameters(simulatorIp, 9010);
    }

    if (al_get_hostname() ==
        "moxi") { // Fullscreen to two monitors on MOXI machine
      std::cout << "On MOXI! --------------------------------------"
                << std::endl;
      int width, height;
      sphere::getFullscreenDimension(&width, &height);
      std::cout << "Setting fullscreen dimensions " << width << "," << height
                << std::endl;
      if (width != 0 && height != 0) {
        decorated(false);
        dimensions(0, 0, width, height);
        //        stereo(true);
        cursorHide(false);
      }
    }

#ifdef AL_EXT_OPENVR
    openVRDomain = OpenVRDomain::enableVR(this);
    openVRDomain->setDrawFunction(
        std::bind(&MyApp::drawVR, this, std::placeholders::_1));

    addWireBox(mControllerMesh, 0.1f);
    mControllerMesh.scale(0.1f, 0.1f, 1);
    mControllerMesh.update();
#endif

    // Information on TINC component
    tinc->print();
  }

  virtual void onAnimate(double /*dt*/) override {
    object_transform.step();

    if (mWindowTitle.processChange()) {
      title(mWindowTitle.get());
    }
    if (isPrimary()) {
      if (showGui) {
        prepareGui();
      }
      navControl().active(!ParameterGUI::usingInput());
    }

    // fbo operations should be done outside onDraw so it does not mess with
    // omni drawing
    for (auto *display : dataDisplays) {
      display->prepare(graphics());
    }
#ifdef AL_EXT_OPENVR
    // Update traking and controller data;
    //    mOpenVR.update();

    auto l = openVRDomain->mOpenVR.LeftController;
    auto r = openVRDomain->mOpenVR.RightController;
    auto ray = r.ray();

    if (l.touchpadPress()) {
      if (fabs(l.touchPos.x) > 0.3f) {
        if (l.touchPos.x > 0) {
          dataDisplays[vdvBundle.currentBundle()]->nextChempot();
        } else {
          dataDisplays[vdvBundle.currentBundle()]->previousChempot();
        }
      }
      if (fabs(l.touchPos.y) > 0.3f) {
        if (l.touchPos.y > 0) {
          dataDisplays[vdvBundle.currentBundle()]->nextTemp();
        } else {
          dataDisplays[vdvBundle.currentBundle()]->previousTemp();
        }
      }
    }

    if (r.touchpadDown()) {
      if (fabs(r.touchPos.x) > 0.4f) {
        dataDisplays[vdvBundle.currentBundle()]->mSlicingPlaneThickness =
            dataDisplays[vdvBundle.currentBundle()]->mSlicingPlaneThickness +
            r.touchPos.x * 0.2;
      }
      if (fabs(r.touchPos.y) > 0.4f) {
        dataDisplays[vdvBundle.currentBundle()]->mPickableManager.event(
            PickEvent(Pick, ray));
        dataDisplays[vdvBundle.currentBundle()]->mPickableManager.event(
            PickEvent(Scale, r.touchPos.y));
        // dataDisplays[vdvBundle.currentBundle()]->perspectivePickable.scale =
        // dataDisplays[vdvBundle.currentBundle()]->perspectivePickable.scale +
        // r.touchPos.y*0.0001;
      }
    }
    if (r.touchpadPress())
      dataDisplays[vdvBundle.currentBundle()]
          ->perspectivePickable.testChildren = true;
    if (r.touchpadRelease())
      dataDisplays[vdvBundle.currentBundle()]
          ->perspectivePickable.testChildren = false;

    for (auto *display : dataDisplays) {
      if (display->mVisible != 0.0f) {
        Pose p = r.pose();
        if (r.triggerPress()) {
          display->mPickableManager.event(PickEvent(Pick, ray));
        } else if (r.gripPress()) {
          display->mPickableManager.event(PickEvent(Pick, ray));
          display->mPickableManager.event(PickEvent(PickPose, p));
        } else if (r.triggerDown()) {
          display->mPickableManager.event(PickEvent(Drag, ray, r.vel));
        } else if (r.gripDown()) {
          display->mPickableManager.event(PickEvent(RotatePose, p));
        } else if (r.triggerRelease() || r.gripRelease()) {
          display->mPickableManager.event(PickEvent(Unpick, ray));
        } else {
          display->mPickableManager.event(PickEvent(al::Point, ray));
        }
      }
    }
#endif
  }

  virtual void onDraw(Graphics &g) override {
    g.clear(backgroundColor);
    drawScene(g);
    if (isPrimary()) {
      processScreenshot();

      if (showGui) {
        imguiDraw();
      }
    }
  }

  virtual bool onKeyDown(const Keyboard &k) override {
    if (isPrimary()) {
      if (!ParameterGUI::usingKeyboard()) {
        for (auto *display : dataDisplays) {
          display->mPickableManager.onKeyDown(k);
          display->perspectivePickable.testChildren = k.shift();
        }

        int index;
        switch (k.key()) {
        case ' ':
          showGui = !showGui;
          break;
        case 'b':
          dataDisplays[vdvBundle.currentBundle()]->mPerspectiveRotY =
              dataDisplays[vdvBundle.currentBundle()]->mPerspectiveRotY + 5;
          break;
        case 'v':
          dataDisplays[vdvBundle.currentBundle()]->mPerspectiveRotY =
              dataDisplays[vdvBundle.currentBundle()]->mPerspectiveRotY - 5;
          break;
        case '-':
          mJumpLayerNeg.trigger(); // This will trigger a change
          break;
        case '=':
          mJumpLayerPos.trigger(); // This will trigger a change
          break;
        case 'o':
          stepTempNeg.trigger(); // This will trigger a change
          break;
        case 'p':
          stepTempPos.trigger(); // This will trigger a change
          break;
        case 'l':
          stepChempotNeg.trigger(); // This will trigger a change
          break;
        case ';':
          stepChempotPos.trigger(); // This will trigger a change
          break;
        case '.':
          stepChempot2Neg.trigger(); // This will trigger a change
          stepTimeNeg.trigger();
          break;
        case '/':
          stepChempot2Pos.trigger(); // This will trigger a change
          stepTimePos.trigger();
          break;
        case 'q':
          index = presetHandler->getCurrentPresetIndex() - 1;
          if (index >= 0) {
            presetHandler->recallPreset(index);
          }
          break;
        case 'w':
          index = presetHandler->getCurrentPresetIndex() + 1;
          presetHandler->recallPreset(index);
          break;

        case 'z':
          mAlignTemperatures.trigger();
          break;
        default:
          break;
        }
      }
    }
    return true;
  }

  virtual bool onKeyUp(Keyboard const &k) override {
    for (auto *display : dataDisplays) {
      display->mPickableManager.onKeyUp(k);
      display->perspectivePickable.testChildren = k.shift();
    }
    return true;
  }

  virtual void onMessage(osc::Message &m) override {
    m.print();
    if ("/point" == m.addressPattern()) {
      float ox, oy, oz, dx, dy, dz;
      int id;
      m >> id >> ox >> oy >> oz >> dx >> dy >> dz;
      Rayd r = Rayd(Vec3f(ox, oy, oz), Vec3f(dx, dy, dz));
      // r = rayTransformAllosphere(r);
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.event(PickEvent(al::Point, r));
        }
      }
    } else if ("/pick" == m.addressPattern()) {
      float ox, oy, oz, dx, dy, dz;
      int id;
      int button;
      m >> id >> button >> ox >> oy >> oz >> dx >> dy >> dz;
      Rayd r = Rayd(Vec3f(ox, oy, oz), Vec3f(dx, dy, dz));
      // r = rayTransformAllosphere(r);
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.event(PickEvent(Pick, r));
        }
      }
    } else if ("/drag" == m.addressPattern()) {
      float ox, oy, oz, dx, dy, dz, x, y, z;
      int id;
      int button;
      m >> id >> button >> ox >> oy >> oz >> dx >> dy >> dz >> x >> y >> z;
      Rayd r = Rayd(Vec3f(ox, oy, oz), Vec3f(dx, dy, dz));
      // r = rayTransformAllosphere(r);
      Vec3f v = Vec3f(x, y, z);
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.event(
              PickEvent(Drag, r, v)); // id, button, v);
        }
      }
    } else if ("/rotate" == m.addressPattern()) {
      float ox, oy, oz, dx, dy, dz, x, y, z;
      int id;
      int button;
      m >> id >> button >> ox >> oy >> oz >> dx >> dy >> dz >> x >> y >> z;
      Rayd r = Rayd(Vec3f(ox, oy, oz), Vec3f(dx, dy, dz));
      // r = rayTransformAllosphere(r);
      //      Vec3f v = Vec3f(x,y,z);
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.event(PickEvent(RotateRay, r));
        }
      }
    } else if ("/scale" == m.addressPattern()) {
      float ox, oy, oz, dx, dy, dz, x, y, z;
      int id;
      int button;
      m >> id >> button >> ox >> oy >> oz >> dx >> dy >> dz >> x >> y >> z;
      Rayd r = Rayd(Vec3f(ox, oy, oz), Vec3f(dx, dy, dz));
      // r = rayTransformAllosphere(r);
      //      Vec3f v = Vec3f(x,y,z);
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.event(PickEvent(Scale, r, y));
        }
      }
    } else if ("/unpick" == m.addressPattern()) {
      float ox, oy, oz, dx, dy, dz;
      int id;
      int button;
      m >> id >> button >> ox >> oy >> oz >> dx >> dy >> dz;
      Rayd r = Rayd(Vec3f(ox, oy, oz), Vec3f(dx, dy, dz));
      // r = rayTransformAllosphere(r);

      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.event(
              PickEvent(Unpick, r)); //, id, button);
        }
      }
    }
  }

  bool onMouseMove(const Mouse &m) override {
    for (auto *display : dataDisplays) {
      if (display->mVisible != 0.0f) {
        if (!ParameterGUI::usingInput()) {
          display->mPickableManager.onMouseMove(graphics(), m, width(),
                                                height());
        } else {
          display->mPickableManager.unhighlightAll();
        }
      }
    }
    return true;
  }

  bool onMouseDown(const Mouse &m) override {
    for (auto *display : dataDisplays) {
      if (display->mVisible != 0.0f) {
        if (!ParameterGUI::usingInput()) {
          display->mPickableManager.onMouseDown(graphics(), m, width(),
                                                height());
        } else {
          display->mPickableManager.onMouseUp(graphics(), m, width(), height());
        }
      }
    }
    return true;
  }

  bool onMouseDrag(const Mouse &m) override {
    if (isPrimary()) {
      if (showGui && ParameterGUI::usingInput())
        return true;
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.onMouseDrag(graphics(), m, width(),
                                                height());
        }
      }
    }
    return true;
  }
  bool onMouseUp(const Mouse &m) override {
    if (isPrimary()) {
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.onMouseUp(graphics(), m, width(), height());
        }
      }
    }
    return true;
  }

  void onSound(AudioIOData &io) override {
    //    static double phase {0};
    //    double phaseIncrement = 440.0 / io.framesPerSecond();

    for (auto *display : dataDisplays) {
      display->synth.render(io);
    }
    //    while(io()){
    //      phase += phaseIncrement;
    //      if(phase > 1) phase -= 1;
    //      float out = 0.05 * cos(phase * M_2PI);
    //      io.out(0) = out;
    //      io.out(1) = out;
    //    }
  }

  virtual void onExit() override {
    if (isPrimary()) {
      storeConfiguration();
      tincServer.stop();
    } else {
      tincClient.stop();
    }
    imguiShutdown();
    //#ifdef AL_EXT_OPENVR
    //    mOpenVR.close();
    //#endif
  }

  // ---- Drawing functions

  void drawScene(Graphics &g) {
    for (auto *display : dataDisplays) {
      display->draw(g);
    }
  }

#ifdef AL_EXT_OPENVR
  void drawVR(Graphics &g) {
    g.clear(backgroundColor);
    drawScene(g);

    // Draw markers for the controllers
    g.pushMatrix();
    g.translate(openVRDomain->mOpenVR.LeftController.pos);
    g.rotate(openVRDomain->mOpenVR.LeftController.quat);
    g.color(0, 1, 1);
    g.draw(mControllerMesh);
    g.popMatrix();

    // right hand
    g.pushMatrix();
    g.translate(openVRDomain->mOpenVR.RightController.pos);
    g.rotate(openVRDomain->mOpenVR.RightController.quat);
    g.color(1, 0, 1);
    g.draw(mControllerMesh);
    g.popMatrix();

    // draw controller rays
    gl::blendTrans();
    auto r1 = openVRDomain->mOpenVR.RightController.ray();
    //    auto r2 = mOpenVR.LeftController.ray();
    Mesh rays;
    // TODO can we make the ray fade out in the distance?
    rays.primitive(Mesh::LINES);
    rays.colors().push_back({1, 0, 1, 1});
    rays.vertices().push_back(r1.o);
    rays.colors().push_back({1, 0, 1, 1});
    rays.vertices().push_back(r1.o + r1.d * 2.5);
    rays.colors().push_back({1, 0, 1, 0.3f});
    rays.vertices().push_back(r1.o + r1.d * 5);
    // rays.vertices().push_back(r2.o);
    // rays.vertices().push_back(r2.o + r2.d*5);
    g.meshColor();
    g.draw(rays);

    //        // Draw "on board" controls
    //        g.pushMatrix();
    //        g.multModelMatrix(mOpenVR.eyePosLeft);
    //        g.translate(0.0, 0.0, 1.0);
    //        g.color(1,1,1);
    //        g.draw(mController);
    //        g.popMatrix();
  }
#endif

  void loadDataset(std::string path, size_t index) {
    assert(index < dataDisplays.size());
    if (!al::File::isDirectory(path)) {
      std::cerr << "ERROR: Can't find directory " << path << std::endl;
      return;
    }

    parameterSpaceProcessor.setRunningDirectory(path);
    parameterSpaceProcessor.setOutputDirectory(path + "/cached_output");

    transfmatExtractor.setDataDirectory(path);

    templateGen.setDataDirectory(path);

    initRootProcessorGraph.setRunningDirectory(path);
    initRootProcessorGraph.process();

    if (index >= 0) {
      if (dataDisplays[index]->mDatasetManager.mCurrentDataset.get() != path) {
        dataDisplays[index]->mDatasetManager.mCurrentDataset.set(path);
      }
      auto previousDatasets = mRecentDatasets.getElements();
      if (std::find(previousDatasets.begin(), previousDatasets.end(), path) ==
          previousDatasets.end()) {
        previousDatasets.insert(previousDatasets.begin(), path);
        mRecentDatasets.setElements(previousDatasets);
        mRecentDatasets.setNoCalls(0);
      }
    } else {
    }
  }

  void prepareMainGui() {
    if (mDatasetSelector) {
      if (dataRoot.size() > 0) {
        ImGui::Text("Data root: %s", dataRoot.c_str());
      }
      ParameterGUI::draw(&mRecentDatasets);
      //      ImGui::Separator();
      // Directory GUI
      if (mDatasetSelector->isActive()) {
        if (mDatasetSelector->drawFileSelector()) { // Selection has been made
          auto selectedItems = mDatasetSelector->getSelection();
          if (selectedItems.count() > 0) {
            auto selection = selectedItems[0];
            mDataset.set(File::conformPathToOS(selection.filepath()));
            mPreviousBrowseDir = selection.path();
          }
          delete mDatasetSelector;
          mDatasetSelector = nullptr;
        }
      } else {
        delete mDatasetSelector;
        mDatasetSelector = nullptr;
      }

    } else {
      static bool loaded = false;
      if (ImGui::Button("Load dataset")) {
        mDatasetSelector = new FileSelector(dataRoot, File::isDirectory);
        mDatasetSelector->start(mPreviousBrowseDir);
        if (this->dataDisplays[vdvBundle.currentBundle()]
                ->mDatasetManager.mParameterSpace.getDimensions()
                .size() > 0) {
          loaded = true;
        }
      }
      if (!loaded && ImGui::Button("Cancel")) {
        loaded = true;
      }
      if (this->dataDisplays[vdvBundle.currentBundle()]->lastError.get() !=
          "") {
        ImGui::Text("%s", this->dataDisplays[vdvBundle.currentBundle()]
                              ->lastError.get()
                              .c_str());
      }
      if (loaded) {
        ParameterGUI::draw(&mComputeCache);
        if (mAutoAdvance == 0.0) {
          ImGui::Separator();
          vis::drawControls(this->dataDisplays[vdvBundle.currentBundle()]
                                ->mDatasetManager.mParameterSpace);
          ImGui::Separator();
        }
        ParameterGUI::drawParameterMeta(
            &this->dataDisplays[vdvBundle.currentBundle()]
                 ->atomrender.mAtomMarkerSize);
        if (this->dataDisplays[vdvBundle.currentBundle()]
                ->mDatasetManager.mParameterSpace.getDimension("time")) {
          ParameterGUI::drawParameterMeta(&mAutoAdvance);
          if (mAutoAdvance == 0.0) {
            ParameterGUI::drawParameterMeta(&mAutoAdvanceFreq);
          } else {
            std::string text = "Incrementing time with freq: " +
                               std::to_string(mAutoAdvanceFreq);
            ImGui::Text("%s", text.c_str());
          }
        }
        ParameterGUI::drawTrigger(&mSaveGraphics);
        ImGui::SameLine();
        ParameterGUI::drawParameterMeta(&ResetSlicing);
        ParameterGUI::draw(
            &this->dataDisplays[vdvBundle.currentBundle()]->mDisplaySlicing);
        ParameterGUI::draw(
            &this->dataDisplays[vdvBundle.currentBundle()]->mDisplayGraphs);

        ImGui::Text("%s", dataDisplays[vdvBundle.currentBundle()]
                              ->mDatasetManager.metaText.c_str());
        ParameterGUI::drawParameterMeta(
            &this->dataDisplays[vdvBundle.currentBundle()]->currentSelection);
        if (this->dataDisplays[vdvBundle.currentBundle()]
                ->mDatasetManager.mShellSiteTypes.getElements()
                .size() > 0) {
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->mDatasetManager.mShellSiteTypes);
        } else {
          ImGui::Text("No Shell site data found.");
        }

        auto percoSize = this->dataDisplays[vdvBundle.currentBundle()]
                             ->mDatasetManager.mPercolationTypes.getElements()
                             .size();
        if (percoSize > 0) {
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->mPercoMarkerScale);
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->mDatasetManager.mPercolationTypes);

          for (size_t i = 0; i < percoSize; i++) {
            ParameterGUI::drawParameterMeta(
                this->dataDisplays[vdvBundle.currentBundle()]
                    ->mPercoColorList[i]
                    .get());
          }
        } else {
          ImGui::Text("No Percolation data found.");
        }

        if (this->dataDisplays[vdvBundle.currentBundle()]
                    ->mDisplaySlicing.get() == 1.0 &&
            ImGui::CollapsingHeader("Slicing")) {
          ImGui::Indent(20.0);
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->atomrender.mAlpha);
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->atomrender.mSlicingPlaneCorner);
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->atomrender.mSlicingPlaneNormal);
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->atomrender.mSlicingPlaneThickness);
          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->atomrender.mSlicingPlaneSize);

          ParameterGUI::drawParameterMeta(
              &this->dataDisplays[vdvBundle.currentBundle()]
                   ->atomrender.mClippedMultiplier);

          ParameterGUI::drawParameterMeta(&mJumpLayerNeg);
          ImGui::SameLine();
          ParameterGUI::drawParameterMeta(&mJumpLayerPos);
          ImGui::SameLine();
          ParameterGUI::drawParameterMeta(&pitchAngleStep);
          ParameterGUI::drawParameterMeta(&stepPitchAngleNeg);
          ImGui::SameLine();
          ParameterGUI::drawParameterMeta(&stepPitchAnglePos);

          ParameterGUI::drawParameterMeta(&rollAngleStep);
          ParameterGUI::drawParameterMeta(&stepRollAngleNeg);
          ImGui::SameLine();
          ParameterGUI::drawParameterMeta(&stepRollAnglePos);
          ImGui::Unindent(20.0);
        }

        if (ImGui::CollapsingHeader("Appearance")) {
          ImGui::Indent(20.0);
          ParameterGUI::drawParameterMeta(&backgroundColor);
          ParameterGUI::drawParameterMeta(&sliceBackground);
          ParameterGUI::drawParameterMeta(&font);
          ImGui::SameLine();
          ParameterGUI::drawParameterMeta(&fontSize);
          ImGui::Unindent(20.0);
        }

        if (ImGui::CollapsingHeader("Elements")) {
          ImGui::Indent(20.0);
          ParameterGUI::drawBundle(
              &dataDisplays[vdvBundle.currentBundle()]->graphPickable.bundle);
          for (auto &graphPick :
               dataDisplays[vdvBundle.currentBundle()]->graphPickables) {
            ParameterGUI::drawBundle(&graphPick->bundle);
          }
          ParameterGUI::drawBundle(&dataDisplays[vdvBundle.currentBundle()]
                                        ->parallelPickable.bundle);
          ParameterGUI::drawBundle(&dataDisplays[vdvBundle.currentBundle()]
                                        ->perspectivePickable.bundle);
          ParameterGUI::draw(
              &dataDisplays[vdvBundle.currentBundle()]->mShowMarkers);
          ParameterGUI::draw(
              &dataDisplays[vdvBundle.currentBundle()]->mShowAxes);
          ImGui::Unindent(20.0);
        }
        ParameterGUI::drawPresetHandler(positionPresets.get());
      } else {
        ImGui::Text("No dataset loaded");
      }
    }
  }

  void prepareGui() {
    if (parameterSpaceProcessor.isRunning()) {
      imguiBeginFrame();
      ParameterGUI::beginPanel("Loading");
      ImGui::Text("Processing parameter space... Please wait.");
      ParameterGUI::endPanel();
      imguiEndFrame();
      return;
    }

    imguiBeginFrame();
    ImGui::SetNextWindowBgAlpha(0.9f);
    ImGui::PushID("CASMViewer");
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(380, 600), ImGuiCond_FirstUseEver);
    ParameterGUI::beginPanel("CASM Viewer");

    static int selected = 0;
    if (!mDatasetSelector) {
      ImGui::Columns(4, nullptr, true);
      if (ImGui::Selectable("Display", selected == 0)) {
        selected = 0;
      }
      ImGui::NextColumn();
      if (ImGui::Selectable("Data", selected == 1)) {
        selected = 1;
      }
      ImGui::NextColumn();
      if (ImGui::Selectable("Presets", selected == 2)) {
        selected = 2;
      }
      ImGui::NextColumn();
      if (ImGui::Selectable("Help", selected == 3)) {
        selected = 3;
      }
      ImGui::NextColumn();
      ImGui::Columns(1);
      ImGui::Separator();
    }

    if (selected == 0) {
      prepareMainGui();
    } else if (selected == 1) {
      ParameterGUI::drawBundleManager(&vdvBundle);
    } else if (selected == 2) {
      ParameterGUI::drawPresetHandler(presetHandler.get(), 10, 4);
      ParameterGUI::draw(&recallPositions);
      int currentSequencerItem = 0;
      ParameterGUI::drawPresetSequencer(sequencer.get(), currentSequencerItem);
      ParameterGUI::drawSequenceRecorder(recorder.get());
    } else if (selected == 3) {
      static std::string helpText;
      if (helpText.size() == 0) {
        std::ifstream t("../readme.md");
        if (t.is_open()) {
          t.seekg(0, std::ios::end);
          helpText.reserve(t.tellg());
          t.seekg(0, std::ios::beg);

          helpText.assign((std::istreambuf_iterator<char>(t)),
                          std::istreambuf_iterator<char>());
        } else {
          std::cerr << "WARNING: Help text file not found. No help text on GUI"
                    << std::endl;
        }
      }

      char buf[512];
      auto currentPython = pythonBinary.get();
      strncpy(buf, currentPython.c_str(), currentPython.size() + 1);
      if (ImGui::InputText("Python executable", buf, 512)) {
        pythonBinary.set(std::string(buf));
      }
      auto currentScriptsPath = pythonScriptsPath.get();
      strncpy(buf, currentScriptsPath.c_str(), currentScriptsPath.size() + 1);
      if (ImGui::InputText("Python Scripts Path", buf, 512)) {
        pythonScriptsPath.set(std::string(buf));
      }
#ifdef AL_WINDOWS
      if (ImGui::Button("Use anaconda")) {
        pythonBinary.set(
            "%userprofile%\\anaconda3\\Scripts\\conda.exe run -n base python");
      }
#endif
      ParameterGUI::draw(&mDebugScripts);

      ImGui::Separator();
      ImGui::Text("%s", helpText.c_str());
    }

    ParameterGUI::endPanel();
    ImGui::PopID();

    auto params = tinc->dimensions();
    ImGui::SetNextWindowPos(ImVec2(400, 10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(300, 400), ImGuiCond_FirstUseEver);
    ImGui::Begin("TINC controls");
    //    if (ImGui::Button("Start KMC analysis script")) {
    //      percoTools.start();
    //    }
    if (isPrimary()) {
      vis::drawTincServerInfo(tincServer, true);
    }
    for (auto *param : params) {
      if (param->getGroup() == "casm") {
        vis::drawControl(param);
      }
    }
    ImGui::End();

    imguiEndFrame();
  }

  void setupParameterServer() {
    // Register parameters with parameter server
    //    for (auto *display : dataDisplays) {
    //    }
    // TODO Allow control of multiple displays remotely
    parameterServer() << dataDisplays[0]->bundle;
    parameterServer() << dataDisplays[0]->graphPickable.bundle;
    for (auto &gp : dataDisplays[0]->graphPickables) {
      parameterServer() << gp->bundle;
    }

    parameterServer() << dataDisplays[0]->parallelPickable.bundle;
    parameterServer() << dataDisplays[0]->perspectivePickable.bundle;
    parameterServer() << dataDisplays[0]->slicePickable.bundle;
    parameterServer() << resetView << stepXpos << stepXneg << stepYpos
                      << stepYneg << stepZpos << stepZneg;
    parameterServer() /*<< CalculateSlicing*/ << ResetSlicing;
    parameterServer() << backgroundColor << sliceBackground;
    parameterServer() << font << fontSize;
    // These parameters are command triggers
    parameterServer() << stepTempPos << stepTempNeg << stepChempotPos
                      << stepChempotNeg << stepChempot2Pos << stepChempot2Neg
                      << mJumpLayerNeg << mJumpLayerPos;

    parameterServer().verbose(true);
  }

  // File IO/ Persistence
  void loadConfiguration() {
    TomlLoader configLoader2;
    configLoader2.setFile("casm_viewer.toml");

    // There is no support for string lists yet in PersistentConfig, so do
    // manually
    auto paths = configLoader2.getVector<string>("openedDatasets");
    mRecentDatasets.setElements(paths);

// Set default values
#ifdef AL_WINDOWS
    font.setNoCalls(al::Font::defaultFont());
#elif defined(AL_OSX)
    font.setNoCalls(al::Font::defaultFont());
#else
    font.setNoCalls(string("/usr/share/fonts/truetype/freefont/FreeMono.ttf"));
#endif
    pythonScriptsPath.setNoCalls(string("./python"));
    pythonBinary.setNoCalls(string("python3"));

    ///
    config.setAppName("casm_viewer");
    config.registerParameter(fontSize);
    config.registerString("previousBrowseDir", &mPreviousBrowseDir);
    config.registerParameter(backgroundColor);
    config.registerParameter(sliceBackground);

    config.registerParameter(pythonScriptsPath);
    config.registerParameter(pythonBinary);
    config.registerParameter(font);

    auto *display = dataDisplays[0];
    for (auto color : display->mPercoColorList) {
      config.registerParameter(*color);
    }
    config.read();
  }

  void storeConfiguration() {
    auto currentDir = File::currentPath();
    std::cout << "Writing config in " << currentDir << std::endl;

    TomlLoader configLoader2;
    configLoader2.setFile("casm_viewer.toml");

    configLoader2.setVector<std::string>("openedDatasets",
                                         mRecentDatasets.getElements());
    configLoader2.writeFile();

    config.write();
  }

  void readElementsIni() {
    // Read elements.ini file
    std::ifstream f("elements.ini");
    if (f.is_open()) {
      std::string line;
      while (getline(f, line)) {
        std::stringstream ss(line);
        std::string field;
        while (field.size() == 0) {
          std::getline(ss, field, ' '); // atom number
        }
        field = "";
        while (field.size() == 0) {
          std::getline(ss, field, ' ');
        }
        string atomName = field;
        field = "";
        while (field.size() == 0) {
          std::getline(ss, field, ' ');
        }
        float radius = stof(field);
        field = "";
        while (field.size() == 0) {
          std::getline(ss, field, ' '); // Ignore VdW radius
        }
        field = "";
        while (field.size() == 0) {
          std::getline(ss, field, ' '); // Ignore Ionic radius
        }
        field = "";
        while (field.size() == 0) {
          std::getline(ss, field, ' ');
        }
        float r = stof(field);
        field = "";
        while (field.size() == 0) {
          std::getline(ss, field, ' ');
        }
        float g = stof(field);
        field = "";
        while (field.size() == 0) {
          std::getline(ss, field, ' ');
        }
        float b = stof(field);
        for (auto *display : dataDisplays) {
          display->elementData[atomName] = {radius, Color(r, g, b)};
        }
      }
    } else {
      std::cout << name() << "Could not open elements.ini - Using defaults."
                << std::endl;
    }
  }

  void processScreenshot() {
    std::unique_lock<std::mutex> lk(mScreenshotMutex);
    if (mScreenshotPrefix.size() > 0) { // Take screenshot
      std::vector<unsigned char> mPixels;
      mPixels.resize(width() * height() * 3);
      unsigned char *pixs = mPixels.data();
      // FIXME this can be done more efficiently
      glReadPixels(1, 1, width(), height(), GL_RGB, GL_UNSIGNED_BYTE, pixs);
      std::string dumpDirectory = File::conformPathToOS(
          dataDisplays[vdvBundle.currentBundle()]->mDatasetManager.mGlobalRoot +
          dataDisplays[vdvBundle.currentBundle()]
              ->mDatasetManager.mCurrentDataset.get() +
          "/graphics/");
      if (!al::File::isDirectory(dataDisplays[vdvBundle.currentBundle()]
                                     ->mDatasetManager.mGlobalRoot +
                                 dataDisplays[vdvBundle.currentBundle()]
                                     ->mDatasetManager.mCurrentDataset.get() +
                                 "/graphics/")) {
        al::Dir::make(dataDisplays[vdvBundle.currentBundle()]
                          ->mDatasetManager.mGlobalRoot +
                      dataDisplays[vdvBundle.currentBundle()]
                          ->mDatasetManager.mCurrentDataset.get() +
                      "/graphics/");
      }
      std::string imagePath = dumpDirectory + mScreenshotPrefix + "_screen.png";

      Image::saveImage(imagePath, pixs, width(), height(), true, 3);

      mScreenshotPrefix = "";
    }
  }

  void initializeComputation() {
    parameterSpaceProcessor.setCommand(pythonBinary);
    parameterSpaceProcessor.setScriptName(pythonScriptsPath.get() +
                                          "/analyze_parameter_space.py");
    parameterSpaceProcessor.setOutputFileNames({"parameter_space.nc"});
    // We can use the simple cache as it doesn't depend on parameter space
    parameterSpaceProcessor.useCache();

    shellSiteFileAnalyzer.setCommand(pythonBinary);
    shellSiteFileAnalyzer.setScriptName(pythonScriptsPath.get() +
                                        "/analyze_shell_site_files.py");

    transfmatExtractor.setCommand(pythonBinary);
    transfmatExtractor.setScriptName(pythonScriptsPath.get() +
                                     "/extract_transfmat.py");
    transfmatExtractor.setOutputFileNames({"transfmat"});
    // We can use the simple cache as it doesn't depend on parameter space
    transfmatExtractor.useCache();

    templateGen.setCommand(pythonBinary);
    templateGen.setScriptName(pythonScriptsPath.get() +
                              "/reassign_occs/template_creator.py");
    templateGen.setOutputFileNames({"cached_output/template.nc"});
    // We can use the simple cache as it doesn't depend on parameter space
    templateGen.useCache();

    bool verbose = false;
    parameterSpaceProcessor.setVerbose(verbose);
    transfmatExtractor.setVerbose(verbose);
    templateGen.setVerbose(verbose);

    templateGen.prepareFunction = [&]() {
      // if transfmat extractor could not generate file, search for alternatives
      if (transfmatExtractor.getOutputFileNames().size() > 0 &&
          File::exists(transfmatExtractor.getOutputDirectory() +
                       transfmatExtractor.getOutputFileNames()[0])) {

        templateGen.setInputFileNames(
            {transfmatExtractor.getOutputFileNames()[0]});
        templateGen.setInputDirectory(transfmatExtractor.getOutputDirectory());
        return true;
      } else if (File::exists(templateGen.getRunningDirectory() +
                              "cached_output/transfmat")) {
        templateGen.setInputFileNames({"cached_output/transfmat"});
      } else if (File::exists(templateGen.getRunningDirectory() +
                              "transfmat")) {
        templateGen.setInputFileNames({"transfmat"});
      } else if (File::exists(templateGen.getRunningDirectory() +
                              "../transfmat")) {
        templateGen.setInputFileNames({"../transfmat"});
      } else {
        std::cerr << "Transformation matrix not found!" << std::endl;
        return false;
      }
      return true;
    };
    initRootProcessorGraph << parameterSpaceProcessor << shellSiteFileAnalyzer
                           << transfmatExtractor << templateGen;

    if (isPrimary()) {
      tinc = &tincServer;
    } else {
      tinc = &tincClient;
    }
    *tinc << initRootProcessorGraph;

    *tinc << dataDisplays[0]->mDatasetManager.sampleProcessorGraph;
    *tinc << dataDisplays[0]->mDatasetManager.mParameterSpace;
    *tinc << dataDisplays[0]->mDatasetManager.mShellSiteTypes;
    *tinc << dataDisplays[0]->mDatasetManager.mPercolationTypes;
    *tinc << dataDisplays[0]->mDatasetManager.dataPool;
    *tinc << dataDisplays[0]->mDatasetManager.trajectoriesPool;
    *tinc << dataDisplays[0]->mDatasetManager.neighborhoodPool;

    *tinc << dataDisplays[0]->imageDiskBuffer;
    for (auto &db : dataDisplays[0]->imageDiskBuffers) {
      *tinc << *db;
    }

    *tinc << dataDisplays[0]->mMarkerColor << dataDisplays[0]->mMarkerScale
          << dataDisplays[0]->mPercoMarkerScale
          << dataDisplays[0]->currentSelection
          << dataDisplays[0]->previousSelection;

    *tinc << mDataset;

    // Expose buffers from renderers to TINC
    dataDisplays[0]->mHistoryRender.trajectoryWidth.max(100);
    dataDisplays[0]->mTrajRender.trajectoryWidth.max(100);

    dataDisplays[0]->mHistoryRender.registerWithTinc(*tinc);
    dataDisplays[0]->mTrajRender.registerWithTinc(*tinc);
    dataDisplays[0]->atomrender.registerWithTinc(*tinc);
    dataDisplays[0]->mTriRender.registerWithTinc(*tinc);

    if (sphere::isSphereMachine()) {
      if (isPrimary()) {
        dataDisplays[0]->mHistoryRender.setRootPath("/Volumes/Data/tmp");
        dataDisplays[0]->mTrajRender.setRootPath("/Volumes/Data/tmp");
        dataDisplays[0]->atomrender.setRootPath("/Volumes/Data/tmp");
        dataDisplays[0]->mTriRender.setRootPath("/Volumes/Data/tmp");
      } else {
        dataDisplays[0]->mHistoryRender.setRootPath("/data/tmp");
        dataDisplays[0]->mTrajRender.setRootPath("/data/tmp");
        dataDisplays[0]->atomrender.setRootPath("/data/tmp");
        dataDisplays[0]->mTriRender.setRootPath("/data/tmp");
      }
    }

    auto primaryHost = getPrimaryHost();
    if (isPrimary()) {
      // Configure TINC server
      tincServer.setVerbose(true);
      if (sphere::isSimulatorMachine()) {
        tincServer.start(34450, "192.168.0.151");
      } else {
        tincServer.start();
      }
      tincServer.setRootMapEntry("/Volumes/Data/tmp", "/data/tmp");
    } else {
      tincClient.setVerbose(true);
      if (sphere::isRendererMachine()) {
        std::cout << "Client on renderer" << std::endl;
        tincClient.start(34450, primaryHost.c_str());
      } else {
        tincClient.start();
      }
      tincClient.setRootMapEntry("/Volumes/Data/tmp", "/data/tmp");
    }
    //    percoTools.init();
  }

  void registerParameterCallbacks() {
    for (auto *display : dataDisplays) {
      display->mVisible.registerChangeCallback(
          [this](float) { updateTitle(); });

      display->mDatasetManager.mCurrentDataset.registerChangeCallback(
          [this, display](std::string value) {
            std::string path = value;
            path = path.substr(path.rfind('/') + 1);
            std::string datasetPath;
            if (display->mDatasetManager.mGlobalRoot.size() > 0 &&
                display->mDatasetManager.mGlobalRoot != "./") {
              datasetPath = display->mDatasetManager.mGlobalRoot;
            }
            datasetPath += value + "/presets";
            presetHandler->setRootPath(datasetPath);
            positionPresets->setRootPath(datasetPath);
            presetHandler->setCurrentPresetMap("default", true);
            std::cout << "Preset Handler sub dir set to " << value << std::endl;
            updateTitle();
            mAutoAdvance = 0.0; // Turn off auto advance
          });
      display->mDatasetManager.currentGraphName.setSynchronousCallbacks();
      display->mDatasetManager.currentGraphName.registerChangeCallback(
          [this, display](std::string value) {
            std::string fullDatasetPath = File::conformPathToOS(
                display->mDatasetManager.getGlobalRootPath() +
                File::conformPathToOS(
                    display->mDatasetManager.mCurrentDataset.get()));

            std::cout << "loading graph " << value << " at " << fullDatasetPath
                      << std::endl;
            // New image module puts origin on top right
            //        if (mGraphTextureLock.try_lock()) {
            //            if (mGraphFilePathToLoad.size() > 0) {
            display->imageDiskBuffer.setRootPath(
                display->mDatasetManager.graphGenerator.getOutputDirectory());
            display->imageDiskBuffer.loadData(value, isPrimary());
          });
      for (size_t i = 0; i < display->graphCount; i++) {
        display->currentGraphNames[i]->setSynchronousCallbacks();
        display->currentGraphNames[i]->registerChangeCallback(
            [&](std::string value) {
              std::string fullDatasetPath = File::conformPathToOS(
                  display->mDatasetManager.getGlobalRootPath() +
                  File::conformPathToOS(
                      display->mDatasetManager.mCurrentDataset.get()));

              std::cout << "loading graph [" << i << "] " << value << " at "
                        << fullDatasetPath << std::endl;
              // New image module puts origin on top right
              //        if (mGraphTextureLock.try_lock()) {
              //            if (mGraphFilePathToLoad.size() > 0) {
              display->imageDiskBuffers[i]->loadData(value, isPrimary());
            });
      }
    }

    // Triggers and callbacks that should only be handled by rank 0
    if (rank == 0) {
      stepXpos.registerChangeCallback([this](float value) {
        ObjectTransformHandler &oth = this->object_transform;
        oth.quat =
            Quatf().fromEuler(Vec3f((1 / 24.0) * 2 * M_PI, 0, 0)) * oth.quat;
      });

      stepXneg.registerChangeCallback([this](float value) {
        ObjectTransformHandler &oth = this->object_transform;
        oth.quat =
            Quatf().fromEuler(-Vec3f((1 / 24.0) * 2 * M_PI, 0, 0)) * oth.quat;
      });

      stepYpos.registerChangeCallback([this](float value) {
        ObjectTransformHandler &oth = this->object_transform;
        oth.quat =
            Quatf().fromEuler(Vec3f(0, (1 / 24.0) * 2 * M_PI, 0)) * oth.quat;
      });

      stepYneg.registerChangeCallback([this](float value) {
        ObjectTransformHandler &oth = this->object_transform;
        oth.quat =
            Quatf().fromEuler(-Vec3f(0, (1 / 24.0) * 2 * M_PI, 0)) * oth.quat;
      });

      stepZpos.registerChangeCallback([this](float value) {
        ObjectTransformHandler &oth = this->object_transform;
        oth.quat =
            Quatf().fromEuler(Vec3f(0, 0, (1 / 24.0) * 2 * M_PI)) * oth.quat;
      });

      stepZneg.registerChangeCallback([this](float value) {
        ObjectTransformHandler &oth = this->object_transform;
        oth.quat =
            Quatf().fromEuler(-Vec3f(0, 0, (1 / 24.0) * 2 * M_PI)) * oth.quat;
      });

      stepTempPos.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->nextTemp();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->nextTemp();
        }
      });

      stepTempNeg.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->previousTemp();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->previousTemp();
        }
      });

      stepChempotPos.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->nextChempot();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->nextChempot();
        }
      });

      stepChempotNeg.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->previousChempot();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->previousChempot();
        }
      });

      stepChempot2Pos.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->nextChempot2();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->nextChempot2();
        }
      });

      stepChempot2Neg.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->previousChempot2();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->previousChempot2();
        }
      });

      stepTimePos.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->nextTime();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->nextTime();
        }
      });

      stepTimeNeg.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->previousTime();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->previousTime();
        }
      });

      ResetSlicing.registerChangeCallback([this](float value) {
        this->dataDisplays[vdvBundle.currentBundle()]->resetSlicing();
      });

      mJumpLayerPos.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->nextLayer();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->nextLayer();
        }
      });

      mJumpLayerNeg.registerChangeCallback([this](float value) {
        if (vdvBundle.bundleGlobal()) {
          for (auto *display : this->dataDisplays) {
            display->previousLayer();
          }
        } else {
          this->dataDisplays[vdvBundle.currentBundle()]->previousLayer();
        }
      });

      mSaveGraphics.registerChangeCallback([this](float value) {
        time_t rawtime;
        struct tm *timeinfo;
        char buffer[80];
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer, sizeof(buffer), "%d-%m-%Y_%H%M%S", timeinfo);
        std::string str(buffer);
        //        int counter = 0;
        //        for (DataDisplay *display : this->dataDisplays) {
        //          if (display->mVisible == 1.0f) {
        //            std::string prefix = str + "_" +
        //            std::to_string(counter++); display->dumpImages(prefix);
        //          }
        //        }
        mScreenshotMutex.lock();
        mScreenshotPrefix = str;
        mScreenshotMutex.unlock();
      });
      if (isPrimary()) {
        // Only connect on primary to avoid feedback
        stepPitchAngleNeg.registerChangeCallback([this](float value) {
          if (vdvBundle.bundleGlobal()) {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationPitch;
            newAngle -= pitchAngleStep * M_2PI / 360.0;
            for (auto *display : dataDisplays) {
              display->atomrender.mSliceRotationPitch = newAngle;
            }

          } else {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationPitch;
            newAngle -= pitchAngleStep * M_2PI / 360.0;
            dataDisplays[vdvBundle.currentBundle()]
                ->atomrender.mSliceRotationPitch = newAngle;
          }
        });
        stepPitchAnglePos.registerChangeCallback([this](float value) {
          if (vdvBundle.bundleGlobal()) {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationPitch;
            newAngle += pitchAngleStep * M_2PI / 360.0;
            newAngle -= pitchAngleStep * M_2PI / 360.0;
            for (auto *display : dataDisplays) {
              display->atomrender.mSliceRotationPitch = newAngle;
            }

          } else {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationPitch;
            newAngle += pitchAngleStep * M_2PI / 360.0;
            dataDisplays[vdvBundle.currentBundle()]
                ->atomrender.mSliceRotationPitch = newAngle;
          }
        });
        stepRollAngleNeg.registerChangeCallback([this](float value) {
          if (vdvBundle.bundleGlobal()) {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationRoll;
            newAngle -= rollAngleStep * M_2PI / 360.0;
            for (auto *display : dataDisplays) {
              display->atomrender.mSliceRotationRoll = newAngle;
            }
          } else {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationRoll;
            newAngle -= rollAngleStep * M_2PI / 360.0;
            dataDisplays[vdvBundle.currentBundle()]
                ->atomrender.mSliceRotationRoll = newAngle;
          }
        });
        stepRollAnglePos.registerChangeCallback([this](float value) {
          if (vdvBundle.bundleGlobal()) {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationRoll;
            newAngle += rollAngleStep * M_2PI / 360.0;
            for (auto *display : dataDisplays) {
              display->atomrender.mSliceRotationRoll = newAngle;
            }
          } else {
            float newAngle = dataDisplays[vdvBundle.currentBundle()]
                                 ->atomrender.mSliceRotationRoll;
            newAngle += rollAngleStep * M_2PI / 360.0;
            dataDisplays[vdvBundle.currentBundle()]
                ->atomrender.mSliceRotationRoll = newAngle;
          }
        });
      }

      mAutoAdvance.registerChangeCallback([this](float value) {
        if (value == 0.0f) {
          mParameterPlayback.stop();
        } else {
          std::chrono::nanoseconds wait =
              std::chrono::nanoseconds(uint64_t(1.0e9 / mAutoAdvanceFreq));
          mParameterPlayback.setWaitTime(wait);
          DataDisplay *d = dataDisplays[vdvBundle.currentBundle()];
          mParameterPlayback.start([d]() {
            d->nextTime();
            return true;
          });
        }
      });

      mDataset.registerChangeCallback([&](std::string value) {
        if (value != mDataset.get()) {
          // Only trigger load if value has changed
          loadDataset(value, vdvBundle.currentBundle());
        }
      });

      // We want to process title changes in the graphics thread.
      mWindowTitle.setSynchronousCallbacks(false);

      mRecentDatasets.registerChangeCallback([&](int index) {
        mDataset.set(mRecentDatasets.getElements()[index]);
        if (mDatasetSelector) {
          mDatasetSelector->cancel();
        }
      });
    }

    resetView.registerChangeCallback(
        [this](float value) { this->object_transform.reset(); });

    Z.registerChangeCallback(
        [this](float value) { this->nav().pos()[2] = value; });
    X.registerChangeCallback(
        [this](float value) { this->nav().pos()[0] = value; });

    pythonScriptsPath.registerChangeCallback([&](std::string value) {
      for (auto *display : dataDisplays) {
        display->mDatasetManager.setPythonScriptPath(value);
      }
    });
    pythonBinary.registerChangeCallback([&](std::string value) {
      for (auto *display : dataDisplays) {
        display->mDatasetManager.setPythonBinary(value);
      }
    });
    sliceBackground.registerChangeCallback([&](Color value) {
      for (auto *display : dataDisplays) {
        display->backgroundColor.set(value);
      }
    });

    font.registerChangeCallback([&](std::string value) {
      for (auto *display : dataDisplays) {
        display->setFont(value, fontSize);
      }
    });
    fontSize.registerChangeCallback([&](float value) {
      for (auto *display : dataDisplays) {
        display->setFont(font.get(), value);
      }
    });
    recallPositions.registerChangeCallback([&](float value) {
      for (auto *display : dataDisplays) {
        for (auto *param : display->graphPickable.bundle.parameters()) {
          std::string addr = display->graphPickable.bundle.bundlePrefix();
          presetHandler->skipParameter(addr + param->getFullAddress(),
                                       value != 1.0f);
        }
        for (auto *param : display->parallelPickable.bundle.parameters()) {
          std::string addr = display->parallelPickable.bundle.bundlePrefix();
          presetHandler->skipParameter(addr + param->getFullAddress(),
                                       value != 1.0f);
        }
        for (auto *param : display->perspectivePickable.bundle.parameters()) {
          std::string addr = display->perspectivePickable.bundle.bundlePrefix();
          presetHandler->skipParameter(addr + param->getFullAddress(),
                                       value != 1.0f);
        }
      }
    });
    mDebugScripts.registerChangeCallback([&](float value) {
      for (auto *proc : initRootProcessorGraph.getProcessors()) {
        proc->setVerbose(value != 0.0);
      }
      for (auto d : dataDisplays) {
        for (auto *proc :
             d->mDatasetManager.sampleProcessorGraph.getProcessors()) {
          proc->setVerbose(value != 0.0);
        }
      }
    });

    mComputeCache.registerChangeCallback([&](float value) {
      if (value != 0.0) {
        dataDisplays[vdvBundle.currentBundle()]
            ->mDatasetManager.mParameterSpace.sweep(
                dataDisplays[vdvBundle.currentBundle()]
                    ->mDatasetManager.sampleProcessorGraph);
      } else {
        dataDisplays[vdvBundle.currentBundle()]
            ->mDatasetManager.mParameterSpace.stopSweep();
      }
    });
  }

  Rayd rayTransformAllosphere(Rayd r) {
    double t = r.intersectAllosphere(); // get t on surface of allosphere screen
    Vec3f pos =
        nav().quat().rotate(r(t)); // rotate point on allosphere to match
                                   // current nav orientation (check this)
    Rayd ray(nav().pos(), pos); // ray from sphere center (camera location) to
                                // intersected location
    return ray;
  }

  void updateTitle() {
    std::string newTitle = "CASM Viewer ";
    for (auto *display : dataDisplays) {
      if (display->mVisible == 1.0f) {
        newTitle += display->mDatasetManager.mCurrentDataset;
        newTitle += " -- ";
      }
    }
    mWindowTitle.set(newTitle);
  }
};

int main() {
  // App is too big for stack!
  std::unique_ptr<MyApp> app = std::make_unique<MyApp>();
  app->start();
  return 0;
}
