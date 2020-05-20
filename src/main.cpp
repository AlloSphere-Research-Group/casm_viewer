
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

#include "tinc/ProcessorServer.hpp"

#undef AL_BUILD_MPI

//#define STB_IMAGE_WRITE_STATIC
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "module/img/stb_image_write.h"

#include "datadisplay.hpp"

#ifdef AL_WINDOWS
#include <direct.h> // for _chdir() and _getcwd()
#define chdir _chdir
#define getcwd _getcwd
#undef rank
#undef CIEXYZ
#endif

using namespace al;
using namespace std;

struct State {
  float zoom = 0.0;
  //    Nav nav;
  Matrix4f transformMatrix;
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
  Parameter pitchAngleStep{"PitchAngleStep", "AngleControl", 15, "", 0.0, 60.0};
  Trigger stepPitchAnglePos{"stepPitchAnglePos", "AngleControl"};
  Trigger stepPitchAngleNeg{"stepPitchAngleNeg", "AngleControl"};
  Parameter rollAngleStep{"RollAngleStep", "AngleControl", 15, "", 0.0, 60.0};
  Trigger stepRollAnglePos{"stepRollAnglePos", "AngleControl"};
  Trigger stepRollAngleNeg{"stepRollAngleNeg", "AngleControl"};
  Trigger CalculateSlicing{"CalculateSlicing"};
  Trigger ResetSlicing{"ResetSlicing"};
  Trigger mSaveGraphics{"SaveSnapshot"};
  Trigger mAlignTemperatures{"AlignTemperatures"};
  Trigger mAlignChempots{"AlignChempots"};
  Trigger mRecomputeSpace{"RecomputeSpace"};

  ParameterBool mAutoAdvance{"autoAdvance"};
  Parameter mAutoAdvanceFreq{"autoAdvanceFreq", "", 2, "", 0.25, 3.0};

  // File selection
  ParameterMenu mDataRootPath{"datarootPath"};
  ParameterMenu mAvailableDatasets{"availableDatasets"};

  // 3D nav parameters
  Parameter Z{"Z", "", 2.8f, "", -15, 15};
  Parameter X{"X", "", 0.0, "", -2, 2};

  // Settings
  ParameterColor backgroundColor{"background", "",
                                 Color(0.0f, 0.0f, 0.0f, 1.0f)};
  ParameterColor sliceBackground{"sliceBackground", "",
                                 Color(0.0f, 0.0f, 0.0f, 1.0f)};
  ParameterString font{"font"};
  Parameter fontSize{"fontSize", "", 32, "", 2, 128};
  ParameterString pythonScriptPath{"pythonScriptPath"};
  ParameterString pythonBinary{"pythonBinary"};

  ParameterBool recallPositions{"recallPositions", "", 1.0};

  std::unique_ptr<PresetHandler> presetHandler;
  std::unique_ptr<PresetSequencer> sequencer;
  std::unique_ptr<SequenceRecorder> recorder;
  std::unique_ptr<PresetServer> presetServer;

  BundleGUIManager vdvBundle;

  std::mutex mScreenshotMutex;
  std::string mScreenshotPrefix;

  PeriodicTask mParameterPlayback;

  // TINC computation chains
  ComputationChain initRootComputationChain{"PrepareDataset"};

  DataScript parameterSpaceProcessor{"ParameterSpaceProcessor"};
  DataScript templateGen{"TemplateGenerator"};

  ProcessorServer processorServer; // To expose processor chains on the network

  // --------
  std::vector<DataDisplay *> dataDisplays;

  ObjectTransformHandler object_transform;
  bool showGui{true};
  bool mMarkUpdateTitle{true};

#ifdef AL_EXT_OPENVR
  std::shared_ptr<OpenVRDomain> openVRDomain;
  VAOMesh mControllerMesh;
#endif

  // ------------- Application callbacks

  void onInit() override {
    presetHandler = std::make_unique<PresetHandler>();
    sequencer = std::make_unique<PresetSequencer>();
    recorder = std::make_unique<SequenceRecorder>();
    presetServer = std::make_unique<PresetServer>();
  }

  virtual void onCreate() override {
    imguiInit();

    // disable nav control mouse drag to look
    navControl().useMouse(false);

    fps(60);
    title("Casm Viewer");
    dimensions(1200, 800);

    // Create and configure displays
    int numDisplays = 2; // 2 is a good number. Usually need one, but some
                         // times need to compare 2

    for (int i = 0; i < numDisplays; i++) {
      dataDisplays.push_back(new DataDisplay);
      dataDisplays.back()->init();
      dataDisplays.back()->mDatasetManager.mRunProcessors = rank == 0;
      dataDisplays.back()->mDatasetManager.mGlobalRoot = dataRoot;
    }

    // Initialize default view.
    dataDisplays[0]->graphPickable.pose.setPos(Vec3d(-2.0, 0.8, -.75));
    dataDisplays[0]->graphPickable.scale = 0.2f;
    dataDisplays[0]->parallelPickable.pose.setPos(Vec3d(-0.68, 1.1, -1.5));
    dataDisplays[0]->parallelPickable.scale = 0.4f;
    dataDisplays[0]->perspectivePickable.pose.setPos(Vec3d(0.55, 0.35, -0.7));
    dataDisplays[0]->perspectivePickable.scale = 0.02f;

    dataDisplays[0]->mShowGraph = true;
    dataDisplays[0]->mShowParallel = true;
    //                dataDisplays[0]->mShowSurface = false;
    dataDisplays[0]->mShowPerspective = true;
    dataDisplays[0]->mBillboarding = false;
    dataDisplays[0]->mSmallLabel = true;
    dataDisplays[0]->mVisible = true;
    for (size_t i = 1; i < dataDisplays.size(); i++) {
      dataDisplays[i]->mVisible = false;
    }

    // ----------
    registerParameterCallbacks();
    loadConfiguration();
    initializeComputation();
    setupParameterServer();
    readElementsIni();

    mDataRootPath.set(0); // Force loading available

    if (isPrimary()) {
      for (auto *display : dataDisplays) {
        *presetHandler << display->bundle;
        *presetHandler << display->graphPickable.bundle
                       << display->parallelPickable.bundle
                       << display->perspectivePickable.bundle;
        vdvBundle << display->bundle;
      }
      // Add parameters that are not part of the bundle
      *presetHandler << Z << X; //

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

    // use object control for model matrix
    append(object_transform);
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
  }

  virtual void onAnimate(double /*dt*/) override {
    object_transform.step();
    if (isPrimary()) {
      if (showGui) {
        prepareGui();
      }
      navControl().active(!ParameterGUI::usingInput());
    }

    // fbo operations should be done outside onDraw so it does not mess with
    // omni drawing
    for (auto *display : dataDisplays) {
      display->prepare(graphics(), state().transformMatrix);
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
          display->mPickableManager.event(PickEvent(Point, ray));
        }
      }
    }
#endif
  }

  virtual void onDraw(Graphics &g) override {
    g.clear(backgroundColor);
    if (mMarkUpdateTitle) {
      std::string newTitle = "CASM Viewer ";
      for (auto *display : dataDisplays) {
        if (display->mVisible == 1.0f) {
          newTitle += display->mDatasetManager.mRootPath.get() + " : ";
          newTitle += display->mDatasetManager.mCurrentDataset;
          newTitle += " -- ";
        }
      }
      title(newTitle);
      mMarkUpdateTitle = false;
    }
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
        case '[':
          dataDisplays[vdvBundle.currentBundle()]->mLayerScaling =
              dataDisplays[vdvBundle.currentBundle()]->mLayerScaling + 0.05f;
          break;
        case ']':
          dataDisplays[vdvBundle.currentBundle()]->mLayerScaling =
              dataDisplays[vdvBundle.currentBundle()]->mLayerScaling - 0.05f;
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
    //      std::cout << name();
    //      m.print();
    if ("/point" == m.addressPattern()) {
      float ox, oy, oz, dx, dy, dz;
      int id;
      m >> id >> ox >> oy >> oz >> dx >> dy >> dz;
      Rayd r = Rayd(Vec3f(ox, oy, oz), Vec3f(dx, dy, dz));
      // r = rayTransformAllosphere(r);
      for (auto *display : dataDisplays) {
        if (display->mVisible != 0.0f) {
          display->mPickableManager.event(PickEvent(Point, r));
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
    }
    imguiShutdown();
    //#ifdef AL_EXT_OPENVR
    //    mOpenVR.close();
    //#endif
  }

  // ---- Drawing functions

  void drawScene(Graphics &g) {
    g.pushMatrix();
    for (auto *display : dataDisplays) {
      display->draw(g);
    }
    g.popMatrix();
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

  // Display parameters GUI. This is called by drawGUI()
  void drawDisplayGui(bool &showLoadingWindow) {
    if (showLoadingWindow) {
      //        if (parameterSpaceProcessor.runningAsync()) {
      //          ImGui::Text("Recalculating parameter space");
      //        }

      if (dataRoot.size() > 0) {
        ImGui::Text("Data root: %s", dataRoot.c_str());
      }
      ParameterGUI::draw(&mDataRootPath);
      ImGui::SameLine();
      if (ImGui::Button("Remove")) {
        auto elements = mDataRootPath.getElements();
        int index = mDataRootPath.get();
        if (index >= 0) {
          elements.erase(elements.begin() + index);
        }
        mDataRootPath.setElements(elements);
        auto current = mDataRootPath.get();
        if (current < 0 && elements.size() > 0) {
          mDataRootPath.set(current); // Force calling callbacks.
        }
      }
      ParameterGUI::draw(&mAvailableDatasets);
      ImGui::SameLine();
      ParameterGUI::draw(&mRecomputeSpace);

      if (ImGui::Button("Load")) {
        std::string conformedPath =
            File::conformDirectory(mDataRootPath.getCurrent());
        std::replace(conformedPath.begin(), conformedPath.end(), '\\', '/');
        int currentBundle = vdvBundle.currentBundle();
        assert(currentBundle < dataDisplays.size());
        if (currentBundle >= 0) {
          dataDisplays[currentBundle]->mDatasetManager.mRootPath.set(
              conformedPath);
          dataDisplays[currentBundle]->mDatasetManager.mCurrentDataset.set(
              mAvailableDatasets.getCurrent());
        }
        showLoadingWindow = false;
      }
      ImGui::SameLine();
      if (ImGui::Button("Cancel##Loading")) {
        showLoadingWindow = false;
      }
      ImGui::Separator();
      // Directory GUI
      static FileSelector *selector;
      if (ImGui::Button("Add directory")) {
        selector = new FileSelector(dataRoot, File::isDirectory);
        selector->start(mDataRootPath.getCurrent());
      }
      if (selector && selector->drawFileSelector()) {
        auto elements = mDataRootPath.getElements();

        auto selectedItems = selector->getSelection();
        if (selectedItems.count() > 0) {
          auto selection = selectedItems[0];
          elements.push_back(selection.filepath());
        }
        delete selector;
        selector = nullptr;
        mDataRootPath.setElements(elements);
      }
    } else {
      if (ImGui::Button("Load dataset")) {
        showLoadingWindow = !showLoadingWindow;
      }
      for (auto &paramSpaces : this->dataDisplays[vdvBundle.currentBundle()]
                                   ->mDatasetManager.mParameterSpaces) {
        ParameterMeta *param = &paramSpaces.second->parameter();
        if (param->getHint("hide") != 1.0f) {
          ParameterGUI::drawParameterMeta(param);
        }
      }
      ParameterGUI::drawParameterMeta(
          &this->dataDisplays[vdvBundle.currentBundle()]
               ->atomrender.mAtomMarkerSize);
      if (this->dataDisplays[vdvBundle.currentBundle()]
              ->mDatasetManager.mParameterSpaces["time"]
              ->size() > 0) {
        ParameterGUI::drawParameterMeta(&mAutoAdvance);
        if (mAutoAdvance == 0.0) {
          ParameterGUI::drawParameterMeta(&mAutoAdvanceFreq);
        } else {
          std::string text = "Incrementing time with freq: " +
                             std::to_string(mAutoAdvanceFreq);
          ImGui::Text("%s", text.c_str());
        }
        ParameterGUI::drawParameterMeta(
            &this->dataDisplays[vdvBundle.currentBundle()]
                 ->mCumulativeTrajectory);
        ParameterGUI::drawParameterMeta(
            &this->dataDisplays[vdvBundle.currentBundle()]
                 ->mIndividualTrajectory);
      }
      //          ImGui::Separator();
      //          ParameterGUI::drawParameterMeta(&mViewMenu);
      ParameterGUI::drawTrigger(&mSaveGraphics);
      ImGui::SameLine();
      ParameterGUI::drawParameterMeta(&ResetSlicing);

      ImGui::Text("%s",
                  dataDisplays[vdvBundle.currentBundle()]->metaText.c_str());

      if (ImGui::CollapsingHeader("Slicing")) {
        ImGui::Indent(20.0);
        ParameterGUI::drawParameterMeta(&dataDisplays[vdvBundle.currentBundle()]
                                             ->atomrender.mSlicingPlaneNormal);
        ParameterGUI::drawParameterMeta(
            &dataDisplays[vdvBundle.currentBundle()]
                 ->atomrender.mSlicingPlaneThickness);

        ParameterGUI::drawParameterMeta(&mJumpLayerNeg);
        ImGui::SameLine();
        ParameterGUI::drawParameterMeta(&mJumpLayerPos);
        ParameterGUI::drawParameterMeta(&CalculateSlicing);
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

      if (ImGui::CollapsingHeader("Graphics")) {
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
        ParameterGUI::drawBundle(
            &dataDisplays[vdvBundle.currentBundle()]->parallelPickable.bundle);
        ParameterGUI::drawBundle(&dataDisplays[vdvBundle.currentBundle()]
                                      ->perspectivePickable.bundle);
        //          ParameterGUI::drawBundle(&dataDisplays[vdvBundle.currentBundle()]->slicePickable.bundle);
        //          // Slice pickable is already exposed enough above

        ImGui::Unindent(20.0);
      }
    }
  }

  void prepareGui() {
    if (parameterSpaceProcessor.runningAsync()) {
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
    ParameterGUI::beginPanel("CASM Viewer");

    ImGui::Columns(4, nullptr, true);
    static int selected = 0;
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

    static bool showLoadingWindow = true;
    if (selected == 0) {
      drawDisplayGui(showLoadingWindow);
    } else if (selected == 1) {
      if (showLoadingWindow) {
        ParameterGUI::drawParameterMeta(&mDataRootPath);
        ParameterGUI::drawParameterMeta(&mAvailableDatasets);

        if (ImGui::Button("Load")) {
          std::string conformedPath =
              File::conformDirectory(mDataRootPath.getCurrent());
          std::replace(conformedPath.begin(), conformedPath.end(), '\\', '/');
          int currentBundle = vdvBundle.currentBundle();
          assert(currentBundle < dataDisplays.size());
          if (currentBundle >= 0) {
            dataDisplays[currentBundle]->mDatasetManager.mRootPath.set(
                conformedPath);
            dataDisplays[currentBundle]->mDatasetManager.mCurrentDataset.set(
                mAvailableDatasets.getCurrent());
          }
          showLoadingWindow = false;
        }
        if (ImGui::Button("Cancel")) {
          showLoadingWindow = false;
        }
      } else {
        if (ImGui::Button("Load dataset")) {
          showLoadingWindow = !showLoadingWindow;
        }
        ParameterGUI::drawBundleManager(&vdvBundle);
      }

    } else if (selected == 2) {
      ParameterGUI::drawPresetHandler(presetHandler.get(), 10, 4);
      ParameterGUI::draw(&recallPositions);
      int currentSequencerItem = 0;
      ParameterGUI::drawPresetSequencer(sequencer.get(), currentSequencerItem);
      ParameterGUI::drawSequenceRecorder(recorder.get());
    } else if (selected == 3) {
      static std::string helpText;
      if (helpText.size() == 0) {
        std::ifstream t("help.txt");
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

      ImGui::Text("%s", helpText.c_str());
    }

    ParameterGUI::endPanel();
    ImGui::PopID();
    imguiEndFrame();
  }

  void setupParameterServer() {
    // Register parameters with parameter server
    for (auto *display : dataDisplays) {
      parameterServer() << display->bundle;
      parameterServer() << display->graphPickable.bundle;
      parameterServer() << display->parallelPickable.bundle;
      parameterServer() << display->perspectivePickable.bundle;
      parameterServer() << display->slicePickable.bundle;
    }
    parameterServer() << resetView << stepXpos << stepXneg << stepYpos
                      << stepYneg << stepZpos << stepZneg;
    parameterServer() << CalculateSlicing << ResetSlicing;
    parameterServer() << backgroundColor << sliceBackground;
    parameterServer() << font << fontSize;
    // These parameters are command triggers
    parameterServer() << stepTempPos << stepTempNeg << stepChempotPos
                      << stepChempotNeg << stepChempot2Pos << stepChempot2Neg
                      << mJumpLayerNeg << mJumpLayerPos;

    parameterServer().verbose(true);
  }

  // File IO/ Persisntece
  void loadConfiguration() {
    TomlLoader configLoader2;
    configLoader2.setFile("casm_viewer.toml");
    // Set default values
#ifdef AL_WINDOWS
    configLoader2.setDefaultValue("font",
                                  string("C:\\Windows\\Fonts\\arial.ttf"));
    configLoader2.setDefaultValueVector(
        "dataRootPaths", std::vector<std::string>({string("e:\\NaxCoO2")}));
    configLoader2.setDefaultValue("pythonScriptsPath",
                                  string("../vdv_group_python"));
#elif defined(AL_OSX)
    configLoader2.setDefaultValue("font", string("/Library/Fonts/arial.ttf"));
    configLoader2.setDefaultValueVector(
        "dataRootPaths",
        std::vector<std::string>({string("/alloshare/vdv group/NaxCoO2")}));
    configLoader2.setDefaultValue("pythonScriptsPath",
                                  string("../vdv_group_python"));
#else
    configLoader2.setDefaultValue(
        "font", string("/usr/share/fonts/truetype/freefont/FreeMono.ttf"));
    //        ["/alloshare/vdv group", "/alloshare/vdv group/TiAlO",
    //        "/alloshare/vdv group/NaxCoO2", "/alloshare/vdv group/proj1"]
    configLoader2.setDefaultValueVector<string>(
        "dataRootPaths",
        std::vector<std::string>({string("/alloshare/vdv group/NaxCoO2"),
                                  string("/alloshare/vdv group")}));
    configLoader2.setDefaultValue("pythonScriptsPath",
                                  string("../vdv_group_python"));
#endif
    configLoader2.setDefaultValue("fontSize", 32);
    configLoader2.setDefaultValue("pythonBinary", string("python3"));

    configLoader2.writeFile();

    std::string scriptPath = configLoader2.gets("pythonScriptsPath");
    if (scriptPath.substr(0, 2) == "..") { // Make path absolute
      char cwd[256];
      getcwd(cwd, sizeof(cwd));
      std::string cwdString = cwd;
      scriptPath =
          cwdString.substr(0, cwdString.rfind('/')) + scriptPath.substr(2);
    }
    vector<double> bgRGB = configLoader2.getVector<double>("background");
    if (bgRGB.size() == 4) {
      backgroundColor = Color(bgRGB[0], bgRGB[1], bgRGB[2], bgRGB[3]);
    }
    bgRGB = configLoader2.getVector<double>("sliceBackground");
    if (bgRGB.size() == 4) {
      sliceBackground = Color(bgRGB[0], bgRGB[1], bgRGB[2], bgRGB[3]);
    }

    std::string pythonBinaryConf = configLoader2.gets("pythonBinary");

    std::replace(scriptPath.begin(), scriptPath.end(), '\\', '/');
    std::replace(pythonBinaryConf.begin(), pythonBinaryConf.end(), '\\', '/');
    pythonScriptPath.set(scriptPath);
    pythonBinary.set(pythonBinaryConf);

    std::string joinedRootPaths;
    auto paths = configLoader2.getVector<string>("dataRootPaths");
    //    if (rank == 0) {
    //      //      for (auto path: paths) {
    //      //        processParameterSpace(path);
    //      //      }
    //    }
    mDataRootPath.setElements(paths);

    std::string fontPath = configLoader2.gets("font");
    std::replace(fontPath.begin(), fontPath.end(), '\\', '/');
    font.set(fontPath);
    fontSize = configLoader2.getd("fontSize");
  }

  void storeConfiguration() {
    char currentDir[512];
    getcwd(currentDir, 512);
    std::cout << "Writing config in " << currentDir << std::endl;

    TomlLoader configLoader2;
    configLoader2.setFile("casm_viewer.toml");

    vector<double> rgba;
    rgba.resize(4);
    Color bg = backgroundColor.get();
    rgba[0] = bg.r;
    rgba[1] = bg.g;
    rgba[2] = bg.b;
    rgba[3] = bg.a;
    configLoader2.setVector<double>("background", rgba);
    bg = sliceBackground.get();
    rgba[0] = bg.r;
    rgba[1] = bg.g;
    rgba[2] = bg.b;
    rgba[3] = bg.a;
    configLoader2.setVector<double>("sliceBackground", rgba);

    configLoader2.set<std::string>("font", font.get());
    configLoader2.set<float>("fontSize", fontSize.get());

    configLoader2.setVector<std::string>("dataRootPaths",
                                         mDataRootPath.getElements());

    configLoader2.writeFile();
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
        //            std::cout << atomName << ":" << radius << " " << r << "  "
        //            << g << " " << b << std::endl;
      }
    } else {
      std::cout << name() << "Could not open elements.ini" << std::endl;
    }
    std::cout << name() << "After Elements.ini" << std::endl;
  }

  void processScreenshot() {
    std::unique_lock<std::mutex> lk(mScreenshotMutex);
    if (mScreenshotPrefix.size() > 0) { // Take screenshot
      std::vector<unsigned char> mPixels;
      mPixels.resize(width() * height() * 3);
      unsigned char *pixs = &mPixels[0];
      // FIXME this can be done more efficiently
      glReadPixels(1, 1, width(), height(), GL_RGB, GL_UNSIGNED_BYTE, pixs);
      std::string dumpDirectory =
          File::conformPathToOS(dataDisplays[vdvBundle.currentBundle()]
                                    ->mDatasetManager.buildRootPath() +
                                dataDisplays[vdvBundle.currentBundle()]
                                    ->mDatasetManager.mCurrentDataset.get() +
                                "/graphics/");
      std::string imagePath = dumpDirectory + mScreenshotPrefix + "_screen.png";

      //      stbi_flip_vertically_on_write(1);
      //      stbi_write_png(imagePath.c_str(), width(), height(), 3, pixs,
      //                     width() * 3);

      mScreenshotPrefix = "";
    }
  }

  void initializeComputation() {
    parameterSpaceProcessor.setCommand(pythonBinary);
    parameterSpaceProcessor.setScriptName(pythonScriptPath.get() +
                                          "/analyze_parameter_space.py");

    parameterSpaceProcessor.setOutputFileNames({"_parameter_space.json"});
    parameterSpaceProcessor.verbose(true);

    templateGen.setCommand(pythonBinary);
    templateGen.setScriptName(pythonScriptPath.get() +
                              "/reassign_occs/template_creator.py");
    templateGen.setOutputFileNames({"template_POSCAR"});
    templateGen.verbose(true);

    initRootComputationChain << parameterSpaceProcessor << templateGen;

    // Only the root chain needs to be
    processorServer << initRootComputationChain
                    << dataDisplays[0]->mDatasetManager.sampleComputationChain;
    processorServer.exposeToNetwork(parameterServer());
  }

  void registerParameterCallbacks() {
    for (auto *display : dataDisplays) {
      display->mVisible.registerChangeCallback(
          [this](float) { updateTitle(); });

      display->mDatasetManager.mCurrentDataset.registerChangeCallback(
          [this, display](std::string value) {
            std::string path = value;
            path = path.substr(path.rfind('/') + 1);
            if (display->mDatasetManager.buildRootPath().size() > 0) {
              presetHandler->setRootPath(
                  display->mDatasetManager.buildRootPath());
              if (!File::exists(display->mDatasetManager.buildRootPath() +
                                value + "/presets")) {
                Dir::make(value + "/presets");
              }
              presetHandler->setSubDirectory(value + "/presets");
              presetHandler->setCurrentPresetMap("default", true);
              std::cout << "Preset Handler sub dir set to " << value
                        << std::endl;
            }
            updateTitle();
            mAutoAdvance = 0.0; // Turn off auto advance
          });
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

      CalculateSlicing.registerChangeCallback([this](float value) {
        this->dataDisplays[vdvBundle.currentBundle()]->computeSlicing();
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
        int counter = 0;
        for (DataDisplay *display : this->dataDisplays) {
          if (display->mVisible == 1.0f) {
            std::string prefix = str + "_" + std::to_string(counter++);
            display->dumpImages(prefix);
          }
        }
        mScreenshotMutex.lock();
        mScreenshotPrefix = str;
        mScreenshotMutex.unlock();
      });

      mAlignTemperatures.registerChangeCallback([this](float value) {
        float baseTemp = dataDisplays[0]
                             ->mDatasetManager.mParameterSpaces["temperature"]
                             ->parameter();
        for (size_t i = 1; i < dataDisplays.size(); i++) {
          dataDisplays[i]
              ->mDatasetManager.mParameterSpaces["temperature"]
              ->parameter() = baseTemp;
          for (size_t j = 0; j < i; j++) {
            dataDisplays[i]->nextTemp();
          }
        }
      });

      mAlignChempots.registerChangeCallback([this](float value) {
        float baseChempot = dataDisplays[0]
                                ->mDatasetManager.mParameterSpaces["chempotA"]
                                ->parameter();
        for (size_t i = 1; i < dataDisplays.size(); i++) {
          dataDisplays[i]
              ->mDatasetManager.mParameterSpaces["chempotA"]
              ->parameter() = baseChempot;
          for (size_t j = 0; j < i; j++) {
            dataDisplays[i]->nextChempot();
          }
        }
      });
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

      mRecomputeSpace.registerChangeCallback([&](float value) {
        if (value == 1.0f) {
          if (mAvailableDatasets.getCurrent() != "") {
            cleanParameterSpace(
                File::conformPathToOS(mDataRootPath.getCurrent()));
            processNewDataRoot(
                File::conformPathToOS(mDataRootPath.getCurrent()));
          }
        }
      });
    }

    resetView.registerChangeCallback(
        [this](float value) { this->object_transform.reset(); });

    Z.registerChangeCallback(
        [this](float value) { this->nav().pos()[2] = value; });
    X.registerChangeCallback(
        [this](float value) { this->nav().pos()[0] = value; });

    pythonScriptPath.registerChangeCallback([&](std::string value) {
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

    mDataRootPath.registerChangeCallback([&](int index) {
      if (this->rank == 0) {
        processNewDataRoot(mDataRootPath.getElements()[index]);
      }
    });
  }

  void cleanParameterSpace(string datasetPath) {
    if (File::isDirectory(dataRoot + datasetPath)) {
      FileList subDirs =
          filterInDir(dataRoot + datasetPath, [&](FilePath const &fp) {
            std::cout << fp.filepath() << std::endl;
            return File::isDirectory(fp.filepath());
          });
      subDirs.add(FilePath(dataRoot + datasetPath));
      for (auto subDir : subDirs) {
        std::string cachedOutputDir = subDir.filepath() + "/cached_output/";
        if (File::exists(cachedOutputDir) &&
            File::isDirectory(cachedOutputDir)) {
          if (!Dir::removeRecursively(cachedOutputDir)) {
            std::cerr << "ERROR cleaning cache at " << cachedOutputDir
                      << std::endl;
          }
        }
        // FIXME This should be removed once no files are generated outside
        // the cache directory
        std::vector<std::string> additionalCacheFiles = {"transfmat",
                                                         "template_POSCAR"};
        for (auto fileToRemove : additionalCacheFiles) {
          if (File::exists(subDir.filepath() + fileToRemove)) {
            if (!File::remove(subDir.filepath() + fileToRemove)) {
              std::cerr << "ERROR deleting cache file: "
                        << subDir.filepath() + fileToRemove << std::endl;
            }
          }
        }
      }
    }
  }

  void processNewDataRoot(string rootPath) {
    auto filelist = itemListInDir(dataRoot + rootPath);
    for (FilePath &element : filelist) {
      if (File::isDirectory(element.filepath())) {
        // Attempt to generate data space for all directories in root dir.

        parameterSpaceProcessor.setRunningDirectory(element.filepath());
        parameterSpaceProcessor.setOutputDirectory(element.filepath() +
                                                   "/cached_output");

        templateGen.setRunningDirectory(element.filepath());
        templateGen.setOutputDirectory(element.filepath() + "/cached_output");

        if (File::exists(templateGen.runningDirectory() +
                         "cached_output/transfmat")) {
          templateGen.configuration["transfmat"] =
              Flag("cached_output/transfmat");
        } else if (File::exists(templateGen.runningDirectory() + "transfmat")) {
          templateGen.configuration["transfmat"] = "transfmat";
        } else {
          templateGen.configuration["transfmat"] = "../transfmat";
        }

        bool ok = initRootComputationChain.process();

        std::cout << "Process parameter space asyc for " << element.filepath()
                  << std::endl;
        if (ok && File::exists(File::conformPathToOS(element.filepath()) +
                               "cached_output/_parameter_space.json")) {
          updateAvailableDatasets(this->dataRoot + rootPath);
        } else {
          //                Dir::removeRecursively(element.filepath() +
          //                "/cached_output");
        }
      }
    }
  }

  void updateAvailableDatasets(std::string datasetRoot) {
    auto entries = filterInDir(datasetRoot, [&](const FilePath &fp) {
      //      std::cout << fp.filepath() <<std::endl;
      return File::isDirectory(fp.filepath());
    });

    std::vector<std::string> directories;
    for (auto entry : entries) {
      auto spaceFileFullPath = File::conformPathToOS(entry.filepath()) +
                               "/cached_output/_parameter_space.json";
      if (File::exists(spaceFileFullPath)) {
        directories.push_back(entry.file());
      }
    }
    mAvailableDatasets.setElements(directories);
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

  //  const char *flowAddress() override { return "interface"; }

  void updateTitle() { mMarkUpdateTitle = true; }
};

int main() {
  std::unique_ptr<MyApp> app = std::make_unique<MyApp>();
  //  app->stereo(true);
  //  app->displayMode(app->displayMode() | Window::STEREO_BUF);
  //  app.print();
  //  app.initAudio();
  app->start();
  return 0;
}
