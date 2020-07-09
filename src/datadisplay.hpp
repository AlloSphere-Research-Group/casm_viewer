#ifndef DATASETDISPLAY_HPP
#define DATASETDISPLAY_HPP

#ifdef AL_WINDOWS
#define NOMINMAX
#endif

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <regex>
#include <string>
#include <thread>
#include <vector>

#undef CIEXYZ

#include "al/graphics/al_Font.hpp"
#include "al/graphics/al_Image.hpp"
#include "al/graphics/al_Light.hpp"
#include "al/io/al_CSVReader.hpp"
#include "al/io/al_File.hpp"
#include "al/io/al_Toml.hpp"
#include "al/math/al_Matrix4.hpp"
#include "al/types/al_Color.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/ui/al_ParameterBundle.hpp"
#include "al/ui/al_Pickable.hpp"
#include "al/ui/al_PickableManager.hpp"
#include "al/ui/al_PickableRotateHandle.hpp"

#include "tinc/AtomRenderer.hpp"
#include "tinc/DataScript.hpp"
#include "tinc/VASPReader.hpp"

#include "datasetmanager.hpp"
#include "slice.hpp"

#undef far
#undef near

using namespace al;
using namespace std;

inline Matrix4f getLookAt(const Vec3f &ux, const Vec3f &uy, const Vec3f &uz,
                          const Vec3f &p) {
  return Matrix4f(ux[0], ux[1], ux[2], -(ux.dot(p)), uy[0], uy[1], uy[2],
                  -(uy.dot(p)), uz[0], uz[1], uz[2], -(uz.dot(p)), 0, 0, 0, 1);
}

inline Matrix4f getLookAt(const Vec3f &eyePos, const Vec3f &at,
                          const Vec3f &up) {
  Vec3f z = (eyePos - at).normalize();
  Vec3f x = cross(up, z).normalize();
  Vec3f y = cross(z, x).normalize();
  return getLookAt(x, y, z, eyePos);
}

inline HSV rgb2hsv(RGB c) {
  Vec4f K = Vec4f(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
  Vec4f p = c.g < c.b ? Vec4f(c.b, c.g, K.w, K.z) : Vec4f(c.g, c.b, K.x, K.y);
  Vec4f q = c.r < p.x ? Vec4f(p.x, p.y, p.w, c.r) : Vec4f(c.r, p.y, p.z, p.x);
  float d = q.x - min(q.w, q.y);
  float e = 1.0e-10;
  return HSV(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

struct ElementData {
  float radius;
  Color color;
};

class DataDisplayParameters {
public:
  ParameterVec3 layerDir{"LayerDir", ""}; // Direction of layers in data

  ParameterChoice mShowAtoms{"ShowAtoms"};
  ParameterBool mCumulativeTrajectory{"ShowCumulativeTrajectory", "", 1.0};
  ParameterBool mIndividualTrajectory{"ShowIndividualTrajectory", "", 1.0};

  ParameterColor backgroundColor{"projectionBackground", "",
                                 Color(0.1f, 0.1f, 0.1f, 0.8f)};
  ParameterBool mVisible{"Visible", "", 1};
  ParameterBool mShowGraph{"ShowGraph", "", 1};
  ParameterBool mShowParallel{"ShowParallel", "", 1};
  ParameterBool mShowPerspective{"ShowPerspective", "", 1};

  ParameterBool mBillboarding{"Billboarding", "", 1};
  ParameterBool mSmallLabel{"SmallLabel", "", 1};
  ParameterBool mDrawLabels{"DrawLabels", "", 1};

  Parameter mLayerScaling{"LayerScaling",
                          "",
                          1.0,
                          "",
                          0,
                          3}; // Increase layer size in projection view

  Parameter mPerspectiveRotY{"perspectiveRotY", "", -75, "", -90, 90};

  ParameterBool mShowGrid{"ShowGrid", "", 0.0};
  Parameter mGridSpacing{"GridSpacing", "", 10.0, "", 1.0, 20.0};
  ParameterMenu mGridType{"GridType"};

  Parameter mGridXOffset{"GridXOffset", "", 0.0, "", 0.0, 5.0};
  Parameter mGridYOffset{"GridYOffset", "", 0.0, "", 0.0, 5.0};

  ParameterBundle bundle{"CASMDataset"};
};

// -----------------------------------------------------------------

class DataDisplay : public DataDisplayParameters {
public:
  std::map<std::string, ElementData> elementData;

  DatasetManager mDatasetManager;

  Mesh mHistoryMesh;
  Mesh mTrajectoryMesh;

  // Pickables
  PickableManager mPickableManager;
  PickableBB graphPickable{"graph"};
  PickableBB parallelPickable{"parallel"};
  PickableBB perspectivePickable{"perspective"};
  PickableBB slicePickable{"slice"};
  // PickableRotateHandle rh;

  std::string mParamText;

  std::mutex mDrawLock;

  SlicingAtomRenderer atomrender;

  std::vector<Color> colorList = {
      Color(0.0, 1.0, 0.0, 0.0), Color(1.0, 0.0, 0.0, 0.0),
      Color(0.0, 1.0, 1.0, 1.0), Color(1.0, 1.0, 0.0, 1.0),
      Color(1.0, 0.0, 1.0, 1.0), Color(0.5, 1.0, 0.5, 1.0),
      Color(1.0, 1.0, 0.0, 1.0), Color(0.0, 0.0, 1.0, 1.0)};

  void init();

  void initDataset();

  // Prepare elements before draw call
  void prepare(Graphics &g, Matrix4f &transformMatrix);
  void prepareHistoryMesh();
  void prepareParallelProjection(Graphics &g, Matrix4f &transformMatrix);

  // Draw function
  void draw(Graphics &g);

  void setFont(std::string name, float size);

  // Functions to change increase/decrease parameters according to internal
  // parameter spaces

  void nextTemp();
  void previousTemp();

  void nextChempot();
  void previousChempot();

  void nextChempot2();
  void previousChempot2();

  void nextTime();
  void previousTime();

  void nextLayer();
  void previousLayer();

  //  void dumpImages(std::string dumpPrefix);

  void computeSlicing() {
    if (mDatasetManager.mRunProcessors) {
      atomrender.mSlicingPlaneThickness =
          findDistanceNormal(mAligned4fData, layerDir.get());
      std::cout << "Data Boundaries:" << std::endl;
      std::cout << "X " << mDataBoundaries.min.x << " -> "
                << mDataBoundaries.max.x << std::endl;
      std::cout << "Y " << mDataBoundaries.min.y << " -> "
                << mDataBoundaries.max.y << std::endl;
      std::cout << "Z " << mDataBoundaries.min.z << " -> "
                << mDataBoundaries.max.z << std::endl;
      Vec3f point = atomrender.mSlicingPlanePoint;
      point.z = mDataBoundaries.min.z;
      atomrender.mSlicingPlanePoint = point;
    }
  }

  void resetSlicing() {
    atomrender.mSlicingPlanePoint.set(slicePickable.bb.min);

    atomrender.mSlicingPlaneThickness =
        mDataBoundaries.max.z - mDataBoundaries.min.z;
    atomrender.mSliceRotationRoll.set(0);
    atomrender.mSliceRotationPitch.set(0);
  }

protected:
  Texture &iso_scene() { return fbo_iso.tex(); }

  // This function should be called whenever there is new atom position data
  void updateDisplayBuffers();

  void drawHistory(Graphics &g);
  void drawPerspective(Graphics &g);
  void drawParallelProjection(Graphics &g);
  void drawGraph(Graphics &g);

private:
  VAOMesh axis;
  VAOMesh orthoMesh;
  VAOMesh graphlinesMesh;
  VAOMesh boxMesh;
  VAOMesh mGridMesh;
  VAOMesh mMarker;

  // Graph
  Texture mGraphTexture;

  Font mLabelFont;

  EasyFBO fbo_iso;

  // Schedule processing for next frame to ensure drawing current
  // frame and processing current window events
  bool mNeedsProcessing{false};

  // For KMC datasets list of atoms added/removed with respect to previous time
  // step
  std::map<std::string, std::vector<DatasetManager::position_t>> atomAdded;
  std::vector<DatasetManager::position_t> atomRemoved;

  std::map<std::string, AtomData> mAtomData;
  std::vector<float> mAligned4fData;
  BoundingBoxData mDataBoundaries;
};

#endif // DATASETDISPLAY_HPP
