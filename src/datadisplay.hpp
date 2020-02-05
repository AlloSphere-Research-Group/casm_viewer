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

#include "al/math/al_Matrix4.hpp"

#include "al/graphics/al_Font.hpp"
#include "al/graphics/al_Light.hpp"
#include "al/io/al_CSVReader.hpp"
#include "al/io/al_File.hpp"
#include "al/io/al_Toml.hpp"
#include "al/types/al_Color.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/ui/al_ParameterBundle.hpp"
#include "al/ui/al_Pickable.hpp"
#include "al/ui/al_PickableManager.hpp"
#include "al/ui/al_PickableRotateHandle.hpp"

#include "al/graphics/al_Image.hpp"

#include "al_DataScript.hpp"
#include "al_VASPReader.hpp"

#include "instanced_mesh.hpp"
#include "processors.hpp"

#include "datasetmanager.hpp"
#include "slice.hpp"

//#include "imgui.h"

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

const std::string instancing_vert =
    R"(
#version 330
uniform mat4 al_ModelViewMatrix;
uniform mat4 al_ProjectionMatrix;
uniform float dataScale;
uniform float markerScale;
uniform float layerSeparation;
uniform float is_line;
uniform float is_omni;
uniform float eye_sep;
uniform float foc_len;

// Region Plane information
uniform vec3 plane_normal = vec3(0, 0, -1);
uniform vec3 plane_point = vec3(0.5, 0.5, 0.5);
uniform float second_plane_distance = 3.0;

//uniform float far_clip;
//uniform float near_clip;
uniform float clipped_mult;
layout (location = 0) in vec4 position;
layout (location = 1) in vec4 offset; // 4th component w is used for color
out vec4 color;

// http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl

vec3 rgb2hsv(vec3 c)
{
vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
vec4 p = c.g < c.b ? vec4(c.bg, K.wz) : vec4(c.gb, K.xy);
vec4 q = c.r < p.x ? vec4(p.xyw, c.r) : vec4(c.r, p.yzx);
float d = q.x - min(q.w, q.y);
float e = 1.0e-10;
return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

vec3 hsv2rgb(vec3 c)
{
vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec4 stereo_displace(vec4 v, float e, float r) {
vec3 OE = vec3(-v.z, 0.0, v.x); // eye direction, orthogonal to vertex vector
OE = normalize(OE);             // but preserving +y up-vector
OE *= e;               // set mag to eye separation
vec3 EV = v.xyz - OE;  // eye to vertex
float ev = length(EV); // save length
EV /= ev;              // normalize
// coefs for polynomial t^2 + 2bt + c = 0
// derived from cosine law r^2 = t^2 + e^2 + 2tecos(theta)
// where theta is angle between OE and EV
// t is distance to sphere surface from eye
float b = -dot(OE, EV);         // multiply -1 to dot product because
// OE needs to be flipped in direction
float c = e * e - r * r;
float t = -b + sqrt(b * b - c); // quadratic formula
v.xyz = OE + t * EV;            // direction from origin to sphere surface
v.xyz = ev * normalize(v.xyz);  // normalize and set mag to eye-to-v distance
return v;
}

bool between_planes(vec3 point) {
vec3 difference = point - plane_point;
float proj = dot(plane_normal, difference);
if (proj >= 0 &&  proj <= second_plane_distance) {
return true;
} else {
return false;
}
}

const float PI_2 = 1.57079632679489661923;

void main()
{
vec4 mesh_center = offset * vec4(1.0, 1.0, 1.0 + layerSeparation, 1.0);
mesh_center.w = 1.0;

float colormult = 1.0;
if (!between_planes(offset.xyz)) {
if (is_line > 0.5) {
colormult = clipped_mult * 0.3;
} else {
colormult = clipped_mult;
}
}

float local_scale = 1.0;

float reflectivity = (0.8 + (0.2 - pow(sin((position.z + 1.0)/2), 0.5)  * 0.2/ sin(1.0)));

if (is_line > 0.5) {
local_scale = 1.03;
color = vec4(hsv2rgb(vec3(offset.w, 1.0, 1.0)), 1.0)* colormult * reflectivity;
} else {
color = vec4(hsv2rgb(vec3(offset.w, 1.0, 0.85)), 1.0)* colormult* reflectivity;
}

vec4 p = vec4(local_scale * markerScale * position.xyz, 0.0) + (mesh_center * dataScale);
if (is_omni > 0.5) {
gl_Position = al_ProjectionMatrix * stereo_displace(al_ModelViewMatrix * p, eye_sep, foc_len);
}
else {
gl_Position = al_ProjectionMatrix * al_ModelViewMatrix * p;
}
}
)";

const string instancing_frag = R"(
                               #version 330
                               //uniform vec4 color;
                               in vec4 color;
                               layout (location = 0) out vec4 frag_out0;
                               void main()
                               {
                               //GLfloat[2] range = glGet(GL_DEPTH_RANGE);
                               frag_out0 = color;
                               }
                               )";

class AtomProperties {
 public:
  string name;
  float drawScale;
  Color color;
  //  gl::PolygonMode polygonMode{Mesh::FILL};
};

struct ElementData {
  float radius;
  Color color;
};

class DataDisplayParameters {
 public:
  ParameterVec3 layerDir{"LayerDir", ""};  // Direction of layers in data

  ParameterMenu mAtomOfInterest{"AtomOfInterest", "", 0};
  ParameterChoice mShowAtoms{"ShowAtoms"};
  Parameter mAtomMarkerSize{"AtomMarkerSize", "", 0.4, "", 0.0, 5.0};
  ParameterBool mShowRadius{"ShowAtomRadius", "", 1};
  ParameterBool mCumulativeTrajectory{"ShowCumulativeTrajectory", "", 1.0};
  ParameterBool mIndividualTrajectory{"ShowIndividualTrajectory", "", 1.0};

  ParameterColor backgroundColor{"projectionBackground", "",
                                 Color(0.1f, 0.1f, 0.1f, 0.8f)};

  ParameterBool mVisible{"Visible", "", 1};
  ParameterBool mShowGraph{"ShowGraph", "", 1};
  ParameterBool mShowParallel{"ShowParallel", "", 1};
  ParameterBool mShowPerspective{"ShowPerspective", "", 1};
  ParameterBool mSingleProjection{"SingleProjection", "", 1};

  ParameterBool mBillboarding{"Billboarding", "", 1};
  ParameterBool mSmallLabel{"SmallLabel", "", 1};
  ParameterBool mDrawLabels{"DrawLabels", "", 1};

  //    ParameterBool mMultiplot {"Multiplot", "", 0.0};
  ParameterMenu mPlotYAxis{"PlotYAxis"};
  ParameterMenu mPlotXAxis{"PlotXAxis"};

  Parameter mLayerSeparation{
      "LayerSeparation",
      "",
      0,
      "",
      0,
      3};  // Increase layer separation (Z- axis scaling) in perspectiveView
  Parameter mLayerScaling{"LayerScaling",
                          "",
                          1.0,
                          "",
                          0,
                          3};  // Increase layer size in projection view

  ParameterVec3 mSlicingPlanePoint{"SlicingPlanePoint", "",
                                   Vec3f(0.0f, 0.0, 0.0)};
  ParameterVec3 mSlicingPlaneNormal{"SliceNormal", "", Vec3f(0.0f, 0.0f, 1.0)};

  Parameter mSliceRotationPitch{
      "SliceRotationPitch", "SliceAngles", 0.0, "", -M_PI, M_PI};
  Parameter mSliceRotationRoll{"SliceRotationRoll", "SliceAngles", 0.0, "",
                               -M_PI / 2.0,         M_PI / 2.0};
  Parameter mSlicingPlaneThickness{
      "SlicingPlaneThickness", "", 3.0, "", 0.0f, 30.0f};

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
  vector<AtomProperties> atomPropertiesProj;

  std::map<std::string, ElementData> elementData;

  std::string metaText;

  DatasetManager mDatasetManager;

  bool mRunComputation{true};

  // #INSTANCED_RENDERING: declare
  ShaderProgram instancing_shader;
  VAOMesh orthoMesh0;
  InstancingMesh instancing_mesh0;
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

  void init();

  void initRootDirectory();

  void updateText();
  // Prepare elements before draw call
  void prepare(Graphics &g, Matrix4f &transformMatrix);

  // Draw function
  void draw(Graphics &g);

  void setFont(std::string name, float size);

  // Functions to change increase/decrease parameters according to internal
  // parameter spaces

  void nextTemp() {
    mDatasetManager.mParameterSpaces["temperature"]->stepIncrement();
  }

  void previousTemp() {
    mDatasetManager.mParameterSpaces["temperature"]->stepDecrease();
  }

  void nextChempot() {
    mDatasetManager.mParameterSpaces["chempotA"]->stepIncrement();
  }

  void previousChempot() {
    mDatasetManager.mParameterSpaces["chempotA"]->stepDecrease();
  }

  void nextChempot2() {
    mDatasetManager.mParameterSpaces["chempotB"]->stepIncrement();
  }

  void previousChempot2() {
    mDatasetManager.mParameterSpaces["chempotB"]->stepDecrease();
  }

  void nextTime() { mDatasetManager.mParameterSpaces["time"]->stepIncrement(); }

  void previousTime() {
    mDatasetManager.mParameterSpaces["time"]->stepDecrease();
  }

  void nextLayer() {
    mSlicingPlanePoint =
        mSlicingPlanePoint.get() +
        mSlicingPlaneNormal.get().normalized() * mSlicingPlaneThickness;
  }

  void previousLayer() {
    mSlicingPlanePoint =
        mSlicingPlanePoint.get() -
        mSlicingPlaneNormal.get().normalized() * mSlicingPlaneThickness;
  }

  //  void requestDataLoad() {
  //    mRequestLoad = true;
  //  }

  //  void requestInitDataset() { mRequestInit = true; }

  void dumpImages(std::string dumpPrefix);

  void computeSlicing() {
    if (mRunComputation) {
      mSlicingPlaneThickness =
          findDistanceNormal(mAligned4fData, layerDir.get());
      std::cout << "Data Boundaries:" << std::endl;
      std::cout << "X " << mDataBoundaries.minx << " -> "
                << mDataBoundaries.maxx << std::endl;
      std::cout << "Y " << mDataBoundaries.miny << " -> "
                << mDataBoundaries.maxy << std::endl;
      std::cout << "Z " << mDataBoundaries.minz << " -> "
                << mDataBoundaries.maxz << std::endl;
      Vec3f point = mSlicingPlanePoint;
      point.z = mDataBoundaries.minz;
      mSlicingPlanePoint = point;
    }
  }

  void resetSlicing() {
    mSlicingPlanePoint.set({0, 0, mDataBoundaries.minz});

    mSlicingPlaneThickness = mDataBoundaries.maxz - mDataBoundaries.minz;
    mSliceRotationRoll.set(0);
    mSliceRotationPitch.set(0);
    //      std::cout << mSlicingPlaneThickness.get() <<std::endl;
  }

 protected:
  Texture &iso_scene() { return fbo_iso.tex(); }

  // This function should be called whenever there is new atom position data
  void updateDisplayBuffers();

  void prepareParallelProjection(Graphics &g, Matrix4f &transformMatrix);

  void drawHistory(Graphics &g) {
    g.pushMatrix();
    g.meshColor();
    gl::polygonFill();
    if (mIndividualTrajectory.get() != 0.0f) {
      g.draw(mHistoryMesh);
    }
    if (mCumulativeTrajectory.get() != 0.0f) {
      g.translate((mDataBoundaries.maxx - mDataBoundaries.minx) / 2.0f,
                  (mDataBoundaries.maxy - mDataBoundaries.miny) / 2.0f,
                  (mDataBoundaries.maxz - mDataBoundaries.minz) / 2.0f);
      g.draw(mTrajectoryMesh);
    }
    g.popMatrix();
  }

  void drawPerspective(Graphics &g);

  void drawParallelProjection(Graphics &g);

  void drawGraph(Graphics &g);

  void updateParameterText() {
    if (mNeedsProcessing.load()) {
      mParamText = "Processing ...";
    } else {
    }
  }

  //    void drawLabel(Graphics &g) {

  //        std::string title = mDatasetManager.currentDataset();
  //        std::string text = getParameterText();

  //        g.depthTesting(false);
  //        g.blending(true);
  //        g.blendModeAdd();
  //        g.pushViewMatrix(Matrix4f::identity());

  //        g.pushMatrix();
  //        if (mSmallLabel) {
  //            g.translate(0.3, 1.45, -6.0);
  //            g.scale(0.2);
  //        } else {
  //            g.translate(-0.6, 1.2, -6.0);
  //            g.scale(0.3);
  //        }

  //        mLabelFont.render(g, title);
  //        g.popMatrix();

  //        g.pushMatrix();
  //        if (mSmallLabel) {
  //            g.translate(0.3, 1.3, -6.0);
  //            g.scale(0.2);
  //        } else {
  //            g.translate(-0.6, 1.0, -6.0);
  //            g.scale(0.3);
  //        }
  //        mLabelFont.render(g, text);
  //        g.popMatrix();

  //        g.blending(false);
  //        g.popViewMatrix();
  //    }

 private:
  VAOMesh axis;
  VAOMesh orthoMesh;
  VAOMesh graphlinesMesh;
  VAOMesh boxMesh;
  VAOMesh mGridMesh;

  // Graph
  //  std::mutex mGraphTextureLock;
  Texture mGraphTexture;
  //  std::string mGraphFilePathToLoad;
  //  std::string mGraphFilePath;

  Font mLabelFont;

  EasyFBO fbo_iso;

  //  bool mRequestLoad {false};
  std::atomic<bool> mNeedsProcessing{false};
  bool mRequestInit{false};

  typedef struct {
    int counts;
    float radius;
    std::string species;
  } AtomData;

  vector<AtomData> mAtomData;
  std::vector<float> mAligned4fData;
  BoundingBox_ mDataBoundaries;

  float mMarkerScale;  // Global marker scaling factor
};

#endif  // DATASETDISPLAY_HPP
