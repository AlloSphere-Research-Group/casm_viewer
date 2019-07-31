#ifndef DATASETDISPLAY_HPP
#define DATASETDISPLAY_HPP

#ifdef AL_WINDOWS
#define NOMINMAX
#endif

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <regex>
#include <algorithm>
#include <memory>
#include <mutex>

#include "al/core.hpp"
#include "al/core/math/al_Matrix4.hpp"

#include "al/util/ui/al_Parameter.hpp"
#include "al/util/ui/al_ParameterBundle.hpp"
#include "al/util/ui/al_Pickable.hpp"
#include "al/util/ui/al_PickableRotateHandle.hpp"
#include "al/util/ui/al_PickableManager.hpp"
#include "al/core/io/al_CSVReader.hpp"
#include "al/core/io/al_File.hpp"
#include "al/core/types/al_Color.hpp"
#include "al/core/graphics/al_Light.hpp"
#include "al/util/al_Toml.hpp"
#include "al/util/al_FontModule.hpp"

#include "module/img/loadImage.hpp"

#include "al_VASPReader.hpp"
#include "al_DataScript.hpp"

#include "processors.hpp"
#include "instanced_mesh.hpp"

#include "slice.hpp"
#include "datasetmanager.hpp"

//#include "imgui.h"

#undef far
#undef near

using namespace al;
using namespace std;

inline Matrix4f getLookAt(const Vec3f& ux, const Vec3f& uy, const Vec3f& uz, const Vec3f& p) {
  return Matrix4f(
        ux[0], ux[1], ux[2], -(ux.dot(p)),
      uy[0], uy[1], uy[2], -(uy.dot(p)),
      uz[0], uz[1], uz[2], -(uz.dot(p)),
      0,     0,     0,       1
      );
}

inline Matrix4f getLookAt(const Vec3f& eyePos, const Vec3f& at, const Vec3f& up) {
  Vec3f z = (eyePos - at).normalize();
  Vec3f x = cross(up, z).normalize();
  Vec3f y = cross(z, x).normalize();
  return getLookAt(x, y, z, eyePos);
}

inline HSV rgb2hsv(RGB c)
{
  Vec4f K = Vec4f(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
  Vec4f p = c.g < c.b ? Vec4f(c.b, c.g, K.w, K.z) : Vec4f(c.g, c.b, K.x, K.y);
  Vec4f q = c.r < p.x ? Vec4f(p.x, p.y, p.w, c.r) : Vec4f(c.r, p.y, p.z, p.x);
  float d = q.x - min(q.w, q.y);
  float e = 1.0e-10;
  return HSV(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

const std::string instancing_vert = R"(
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
  Graphics::PolygonMode polygonMode {Graphics::FILL};
};


struct ElementData {
  float radius;
  Color color;
};


class DataDisplayParameters {
public:

  ParameterVec3 layerDir {"LayerDir", ""}; // Direction of layers in data

  ParameterMenu mAtomOfInterest{ "AtomOfInterest", "", 0 };
  ParameterChoice mShowAtoms{ "ShowAtoms" };
  Parameter mAtomMarkerSize {"AtomMarkerSize", "", 0.4, "", 0.0, 5.0};
  ParameterBool mShowRadius {"ShowAtomRadius", "", 1};
  ParameterBool mCumulativeTrajectory {"ShowCumulativeTrajectory", "", 1.0};
  ParameterBool mIndividualTrajectory {"ShowIndividualTrajectory", "", 1.0};

  ParameterColor backgroundColor{"projectionBackground", "", Color(0.1f, 0.1f, 0.1f, 0.8f) };

  ParameterBool mVisible {"Visible", "", 1};
  ParameterBool mShowGraph {"ShowGraph", "", 1};
  ParameterBool mShowParallel {"ShowParallel", "", 1};
  ParameterBool mShowPerspective {"ShowPerspective", "", 1};
  ParameterBool mSingleProjection {"SingleProjection", "", 1};

  ParameterBool mBillboarding {"Billboarding", "", 1};
  ParameterBool mSmallLabel {"SmallLabel", "", 1};
  ParameterBool mDrawLabels {"DrawLabels", "", 1};

  //    ParameterBool mMultiplot {"Multiplot", "", 0.0};
  ParameterMenu mPlotYAxis{ "PlotYAxis" };
  ParameterMenu mPlotXAxis{ "PlotXAxis" };

  Parameter mLayerSeparation{"LayerSeparation", "", 0, "", 0, 3}; // Increase layer separation (Z- axis scaling) in perspectiveView
  Parameter mLayerScaling{"LayerScaling", "", 1.0, "", 0, 3}; // Increase layer size in projection view

  ParameterVec3 mSlicingPlanePoint{ "SlicingPlanePoint", "", Vec3f(0.0f, 0.0, 0.0)};
  ParameterVec3 mSlicingPlaneNormal {"SliceNormal", "", Vec3f(0.0f, 0.0f, 1.0)};

  Parameter mSliceRotationPitch{"SliceRotationPitch", "SliceAngles", 0.0, "", - M_PI, M_PI};
  Parameter mSliceRotationRoll{"SliceRotationRoll", "SliceAngles", 0.0, "", - M_PI/2.0, M_PI/2.0};
  Parameter mSlicingPlaneDistance{ "SlicingPlaneThickness", "", 3.0, "", 0.0f, 30.0f };

  Parameter mPerspectiveRotY {"perspectiveRotY", "", -75, "", -90, 90};

  ParameterBool mShowGrid {"ShowGrid", "", 0.0};
  Parameter mGridSpacing {"GridSpacing", "", 10.0, "", 1.0, 20.0};
  ParameterMenu mGridType {"GridType"};

  Parameter mGridXOffset {"GridXOffset", "", 0.0, "", 0.0, 5.0};
  Parameter mGridYOffset {"GridYOffset", "", 0.0, "", 0.0, 5.0};

  ParameterBundle bundle{"CASMDataset"};
};

// -----------------------------------------------------------------

class DataDisplay : public DataDisplayParameters{
public:

  vector<AtomProperties> atomPropertiesProj;

  std::map<std::string, ElementData> elementData;

  DatasetManager mDatasetManager;

  bool mRunComputation {true};

  // #INSTANCED_RENDERING: declare
  ShaderProgram instancing_shader;
  VAOMesh orthoMesh0;
  InstancingMesh instancing_mesh0;
  Mesh mHistoryMesh;
  Mesh mTrajectoryMesh;

  // Pickables
  PickableManager mPickableManager;
  PickableBB graphPickable {"graph"};
  PickableBB parallelPickable{"parallel"};
  PickableBB perspectivePickable{"perspective"};
  PickableBB slicePickable{"slice"};
  // PickableRotateHandle rh;

  std::string mParamText;

  std::mutex mDrawLock;

  void init();

  void initRootDirectory();

  // Prepare elements before draw call
  void prepare(Graphics &g, Matrix4f &transformMatrix);

  // Draw function
  void draw(Graphics &g);

  void setFont(std::string name, float size);

  // Functions to change increase/decrease parameters according to internal parameter spaces

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

  void nextTime() {
    mDatasetManager.mParameterSpaces["time"]->stepIncrement();
  }

  void previousTime() {
    mDatasetManager.mParameterSpaces["time"]->stepDecrease();
  }

  void nextLayer() {
    mSlicingPlanePoint = mSlicingPlanePoint.get() + mSlicingPlaneNormal.get().normalized() * mSlicingPlaneDistance;
  }

  void previousLayer() {
    mSlicingPlanePoint = mSlicingPlanePoint.get() - mSlicingPlaneNormal.get().normalized() * mSlicingPlaneDistance;
  }

//  void requestDataLoad() {
//    mRequestLoad = true;
//  }

  void requestInitDataset() {
    mRequestInit = true;
  }

  void dumpImages(std::string dumpPrefix);

  void computeSlicing() {

    if (mRunComputation) {
      mSlicingPlaneDistance = findDistanceNormal(mAligned4fData, layerDir.get());
      std::cout << "Data Boundaries:" << std::endl;
      std::cout << "X " << mDataBoundaries.minx << " -> " << mDataBoundaries.maxx << std::endl;
      std::cout << "Y " << mDataBoundaries.miny << " -> " << mDataBoundaries.maxy << std::endl;
      std::cout << "Z " << mDataBoundaries.minz << " -> " << mDataBoundaries.maxz << std::endl;
      Vec3f point = mSlicingPlanePoint;
      point.z = mDataBoundaries.minz;
      mSlicingPlanePoint = point;
    }
  }

  void resetSlicing() {
    mSlicingPlanePoint.set({mDataBoundaries.minx, mDataBoundaries.miny,
                            mDataBoundaries.minz});

    mSlicingPlaneDistance = mDataBoundaries.maxz - mDataBoundaries.minz;
    mSliceRotationRoll.set(0);
    mSliceRotationPitch.set(0);
    //      std::cout << mSlicingPlaneDistance.get() <<std::endl;
  }

protected:

  Texture& iso_scene() { return fbo_iso.tex(); }


  // This function should be called whenever there is new atom position data
  void updateDisplayBuffers() {
    if (mDatasetManager.positionBuffers.newDataAvailable()) {

      std::map<string, int> elementCounts;

      auto allPositions = mDatasetManager.positionBuffers.get();

      mDatasetManager.mCurrentLoadedIndeces.clear();
      for (auto &paramSpace: mDatasetManager.mParameterSpaces) {
        mDatasetManager.mCurrentLoadedIndeces[paramSpace.first] = paramSpace.second->getCurrentIndex();
      }

      vector<vector<float> *> elemPositions;

      // Now that the layer direction has been computed from the atom of interest, we need to
      // generate aligned data from the visible atoms.
      elemPositions.clear();
      auto visibleAtoms = mShowAtoms.getSelectedElements();
      for (auto &elementData: *allPositions) {
        if (std::find(visibleAtoms.begin(), visibleAtoms.end(), elementData.first) != visibleAtoms.end()) {
          elemPositions.push_back(&(elementData.second));
        }
      }
      // also prepare layer normal direction aligned data
      size_t totalSize = 0;
      for(auto *elems:elemPositions) {
        totalSize += elems->size();
      }

      mDataBoundaries.reset();
      mAligned4fData.resize(totalSize*4);
      auto outit = mAligned4fData.begin();
      // now fill the
      for(auto *elems:elemPositions) {
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

      if (mAligned4fData.size() > 0 ) {
        auto &b = mDataBoundaries;
        perspectivePickable.bb.set(Vec3f(b.minx,b.miny,b.minz),Vec3f(b.maxx,b.maxy,b.maxz));
        slicePickable.bb.set(Vec3f(b.minx,b.miny,b.minz),Vec3f(b.maxx,b.maxy,(b.maxz-b.minz)*0.5f));
        // rh.pose.pos().set(perspectivePickable.bb.cen);
        mSlicingPlanePoint.setHint("maxx", b.maxx);
        mSlicingPlanePoint.setHint("minx", b.minx- (b.maxx));
        mSlicingPlanePoint.setHint("maxy", b.maxy);
        mSlicingPlanePoint.setHint("miny", b.miny- (b.maxy));
        mSlicingPlanePoint.setHint("maxz", b.maxz);
        mSlicingPlanePoint.setHint("minz", b.minz- (b.maxz));
        mSlicingPlaneDistance.min(mDataBoundaries.minz);
        mSlicingPlaneDistance.max(mDataBoundaries.maxz);
      }

      // Set active atoms and colors
      atomPropertiesProj.clear();
      vector<AtomProperties> atomPropertiesPersp;

      // TODO these colors should be exposed as a preference
      vector<Color> colorList = {
        Color(0.0, 1.0, 1.0, 1.0), Color(1.0, 1.0, 0.0, 1.0), Color(1.0, 0.0, 1.0, 1.0),
        Color(0.0, 1.0, 1.0, 1.0), Color(1.0, 1.0, 0.0, 1.0), Color(1.0, 0.0, 1.0, 1.0)
      };

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
      for (auto atom: mShowAtoms.getElements()) {
        if (std::find(selectedElements.begin(), selectedElements.end(), atom) != selectedElements.end()) {

          if (elementData.find(atom) != elementData.end()) {
            //                    std::cout << "Color for: " << atom << ":" << elementData[atom].color.r << " " << elementData[atom].color.g << " "<< elementData[atom].color.b << " "  <<std::endl ;
            // Atom was matched in elements.ini file, so use those colors

            atomPropertiesProj.push_back(AtomProperties{ atom, elementData[atom].radius, elementData[atom].color});

            atomPropertiesPersp.push_back(AtomProperties{ atom, elementData[atom].radius, elementData[atom].color, Graphics::FILL });
            atomPropertiesPersp.push_back(AtomProperties{ atom, elementData[atom].radius, elementData[atom].color, Graphics::LINE });
          } else { // Use defaults
            atomPropertiesProj.push_back(AtomProperties{ atom, 1.0f , *colorListIt});

            atomPropertiesPersp.push_back(AtomProperties{ atom, 1.0f, *colorList2It, Graphics::FILL });
            atomPropertiesPersp.push_back(AtomProperties{ atom, 1.0f, *colorList2It, Graphics::LINE });
          }
        }
        colorListIt++;
        colorList2It++;
        colorList2It++;
      }

      // Apply colors to aligned data -----------
      list<Color> colors;
      mAtomData.clear();
      for (auto elem: *allPositions) {
        elementCounts[elem.first] = elem.second.size()/4;
      }
      for (auto atomProps: atomPropertiesProj) {
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
        for (size_t i = 0; i < mAligned4fData.size() /4 ; i++) {
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
        mDatasetManager.generateGraph(xData, yData,  mDatasetManager.mCurrentDataset.get(), false);
      }
    }
  }

  void prepareParallelProjection(Graphics &g, Matrix4f &transformMatrix) {
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
    float centerX = (mDataBoundaries.maxx + mDataBoundaries.minx)/ 2.0f;
    float centerY = (mDataBoundaries.maxy + mDataBoundaries.miny)/ 2.0f;
    float rangeSizeX = mDataBoundaries.maxx - mDataBoundaries.minx;
    float rangeSizeY = mDataBoundaries.maxy - mDataBoundaries.miny;
    float maxrange = std::max(rangeSizeX, rangeSizeY);
    float padding = maxrange*0.05f;
    float left = centerX- 0.5f*maxrange* mLayerScaling - padding;
    float right = centerX+ 0.5f*maxrange* mLayerScaling + padding;
    float top = centerY+ 0.5f*maxrange* mLayerScaling + padding;
    float bottom = centerY- 0.5f*maxrange* mLayerScaling - padding;
    //        const float near = mDataBoundaries.minz + (mNearClip * rangeSizeZ);
    //        const float farClip = mDataBoundaries.minz + (mFarClip * rangeSizeZ);
    //        float minxy = std::min(mDataBoundaries.minx, mDataBoundaries.miny);
    g.projMatrix(Matrix4f::ortho(ar *left, ar* right,
                                 bottom, top,
                                 mDataBoundaries.minz - 100 , mDataBoundaries.maxz + 100));

    bool mAlignData = true;
    double scalingFactor = 1.0/(mDataBoundaries.maxy - mDataBoundaries.miny);

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
                             mSlicingPlanePoint.get() - mSlicingPlaneNormal.get().normalized(),
      { 0.0f, 1.0f, 0.0f }));

      g.blending(false);
      g.depthTesting(true);
      g.meshColor();
      g.draw(axis);


      if (mShowGrid == 1.0f) {
        mGridMesh.reset();
        addRect(mGridMesh,0,-2.0, 0.1f/(maxrange* mLayerScaling), 4.0f);
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
            g.scale(maxrange* mLayerScaling);
            g.rotate(90, 0,0,1);
            g.draw(mGridMesh);
            g.popMatrix();
            g.pushMatrix();
            g.translate(mGridXOffset + mGridSpacing * i, mGridYOffset, 0);
            g.scale(maxrange* mLayerScaling);
            g.draw(mGridMesh);
            g.popMatrix();
          } else if (mGridType.getCurrent() == "triangle") {
            float spacing = 0.86602540378f * mGridSpacing; // sin(60)
            g.pushMatrix();
            g.translate(mGridXOffset, mGridYOffset + spacing * i, 0);
            g.scale(maxrange* mLayerScaling);
            g.rotate(90, 0,0,1);
            g.draw(mGridMesh);
            g.popMatrix();
            g.pushMatrix();
            g.translate(mGridXOffset, mGridYOffset, 0);
            g.rotate(30, 0,0,1);
            g.translate(spacing * i, 0, 0);
            g.scale(maxrange* mLayerScaling);
            g.draw(mGridMesh);
            g.popMatrix();
            g.pushMatrix();
            g.translate(mGridXOffset, mGridYOffset, 0);
            g.rotate(-30, 0,0,1);
            g.translate(spacing * i, 0, 0);
            g.scale(maxrange* mLayerScaling);
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
        for(auto &data: mAtomData) {

          if (mAlignData) {
            //            instancing_mesh0.attrib_data(
            //              mAligned4fData.size() * sizeof(float),
            //              mAligned4fData.data(),
            //              mAligned4fData.size()/4
            //            );
            int count = data.counts;
            assert((int) mAligned4fData.size() >= (cumulativeCount + count) * 4 );
            instancing_mesh0.attrib_data(
                  count * 4 * sizeof(float),
                  mAligned4fData.data() + (cumulativeCount * 4),
                  count
                  );
            cumulativeCount += count;
            //                        std::cout << "Drawing " << counts << " of " << cumulativeCount << std::endl;
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
            g.shader().uniform("markerScale", data.radius * mAtomMarkerSize * mMarkerScale);
          } else {
            g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale);
          }
          g.shader().uniform("is_line", 0.0f);
          g.shader().uniform("is_omni", 0.0f);

          g.shader().uniform("plane_point", mSlicingPlanePoint.get());
          g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
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
    }
    else {
      // (from x, y up), (right bottom)
      {
        g.pushViewport(w, 0, w, h);
        g.scissor(w, 0, w, h);
        g.clear(backgroundColor);

        g.pushViewMatrix();
        g.viewMatrix(getLookAt({ 4, 0, 0 }, { 0, 0, 0 }, { 0, 1, 0 }));

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
          g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
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
        g.viewMatrix(getLookAt({ 0, 4, 0 }, { 0, 0, 0 }, { 0, 0, -1 }));

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
          g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
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
        g.viewMatrix(getLookAt({ 0, 0, 4 }, { 0, 0, 0 }, { 0, 1, 0 }));

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
          g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
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
        g.viewMatrix(getLookAt({ 3, 3, 3 }, { 0, 0, 0 }, { 0, 1, 0 }));

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
          g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
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

  void drawHistory(Graphics &g) {
    g.pushMatrix();
    g.meshColor();
    g.polygonMode(Graphics::FILL);
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

  void drawPerspective(Graphics &g) {
    //        Vec3d normVector = reader.getNormalizingVector();
    //        Vec3d centeringVector = reader.getCenteringVector();
    if (mAligned4fData.size() == 0) {
      return; // No data has been loaded
    }

    g.lens().far(2000);
    perspectivePickable.pushMatrix(g);

    // move to center
    // g.translate(-(mDataBoundaries.maxx - mDataBoundaries.minx)/2,-(mDataBoundaries.maxy - mDataBoundaries.miny)/2,-(mDataBoundaries.maxz - mDataBoundaries.minz)/2);
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

    //        float near = mDataBoundaries.minz + (mNearClip * (mDataBoundaries.maxz- mDataBoundaries.minz));
    //        float farClip = mDataBoundaries.minz + (mFarClip * (mDataBoundaries.maxz- mDataBoundaries.minz));

    int cumulativeCount = 0;
    // now draw data with custom shaderg.shader(instancing_mesh0.shader);
    g.shader(instancing_mesh0.shader);
    g.shader().uniform("dataScale", 1.0f/((mDataBoundaries.maxy - mDataBoundaries.miny)* perspectivePickable.scale));
    g.shader().uniform("layerSeparation", mLayerSeparation);
    g.shader().uniform("is_omni", 1.0f);
    g.shader().uniform("eye_sep", perspectivePickable.scale * g.lens().eyeSep() * g.eye() / 2.0f);
    // g.shader().uniform("eye_sep", g.lens().eyeSep() * g.eye() / 2.0f);
    g.shader().uniform("foc_len", g.lens().focalLength());

    g.shader().uniform("plane_point", mSlicingPlanePoint.get());
    g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
    g.shader().uniform("second_plane_distance", mSlicingPlaneDistance);
    //        g.shader().uniform("near_clip", near);
    //        g.shader().uniform("far_clip", farClip);
    g.shader().uniform("clipped_mult", 0.45);
    g.update();

    for(auto data: mAtomData) {
      if (mShowRadius == 1.0f) {
        g.shader().uniform("markerScale", data.radius * mAtomMarkerSize * mMarkerScale / perspectivePickable.scale);
        //                std::cout << data.radius << std::endl;
      } else {
        g.shader().uniform("markerScale", mAtomMarkerSize * mMarkerScale / perspectivePickable.scale);
      }
      int count = data.counts;
      assert((int) mAligned4fData.size() >= (cumulativeCount + count) * 4 );
      instancing_mesh0.attrib_data(
            count * 4 * sizeof(float),
            mAligned4fData.data() + (cumulativeCount * 4),
            count
            );
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
    //const float z[2] = {mNearClip.get(), mFarClip.get()};
    // double scalingFactor = 1.0/(mDataBoundaries.maxz - mDataBoundaries.minz);

    // Plane equation Ax + By + Cz - D = 0
    //        float D =  -mSlicingPlanePoint.get().dot(mSlicingPlaneNormal);
    // Line Equation

    const float x[2] = { mDataBoundaries.maxx, mDataBoundaries.minx };
    const float y[2] = { mDataBoundaries.maxy, mDataBoundaries.miny };

    const float z[2] = {
      0.0,
      (mSlicingPlaneDistance)
    };

    for(int k=0; k<=1; k++){
      for(int j=0; j<=1; j++){
        for(int i=0; i<=1; i++){
          boxMesh.vertex(x[i], y[j], z[k]) ;
        }}}

    static const int I[] = {
      0,1, 2,3, 4,5, 6,7,
      0,2, 1,3, 4,6, 5,7,
      0,4, 1,5, 2,6, 3,7
    };
    boxMesh.index(I, sizeof(I)/sizeof(*I), 0);
    boxMesh.update();
    g.shader().uniform("eye_sep", perspectivePickable.scale * g.lens().eyeSep() * g.eye() / 2.0f);

    g.shader().uniform("eye_sep", /*perspectivePickable.scale * */g.lens().eyeSep() * g.eye() / 2.0f);

    glLineWidth(5);
    g.polygonMode(Graphics::LINE);
    g.color(0.8f, 0.8f, 1.0f, 0.9f);
    g.scale(1.0, 1.0, 1.0f + mLayerSeparation);
    g.translate(mSlicingPlanePoint.get());
    g.rotate(mSliceRotationPitch * 360.0f/(M_2PI), 1.0, 0.0, 0.0);
    g.rotate(mSliceRotationRoll * 360.0f/(M_2PI), 0.0, -1.0, 0.0);
    g.draw(boxMesh);
    g.polygonMode(Graphics::FILL);
    glLineWidth(1);

    //  -----------

    g.popMatrix();

    drawHistory(g);

    g.popMatrix(); //pickable

    // Draw label

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
      mLabelFont.render(g, mDatasetManager.currentDataset() + " " + mParamText);
    }
    g.popMatrix();

    perspectivePickable.drawBB(g);
    //        perspectivePickable.drawChildren(g);
    // if(perspectivePickable.hover || rh.hover[0] || rh.hover[1] || rh.hover[2] ){
    //     g.depthTesting(false);
    //     perspectivePickable.drawChildren(g);
    //     g.depthTesting(true);
    // }

  }

  void drawParallelProjection(Graphics &g) {
    float pixel_w = float(fbo_iso.width());
    float pixel_h = float(fbo_iso.height());
    float w = 2;
    float h = w * pixel_h / pixel_w;
    Mesh m;
    addTexQuad(m, w, h);
    iso_scene().filter(Texture::LINEAR_MIPMAP_LINEAR);
    iso_scene().generateMipmap(); //XXX this works.. confusing api needs help..
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
      mLabelFont.render(g, mDatasetManager.currentDataset() + " " + mParamText);
    }
    g.popMatrix();
    parallelPickable.drawBB(g);
  }

  void drawGraph(Graphics &g) {
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

        auto imageData = imgModule::loadImage(mGraphFilePathToLoad.c_str());
        if (imageData.data.size() == 0) {
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
        }
        else {
          //                         cout << "loaded image size: " << imageData.width << ", " << imageData.height << endl;

          mGraphTexture.resize(imageData.width, imageData.height);
          mGraphTexture.submit(imageData.data.data(), GL_RGBA, GL_UNSIGNED_BYTE);

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
      mLabelFont.render(g, mDatasetManager.currentDataset() + " " + mParamText);
    }
    g.popMatrix();
    graphPickable.drawBB(g);
  }

  void loadGraphTexture(std::string graphFilePath) {
    std::unique_lock<std::mutex> lk(mGraphTextureLock);
    mGraphFilePathToLoad = graphFilePath;
  }

  void updateParameterText() {
    mParamText = "";
    if (mProcessing.load()) {
      mParamText = "Processing ...";
    } else {
      auto subDir = mDatasetManager.getSubDir();
      auto temperatureId = mDatasetManager.mParameterSpaces["temperature"]->getCurrentId();
      std::string timeId;
      if (mDatasetManager.mParameterSpaces["time"]->size() > 0) {
        timeId = mDatasetManager.mParameterSpaces["time"]->getCurrentId();
      }
      if (subDir.size() > 0) {
        mParamText = subDir;
      }
      if (temperatureId.size()  > 0) {
        mParamText += " temp:" + temperatureId + " ";
      }
      if (timeId.size()  > 0) {
        mParamText += " time:" + timeId;
      }
      if (mParamText.size() == 0) {
        mParamText = "Dataset unavailable";
      }
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
  std::mutex mGraphTextureLock;
  Texture mGraphTexture;
  std::string mGraphFilePathToLoad;
  std::string mGraphFilePath;

  FontModule mLabelFont;

  EasyFBO fbo_iso;

//  bool mRequestLoad {false};
  std::atomic<bool> mProcessing{ false };
  bool mRequestInit {false};

  typedef struct {
    int counts;
    float radius;
  } AtomData;

  vector<AtomData> mAtomData;
  std::vector<float> mAligned4fData;
  BoundingBox_ mDataBoundaries;

  float mMarkerScale; // Global marker scaling factor

};

#endif // DATASETDISPLAY_HPP
