

#ifndef AL_VASPRENDER
#define AL_VASPRENDER
#include "al/graphics/al_Graphics.hpp"

#undef CIEXYZ
#include "al/graphics/al_Shader.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/graphics/al_VAOMesh.hpp"
#include "al/ui/al_BoundingBox.hpp"
#include "al/ui/al_Parameter.hpp"

namespace al {

// When setting per instance attribute position,
// remember that in VAOMesh
// position (vertices) uses location 0,
// color 1, texcoord 2, normal 3
//
// Avoid attribute indices that are going to be used.
// ex) if mesh to be instanced will have vertices and color,
//     set instancing attribute index to be 2
//
// Also note that VAOMesh's update method will disable attribute
// if there's no data for that index, so it is better to just set
// instance attribute's index to 4...
struct InstancingMesh {
  VAOMesh mesh;
  BufferObject buffer;
  ShaderProgram shader;
  size_t num_instances = 0;

  // attrib_loc: per instance attribute location, should be set
  //             with "layout (location=#)" in vert shader
  // attrib_num_elems: vec4? vec3?
  // attrib_type: GL_FLOAT? GL_UNSIGNED_BYTE?
  void init(const std::string &vert_str, const std::string &frag_str,
            GLuint attrib_loc, GLint attrib_num_elems, GLenum attrib_type) {
    shader.compile(vert_str, frag_str);
    buffer.bufferType(GL_ARRAY_BUFFER);
    buffer.usage(GL_DYNAMIC_DRAW);  // assumes buffer will change every frame
                                    // and will be used for drawing
    buffer.create();

    auto &v = mesh.vao();
    v.bind();
    v.enableAttrib(attrib_loc);
    // for normalizing, this code only considers GL_FLOAT AND GL_UNSIGNED_BYTE,
    // (does not normalize floats and normalizes unsigned bytes)
    v.attribPointer(
        attrib_loc, buffer, attrib_num_elems, attrib_type,
        (attrib_type == GL_FLOAT) ? GL_FALSE : GL_TRUE,  // normalize?
        0,                                               // stride
        0);                                              // offset
    glVertexAttribDivisor(attrib_loc, 1);  // step attribute once per instance
  }

  // size: size (in bytes) of whole data. if data is 10 vec4,
  //       it should be 10 * 4 * sizeof(float)
  // count: number of instances to draw with this data
  void attrib_data(size_t size, const void *data, size_t count) {
    buffer.bind();
    buffer.data(size, data);
    num_instances = count;
  }

  // must set shader before calling draw
  // g.shader(instancing_mesh.shader);
  // g.shader().uniform("my_uniform", my_uniform_data);
  // g.update();
  // instancing_mesh.draw(instance_count);
  void draw() {
    mesh.vao().bind();
    if (mesh.indices().size()) {
      mesh.indexBuffer().bind();
      glDrawElementsInstanced(mesh.vaoWrapper->GLPrimMode,
                              mesh.indices().size(), GL_UNSIGNED_INT, 0,
                              num_instances);
    } else {
      glDrawArraysInstanced(mesh.vaoWrapper->GLPrimMode, 0,
                            mesh.vertices().size(), num_instances);
    }
  }
};

typedef struct {
  int counts;
  float radius;
  std::string species;
} AtomData;

class VASPRender {
 public:
  ParameterVec3 mSlicingPlanePoint{"SlicingPlanePoint", "",
                                   Vec3f(0.0f, 0.0, 0.0)};
  ParameterVec3 mSlicingPlaneNormal{"SliceNormal", "", Vec3f(0.0f, 0.0f, 1.0)};
  Parameter mSlicingPlaneThickness{
      "SlicingPlaneThickness", "", 3.0, "", 0.0f, 30.0f};

  Parameter mSliceRotationPitch{
      "SliceRotationPitch", "SliceAngles", 0.0, "", -M_PI, M_PI};
  Parameter mSliceRotationRoll{"SliceRotationRoll", "SliceAngles", 0.0, "",
                               -M_PI / 2.0,         M_PI / 2.0};

  Parameter mAtomMarkerSize{"AtomMarkerSize", "", 0.4, "", 0.0, 5.0};
  ParameterBool mShowRadius{"ShowAtomRadius", "", 1};

  // Increase layer separation (Z- axis scaling) in perspectiveView
  Parameter mLayerSeparation{"LayerSeparation", "", 0, "", 0, 3};

  ShaderProgram instancing_shader;
  InstancingMesh instancing_mesh0;
  float mMarkerScale;  // Global marker scaling factor

  void init() {
    addSphere(instancing_mesh0.mesh, 1, 12, 6);
    instancing_mesh0.mesh.update();
    instancing_mesh0.init(instancing_vert, instancing_frag,
                          1,          // location
                          4,          // num elements
                          GL_FLOAT);  // type

    instancing_shader.compile(instancing_vert, instancing_frag);

    mSliceRotationPitch.registerChangeCallback([this](float value) {
      mSlicingPlaneNormal.setNoCalls(Vec3f(sin(mSliceRotationRoll),
                                           cos(mSliceRotationRoll) * sin(value),
                                           cos(value))
                                         .normalize());
    });

    mSliceRotationRoll.registerChangeCallback([this](float value) {
      mSlicingPlaneNormal.setNoCalls(
          Vec3f(sin(value), cos(value) * sin(mSliceRotationPitch),
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

    mMarkerScale = 0.01f;
  }

  void setDataBoundaries(BoundingBoxData &b) {
    mSlicingPlanePoint.setHint("maxx", b.max.x);
    mSlicingPlanePoint.setHint("minx", b.min.x - (b.max.x));
    mSlicingPlanePoint.setHint("maxy", b.max.y);
    mSlicingPlanePoint.setHint("miny", b.min.y - (b.max.y));
    mSlicingPlanePoint.setHint("maxz", b.max.z);
    mSlicingPlanePoint.setHint("minz", b.min.z - (b.max.z));
    mSlicingPlaneThickness.min(0.0);
    mSlicingPlaneThickness.max(b.max.z - b.min.z);
  }

  void drawVASP(Graphics &g, float scale, std::vector<AtomData> &mAtomData,
                std::vector<float> &mAligned4fData) {
    int cumulativeCount = 0;
    // now draw data with custom shaderg.shader(instancing_mesh0.shader);
    g.shader(instancing_mesh0.shader);
    g.shader().uniform("dataScale",
                       1.0f / ((mSlicingPlanePoint.getHint("maxy") -
                                mSlicingPlanePoint.getHint("miny")) *
                               scale));
    g.shader().uniform("layerSeparation", mLayerSeparation);
    g.shader().uniform("is_omni", 1.0f);
    g.shader().uniform("eye_sep", scale * g.lens().eyeSep() * g.eye() / 2.0f);
    // g.shader().uniform("eye_sep", g.lens().eyeSep() * g.eye() / 2.0f);
    g.shader().uniform("foc_len", g.lens().focalLength());

    g.shader().uniform("plane_point", mSlicingPlanePoint.get());
    g.shader().uniform("plane_normal", mSlicingPlaneNormal.get().normalized());
    g.shader().uniform("second_plane_distance", mSlicingPlaneThickness);
    //        g.shader().uniform("near_clip", near);
    //        g.shader().uniform("far_clip", farClip);
    g.shader().uniform("clipped_mult", 0.45);
    g.update();

    for (auto data : mAtomData) {
      if (mShowRadius == 1.0f) {
        g.shader().uniform("markerScale", data.radius * mAtomMarkerSize *
                                              mMarkerScale / scale);
        //                std::cout << data.radius << std::endl;
      } else {
        g.shader().uniform("markerScale",
                           mAtomMarkerSize * mMarkerScale / scale);
      }
      int count = data.counts;
      assert((int)mAligned4fData.size() >= (cumulativeCount + count) * 4);
      instancing_mesh0.attrib_data(
          count * 4 * sizeof(float),
          mAligned4fData.data() + (cumulativeCount * 4), count);
      cumulativeCount += count;

      gl::polygonFill();
      g.shader().uniform("is_line", 0.0f);
      instancing_mesh0.draw();

      //            g.shader().uniform("is_line", 1.0f);
      //            g.polygonMode(Graphics::LINE);
      //            instancing_mesh0.draw();
      //            g.polygonMode(Graphics::FILL);
    }
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

  const std::string instancing_frag = R"(
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
};

}  // namespace al

#endif AL_VASPRENDER
