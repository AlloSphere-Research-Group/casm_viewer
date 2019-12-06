#ifndef INCLUDE_INSTANCED_MESH_HPP
#define INCLUDE_INSTANCED_MESH_HPP

#undef CIEXYZ
#include "al/graphics/al_Shader.hpp"
#include "al/graphics/al_VAOMesh.hpp"

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
  al::VAOMesh mesh;
  al::BufferObject buffer;
  al::ShaderProgram shader;
  size_t num_instances = 0;

  // attrib_loc: per instance attribute location, should be set
  //             with "layout (location=#)" in vert shader
  // attrib_num_elems: vec4? vec3?
  // attrib_type: GL_FLOAT? GL_UNSIGNED_BYTE?
  void init(const std::string &vert_str, const std::string &frag_str,
            GLuint attrib_loc, GLint attrib_num_elems, GLenum attrib_type) {
    shader.compile(vert_str, frag_str);
    buffer.bufferType(GL_ARRAY_BUFFER);
    buffer.usage(GL_DYNAMIC_DRAW); // assumes buffer will change every frame
                                   // and will be used for drawing
    buffer.create();

    auto &v = mesh.vao();
    v.bind();
    v.enableAttrib(attrib_loc);
    // for normalizing, this code only considers GL_FLOAT AND GL_UNSIGNED_BYTE,
    // (does not normalize floats and normalizes unsigned bytes)
    v.attribPointer(attrib_loc, buffer, attrib_num_elems, attrib_type,
                    (attrib_type == GL_FLOAT) ? GL_FALSE
                                              : GL_TRUE, // normalize?
                    0,                                   // stride
                    0);                                  // offset
    glVertexAttribDivisor(attrib_loc, 1); // step attribute once per instance
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

#endif
