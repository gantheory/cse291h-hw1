#include "cuboid.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <random>

Cuboid::Cuboid(glm::vec3 cuboidMin, glm::vec3 cuboidMax) {
  model = glm::mat4(1.0f);

  // Set color.
  color = glm::vec3(1.f, 1.f, 1.f);

  // Set positions.
  glm::vec3 one(1.f);
  numOfParticles = glm::vec3(kNBlocks) * (cuboidMax - cuboidMin) + one;

  indexOfParticles.resize(
      numOfParticles.x,
      std::vector<std::vector<int>>(numOfParticles.y,
                                    std::vector<int>(numOfParticles.z)));

  positions.resize(numOfParticles.x * numOfParticles.y * numOfParticles.z);

  int index = 0;
  for (size_t i = 0; i < indexOfParticles.size(); ++i)
    for (size_t j = 0; j < indexOfParticles[i].size(); ++j)
      for (size_t k = 0; k < indexOfParticles[i][j].size(); ++k) {
        indexOfParticles[i][j][k] = index;
        positions[index] = cuboidMin + glm::vec3(i, j, k) *
                                           (cuboidMax - cuboidMin) /
                                           (numOfParticles - one);
        ++index;
      }

  // TODO: Set mass.
  mass.resize(positions.size());

  // TODO: Set velocities.
  velocities.resize(positions.size());

  numOfSquares = {// Front
                  ((int)numOfParticles.x - 1) * ((int)numOfParticles.y - 1),
                  // Back
                  ((int)numOfParticles.x - 1) * ((int)numOfParticles.y - 1),
                  // Top
                  ((int)numOfParticles.x - 1) * ((int)numOfParticles.z - 1),
                  // Bottom
                  ((int)numOfParticles.x - 1) * ((int)numOfParticles.z - 1),
                  // Left
                  ((int)numOfParticles.y - 1) * ((int)numOfParticles.z - 1),
                  // Right
                  ((int)numOfParticles.y - 1) * ((int)numOfParticles.z - 1)};
  totalNumOfSquares =
      std::accumulate(numOfSquares.begin(), numOfSquares.end(), 0);

  // Generate a vertex array (VAO) and two vertex buffer objects (VBO).
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO_positions);
  glGenBuffers(1, &VBO_normals);

  update();
}

Cuboid::~Cuboid() {}

void Cuboid::draw(const glm::mat4& viewProjMtx, GLuint shader) {
  // actiavte the shader program
  glUseProgram(shader);

  // get the locations and send the uniforms to the shader
  glUniformMatrix4fv(glGetUniformLocation(shader, "viewProj"), 1, false,
                     (float*)&viewProjMtx);
  glUniformMatrix4fv(glGetUniformLocation(shader, "model"), 1, GL_FALSE,
                     (float*)&model);
  glUniform3fv(glGetUniformLocation(shader, "DiffuseColor"), 1, &color[0]);

  // Bind the VAO
  glBindVertexArray(VAO);

  // draw the points using triangles, indexed with the EBO
  glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

  // Unbind the VAO and shader program
  glBindVertexArray(0);
  glUseProgram(0);
}

glm::vec3 Cuboid::OuterProduct(glm::vec3 a, glm::vec3 b) {
  glm::vec3 ans(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
  float norm = sqrt(pow(ans.x, 2) + pow(ans.y, 2) + pow(ans.z, 2));
  ans.x /= norm, ans.y /= norm, ans.z /= norm;
  return ans;
}

// SetSurface {{{
void Cuboid::SetSurface() {
  if (surface_positions.empty())
    surface_positions.resize(6 * totalNumOfSquares);
  if (surface_normals.empty()) surface_normals.resize(6 * totalNumOfSquares);
  if (indices.empty()) indices.resize(6 * totalNumOfSquares);
  iota(indices.begin(), indices.end(), 0);

  int index = 0;

  // Front
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int j = 0; j < numOfParticles.y - 1; ++j) {
      int z_axis = numOfParticles.z - 1;

      glm::vec3& p0 = surface_positions[index];
      glm::vec3& p1 = surface_positions[index + 1];
      glm::vec3& p2 = surface_positions[index + 2];
      glm::vec3& p3 = surface_positions[index + 3];
      glm::vec3& p4 = surface_positions[index + 4];
      glm::vec3& p5 = surface_positions[index + 5];

      p0 = positions[indexOfParticles[i][j][z_axis]];
      p1 = positions[indexOfParticles[i + 1][j][z_axis]];
      p2 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p3 = positions[indexOfParticles[i][j][z_axis]];
      p4 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p5 = positions[indexOfParticles[i][j + 1][z_axis]];

      glm::vec3& n0 = surface_normals[index];
      glm::vec3& n1 = surface_normals[index + 1];
      glm::vec3& n2 = surface_normals[index + 2];
      glm::vec3& n3 = surface_normals[index + 3];
      glm::vec3& n4 = surface_normals[index + 4];
      glm::vec3& n5 = surface_normals[index + 5];

      n0 = n1 = n2 = OuterProduct(p1 - p0, p2 - p0);
      n3 = n4 = n5 = OuterProduct(p4 - p3, p5 - p3);

      index += 6;
    }

  // Back
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int j = 0; j < numOfParticles.y - 1; ++j) {
      int z_axis = 0;

      glm::vec3& p0 = surface_positions[index];
      glm::vec3& p1 = surface_positions[index + 1];
      glm::vec3& p2 = surface_positions[index + 2];
      glm::vec3& p3 = surface_positions[index + 3];
      glm::vec3& p4 = surface_positions[index + 4];
      glm::vec3& p5 = surface_positions[index + 5];

      p0 = positions[indexOfParticles[i][j][z_axis]];
      p1 = positions[indexOfParticles[i + 1][j][z_axis]];
      p2 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p3 = positions[indexOfParticles[i][j][z_axis]];
      p4 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p5 = positions[indexOfParticles[i][j + 1][z_axis]];

      glm::vec3& n0 = surface_normals[index];
      glm::vec3& n1 = surface_normals[index + 1];
      glm::vec3& n2 = surface_normals[index + 2];
      glm::vec3& n3 = surface_normals[index + 3];
      glm::vec3& n4 = surface_normals[index + 4];
      glm::vec3& n5 = surface_normals[index + 5];

      n0 = n1 = n2 = -OuterProduct(p1 - p0, p2 - p0);
      n3 = n4 = n5 = -OuterProduct(p4 - p3, p5 - p3);

      index += 6;
    }

  // Top
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int y_axis = numOfParticles.y - 1;

      glm::vec3& p0 = surface_positions[index];
      glm::vec3& p1 = surface_positions[index + 1];
      glm::vec3& p2 = surface_positions[index + 2];
      glm::vec3& p3 = surface_positions[index + 3];
      glm::vec3& p4 = surface_positions[index + 4];
      glm::vec3& p5 = surface_positions[index + 5];

      p0 = positions[indexOfParticles[i][y_axis][k]];
      p1 = positions[indexOfParticles[i][y_axis][k + 1]];
      p2 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p3 = positions[indexOfParticles[i][y_axis][k]];
      p4 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p5 = positions[indexOfParticles[i + 1][y_axis][k]];

      glm::vec3& n0 = surface_normals[index];
      glm::vec3& n1 = surface_normals[index + 1];
      glm::vec3& n2 = surface_normals[index + 2];
      glm::vec3& n3 = surface_normals[index + 3];
      glm::vec3& n4 = surface_normals[index + 4];
      glm::vec3& n5 = surface_normals[index + 5];

      n0 = n1 = n2 = OuterProduct(p1 - p0, p2 - p0);
      n3 = n4 = n5 = OuterProduct(p4 - p3, p5 - p3);

      index += 6;
    }

  // Back
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int y_axis = 0;

      glm::vec3& p0 = surface_positions[index];
      glm::vec3& p1 = surface_positions[index + 1];
      glm::vec3& p2 = surface_positions[index + 2];
      glm::vec3& p3 = surface_positions[index + 3];
      glm::vec3& p4 = surface_positions[index + 4];
      glm::vec3& p5 = surface_positions[index + 5];

      p0 = positions[indexOfParticles[i][y_axis][k]];
      p1 = positions[indexOfParticles[i][y_axis][k + 1]];
      p2 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p3 = positions[indexOfParticles[i][y_axis][k]];
      p4 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p5 = positions[indexOfParticles[i + 1][y_axis][k]];

      glm::vec3& n0 = surface_normals[index];
      glm::vec3& n1 = surface_normals[index + 1];
      glm::vec3& n2 = surface_normals[index + 2];
      glm::vec3& n3 = surface_normals[index + 3];
      glm::vec3& n4 = surface_normals[index + 4];
      glm::vec3& n5 = surface_normals[index + 5];

      n0 = n1 = n2 = -OuterProduct(p1 - p0, p2 - p0);
      n3 = n4 = n5 = -OuterProduct(p4 - p3, p5 - p3);

      index += 6;
    }

  // Left
  for (int j = 0; j < numOfParticles.y - 1; ++j)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int x_axis = 0;

      glm::vec3& p0 = surface_positions[index];
      glm::vec3& p1 = surface_positions[index + 1];
      glm::vec3& p2 = surface_positions[index + 2];
      glm::vec3& p3 = surface_positions[index + 3];
      glm::vec3& p4 = surface_positions[index + 4];
      glm::vec3& p5 = surface_positions[index + 5];

      p0 = positions[indexOfParticles[x_axis][j][k]];
      p1 = positions[indexOfParticles[x_axis][j + 1][k]];
      p2 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p3 = positions[indexOfParticles[x_axis][j][k]];
      p4 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p5 = positions[indexOfParticles[x_axis][j][k + 1]];

      glm::vec3& n0 = surface_normals[index];
      glm::vec3& n1 = surface_normals[index + 1];
      glm::vec3& n2 = surface_normals[index + 2];
      glm::vec3& n3 = surface_normals[index + 3];
      glm::vec3& n4 = surface_normals[index + 4];
      glm::vec3& n5 = surface_normals[index + 5];

      n0 = n1 = n2 = -OuterProduct(p1 - p0, p2 - p0);
      n3 = n4 = n5 = -OuterProduct(p4 - p3, p5 - p3);

      index += 6;
    }

  // Right
  for (int j = 0; j < numOfParticles.y - 1; ++j)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int x_axis = numOfParticles.x - 1;

      glm::vec3& p0 = surface_positions[index];
      glm::vec3& p1 = surface_positions[index + 1];
      glm::vec3& p2 = surface_positions[index + 2];
      glm::vec3& p3 = surface_positions[index + 3];
      glm::vec3& p4 = surface_positions[index + 4];
      glm::vec3& p5 = surface_positions[index + 5];

      p0 = positions[indexOfParticles[x_axis][j][k]];
      p1 = positions[indexOfParticles[x_axis][j + 1][k]];
      p2 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p3 = positions[indexOfParticles[x_axis][j][k]];
      p4 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p5 = positions[indexOfParticles[x_axis][j][k + 1]];

      glm::vec3& n0 = surface_normals[index];
      glm::vec3& n1 = surface_normals[index + 1];
      glm::vec3& n2 = surface_normals[index + 2];
      glm::vec3& n3 = surface_normals[index + 3];
      glm::vec3& n4 = surface_normals[index + 4];
      glm::vec3& n5 = surface_normals[index + 5];

      n0 = n1 = n2 = OuterProduct(p1 - p0, p2 - p0);
      n3 = n4 = n5 = OuterProduct(p4 - p3, p5 - p3);

      index += 6;
    }
}
// }}}

void Cuboid::update() {
  SetSurface();

  // Bind to the VAO.
  glBindVertexArray(VAO);

  // Bind to the first VBO - We will use it to store the vertices
  glBindBuffer(GL_ARRAY_BUFFER, VBO_positions);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * surface_positions.size(),
               surface_positions.data(), GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  // Bind to the second VBO - We will use it to store the normals
  glBindBuffer(GL_ARRAY_BUFFER, VBO_normals);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * surface_normals.size(),
               surface_normals.data(), GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  // Generate EBO, bind the EBO to the bound VAO and send the data
  glGenBuffers(1, &EBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(),
               indices.data(), GL_DYNAMIC_DRAW);

  // Unbind the VBOs.
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}
