#include "cuboid.h"
#ifdef __APPLE__
#define GLFW_INCLUDE_GLCOREARB
#include <OpenGL/gl3.h>
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <random>

Cuboid::Cuboid(glm::vec3 cuboidMin, glm::vec3 cuboidMax) {
  model = glm::mat4(1.0f);

  // Set color.
  color = glm::vec3(1.f, 1.f, 1.f);

  volume = (cuboidMax.x - cuboidMin.x) * (cuboidMax.y - cuboidMin.y) *
           (cuboidMax.z - cuboidMin.z);
  totalMass = volume * kDensity;

  lambda = kE * kNu / (1 + kNu) / (1 - 2 * kNu);
  mu = kE / 2 / (1 + kNu);
  std::cout << "lambda: " << lambda << ", mu: " << mu << std::endl;

  // Set positions.
  glm::vec3 one(1.f);
  numOfParticles = glm::vec3(kNBlocks) * (cuboidMax - cuboidMin) + one;

  indexOfParticles.resize(
      numOfParticles.x,
      std::vector<std::vector<int>>(numOfParticles.y,
                                    std::vector<int>(numOfParticles.z)));

  positions.resize(numOfParticles.x * numOfParticles.y * numOfParticles.z);

  glm::mat3 kRotation = {{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, {0.f, 0.f, 1.f}};
  glm::vec3 kTranslation = {0.f, 0.f, 0.f};

  int index = 0;
  for (size_t i = 0; i < indexOfParticles.size(); ++i)
    for (size_t j = 0; j < indexOfParticles[i].size(); ++j)
      for (size_t k = 0; k < indexOfParticles[i][j].size(); ++k) {
        indexOfParticles[i][j][k] = index;
        positions[index] = cuboidMin + glm::vec3(i, j, k) *
                                           (cuboidMax - cuboidMin) /
                                           (numOfParticles - one);
        positions[index] = kRotation * positions[index];
        positions[index] += kTranslation;
        ++index;
      }

  tetrahedrons.resize(5 * (numOfParticles.x - 1) * (numOfParticles.y - 1) *
                          (numOfParticles.z - 1),
                      std::vector<int>(4));
  index = 0;
  for (size_t i = 0; i < indexOfParticles.size() - 1; ++i)
    for (size_t j = 0; j < indexOfParticles[i].size() - 1; ++j)
      for (size_t k = 0; k < indexOfParticles[i][j].size() - 1; ++k)
        for (size_t l = 0; l < kCubeToTetrahedron.size(); ++l) {
          for (size_t m = 0; m < 4; ++m) {
            int di = kCubeToTetrahedron[l][m].x;
            int dj = kCubeToTetrahedron[l][m].y;
            int dk = kCubeToTetrahedron[l][m].z;
            tetrahedrons[index][m] = indexOfParticles[i + di][j + dj][k + dk];
          }
          ++index;
        }

  mass.resize(positions.size());
  for (std::vector<int>& tetrahedron : tetrahedrons) {
    int& i0 = tetrahedron[0];
    int& i1 = tetrahedron[1];
    int& i2 = tetrahedron[2];
    int& i3 = tetrahedron[3];
    glm::vec3& p0 = positions[i0];
    glm::vec3& p1 = positions[i1];
    glm::vec3& p2 = positions[i2];
    glm::vec3& p3 = positions[i3];
    float vol = glm::dot(glm::cross(p1 - p0, p2 - p0), p3 - p0) / 6.f;
    float m = abs(vol) * kDensity * 0.25f;
    mass[i0] += m;
    mass[i1] += m;
    mass[i2] += m;
    mass[i3] += m;
  }

  force.resize(positions.size());

  velocities.resize(positions.size());

  // Calculate rest tetrahedral frame (R^-1).
  inverseRs.resize(tetrahedrons.size(),
                   std::vector<glm::mat3>(kSurfaceOfTetrahedron.size()));
  n_stars.resize(tetrahedrons.size(),
                 std::vector<glm::vec3>(kSurfaceOfTetrahedron.size()));
  glm::vec3 half(0.5f);
  for (size_t i = 0; i < tetrahedrons.size(); ++i) {
    std::vector<int>& tetrahedron = tetrahedrons[i];
    for (size_t j = 0; j < kSurfaceOfTetrahedron.size(); ++j) {
      const std::vector<int>& order = kSurfaceOfTetrahedron[j];
      int& i0 = tetrahedron[order[0]];
      int& i1 = tetrahedron[order[1]];
      int& i2 = tetrahedron[order[2]];
      int& i3 = tetrahedron[order[3]];
      glm::vec3& r0 = positions[i0];
      glm::vec3& r1 = positions[i1];
      glm::vec3& r2 = positions[i2];
      glm::vec3& r3 = positions[i3];
      inverseRs[i][j] = inverse(glm::mat3(r1 - r0, r2 - r0, r3 - r0));
      n_stars[i][j] = cross(r2 - r1, r3 - r1) * half;
    }
  }

  totalNumOfSquares =
      2 * ((((int)numOfParticles.x - 1) * ((int)numOfParticles.y - 1)) +
           (((int)numOfParticles.x - 1) * ((int)numOfParticles.z - 1)) +
           (((int)numOfParticles.y - 1) * ((int)numOfParticles.z - 1)));

  // Generate a vertex array (VAO) and two vertex buffer objects (VBO).
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO_positions);
  glGenBuffers(1, &VBO_normals);

  Update();
}

Cuboid::~Cuboid() {
  // Delete the VBOs and the VAO.
  glDeleteBuffers(1, &VBO_positions);
  glDeleteBuffers(1, &VBO_normals);
  glDeleteBuffers(1, &EBO);
  glDeleteVertexArrays(1, &VAO);
}

void Cuboid::Draw(const glm::mat4& viewProjMtx, GLuint shader) {
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

void Cuboid::Update() {
  CalculateForce();

  ApplyForce();

  SetSurface();

  // Bind to the VAO.
  glBindVertexArray(VAO);

  // Bind to the first VBO - We will use it to store the vertices
  glBindBuffer(GL_ARRAY_BUFFER, VBO_positions);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * surfacePositions.size(),
               surfacePositions.data(), GL_STATIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  // Bind to the second VBO - We will use it to store the normals
  glBindBuffer(GL_ARRAY_BUFFER, VBO_normals);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * surfaceNormals.size(),
               surfaceNormals.data(), GL_STATIC_DRAW);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

  // Generate EBO, bind the EBO to the bound VAO and send the data
  glGenBuffers(1, &EBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(),
               indices.data(), GL_STATIC_DRAW);

  // Unbind the VBOs.
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}
std::pair<glm::mat3, glm::mat3> Cuboid::ComputeStrain(glm::vec3& r0,
                                                      glm::vec3& r1,
                                                      glm::vec3& r2,
                                                      glm::vec3& r3,
                                                      glm::mat3& inverseR) {
  glm::vec3 e1 = r1 - r0;
  glm::vec3 e2 = r2 - r0;
  glm::vec3 e3 = r3 - r0;
  glm::mat3 T(e1, e2, e3);
  glm::mat3 F = T * inverseR;
  glm::mat3 epsilon = 0.5f * (transpose(F) * F - glm::mat3(1.f));
  return {F, epsilon};
}

glm::mat3 Cuboid::ComputeStress(glm::mat3& epsilon) {
  glm::mat3 mu_mat(glm::vec3(2.f * mu), glm::vec3(2.f * mu),
                   glm::vec3(2.f * mu));
  float trace = epsilon[0][0] + epsilon[1][1] + epsilon[2][2];
  return matrixCompMult(mu_mat, epsilon) + glm::mat3(trace * lambda);
}

void Cuboid::CalculateForce() {
  fill(force.begin(), force.end(), glm::vec3(0.f));

  // Gravity.
  for (size_t i = 0; i < positions.size(); ++i) {
    force[i].y -= mass[i] * kG;
  }

  // Linear elasicity.
  for (size_t i = 0; i < tetrahedrons.size(); ++i) {
    std::vector<int>& tetrahedron = tetrahedrons[i];
    for (size_t j = 0; j < kSurfaceOfTetrahedron.size(); ++j) {
      const std::vector<int>& order = kSurfaceOfTetrahedron[j];
      int& i0 = tetrahedron[order[0]];
      int& i1 = tetrahedron[order[1]];
      int& i2 = tetrahedron[order[2]];
      int& i3 = tetrahedron[order[3]];
      glm::vec3& r0 = positions[i0];
      glm::vec3& r1 = positions[i1];
      glm::vec3& r2 = positions[i2];
      glm::vec3& r3 = positions[i3];

      if (j == 0) {
        float vol = glm::dot(glm::cross(r1 - r0, r2 - r0), r3 - r0) / 6.f;
        if (vol < kEpsilon) std::cerr << "vol: " << vol << std::endl;
        assert(vol > kEpsilon);
      }

      auto [F, epsilon] = ComputeStrain(r0, r1, r2, r3, inverseRs[i][j]);

      glm::mat3 sigma = ComputeStress(epsilon);

      force[i0] += F * sigma * n_stars[i][j];
      if (isnan(force[i0].x) || isnan(force[i0].y) || isnan(force[i0].z)) {
        std::cerr << to_string(force[i0]) << " should not be nan!" << std::endl;
        exit(0);
      }
    }
  }
}

void Cuboid::ApplyForce() {
  for (size_t i = 0; i < positions.size(); ++i) {
    glm::vec3 acceleration = force[i] / mass[i];
    velocities[i] += acceleration * kdt;
    positions[i] += velocities[i] * kdt;
  }

  // Ground detection.
  for (size_t i = 0; i < positions.size(); ++i)
    if (positions[i].y < kEpsilon) {
      positions[i].y = 0, velocities[i].y *= -1;
    }

  // Velocity decay.
  for (size_t i = 0; i < positions.size(); ++i) velocities[i] *= kVelocityDecay;
}

void Cuboid::SetSurface() {
  if (surfacePositions.empty()) surfacePositions.resize(6 * totalNumOfSquares);
  if (surfaceNormals.empty()) surfaceNormals.resize(6 * totalNumOfSquares);
  if (indices.empty()) indices.resize(6 * totalNumOfSquares);
  iota(indices.begin(), indices.end(), 0);

  int index = 0;

  // Front
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int j = 0; j < numOfParticles.y - 1; ++j) {
      int z_axis = numOfParticles.z - 1;

      glm::vec3& p0 = surfacePositions[index];
      glm::vec3& p1 = surfacePositions[index + 1];
      glm::vec3& p2 = surfacePositions[index + 2];
      glm::vec3& p3 = surfacePositions[index + 3];
      glm::vec3& p4 = surfacePositions[index + 4];
      glm::vec3& p5 = surfacePositions[index + 5];

      p0 = positions[indexOfParticles[i][j][z_axis]];
      p1 = positions[indexOfParticles[i + 1][j][z_axis]];
      p2 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p3 = positions[indexOfParticles[i][j][z_axis]];
      p4 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p5 = positions[indexOfParticles[i][j + 1][z_axis]];

      glm::vec3& n0 = surfaceNormals[index];
      glm::vec3& n1 = surfaceNormals[index + 1];
      glm::vec3& n2 = surfaceNormals[index + 2];
      glm::vec3& n3 = surfaceNormals[index + 3];
      glm::vec3& n4 = surfaceNormals[index + 4];
      glm::vec3& n5 = surfaceNormals[index + 5];

      n0 = n1 = n2 = glm::normalize(glm::cross(p1 - p0, p2 - p0));
      n3 = n4 = n5 = glm::normalize(glm::cross(p4 - p3, p5 - p3));

      index += 6;
    }

  // Back
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int j = 0; j < numOfParticles.y - 1; ++j) {
      int z_axis = 0;

      glm::vec3& p0 = surfacePositions[index];
      glm::vec3& p1 = surfacePositions[index + 1];
      glm::vec3& p2 = surfacePositions[index + 2];
      glm::vec3& p3 = surfacePositions[index + 3];
      glm::vec3& p4 = surfacePositions[index + 4];
      glm::vec3& p5 = surfacePositions[index + 5];

      p0 = positions[indexOfParticles[i][j][z_axis]];
      p1 = positions[indexOfParticles[i + 1][j][z_axis]];
      p2 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p3 = positions[indexOfParticles[i][j][z_axis]];
      p4 = positions[indexOfParticles[i + 1][j + 1][z_axis]];
      p5 = positions[indexOfParticles[i][j + 1][z_axis]];

      glm::vec3& n0 = surfaceNormals[index];
      glm::vec3& n1 = surfaceNormals[index + 1];
      glm::vec3& n2 = surfaceNormals[index + 2];
      glm::vec3& n3 = surfaceNormals[index + 3];
      glm::vec3& n4 = surfaceNormals[index + 4];
      glm::vec3& n5 = surfaceNormals[index + 5];

      n0 = n1 = n2 = -glm::normalize(glm::cross(p1 - p0, p2 - p0));
      n3 = n4 = n5 = -glm::normalize(glm::cross(p4 - p3, p5 - p3));

      index += 6;
    }

  // Top
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int y_axis = numOfParticles.y - 1;

      glm::vec3& p0 = surfacePositions[index];
      glm::vec3& p1 = surfacePositions[index + 1];
      glm::vec3& p2 = surfacePositions[index + 2];
      glm::vec3& p3 = surfacePositions[index + 3];
      glm::vec3& p4 = surfacePositions[index + 4];
      glm::vec3& p5 = surfacePositions[index + 5];

      p0 = positions[indexOfParticles[i][y_axis][k]];
      p1 = positions[indexOfParticles[i][y_axis][k + 1]];
      p2 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p3 = positions[indexOfParticles[i][y_axis][k]];
      p4 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p5 = positions[indexOfParticles[i + 1][y_axis][k]];

      glm::vec3& n0 = surfaceNormals[index];
      glm::vec3& n1 = surfaceNormals[index + 1];
      glm::vec3& n2 = surfaceNormals[index + 2];
      glm::vec3& n3 = surfaceNormals[index + 3];
      glm::vec3& n4 = surfaceNormals[index + 4];
      glm::vec3& n5 = surfaceNormals[index + 5];

      n0 = n1 = n2 = glm::normalize(glm::cross(p1 - p0, p2 - p0));
      n3 = n4 = n5 = glm::normalize(glm::cross(p4 - p3, p5 - p3));

      index += 6;
    }

  // Back
  for (int i = 0; i < numOfParticles.x - 1; ++i)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int y_axis = 0;

      glm::vec3& p0 = surfacePositions[index];
      glm::vec3& p1 = surfacePositions[index + 1];
      glm::vec3& p2 = surfacePositions[index + 2];
      glm::vec3& p3 = surfacePositions[index + 3];
      glm::vec3& p4 = surfacePositions[index + 4];
      glm::vec3& p5 = surfacePositions[index + 5];

      p0 = positions[indexOfParticles[i][y_axis][k]];
      p1 = positions[indexOfParticles[i][y_axis][k + 1]];
      p2 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p3 = positions[indexOfParticles[i][y_axis][k]];
      p4 = positions[indexOfParticles[i + 1][y_axis][k + 1]];
      p5 = positions[indexOfParticles[i + 1][y_axis][k]];

      glm::vec3& n0 = surfaceNormals[index];
      glm::vec3& n1 = surfaceNormals[index + 1];
      glm::vec3& n2 = surfaceNormals[index + 2];
      glm::vec3& n3 = surfaceNormals[index + 3];
      glm::vec3& n4 = surfaceNormals[index + 4];
      glm::vec3& n5 = surfaceNormals[index + 5];

      n0 = n1 = n2 = -glm::normalize(glm::cross(p1 - p0, p2 - p0));
      n3 = n4 = n5 = -glm::normalize(glm::cross(p4 - p3, p5 - p3));

      index += 6;
    }

  // Left
  for (int j = 0; j < numOfParticles.y - 1; ++j)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int x_axis = 0;

      glm::vec3& p0 = surfacePositions[index];
      glm::vec3& p1 = surfacePositions[index + 1];
      glm::vec3& p2 = surfacePositions[index + 2];
      glm::vec3& p3 = surfacePositions[index + 3];
      glm::vec3& p4 = surfacePositions[index + 4];
      glm::vec3& p5 = surfacePositions[index + 5];

      p0 = positions[indexOfParticles[x_axis][j][k]];
      p1 = positions[indexOfParticles[x_axis][j + 1][k]];
      p2 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p3 = positions[indexOfParticles[x_axis][j][k]];
      p4 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p5 = positions[indexOfParticles[x_axis][j][k + 1]];

      glm::vec3& n0 = surfaceNormals[index];
      glm::vec3& n1 = surfaceNormals[index + 1];
      glm::vec3& n2 = surfaceNormals[index + 2];
      glm::vec3& n3 = surfaceNormals[index + 3];
      glm::vec3& n4 = surfaceNormals[index + 4];
      glm::vec3& n5 = surfaceNormals[index + 5];

      n0 = n1 = n2 = -glm::normalize(glm::cross(p1 - p0, p2 - p0));
      n3 = n4 = n5 = -glm::normalize(glm::cross(p4 - p3, p5 - p3));

      index += 6;
    }

  // Right
  for (int j = 0; j < numOfParticles.y - 1; ++j)
    for (int k = 0; k < numOfParticles.z - 1; ++k) {
      int x_axis = numOfParticles.x - 1;

      glm::vec3& p0 = surfacePositions[index];
      glm::vec3& p1 = surfacePositions[index + 1];
      glm::vec3& p2 = surfacePositions[index + 2];
      glm::vec3& p3 = surfacePositions[index + 3];
      glm::vec3& p4 = surfacePositions[index + 4];
      glm::vec3& p5 = surfacePositions[index + 5];

      p0 = positions[indexOfParticles[x_axis][j][k]];
      p1 = positions[indexOfParticles[x_axis][j + 1][k]];
      p2 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p3 = positions[indexOfParticles[x_axis][j][k]];
      p4 = positions[indexOfParticles[x_axis][j + 1][k + 1]];
      p5 = positions[indexOfParticles[x_axis][j][k + 1]];

      glm::vec3& n0 = surfaceNormals[index];
      glm::vec3& n1 = surfaceNormals[index + 1];
      glm::vec3& n2 = surfaceNormals[index + 2];
      glm::vec3& n3 = surfaceNormals[index + 3];
      glm::vec3& n4 = surfaceNormals[index + 4];
      glm::vec3& n5 = surfaceNormals[index + 5];

      n0 = n1 = n2 = glm::normalize(glm::cross(p1 - p0, p2 - p0));
      n3 = n4 = n5 = glm::normalize(glm::cross(p4 - p3, p5 - p3));

      index += 6;
    }
}
