#ifndef _CUBOID_H_
#define _CUBOID_H_

#include "core.h"

class Cuboid {
 private:
  const int kNBlocks = 1;         // Should be at least 1.
  const float kDensity = 1.0;     // kg / m^3
  const float kG = 0.0003 * 9.8;  // N / kg
  const float kEpsilon = 1e-6;
  const float kVelocityDecay = 0.9995;

  GLuint VAO;
  GLuint VBO_positions, VBO_normals, EBO;

  glm::mat4 model;
  glm::vec3 color;
  float totalMass, volume;
  GLfloat lastFrame;

  // Properties for linear elasticity.
  std::vector<float> mass;
  std::vector<glm::vec3> positions, velocities, force;
  std::vector<std::vector<std::vector<int>>> indexOfParticles;
  std::vector<std::vector<int>> tetrahedrons;

  // For rendering.
  std::vector<glm::vec3> surfacePositions, surfaceNormals;
  std::vector<unsigned int> indices;

  // Helpful member variables.
  glm::vec3 numOfParticles;
  int totalNumOfSquares;
  std::vector<int> numOfSquares;

  void SetSurface();

  void CalculateForce();
  void ApplyForce();

 public:
  Cuboid(glm::vec3 cuboidMin = glm::vec3(-1, -1, -1),
         glm::vec3 cuboidMax = glm::vec3(1, 1, 1));
  ~Cuboid();

  void draw(const glm::mat4& viewProjMtx, GLuint shader);
  void update();
};

#endif
