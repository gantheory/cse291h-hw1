#ifndef _CUBOID_H_
#define _CUBOID_H_

#include "core.h"

class Cuboid {
 private:
  // const int kNBlocks = 3;         // Should be at least 1.
  // const float kDensity = 1.0;     // kg / m^3
  // const float kG = 0.0003 * 9.8;  // N / kg
  // const float kEpsilon = 1e-6;
  // const float kVelocityDecay = 1.0;  // 0.9995;
  // const float E = 10.f;
  // const float nu = 0.25f;  // [0.0, 0.5]

  GLuint VAO;
  GLuint VBO_positions, VBO_normals, EBO;

  glm::mat4 model;
  glm::vec3 color;
  float totalMass, volume;
  // GLfloat lastFrame;
  float lambda, mu;

  // Properties for linear elasticity.
  std::vector<float> mass;
  std::vector<glm::vec3> positions, velocities, force;
  std::vector<std::vector<std::vector<int>>> indexOfParticles;
  std::vector<std::vector<int>> tetrahedrons;
  std::vector<std::vector<glm::mat3>> inverseRs;
  std::vector<std::vector<glm::vec3>> n_stars;

  // For rendering.
  std::vector<glm::vec3> surfacePositions, surfaceNormals;
  std::vector<unsigned int> indices;

  // Helpful member variables.
  glm::vec3 numOfParticles;
  int totalNumOfSquares;
  std::vector<int> numOfSquares;
  const std::vector<std::vector<int>> kSurfaceOfTetrahedron = {
      {2, 0, 1, 3}, {1, 0, 3, 2}, {3, 0, 2, 1}, {0, 1, 2, 3}};

  void SetSurface();

  void CalculateForce();
  void ApplyForce();

  std::pair<glm::mat3, glm::mat3> ComputeStrain(glm::vec3& r0, glm::vec3& r1,
                                                glm::vec3& r2, glm::vec3& r3,
                                                glm::mat3& inverseR);

  glm::mat3 ComputeStress(glm::mat3& epsilon);

 public:
  Cuboid(glm::vec3 cuboidMin = glm::vec3(-1, -1, -1),
         glm::vec3 cuboidMax = glm::vec3(1, 1, 1));
  ~Cuboid();

  void draw(const glm::mat4& viewProjMtx, GLuint shader);
  void update();
};

#endif
