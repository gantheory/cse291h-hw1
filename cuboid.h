#ifndef _CUBOID_H_
#define _CUBOID_H_

#include "core.h"

class Cuboid {
 private:
  const int kNBlocks = 2;
  const float kDensity = 5.000;
  const float kG = 10.0;
  const float kEpsilon = 1e-6;
  const float kVelocityDecay = 0.999;
  const float kE = 680.0;
  const float kNu = 0.487;
  const float kdt = 0.001;

  const std::vector<std::vector<int>> kSurfaceOfTetrahedron = {
      {2, 0, 1, 3}, {1, 0, 3, 2}, {3, 0, 2, 1}, {0, 1, 2, 3}};
  const std::vector<std::vector<glm::vec3>> kCubeToTetrahedron = {
      {
          glm::vec3(1, 1, 1),
          glm::vec3(1, 0, 0),
          glm::vec3(1, 0, 1),
          glm::vec3(0, 0, 1),
      },
      {
          glm::vec3(0, 0, 1),
          glm::vec3(1, 1, 1),
          glm::vec3(1, 0, 0),
          glm::vec3(0, 1, 0),
      },
      {
          glm::vec3(0, 0, 1),
          glm::vec3(0, 1, 1),
          glm::vec3(1, 1, 1),
          glm::vec3(0, 1, 0),
      },
      {
          glm::vec3(1, 0, 0),
          glm::vec3(0, 1, 0),
          glm::vec3(1, 1, 1),
          glm::vec3(1, 1, 0),
      },
      {
          glm::vec3(0, 0, 0),
          glm::vec3(0, 1, 0),
          glm::vec3(0, 0, 1),
          glm::vec3(1, 0, 0),
      },
  };

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

  std::pair<glm::mat3, glm::mat3> ComputeStrain(glm::vec3& r0, glm::vec3& r1,
                                                glm::vec3& r2, glm::vec3& r3,
                                                glm::mat3& inverseR);

  glm::mat3 ComputeStress(glm::mat3& epsilon);

  void CalculateForce();
  void ApplyForce();

  void SetSurface();

 public:
  Cuboid(glm::vec3 cuboidSide);
  ~Cuboid();

  void Draw(const glm::mat4& viewProjMtx, GLuint shader);
  void Update();
};

#endif
