#ifndef _CUBOID_H_
#define _CUBOID_H_

#include "core.h"

class Cuboid {
 private:
  const int kNBlocks = 2;  // Should be at least 1.

  GLuint VAO;
  GLuint VBO_positions, VBO_normals, EBO;

  glm::mat4 model;
  glm::vec3 color;

  std::vector<float> mass;
  std::vector<glm::vec3> positions, velocities;
  std::vector<std::vector<std::vector<int>>> indexOfParticles;

  std::vector<glm::vec3> surface_positions, surface_normals;
  std::vector<unsigned int> indices;

  glm::vec3 numOfParticles;
  int totalNumOfSquares;
  std::vector<int> numOfSquares;

  void SetSurface();

  glm::vec3 OuterProduct(glm::vec3 a, glm::vec3 b);

 public:
  Cuboid(glm::vec3 cuboidMin = glm::vec3(-1, -1, -1),
         glm::vec3 cuboidMax = glm::vec3(1, 1, 1));
  ~Cuboid();

  void draw(const glm::mat4& viewProjMtx, GLuint shader);
  void update();
};

#endif
