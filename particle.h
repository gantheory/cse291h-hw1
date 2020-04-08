#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "core.h"

class Particle {
 private:
  float mass;

  glm::vec3 position, velocity, color;

 public:
  Particle();

  Particle(glm::vec3 _position);

  ~Particle();

  void SetPosition(glm::vec3 _position);
};

#endif
