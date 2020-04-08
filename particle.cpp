#include "particle.h"

Particle::Particle()
    : mass(0.0),
      position(0.f, 0.f, 0.f),
      velocity(0.f, 0.f, 0.f),
      color(1.f, 1.f, 1.f) {}

Particle::Particle(glm::vec3 _position)
    : mass(0.0),
      position(_position),
      velocity(0.f, 0.f, 0.f),
      color(1.f, 1.f, 1.f) {}

Particle::~Particle() {}

void Particle::SetPosition(glm::vec3 _position) { position = _position; }
