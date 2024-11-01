//
// Created by jkr on 10/18/24.
//
#include "ParticleContainer.h"

ParticleContainer::ParticleContainer() = default;
ParticleContainer::~ParticleContainer() = default;

ParticleContainer::ParticleContainer(const std::list<Particle>& particles) {
  for (const auto& particle : particles) {
    this->particles.push_back(particle);
  }
}

ParticleContainer::ParticleContainer(const std::vector<Particle>& particles) {
  for (const auto& particle : particles) {
    this->particles.push_back(particle);
  }
}

void ParticleContainer::addParticle(const Particle& p) {
  particles.push_back(p);
}

auto ParticleContainer::begin() { return particles.begin(); }

auto ParticleContainer::end() { return particles.end(); }

Particle ParticleContainer::operator[](const int n) { return particles[n]; }

std::vector<Particle> ParticleContainer::getParticles() { return particles; }

std::vector<Particle>& ParticleContainer::getParticlesReference() { return particles; }

std::size_t ParticleContainer::size() const { return particles.size(); }

void ParticleContainer::setParticles(const std::vector<Particle>& Particles) { particles = Particles; }