#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
    
    for (int particleIdx = 0; particleIdx < jelly.getParticleNum(); particleIdx++) {
        Particle& P = jelly.getParticle(particleIdx);

        Eigen::Vector3f ParticlePos = P.getPosition();
        Eigen::Vector3f ParticleVel = P.getVelocity();
        Eigen::Vector3f ParticleFor = P.getForce();

        if (this->normal.dot(ParticlePos - this->position) < eEPSILON && this->normal.dot(ParticleVel) < 0 &&
            (ParticlePos - this->hole_position).norm() > this->hole_radius) {
            Eigen::Vector3f VelN = ParticleVel.dot(this->normal) * this->normal;
            Eigen::Vector3f VelT = ParticleVel - VelN;
            P.setVelocity(-coefResist * VelN + VelT);
            Eigen::Vector3f ContactForce = Eigen::Vector3f::Zero();
            Eigen::Vector3f FrictionForce = Eigen::Vector3f::Zero();
            if (this->normal.dot(ParticleFor) < 0) {
                ContactForce = (-1) * (this->normal.dot(ParticleFor)) * this->normal;            
                FrictionForce = (-1) * coefFriction * ((-1)*this->normal.dot(ParticleFor)) * VelT;
            }
            P.addForce(FrictionForce);
            P.addForce(ContactForce);     
        }
    }
}
        // BowlTerrain //

BowlTerrain::BowlTerrain() {
    modelMatrix = util::translate(position) * util::rotateDegree(-90, 0, 0) * util::scale(radius, radius, radius);
}

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-2: Handle collision when a particle collide with the sphere terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data. 
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
    for (int particleIdx = 0; particleIdx < jelly.getParticleNum(); particleIdx++) {
        Particle& P = jelly.getParticle(particleIdx);
        Eigen::Vector3f ParticlePos = P.getPosition();
        Eigen::Vector3f ParticleVel = P.getVelocity();
        Eigen::Vector3f ParticleFor = P.getForce();

        Eigen::Vector3f normal = (this->position - ParticlePos).normalized();
        Eigen::Vector3f pos = ParticlePos;  
        if (normal.dot(ParticlePos - this->position) < eEPSILON && normal.dot(ParticleVel) < 0 &&
            (ParticlePos - this->position).norm() > this->radius &&ParticlePos[1] <= 0) {
            Eigen::Vector3f VelN = ParticleVel.dot(normal) * normal;
            Eigen::Vector3f VelT = ParticleVel - VelN;
            P.setVelocity(VelN * (P.getMass() - this->mass) / (P.getMass() + this->mass) + VelT);
            Eigen::Vector3f ContactForce = Eigen::Vector3f::Zero();
            Eigen::Vector3f FrictionForce = Eigen::Vector3f::Zero();
            if (normal.dot(ParticleFor) < 0) {
                ContactForce = (-1) * (normal.dot(ParticleFor)) * normal;
                FrictionForce = (-1) * coefFriction * ((-1) * normal.dot(ParticleFor)) * VelT;
            }
            P.addForce(FrictionForce);
            P.addForce(ContactForce);  
        }

    }
}
}  // namespace simulation
