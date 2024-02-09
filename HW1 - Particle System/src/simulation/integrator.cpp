#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {

    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first. Then you can check whether your collision is correct or not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.15 - p.16
    for (int jellyIdx = 0; jellyIdx < particleSystem.jellyCount; jellyIdx++) {
        Jelly* J = particleSystem.getJellyPointer(jellyIdx);
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            P.addPosition(Vel * particleSystem.deltaTime);
            P.addVelocity(Acc * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
    }

}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1
    //   2. Review ¡§ODE_basics.pptx¡¨ from p.18 - p.19
    Jelly tmp;
   
    for (int jellyIdx = 0; jellyIdx < particleSystem.jellyCount; jellyIdx++) {
        Jelly* J = particleSystem.getJellyPointer(jellyIdx);
        tmp = *J;
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            //Eigen::Vector3f Vel = P.getVelocity();
            //Eigen::Vector3f Acc = P.getAcceleration();
            //P.addPosition(Vel * particleSystem.deltaTime);
            //P.addVelocity(Acc * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*J);
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            P.setPosition(tmp.getParticle(particleIdx).getPosition() + Vel * particleSystem.deltaTime);
            P.setVelocity(tmp.getParticle(particleIdx).getVelocity() + Acc * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
    }   
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. Review ¡§ODE_basics.pptx¡¨ from p .18 - p .20
    Jelly tmp;

    for (int jellyIdx = 0; jellyIdx < particleSystem.jellyCount; jellyIdx++) {
        Jelly* J = particleSystem.getJellyPointer(jellyIdx);
        tmp = *J;
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            P.addPosition(Vel * 0.5*particleSystem.deltaTime);
            P.addVelocity(Acc * 0.5*particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*J);
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            P.setPosition(tmp.getParticle(particleIdx).getPosition() + Vel * particleSystem.deltaTime);
            P.setVelocity(tmp.getParticle(particleIdx).getVelocity() + Acc * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
    }
    
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. StateStep struct is just a hint, you can use whatever you want.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.21
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    std::vector<StateStep> steps_1, steps_2, steps_3, steps_4;
    Jelly tmp;
    for (int jellyIdx = 0; jellyIdx < particleSystem.jellyCount; jellyIdx++) {
        Jelly* J = particleSystem.getJellyPointer(jellyIdx);
        tmp = *J;
        StateStep step;
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            Eigen::Vector3f Pos = P.getPosition();
            step.deltaPos = Pos;
            step.deltaVel = Vel;
            steps_1.push_back(step);
            P.addPosition(Vel * particleSystem.deltaTime);
            P.addVelocity(Acc * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*J);
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            Eigen::Vector3f Pos = P.getPosition();
            step.deltaPos = Pos;
            step.deltaVel = Vel;
            steps_2.push_back(step);
            P.setPosition(steps_1[particleIdx].deltaPos + Vel * 0.5 * particleSystem.deltaTime);
            P.setVelocity(steps_1[particleIdx].deltaVel + Acc * 0.5 * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*J);
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            Eigen::Vector3f Pos = P.getPosition();
            step.deltaPos = Pos;
            step.deltaVel = Vel;
            steps_3.push_back(step);
            P.setPosition(steps_2[particleIdx].deltaPos + Vel * 0.5 * particleSystem.deltaTime);
            P.setVelocity(steps_2[particleIdx].deltaVel + Acc * 0.5 * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*J);
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            Eigen::Vector3f Pos = P.getPosition();
            step.deltaPos = Pos;
            step.deltaVel = Vel;
            steps_4.push_back(step);
            P.setPosition(steps_3[particleIdx].deltaPos + Vel * particleSystem.deltaTime);
            P.setVelocity(steps_3[particleIdx].deltaVel + Acc * particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
        particleSystem.computeJellyForce(*J);
        for (int particleIdx = 0; particleIdx < J->getParticleNum(); particleIdx++) {
            Particle& P = J->getParticle(particleIdx);
            Eigen::Vector3f Vel = P.getVelocity();
            Eigen::Vector3f Acc = P.getAcceleration();
            P.setPosition(tmp.getParticle(particleIdx).getPosition() +
                          1 / 6 *
                              (steps_1[particleIdx].deltaPos + steps_2[particleIdx].deltaPos +
                               steps_3[particleIdx].deltaPos + steps_4[particleIdx].deltaPos) +
                          Vel * particleSystem.deltaTime);
            P.setVelocity(tmp.getParticle(particleIdx).getVelocity() +
                          1 / 6 *
                              (steps_1[particleIdx].deltaVel + 2*steps_2[particleIdx].deltaVel +
                               2 * steps_3[particleIdx].deltaVel + steps_4[particleIdx].deltaVel) + Acc*
                              particleSystem.deltaTime);
            P.setForce(Eigen::Vector3f::Zero());
        }
    }
    }
}  // namespace simulation
