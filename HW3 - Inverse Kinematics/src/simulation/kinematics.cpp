#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>
#include "acclaim/bone.h"
#include "util/helper.h"
#include <iostream>
namespace kinematics {

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO (FK)
    // Same as HW2
    // Hint:
    //   1. If you don't use `axis` in this function, you can copy-paste your code

    bone->start_position = Eigen::Vector4d::Zero();
    bone->end_position = Eigen::Vector4d::Zero();
    bone->rotation = Eigen::Matrix4d::Zero();

    Eigen::Quaterniond R = util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);

    if (bone->name == "root") {
        bone->rotation = R.toRotationMatrix();
        bone->start_position = posture.bone_translations[bone->idx];
        bone->end_position = bone->start_position;
    } else {
        bone->rotation = bone->parent->rotation * bone->rot_parent_current * R.toRotationMatrix();
        bone->start_position = bone->parent->end_position + posture.bone_translations[bone->idx];
        bone->end_position = bone->start_position + bone->rotation * bone->dir * bone->length;
    }

    if (bone->child) forwardSolver(posture, bone->child);
    if (bone->sibling) forwardSolver(posture, bone->sibling);
}

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    // TODO (find x which min(| jacobian * x - target |))
    // Hint:
    //   1. Linear algebra - least squares solution
    //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
    // Note:
    //   1. SVD or other pseudo-inverse method is useful
    //   2. Some of them have some limitation, if you use that method you should check it.
    Eigen::VectorXd deltatheta(Jacobian.cols());
    deltatheta.setZero();
    deltatheta = Jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(target);
    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW3 bonus)
 */

bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
                             acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;

    // TODO
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    size_t bone_num = 0;
    std::vector<acclaim::Bone*> boneList;
    // TODO
    // Calculate number of bones need to move to perform IK, store in `bone_num` 
    // a.k.a. how may bones from end_bone to its parent then to start_bone (include both start_bone and end_bone)
    // Store the bones need to move to perform IK into boneList
    // Hint:
    //   1. Traverse from end_bone to start_bone is easier than start to end (since there is only 1 parent)
    //   2. If start bone is not reachable from end. Go to root first.
    // Note:
    //   1. Both start_bone and end_bone should be in the list
    acclaim::Bone* current = end_bone;
    bool visited = false;
    while (current != root_bone) {
        boneList.push_back(current);
        if (current == start_bone) {
            visited = true;
            break;
        }
        current = current->parent;
    }
    if (!visited) {
        while (current != root_bone) {
            boneList.push_back(current);
            current = current->parent;
        }
    }
    bone_num = boneList.size();

    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    bool stable;
    for (int iter = 0; iter < max_iteration; ++iter) {
        stable = true;
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
        
        if (desiredVector.norm() < epsilon) {
            break;
        }

        // TODO (compute jacobian)
        //   1. Compute arm vectors
        //   2. Compute jacobian columns, store in `Jacobian`
        // Hint:
        //   1. You should not put rotation in jacobian if it doesn't have that DoF.
        //   2. jacobian.col(/* some column index */) = /* jacobian column */
        Jacobian.setZero();
        for (long long i = 0; i < bone_num; i++) {
            Eigen::Vector3d arm = end_bone->end_position.head<3>() - boneList[i]->start_position.head<3>();
            Eigen::Affine3d rotation = boneList[i]->rotation;
            if (boneList[i]->dofrx) {
                Eigen::Vector3d unit_rotation = rotation.matrix().col(0).head<3>();  // x
                Eigen::Vector3d J = unit_rotation.cross(arm);
                Jacobian.col(3 * i) = Eigen::Vector4d(J[0], J[1], J[2], 0);
            }
            if (boneList[i]->dofry) {
                Eigen::Vector3d unit_rotation = rotation.matrix().col(1).head<3>();  // y
                Eigen::Vector3d J = unit_rotation.cross(arm);
                Jacobian.col(3 * i + 1) = Eigen::Vector4d(J[0], J[1], J[2], 0);
            }
            if (boneList[i]->dofrz) {
                Eigen::Vector3d unit_rotation = rotation.matrix().col(2).head<3>();  // z
                Eigen::Vector3d J = unit_rotation.cross(arm);
                Jacobian.col(3 * i + 2) = Eigen::Vector4d(J[0], J[1], J[2], 0);
            } 
        }

        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);

        // TODO (update rotation)
        //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
        // Hint:
        //   1. You can ignore rotation limit of the bone.
        // Bonus:
        //   1. You cannot ignore rotation limit of the bone.
        for (long long i = 0; i < bone_num; i++) {
            acclaim::Bone curr = *boneList[i];
            Eigen::Vector3d delta = deltatheta.segment(i * 3, 3);
            posture.bone_rotations[curr.idx] += util::toDegree(Eigen::Vector4d(delta[0], delta[1], delta[2], 0));
            
            //-----Bonus-----//
            if (posture.bone_rotations[curr.idx][0] < curr.rxmin) {
                stable = false;
                posture.bone_rotations[curr.idx][0] = curr.rxmin;
            } else if (posture.bone_rotations[curr.idx][0] > curr.rxmax) {
                stable = false;
                posture.bone_rotations[curr.idx][0] = curr.rxmax;
            }
            if (posture.bone_rotations[curr.idx][1] < curr.rymin) {
                stable = false;
                posture.bone_rotations[curr.idx][1] = curr.rymin;
            } else if (posture.bone_rotations[curr.idx][1] > curr.rymax) {
                stable = false;
                posture.bone_rotations[curr.idx][1] = curr.rymax;
            } 

            if (posture.bone_rotations[curr.idx][2] < curr.rzmin) {
                stable = false;
                posture.bone_rotations[curr.idx][2] = curr.rzmin;
            } else if (posture.bone_rotations[curr.idx][2] > curr.rzmax) {
                stable = false;
                posture.bone_rotations[curr.idx][2] = curr.rzmax;
            } 
            //---------------//
        }
    }
    // TODO (Bonus)
    // Return whether IK is stable (i.e. whether the ball is reachable) and let the skeleton not swing its hand in the air
    
    if (!stable) {
        posture = original_posture;
        forwardSolver(posture, root_bone);
    }
    return stable;
}
}  // namespace kinematics
