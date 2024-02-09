#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>
#include "acclaim/bone.h"
#include "util/helper.h"
#include <iostream>

namespace kinematics {
void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO#1 (FK)
    // You should set these variables:
    //     bone->start_position = Eigen::Vector4d::Zero();
    //     bone->end_position = Eigen::Vector4d::Zero();
    //     bone->rotation = Eigen::Matrix4d::Zero();
    // The sample above just set everything to zero
    // Hint:
    //   1. posture.bone_translations, posture.bone_rotations
    // Note:
    //   1. This function will be called with bone == root bone of the skeleton
    //   2. we use 4D vector to represent 3D vector, so keep the last dimension as "0"
    //   3. util::rotate{Degree | Radian} {XYZ | ZYX}
    //      e.g. rotateDegreeXYZ(x, y, z) means:
    //      x, y and z are presented in degree rotate z degrees along z - axis first, then y degrees along y - axis, and
    //      then x degrees along x - axis

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

std::vector<acclaim::Posture> timeWarper(const std::vector<acclaim::Posture>& postures, int allframe_old, int allframe_new) {

    int total_frames = static_cast<int>(postures.size());
    int total_bones = static_cast<int>(postures[0].bone_rotations.size());
    std::vector<acclaim::Posture> new_postures;
    for (int i = 0; i <= allframe_new; ++i) {
        acclaim::Posture new_poseture(total_bones);
        for (int j = 0; j < total_bones; ++j) {

            // TODO#2 (Time warping)
            // original: |--------------|
            // new     : |----------------------|
            // OR
            // original: |--------------|
            // new     : |-------|
            // You should set these variables:
            //     new_postures[i].bone_translations[j] = Eigen::Vector4d::Zero();
            //     new_postures[i].bone_rotations[j] = Eigen::Vector4d::Zero();
            // The sample above just set everything to zero
            // Hint:
            //   1. Scale the frames.
            //   2. You can use linear interpolation with translations.
            //   3. You should use spherical linear interpolation for rotations.

            new_poseture.bone_translations[j] = Eigen::Vector4d::Zero();
            new_poseture.bone_rotations[j] = Eigen::Vector4d::Zero();

            double scale_factor = static_cast<double>(allframe_old) / allframe_new;
            double old_frame_idx = i * scale_factor;
            int left_frame_idx = static_cast<int>(old_frame_idx);
            int right_frame_idx = std::min(left_frame_idx + 1, total_frames - 1);      
            double alpha = old_frame_idx - left_frame_idx; //right_frame_idx - left_frame_idx = 1

            Eigen::Vector4d left_translation = postures[left_frame_idx].bone_translations[j];
            Eigen::Vector4d right_translation = postures[right_frame_idx].bone_translations[j];
            Eigen::Vector4d new_translation = (1 - alpha) * left_translation + alpha * right_translation;
            new_poseture.bone_translations[j] = new_translation;

            Eigen::Quaterniond left_rotation(postures[left_frame_idx].bone_rotations[j]);
            Eigen::Quaterniond right_rotation(postures[right_frame_idx].bone_rotations[j]);
            Eigen::Quaterniond new_rotation = left_rotation.slerp(alpha, right_rotation);
            new_poseture.bone_rotations[j] = new_rotation.coeffs();
        }
        new_postures.push_back(new_poseture);
    }
    return new_postures;
}
}  // namespace 
