#include "SdfVolume.h"

namespace chrono_sdf_contact {

SdfVolume::Sample SdfVolume::SampleLocal(const chrono::ChVector3d& p_local) const {
    (void)p_local;
    return {1.0, chrono::ChVector3d(1.0, 0.0, 0.0)};
}

SdfVolume::Sample SdfVolume::SampleBodyPose(const chrono::ChVector3d& body_pos,
                                            const chrono::ChQuaternion<>& body_rot,
                                            const chrono::ChVector3d& p_world) const {
    const auto p_local = body_rot.RotateBack(p_world - body_pos);
    auto sample = SampleLocal(p_local);
    sample.grad = body_rot.Rotate(sample.grad);
    if (sample.grad.Length2() > 0.0) {
        sample.grad.Normalize();
    }
    return sample;
}

SdfVolume::Sample SdfVolume::SampleWorld(const chrono::ChFrame<>& body_frame, const chrono::ChVector3d& p_world) const {
    return SampleBodyPose(body_frame.GetPos(), body_frame.GetRot(), p_world);
}

}  // namespace chrono_sdf_contact
