#pragma once

#include "chrono/core/ChFrame.h"
#include "chrono/core/ChQuaternion.h"
#include "chrono/core/ChVector3.h"

namespace chrono_sdf_contact {

class SdfVolume {
  public:
    struct Sample {
        double phi;
        chrono::ChVector3d grad;
    };

    virtual ~SdfVolume() = default;

    virtual Sample SampleLocal(const chrono::ChVector3d& p_local) const;
    virtual Sample SampleBodyPose(const chrono::ChVector3d& body_pos,
                                  const chrono::ChQuaternion<>& body_rot,
                                  const chrono::ChVector3d& p_world) const;
    virtual Sample SampleWorld(const chrono::ChFrame<>& body_frame, const chrono::ChVector3d& p_world) const;
};

}  // namespace chrono_sdf_contact
