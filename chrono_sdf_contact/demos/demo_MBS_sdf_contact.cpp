#include <algorithm>
#include <iostream>
#include <memory>

#include "chrono/core/ChTypes.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactMaterialSMC.h"
#include "chrono/physics/ChSystemSMC.h"

#include "SdfClustering.h"
#include "SdfCustomCollisionCallback.h"
#include "SdfNarrowphase.h"
#include "SdfNarrowphaseCallback.h"
#include "SdfPairCache.h"
#include "SdfRegistry.h"
#include "SdfVolume.h"

namespace {

class PlaneSdfVolume : public chrono_sdf_contact::SdfVolume {
  public:
    explicit PlaneSdfVolume(double plane_y_local) : plane_y_local_(plane_y_local) {}

    Sample SampleLocal(const chrono::ChVector3d& p_local) const override {
        return {p_local.y() - plane_y_local_, chrono::ChVector3d(0.0, 1.0, 0.0)};
    }

  private:
    double plane_y_local_;
};

}  // namespace

int main() {
    using namespace chrono;
    using namespace chrono_sdf_contact;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0.0, -9.81, 0.0));

    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetFriction(0.4f);
    material->SetRestitution(0.0f);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(2.0, 0.2, 2.0, 1000.0, true, true, material);
    ground->SetFixed(true);
    ground->SetPos(ChVector3d(0.0, -0.1, 0.0));
    sys.AddBody(ground);

    auto sphere = chrono_types::make_shared<ChBodyEasySphere>(0.1, 1000.0, true, true, material);
    sphere->SetPos(ChVector3d(0.0, 0.6, 0.0));
    sys.AddBody(sphere);

    auto box = chrono_types::make_shared<ChBodyEasyBox>(0.16, 0.16, 0.16, 1000.0, true, true, material);
    box->SetPos(ChVector3d(0.35, 0.65, 0.0));
    sys.AddBody(box);

    auto registry = std::make_shared<SdfRegistry>();
    auto pair_cache = std::make_shared<SdfPairCache>();
    auto narrowphase = std::make_shared<SdfNarrowphase>();
    auto clustering = std::make_shared<SdfClustering>(32);

    auto narrow_cb = std::make_shared<SdfNarrowphaseCallback>(registry, pair_cache);
    auto custom_cb = std::make_shared<SdfCustomCollisionCallback>(registry, pair_cache, narrowphase, clustering);

    // Ground top surface in local coordinates is y = +0.1 (box half-height).
    SdfProxy ground_proxy;
    ground_proxy.sdf = std::make_shared<PlaneSdfVolume>(0.1);
    ground_proxy.material = material;
    ground_proxy.envelope = 0.0;
    registry->Register(ground->GetCollisionModel().get(), ground_proxy);

    sys.GetCollisionSystem()->RegisterNarrowphaseCallback(narrow_cb);
    sys.RegisterCustomCollisionCallback(custom_cb);

    constexpr int kNumSteps = 1200;
    constexpr double kStepSize = 1e-3;
    double min_sphere_height = sphere->GetPos().y();
    double min_box_height = box->GetPos().y();
    unsigned int max_contacts = 0;
    for (int step = 0; step < kNumSteps; ++step) {
        pair_cache->BeginStep();
        sys.DoStepDynamics(kStepSize);
        min_sphere_height = std::min(min_sphere_height, sphere->GetPos().y());
        min_box_height = std::min(min_box_height, box->GetPos().y());
        max_contacts = std::max(max_contacts, sys.GetNumContacts());
    }

    std::cout << "SDF contact skeleton demo finished." << std::endl;
    std::cout << "Final contact count: " << sys.GetNumContacts() << std::endl;
    std::cout << "Max contact count: " << max_contacts << std::endl;
    std::cout << "Min sphere height: " << min_sphere_height << std::endl;
    std::cout << "Final sphere height: " << sphere->GetPos().y() << std::endl;
    std::cout << "Min box center height: " << min_box_height << std::endl;
    std::cout << "Final box center height: " << box->GetPos().y() << std::endl;
    return 0;
}
