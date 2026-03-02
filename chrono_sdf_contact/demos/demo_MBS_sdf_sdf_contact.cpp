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

class SphereSdfVolume : public chrono_sdf_contact::SdfVolume {
  public:
    explicit SphereSdfVolume(double radius) : radius_(radius) {}

    Sample SampleLocal(const chrono::ChVector3d& p_local) const override {
        const double len = p_local.Length();
        chrono::ChVector3d grad(1.0, 0.0, 0.0);
        if (len > 1e-12) {
            grad = p_local / len;
        }
        return {len - radius_, grad};
    }

  private:
    double radius_;
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

    constexpr double kRadius = 0.15;

    auto fixed_sphere = chrono_types::make_shared<ChBodyEasySphere>(kRadius, 1500.0, true, true, material);
    fixed_sphere->SetFixed(true);
    fixed_sphere->SetPos(ChVector3d(0.0, 0.0, 0.0));
    sys.AddBody(fixed_sphere);

    auto moving_sphere = chrono_types::make_shared<ChBodyEasySphere>(kRadius, 1500.0, true, true, material);
    moving_sphere->SetPos(ChVector3d(0.0, 0.45, 0.0));
    sys.AddBody(moving_sphere);

    auto registry = std::make_shared<SdfRegistry>();
    auto pair_cache = std::make_shared<SdfPairCache>();
    auto narrowphase = std::make_shared<SdfNarrowphase>();
    auto clustering = std::make_shared<SdfClustering>(32, 0.04);

    auto narrow_cb = std::make_shared<SdfNarrowphaseCallback>(registry, pair_cache);
    auto custom_cb = std::make_shared<SdfCustomCollisionCallback>(registry, pair_cache, narrowphase, clustering);

    SdfProxy fixed_proxy;
    fixed_proxy.sdf = std::make_shared<SphereSdfVolume>(kRadius);
    fixed_proxy.material = material;
    registry->Register(fixed_sphere->GetCollisionModel().get(), fixed_proxy);

    SdfProxy moving_proxy;
    moving_proxy.sdf = std::make_shared<SphereSdfVolume>(kRadius);
    moving_proxy.material = material;
    registry->Register(moving_sphere->GetCollisionModel().get(), moving_proxy);

    sys.GetCollisionSystem()->RegisterNarrowphaseCallback(narrow_cb);
    sys.RegisterCustomCollisionCallback(custom_cb);

    constexpr int kNumSteps = 1800;
    constexpr double kStepSize = 1e-3;
    unsigned int max_contacts = 0;
    double min_height = moving_sphere->GetPos().y();
    for (int step = 0; step < kNumSteps; ++step) {
        pair_cache->BeginStep();
        sys.DoStepDynamics(kStepSize);
        max_contacts = std::max(max_contacts, sys.GetNumContacts());
        min_height = std::min(min_height, moving_sphere->GetPos().y());
    }

    std::cout << "SDF vs SDF demo finished." << std::endl;
    std::cout << "Final contact count: " << sys.GetNumContacts() << std::endl;
    std::cout << "Max contact count: " << max_contacts << std::endl;
    std::cout << "Min moving sphere center height: " << min_height << std::endl;
    std::cout << "Final moving sphere center height: " << moving_sphere->GetPos().y() << std::endl;

    if (max_contacts == 0) {
        std::cerr << "Validation failed: no SDF-SDF contacts were generated." << std::endl;
        return 3;
    }

    return 0;
}
