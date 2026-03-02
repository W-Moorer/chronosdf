#include <algorithm>
#include <iostream>
#include <memory>

#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/core/ChDataPath.h"
#include "chrono/core/ChTypes.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBody.h"
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
    SetChronoDataPath(CHRONO_DATA_DIR);
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0.0, -9.81, 0.0));

    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetFriction(0.4f);
    material->SetRestitution(0.0f);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(3.0, 0.2, 3.0, 1000.0, true, true, material);
    ground->SetFixed(true);
    ground->SetPos(ChVector3d(0.0, -0.1, 0.0));
    sys.AddBody(ground);

    auto mesh_body = chrono_types::make_shared<ChBody>();
    mesh_body->SetMass(20.0);
    mesh_body->SetInertiaXX(ChVector3d(0.2, 0.2, 0.2));
    mesh_body->SetPos(ChVector3d(0.0, 0.9, 0.0));
    mesh_body->EnableCollision(true);

    auto cube_mesh = ChTriangleMeshConnected::CreateFromWavefrontFile(GetChronoDataFile("models/cube.obj"), true, true);
    if (!cube_mesh) {
        std::cerr << "Failed to load mesh: models/cube.obj" << std::endl;
        return 2;
    }

    auto tri_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(material, cube_mesh, false, false, 0.0);
    mesh_body->AddCollisionShape(tri_shape);
    sys.AddBody(mesh_body);

    auto registry = std::make_shared<SdfRegistry>();
    auto pair_cache = std::make_shared<SdfPairCache>();
    auto narrowphase = std::make_shared<SdfNarrowphase>();
    auto clustering = std::make_shared<SdfClustering>(64, 0.06);

    auto narrow_cb = std::make_shared<SdfNarrowphaseCallback>(registry, pair_cache);
    auto custom_cb = std::make_shared<SdfCustomCollisionCallback>(registry, pair_cache, narrowphase, clustering);

    SdfProxy ground_proxy;
    ground_proxy.sdf = std::make_shared<PlaneSdfVolume>(0.1);
    ground_proxy.material = material;
    ground_proxy.envelope = 0.0;
    registry->Register(ground->GetCollisionModel().get(), ground_proxy);

    sys.GetCollisionSystem()->RegisterNarrowphaseCallback(narrow_cb);
    sys.RegisterCustomCollisionCallback(custom_cb);

    constexpr int kNumSteps = 1600;
    constexpr double kStepSize = 1e-3;
    unsigned int max_contacts = 0;
    double min_mesh_height = mesh_body->GetPos().y();
    for (int step = 0; step < kNumSteps; ++step) {
        pair_cache->BeginStep();
        sys.DoStepDynamics(kStepSize);
        max_contacts = std::max(max_contacts, sys.GetNumContacts());
        min_mesh_height = std::min(min_mesh_height, mesh_body->GetPos().y());
    }

    std::cout << "SDF vs trianglemesh demo finished." << std::endl;
    std::cout << "Final contact count: " << sys.GetNumContacts() << std::endl;
    std::cout << "Max contact count: " << max_contacts << std::endl;
    std::cout << "Min mesh center height: " << min_mesh_height << std::endl;
    std::cout << "Final mesh center height: " << mesh_body->GetPos().y() << std::endl;

    if (max_contacts == 0) {
        std::cerr << "Validation failed: no SDF-trianglemesh contacts were generated." << std::endl;
        return 3;
    }

    return 0;
}
