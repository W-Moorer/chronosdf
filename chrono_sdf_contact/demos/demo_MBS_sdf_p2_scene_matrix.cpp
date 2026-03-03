#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "chrono/collision/ChCollisionShapeBox.h"
#include "chrono/core/ChFrame.h"
#include "chrono/core/ChTypes.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
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

class BoxSdfVolume : public chrono_sdf_contact::SdfVolume {
  public:
    explicit BoxSdfVolume(const chrono::ChVector3d& half_lengths) : h_(half_lengths) {}

    Sample SampleLocal(const chrono::ChVector3d& p_local) const override {
        return SampleOffsetBox(p_local, chrono::VNULL, h_);
    }

  protected:
    static Sample SampleOffsetBox(const chrono::ChVector3d& p_local,
                                  const chrono::ChVector3d& center,
                                  const chrono::ChVector3d& half_lengths) {
        const auto p = p_local - center;
        const chrono::ChVector3d q(std::abs(p.x()) - half_lengths.x(), std::abs(p.y()) - half_lengths.y(),
                                   std::abs(p.z()) - half_lengths.z());

        const chrono::ChVector3d q_pos(std::max(q.x(), 0.0), std::max(q.y(), 0.0), std::max(q.z(), 0.0));
        const double outside_len = q_pos.Length();
        const double inside_term = std::min(std::max(q.x(), std::max(q.y(), q.z())), 0.0);
        const double phi = outside_len + inside_term;

        chrono::ChVector3d grad_local(1.0, 0.0, 0.0);
        if (outside_len > 1e-12) {
            grad_local = chrono::ChVector3d(Sign(p.x()) * q_pos.x(), Sign(p.y()) * q_pos.y(), Sign(p.z()) * q_pos.z()) /
                         outside_len;
        } else {
            if (q.y() >= q.x() && q.y() >= q.z()) {
                grad_local = chrono::ChVector3d(0.0, Sign(p.y()), 0.0);
            } else if (q.z() >= q.x() && q.z() >= q.y()) {
                grad_local = chrono::ChVector3d(0.0, 0.0, Sign(p.z()));
            } else {
                grad_local = chrono::ChVector3d(Sign(p.x()), 0.0, 0.0);
            }
        }
        return {phi, grad_local};
    }

  private:
    static double Sign(double v) {
        if (v > 0.0) {
            return 1.0;
        }
        if (v < 0.0) {
            return -1.0;
        }
        return 0.0;
    }

    chrono::ChVector3d h_;
};

class NonconvexTSdfVolume : public chrono_sdf_contact::SdfVolume {
  public:
    NonconvexTSdfVolume(const chrono::ChVector3d& stem_center,
                        const chrono::ChVector3d& stem_half,
                        const chrono::ChVector3d& cap_center,
                        const chrono::ChVector3d& cap_half)
        : stem_center_(stem_center), stem_half_(stem_half), cap_center_(cap_center), cap_half_(cap_half) {}

    Sample SampleLocal(const chrono::ChVector3d& p_local) const override {
        const auto stem = SampleBox(p_local, stem_center_, stem_half_);
        const auto cap = SampleBox(p_local, cap_center_, cap_half_);
        return (stem.phi <= cap.phi) ? stem : cap;
    }

  private:
    static double Sign(double v) {
        if (v > 0.0) {
            return 1.0;
        }
        if (v < 0.0) {
            return -1.0;
        }
        return 0.0;
    }

    static Sample SampleBox(const chrono::ChVector3d& p_local,
                            const chrono::ChVector3d& center,
                            const chrono::ChVector3d& half_lengths) {
        const auto p = p_local - center;
        const chrono::ChVector3d q(std::abs(p.x()) - half_lengths.x(), std::abs(p.y()) - half_lengths.y(),
                                   std::abs(p.z()) - half_lengths.z());
        const chrono::ChVector3d q_pos(std::max(q.x(), 0.0), std::max(q.y(), 0.0), std::max(q.z(), 0.0));
        const double outside_len = q_pos.Length();
        const double inside_term = std::min(std::max(q.x(), std::max(q.y(), q.z())), 0.0);
        const double phi = outside_len + inside_term;

        chrono::ChVector3d grad_local(1.0, 0.0, 0.0);
        if (outside_len > 1e-12) {
            grad_local = chrono::ChVector3d(Sign(p.x()) * q_pos.x(), Sign(p.y()) * q_pos.y(), Sign(p.z()) * q_pos.z()) /
                         outside_len;
        } else {
            if (q.y() >= q.x() && q.y() >= q.z()) {
                grad_local = chrono::ChVector3d(0.0, Sign(p.y()), 0.0);
            } else if (q.z() >= q.x() && q.z() >= q.y()) {
                grad_local = chrono::ChVector3d(0.0, 0.0, Sign(p.z()));
            } else {
                grad_local = chrono::ChVector3d(Sign(p.x()), 0.0, 0.0);
            }
        }
        return {phi, grad_local};
    }

    chrono::ChVector3d stem_center_;
    chrono::ChVector3d stem_half_;
    chrono::ChVector3d cap_center_;
    chrono::ChVector3d cap_half_;
};

class MinDistanceReporter : public chrono::ChContactContainer::ReportContactCallback {
  public:
    void Reset() {
        min_distance_ = std::numeric_limits<double>::infinity();
        has_contact_ = false;
    }

    bool OnReportContact(const chrono::ChVector3d& pA,
                         const chrono::ChVector3d& pB,
                         const chrono::ChMatrix33<>& plane_coord,
                         double distance,
                         double eff_radius,
                         const chrono::ChVector3d& react_forces,
                         const chrono::ChVector3d& react_torques,
                         chrono::ChContactable* contactobjA,
                         chrono::ChContactable* contactobjB,
                         int constraint_offset) override {
        (void)pA;
        (void)pB;
        (void)plane_coord;
        (void)eff_radius;
        (void)react_forces;
        (void)react_torques;
        (void)contactobjA;
        (void)contactobjB;
        (void)constraint_offset;
        has_contact_ = true;
        min_distance_ = std::min(min_distance_, distance);
        return true;
    }

    bool HasContact() const { return has_contact_; }
    double GetMinDistance() const { return min_distance_; }

  private:
    double min_distance_ = std::numeric_limits<double>::infinity();
    bool has_contact_ = false;
};

enum class P2SceneId {
    BoxRest = 0,
    BoxImpact = 1,
    NonconvexDrop = 2,
    Stack3 = 3,
};

const char* SceneName(P2SceneId scene) {
    switch (scene) {
        case P2SceneId::BoxRest:
            return "box_plane_rest";
        case P2SceneId::BoxImpact:
            return "box_plane_impact";
        case P2SceneId::NonconvexDrop:
            return "nonconvex_t_drop";
        case P2SceneId::Stack3:
            return "stack3_boxes";
        default:
            return "unknown_scene";
    }
}

std::array<P2SceneId, 4> AllScenes() {
    return {P2SceneId::BoxRest, P2SceneId::BoxImpact, P2SceneId::NonconvexDrop, P2SceneId::Stack3};
}

struct RunResult {
    std::string scene;
    std::string mode;
    double step_size = 0.0;
    int steps = 0;
    int dynamic_body_count = 0;
    int clustering_enabled = 1;
    int cache_enabled = 1;
    int bound_cull_enabled = 1;
    double wall_time_s = 0.0;
    double avg_contacts = 0.0;
    unsigned int max_contacts = 0;
    double max_penetration = 0.0;
    double min_height = 0.0;
    double final_height = 0.0;
    double max_relative_energy_drift = 0.0;
    bool exploded = false;
};

struct TrajectoryRow {
    std::string scene;
    std::string mode;
    double step_size = 0.0;
    int step = 0;
    double time = 0.0;
    int body_id = 0;
    double com_x = 0.0;
    double com_y = 0.0;
    double com_z = 0.0;
};

struct SceneMatrixConfig {
    std::string csv_path = "sdf_p2_scene_matrix.csv";
    std::string trajectory_csv_path;
    double duration = 1.2;
    std::vector<double> step_sizes = {1e-3, 2e-3, 3e-3};
    std::size_t kmax = 64;
    double cell_size = 0.05;
    bool clustering_enabled = true;
    bool cache_enabled = true;
    bool bound_cull_enabled = true;
    std::uint64_t seed = 42;
};

double ComputeBodyMechanicalEnergy(chrono::ChBody* body, const chrono::ChVector3d& gravity_acc) {
    if (!body) {
        return 0.0;
    }
    const double m = body->GetMass();
    const auto v = body->GetPosDt();
    const auto w_loc = body->GetAngVelLocal();
    const auto I_loc = body->GetInertia();
    const double translational_ke = 0.5 * m * v.Length2();
    const double rotational_ke = 0.5 * w_loc.Dot(I_loc * w_loc);
    const double potential = -m * gravity_acc.Dot(body->GetPos());
    return translational_ke + rotational_ke + potential;
}

double ComputeBodiesMechanicalEnergy(const std::vector<std::shared_ptr<chrono::ChBody>>& bodies,
                                     const chrono::ChVector3d& gravity_acc) {
    double total = 0.0;
    for (const auto& b : bodies) {
        total += ComputeBodyMechanicalEnergy(b.get(), gravity_acc);
    }
    return total;
}

bool ParseDouble(const std::string& text, double& out) {
    try {
        out = std::stod(text);
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseSize(const std::string& text, std::size_t& out) {
    try {
        out = static_cast<std::size_t>(std::stoull(text));
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseUInt64(const std::string& text, std::uint64_t& out) {
    try {
        out = static_cast<std::uint64_t>(std::stoull(text));
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseStepSizes(const std::string& text, std::vector<double>& out_steps) {
    out_steps.clear();
    std::size_t start = 0;
    while (start < text.size()) {
        std::size_t end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        auto token = text.substr(start, end - start);
        double value = 0.0;
        if (!ParseDouble(token, value) || value <= 0.0) {
            return false;
        }
        out_steps.push_back(value);
        start = end + 1;
    }
    return !out_steps.empty();
}

void PrintUsage(const char* exe_name) {
    std::cout << "Usage: " << exe_name
              << " [--csv <path>] [--trajectory-csv <path>] [--duration <seconds>] [--dts <h1,h2,...>]"
                 " [--kmax <int>] [--cell-size <meters>]"
                 " [--seed <u64>] [--no-cluster] [--no-cache] [--no-bound-cull]\n";
}

bool ParseArgs(int argc, char* argv[], SceneMatrixConfig& cfg) {
    bool positional_csv_set = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return false;
        }

        if (!arg.empty() && arg[0] != '-') {
            if (!positional_csv_set) {
                cfg.csv_path = arg;
                positional_csv_set = true;
                continue;
            }
            std::cerr << "Unexpected positional argument: " << arg << std::endl;
            return false;
        }

        auto need_value = [&](const char* opt_name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << opt_name << std::endl;
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--csv") {
            const auto* val = need_value("--csv");
            if (!val) {
                return false;
            }
            cfg.csv_path = val;
            continue;
        }

        if (arg == "--trajectory-csv") {
            const auto* val = need_value("--trajectory-csv");
            if (!val) {
                return false;
            }
            cfg.trajectory_csv_path = val;
            continue;
        }

        if (arg == "--duration") {
            const auto* val = need_value("--duration");
            if (!val || !ParseDouble(val, cfg.duration)) {
                std::cerr << "Invalid --duration value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--dts") {
            const auto* val = need_value("--dts");
            if (!val || !ParseStepSizes(val, cfg.step_sizes)) {
                std::cerr << "Invalid --dts value. Example: --dts 0.001,0.002,0.003" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--kmax") {
            const auto* val = need_value("--kmax");
            if (!val || !ParseSize(val, cfg.kmax)) {
                std::cerr << "Invalid --kmax value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--cell-size") {
            const auto* val = need_value("--cell-size");
            if (!val || !ParseDouble(val, cfg.cell_size)) {
                std::cerr << "Invalid --cell-size value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--seed") {
            const auto* val = need_value("--seed");
            if (!val || !ParseUInt64(val, cfg.seed)) {
                std::cerr << "Invalid --seed value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--no-cluster") {
            cfg.clustering_enabled = false;
            continue;
        }

        if (arg == "--no-cache") {
            cfg.cache_enabled = false;
            continue;
        }

        if (arg == "--no-bound-cull") {
            cfg.bound_cull_enabled = false;
            continue;
        }

        std::cerr << "Unknown option: " << arg << std::endl;
        return false;
    }

    if (cfg.duration <= 0.0 || cfg.step_sizes.empty() || cfg.kmax == 0 || cfg.cell_size <= 0.0) {
        std::cerr << "Invalid configuration: duration/dts/kmax/cell-size must be positive." << std::endl;
        return false;
    }
    return true;
}

struct DynamicBodySpec {
    std::shared_ptr<chrono::ChBody> body;
    std::shared_ptr<chrono_sdf_contact::SdfVolume> sdf;
};

std::uint64_t BuildSceneSeed(std::uint64_t base_seed, P2SceneId scene, double step_size) {
    const auto scene_code = static_cast<std::uint64_t>(static_cast<int>(scene) + 1);
    const auto dt_code = static_cast<std::uint64_t>(std::llround(step_size * 1e9));
    return base_seed ^ (scene_code * 0x9e3779b97f4a7c15ULL) ^ (dt_code * 0x94d049bb133111ebULL);
}

void BuildSceneBodies(P2SceneId scene,
                      std::uint64_t seed,
                      std::shared_ptr<chrono::ChContactMaterialSMC> material,
                      chrono::ChSystemSMC& sys,
                      std::vector<DynamicBodySpec>& out_dynamic) {
    using namespace chrono;
    using namespace chrono_sdf_contact;

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> jitter(-0.03, 0.03);

    out_dynamic.clear();

    if (scene == P2SceneId::BoxRest || scene == P2SceneId::BoxImpact) {
        const ChVector3d half(0.25, 0.15, 0.20);
        auto box = chrono_types::make_shared<ChBodyEasyBox>(2.0 * half.x(), 2.0 * half.y(), 2.0 * half.z(), 1200.0, true,
                                                            true, material);
        const double y0 = (scene == P2SceneId::BoxImpact) ? 1.05 : 0.56;
        box->SetPos(ChVector3d(0.04 + jitter(rng), y0, 0.02 + 0.5 * jitter(rng)));
        box->SetRot(QuatFromAngleZ(0.12 + 0.25 * jitter(rng)) * QuatFromAngleX(0.05 + 0.2 * jitter(rng)));
        if (scene == P2SceneId::BoxImpact) {
            box->SetPosDt(ChVector3d(0.8, -3.5, 0.0));
        }
        sys.AddBody(box);
        out_dynamic.push_back({box, std::make_shared<BoxSdfVolume>(half)});
        return;
    }

    if (scene == P2SceneId::NonconvexDrop) {
        const ChVector3d stem_half(0.08, 0.20, 0.08);
        const ChVector3d cap_half(0.24, 0.06, 0.08);
        const ChVector3d stem_center(0.0, -0.04, 0.0);
        const ChVector3d cap_center(0.0, 0.16, 0.0);

        auto body = chrono_types::make_shared<ChBody>();
        body->SetMass(35.0);
        body->SetInertiaXX(ChVector3d(0.55, 0.70, 0.55));
        body->SetPos(ChVector3d(0.0 + 0.5 * jitter(rng), 0.92, 0.0 + 0.5 * jitter(rng)));
        body->SetRot(QuatFromAngleY(0.22 + 0.25 * jitter(rng)) * QuatFromAngleX(0.14 + 0.20 * jitter(rng)));
        body->SetPosDt(ChVector3d(0.5, -1.8, 0.0));
        body->SetAngVelLocal(ChVector3d(0.0, 2.0, 0.0));
        body->EnableCollision(true);

        auto stem_shape = chrono_types::make_shared<ChCollisionShapeBox>(material, stem_half);
        auto cap_shape = chrono_types::make_shared<ChCollisionShapeBox>(material, cap_half);
        body->AddCollisionShape(stem_shape, ChFrame<>(stem_center, QUNIT));
        body->AddCollisionShape(cap_shape, ChFrame<>(cap_center, QUNIT));
        sys.AddBody(body);

        out_dynamic.push_back(
            {body, std::make_shared<NonconvexTSdfVolume>(stem_center, stem_half, cap_center, cap_half)});
        return;
    }

    if (scene == P2SceneId::Stack3) {
        const ChVector3d half(0.22, 0.12, 0.22);
        for (int i = 0; i < 3; ++i) {
            auto box = chrono_types::make_shared<ChBodyEasyBox>(2.0 * half.x(), 2.0 * half.y(), 2.0 * half.z(), 1000.0,
                                                                true, true, material);
            const double y = half.y() + 0.01 + i * (2.0 * half.y() + 0.008);
            const double x = ((i % 2 == 0) ? 0.0 : 0.03) + 0.4 * jitter(rng);
            box->SetPos(ChVector3d(x, y, 0.0 + 0.4 * jitter(rng)));
            box->SetRot(QuatFromAngleZ(0.06 * (i + 1) + 0.2 * jitter(rng)));
            if (i == 2) {
                box->SetPosDt(ChVector3d(0.6, 0.0, 0.0));
            }
            sys.AddBody(box);
            out_dynamic.push_back({box, std::make_shared<BoxSdfVolume>(half)});
        }
    }
}

RunResult RunCase(P2SceneId scene,
                  bool use_sdf_contact,
                  double step_size,
                  const SceneMatrixConfig& cfg,
                  std::vector<TrajectoryRow>* trajectory_rows) {
    using namespace chrono;
    using namespace chrono_sdf_contact;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0.0, -9.81, 0.0));

    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetFriction(0.45f);
    material->SetRestitution(0.0f);
    material->SetYoungModulus(3.0e6f);
    material->SetKn(6.0e4f);
    material->SetGn(240.0f);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(6.0, 0.2, 6.0, 1000.0, true, true, material);
    ground->SetFixed(true);
    ground->SetPos(ChVector3d(0.0, -0.1, 0.0));
    sys.AddBody(ground);

    std::vector<DynamicBodySpec> dynamic_specs;
    BuildSceneBodies(scene, BuildSceneSeed(cfg.seed, scene, step_size), material, sys, dynamic_specs);

    std::shared_ptr<SdfPairCache> pair_cache;
    if (use_sdf_contact) {
        auto registry = std::make_shared<SdfRegistry>();
        pair_cache = std::make_shared<SdfPairCache>();
        auto narrowphase = std::make_shared<SdfNarrowphase>();
        narrowphase->SetCacheEnabled(cfg.cache_enabled);
        narrowphase->SetBoundCullEnabled(cfg.bound_cull_enabled);
        auto clustering = std::make_shared<SdfClustering>(cfg.kmax, cfg.cell_size);
        auto narrow_cb = std::make_shared<SdfNarrowphaseCallback>(registry, pair_cache);
        auto custom_cb = std::make_shared<SdfCustomCollisionCallback>(registry, pair_cache, narrowphase, clustering);
        custom_cb->SetClusteringEnabled(cfg.clustering_enabled);

        SdfProxy ground_proxy;
        ground_proxy.sdf = std::make_shared<PlaneSdfVolume>(0.1);
        ground_proxy.material = material;
        registry->Register(ground->GetCollisionModel().get(), ground_proxy);

        for (const auto& spec : dynamic_specs) {
            if (!spec.body || !spec.sdf) {
                continue;
            }
            SdfProxy body_proxy;
            body_proxy.sdf = spec.sdf;
            body_proxy.material = material;
            registry->Register(spec.body->GetCollisionModel().get(), body_proxy);
        }

        sys.GetCollisionSystem()->RegisterNarrowphaseCallback(narrow_cb);
        sys.RegisterCustomCollisionCallback(custom_cb);
    }

    std::vector<std::shared_ptr<ChBody>> dynamic_bodies;
    dynamic_bodies.reserve(dynamic_specs.size());
    for (const auto& spec : dynamic_specs) {
        dynamic_bodies.push_back(spec.body);
    }

    const std::string scene_name = SceneName(scene);
    const std::string mode_name = use_sdf_contact ? "sdf_contact" : "chrono_baseline";
    auto reporter = std::make_shared<MinDistanceReporter>();
    double min_dist = std::numeric_limits<double>::infinity();
    double min_height = std::numeric_limits<double>::infinity();
    unsigned int max_contacts = 0;
    unsigned long long total_contacts = 0;
    bool exploded = false;

    const auto gravity = sys.GetGravitationalAcceleration();
    const double energy0 = ComputeBodiesMechanicalEnergy(dynamic_bodies, gravity);
    double max_relative_energy_drift = 0.0;
    if (trajectory_rows) {
        for (std::size_t i = 0; i < dynamic_bodies.size(); ++i) {
            const auto p = dynamic_bodies[i]->GetPos();
            trajectory_rows->push_back({scene_name, mode_name, step_size, 0, 0.0, static_cast<int>(i), p.x(), p.y(), p.z()});
        }
    }

    const int num_steps = std::max(1, static_cast<int>(cfg.duration / step_size));
    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; ++step) {
        if (pair_cache) {
            pair_cache->BeginStep();
        }
        sys.DoStepDynamics(step_size);

        for (const auto& body : dynamic_bodies) {
            const auto p = body->GetPos();
            if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()) || std::abs(p.y()) > 5.0) {
                exploded = true;
                break;
            }
            min_height = std::min(min_height, p.y());
        }
        if (exploded) {
            break;
        }

        if (trajectory_rows) {
            for (std::size_t i = 0; i < dynamic_bodies.size(); ++i) {
                const auto p = dynamic_bodies[i]->GetPos();
                trajectory_rows->push_back(
                    {scene_name, mode_name, step_size, step + 1, (step + 1) * step_size, static_cast<int>(i), p.x(),
                     p.y(), p.z()});
            }
        }

        const auto ncontacts = sys.GetNumContacts();
        total_contacts += ncontacts;
        max_contacts = std::max(max_contacts, ncontacts);

        const double energy = ComputeBodiesMechanicalEnergy(dynamic_bodies, gravity);
        const double rel_drift = std::abs((energy - energy0) / (std::abs(energy0) + 1e-12));
        max_relative_energy_drift = std::max(max_relative_energy_drift, rel_drift);

        if (ncontacts > 0) {
            reporter->Reset();
            sys.GetContactContainer()->ReportAllContacts(reporter);
            if (reporter->HasContact()) {
                min_dist = std::min(min_dist, reporter->GetMinDistance());
            }
        }
    }
    const auto t1 = std::chrono::high_resolution_clock::now();

    double final_height = std::numeric_limits<double>::infinity();
    for (const auto& body : dynamic_bodies) {
        final_height = std::min(final_height, body->GetPos().y());
    }
    if (!std::isfinite(min_height)) {
        min_height = final_height;
    }

    RunResult out;
    out.scene = scene_name;
    out.mode = mode_name;
    out.step_size = step_size;
    out.steps = num_steps;
    out.dynamic_body_count = static_cast<int>(dynamic_bodies.size());
    out.clustering_enabled = cfg.clustering_enabled ? 1 : 0;
    out.cache_enabled = cfg.cache_enabled ? 1 : 0;
    out.bound_cull_enabled = cfg.bound_cull_enabled ? 1 : 0;
    out.wall_time_s = std::chrono::duration<double>(t1 - t0).count();
    out.avg_contacts = (num_steps > 0) ? static_cast<double>(total_contacts) / static_cast<double>(num_steps) : 0.0;
    out.max_contacts = max_contacts;
    out.max_penetration = (min_dist < 0.0) ? -min_dist : 0.0;
    out.min_height = std::isfinite(min_height) ? min_height : 0.0;
    out.final_height = std::isfinite(final_height) ? final_height : 0.0;
    out.max_relative_energy_drift = max_relative_energy_drift;
    out.exploded = exploded;
    return out;
}

void WriteTrajectoryCsv(const std::string& path, const std::vector<TrajectoryRow>& rows) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open trajectory CSV file: " + path);
    }
    out << "scene,mode,step_size,step,time,body_id,com_x,com_y,com_z\n";
    for (const auto& r : rows) {
        out << r.scene << "," << r.mode << "," << r.step_size << "," << r.step << "," << r.time << "," << r.body_id
            << "," << r.com_x << "," << r.com_y << "," << r.com_z << "\n";
    }
}

void PrintResult(const RunResult& r) {
    std::cout << "[" << r.scene << ", " << r.mode << ", dt=" << r.step_size << "] "
              << "wall_time_s=" << std::fixed << std::setprecision(6) << r.wall_time_s
              << ", avg_contacts=" << r.avg_contacts << ", max_contacts=" << r.max_contacts
              << ", max_penetration=" << r.max_penetration << ", min_height=" << r.min_height
              << ", final_height=" << r.final_height
              << ", max_relative_energy_drift=" << r.max_relative_energy_drift
              << ", exploded=" << (r.exploded ? 1 : 0)
              << ", clustering=" << r.clustering_enabled
              << ", cache=" << r.cache_enabled
              << ", bound_cull=" << r.bound_cull_enabled << std::endl;
}

void WriteCsv(const std::string& path, const std::vector<RunResult>& rows) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open CSV file: " + path);
    }
    out << "scene,mode,step_size,steps,dynamic_body_count,clustering_enabled,cache_enabled,bound_cull_enabled,wall_time_s,avg_contacts,max_contacts,max_penetration,min_height,final_height,max_relative_energy_drift,exploded\n";
    for (const auto& r : rows) {
        out << r.scene << "," << r.mode << "," << r.step_size << "," << r.steps << "," << r.dynamic_body_count
            << "," << r.clustering_enabled << "," << r.cache_enabled << "," << r.bound_cull_enabled << ","
            << r.wall_time_s << "," << r.avg_contacts << "," << r.max_contacts << "," << r.max_penetration << ","
            << r.min_height << "," << r.final_height << ","
            << r.max_relative_energy_drift << ","
            << (r.exploded ? 1 : 0) << "\n";
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    SceneMatrixConfig cfg;
    if (!ParseArgs(argc, argv, cfg)) {
        return 1;
    }

    try {
        std::vector<RunResult> results;
        const auto scenes = AllScenes();
        results.reserve(cfg.step_sizes.size() * scenes.size() * 2);
        std::vector<TrajectoryRow> trajectories;

        for (double dt : cfg.step_sizes) {
            for (auto scene : scenes) {
                results.push_back(
                    RunCase(scene, false, dt, cfg, cfg.trajectory_csv_path.empty() ? nullptr : &trajectories));
                results.push_back(
                    RunCase(scene, true, dt, cfg, cfg.trajectory_csv_path.empty() ? nullptr : &trajectories));
            }
        }

        for (const auto& r : results) {
            PrintResult(r);
        }
        WriteCsv(cfg.csv_path, results);
        std::cout << "CSV written to: " << cfg.csv_path << std::endl;
        if (!cfg.trajectory_csv_path.empty()) {
            WriteTrajectoryCsv(cfg.trajectory_csv_path, trajectories);
            std::cout << "Trajectory CSV written to: " << cfg.trajectory_csv_path << std::endl;
        }

        for (double dt : cfg.step_sizes) {
            for (auto scene : scenes) {
                bool ok = false;
                for (const auto& r : results) {
                    if (r.scene == SceneName(scene) && std::abs(r.step_size - dt) < 1e-12 && r.mode == "sdf_contact" &&
                        !r.exploded && r.max_contacts > 0) {
                        ok = true;
                        break;
                    }
                }
                if (!ok) {
                    std::cerr << "Validation failed: no stable SDF contacts for scene=" << SceneName(scene)
                              << ", dt=" << dt << std::endl;
                    return 3;
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    return 0;
}
