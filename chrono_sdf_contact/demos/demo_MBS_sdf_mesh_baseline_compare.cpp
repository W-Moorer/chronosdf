#include <algorithm>
#include <cstdint>
#include <cstddef>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>

#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/core/ChDataPath.h"
#include "chrono/core/ChTypes.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/physics/ChContactMaterialSMC.h"
#include "chrono/physics/ChMassProperties.h"
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

struct RunResult {
    std::string mode;
    double wall_time_s = 0.0;
    double avg_contacts = 0.0;
    unsigned int max_contacts = 0;
    double max_penetration = 0.0;
    double min_mesh_height = 0.0;
    double final_mesh_height = 0.0;
    double max_relative_energy_drift = 0.0;
};

struct TrajectoryRow {
    std::string mode;
    int step = 0;
    double time = 0.0;
    int body_id = 0;
    double com_x = 0.0;
    double com_y = 0.0;
    double com_z = 0.0;
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

struct MeshCompareConfig {
    std::string csv_path = "sdf_mesh_baseline_compare.csv";
    std::string trajectory_csv_path;
    int num_steps = 1200;
    double step_size = 1e-3;
    std::size_t kmax = 64;
    double cell_size = 0.06;
    std::uint64_t seed = 42;
};

bool ParseInt(const std::string& text, int& out) {
    try {
        const auto value = std::stoi(text);
        out = value;
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseSize(const std::string& text, std::size_t& out) {
    try {
        const auto value = std::stoull(text);
        out = static_cast<std::size_t>(value);
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseDouble(const std::string& text, double& out) {
    try {
        const auto value = std::stod(text);
        out = value;
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

void BuildMeshInitialPose(std::uint64_t seed, chrono::ChVector3d& pos, chrono::ChQuaternion<>& rot) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> xz_jitter(-0.025, 0.025);
    std::uniform_real_distribution<double> yaw_jitter(-0.12, 0.12);
    std::uniform_real_distribution<double> pitch_jitter(-0.05, 0.05);
    std::uniform_real_distribution<double> roll_jitter(-0.05, 0.05);

    pos = chrono::ChVector3d(xz_jitter(rng), 0.7, xz_jitter(rng));
    rot = chrono::QuatFromAngleY(yaw_jitter(rng)) * chrono::QuatFromAngleX(pitch_jitter(rng)) *
          chrono::QuatFromAngleZ(roll_jitter(rng));
}

void PrintUsage(const char* exe_name) {
    std::cout << "Usage: " << exe_name
              << " [--csv <path>] [--trajectory-csv <path>] [--steps <int>] [--dt <seconds>] [--kmax <int>]"
                 " [--cell-size <meters>] [--seed <u64>]\n";
}

bool ParseArgs(int argc, char* argv[], MeshCompareConfig& cfg) {
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

        if (arg == "--steps") {
            const auto* val = need_value("--steps");
            if (!val || !ParseInt(val, cfg.num_steps)) {
                std::cerr << "Invalid --steps value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--dt") {
            const auto* val = need_value("--dt");
            if (!val || !ParseDouble(val, cfg.step_size)) {
                std::cerr << "Invalid --dt value" << std::endl;
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

        std::cerr << "Unknown option: " << arg << std::endl;
        return false;
    }

    if (cfg.num_steps <= 0 || cfg.step_size <= 0.0 || cfg.kmax == 0 || cfg.cell_size <= 0.0) {
        std::cerr << "Invalid configuration: steps/dt/kmax/cell-size must be positive." << std::endl;
        return false;
    }
    return true;
}

RunResult RunScenario(bool use_sdf_contact, const MeshCompareConfig& cfg, std::vector<TrajectoryRow>* trajectory_rows) {
    using namespace chrono;
    using namespace chrono_sdf_contact;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0.0, -9.81, 0.0));

    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetFriction(0.4f);
    material->SetRestitution(0.0f);
    material->SetYoungModulus(2e6f);
    material->SetKn(4e4f);
    material->SetGn(150.0f);

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(3.0, 0.2, 3.0, 1000.0, true, true, material);
    ground->SetFixed(true);
    ground->SetPos(ChVector3d(0.0, -0.1, 0.0));
    sys.AddBody(ground);

    ChVector3d mesh_init_pos;
    ChQuaternion<> mesh_init_rot;
    BuildMeshInitialPose(cfg.seed, mesh_init_pos, mesh_init_rot);

    auto mesh_body = chrono_types::make_shared<ChBody>();
    mesh_body->SetPos(mesh_init_pos);
    mesh_body->SetRot(mesh_init_rot);
    mesh_body->EnableCollision(true);

    auto cube_mesh = ChTriangleMeshConnected::CreateFromWavefrontFile(GetChronoDataFile("models/cube.obj"), true, true);
    if (!cube_mesh) {
        throw std::runtime_error("Failed to load mesh: models/cube.obj");
    }

    constexpr double kDensity = 800.0;
    double mesh_mass = 0.0;
    ChVector3d mesh_cog = VNULL;
    ChMatrix33<> mesh_inertia;
    cube_mesh->ComputeMassProperties(true, mesh_mass, mesh_cog, mesh_inertia);
    mesh_body->SetMass(std::max(1e-6, mesh_mass * kDensity));
    mesh_body->SetInertiaXX(
        kDensity * ChVector3d(mesh_inertia(0, 0), mesh_inertia(1, 1), mesh_inertia(2, 2)));
    mesh_body->SetInertiaXY(
        kDensity * ChVector3d(mesh_inertia(0, 1), mesh_inertia(0, 2), mesh_inertia(1, 2)));

    auto tri_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(material, cube_mesh, false, false, 0.0);
    mesh_body->AddCollisionShape(tri_shape);
    sys.AddBody(mesh_body);

    std::shared_ptr<SdfPairCache> pair_cache;
    if (use_sdf_contact) {
        auto registry = std::make_shared<SdfRegistry>();
        pair_cache = std::make_shared<SdfPairCache>();
        auto narrowphase = std::make_shared<SdfNarrowphase>();
        auto clustering = std::make_shared<SdfClustering>(cfg.kmax, cfg.cell_size);
        auto narrow_cb = std::make_shared<SdfNarrowphaseCallback>(registry, pair_cache);
        auto custom_cb = std::make_shared<SdfCustomCollisionCallback>(registry, pair_cache, narrowphase, clustering);

        SdfProxy ground_proxy;
        ground_proxy.sdf = std::make_shared<PlaneSdfVolume>(0.1);
        ground_proxy.material = material;
        registry->Register(ground->GetCollisionModel().get(), ground_proxy);

        sys.GetCollisionSystem()->RegisterNarrowphaseCallback(narrow_cb);
        sys.RegisterCustomCollisionCallback(custom_cb);
    }

    const std::string mode_name = use_sdf_contact ? "sdf_contact" : "chrono_baseline";
    auto reporter = std::make_shared<MinDistanceReporter>();
    double min_distance = std::numeric_limits<double>::infinity();
    double min_mesh_height = mesh_body->GetPos().y();
    unsigned int max_contacts = 0;
    unsigned long long total_contacts = 0;
    const auto gravity = sys.GetGravitationalAcceleration();
    const double energy0 = ComputeBodyMechanicalEnergy(mesh_body.get(), gravity);
    double max_relative_energy_drift = 0.0;
    if (trajectory_rows) {
        const auto p = mesh_body->GetPos();
        trajectory_rows->push_back({mode_name, 0, 0.0, 0, p.x(), p.y(), p.z()});
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < cfg.num_steps; ++step) {
        if (pair_cache) {
            pair_cache->BeginStep();
        }

        sys.DoStepDynamics(cfg.step_size);

        const auto num_contacts = sys.GetNumContacts();
        total_contacts += num_contacts;
        max_contacts = std::max(max_contacts, num_contacts);
        min_mesh_height = std::min(min_mesh_height, mesh_body->GetPos().y());
        const double energy = ComputeBodyMechanicalEnergy(mesh_body.get(), gravity);
        const double rel_drift = std::abs((energy - energy0) / (std::abs(energy0) + 1e-12));
        max_relative_energy_drift = std::max(max_relative_energy_drift, rel_drift);

        if (num_contacts > 0) {
            reporter->Reset();
            sys.GetContactContainer()->ReportAllContacts(reporter);
            if (reporter->HasContact()) {
                min_distance = std::min(min_distance, reporter->GetMinDistance());
            }
        }

        if (trajectory_rows) {
            const auto p = mesh_body->GetPos();
            trajectory_rows->push_back({mode_name, step + 1, (step + 1) * cfg.step_size, 0, p.x(), p.y(), p.z()});
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    RunResult result;
    result.mode = mode_name;
    result.wall_time_s = std::chrono::duration<double>(t1 - t0).count();
    result.avg_contacts = static_cast<double>(total_contacts) / static_cast<double>(cfg.num_steps);
    result.max_contacts = max_contacts;
    result.max_penetration = (min_distance < 0.0) ? -min_distance : 0.0;
    result.min_mesh_height = min_mesh_height;
    result.final_mesh_height = mesh_body->GetPos().y();
    result.max_relative_energy_drift = max_relative_energy_drift;
    return result;
}

void WriteTrajectoryCsv(const std::string& path, const std::vector<TrajectoryRow>& rows) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open trajectory CSV file: " + path);
    }
    out << "mode,step,time,body_id,com_x,com_y,com_z\n";
    for (const auto& r : rows) {
        out << r.mode << "," << r.step << "," << r.time << "," << r.body_id << "," << r.com_x << "," << r.com_y << ","
            << r.com_z << "\n";
    }
}

void PrintResult(const RunResult& result) {
    std::cout << "[" << result.mode << "] wall_time_s=" << std::fixed << std::setprecision(6) << result.wall_time_s
              << ", avg_contacts=" << result.avg_contacts << ", max_contacts=" << result.max_contacts
              << ", max_penetration=" << result.max_penetration << ", min_mesh_height=" << result.min_mesh_height
              << ", final_mesh_height=" << result.final_mesh_height
              << ", max_relative_energy_drift=" << result.max_relative_energy_drift << std::endl;
}

void WriteCsv(const std::string& path, int num_steps, double step_size, const RunResult& baseline, const RunResult& sdf) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open CSV file: " + path);
    }

    out << "mode,steps,step_size,wall_time_s,avg_contacts,max_contacts,max_penetration,min_mesh_height,final_mesh_height,max_relative_energy_drift\n";
    auto write_row = [&](const RunResult& r) {
        out << r.mode << "," << num_steps << "," << step_size << "," << r.wall_time_s << "," << r.avg_contacts << ","
            << r.max_contacts << "," << r.max_penetration << "," << r.min_mesh_height << "," << r.final_mesh_height
            << "," << r.max_relative_energy_drift << "\n";
    };

    write_row(baseline);
    write_row(sdf);
}

}  // namespace

int main(int argc, char* argv[]) {
    using namespace chrono;

    SetChronoDataPath(CHRONO_DATA_DIR);

    MeshCompareConfig cfg;
    if (!ParseArgs(argc, argv, cfg)) {
        return 1;
    }

    try {
        std::vector<TrajectoryRow> trajectories;
        if (!cfg.trajectory_csv_path.empty()) {
            trajectories.reserve(static_cast<std::size_t>(2 * (cfg.num_steps + 1)));
        }

        auto baseline = RunScenario(false, cfg, cfg.trajectory_csv_path.empty() ? nullptr : &trajectories);
        auto sdf = RunScenario(true, cfg, cfg.trajectory_csv_path.empty() ? nullptr : &trajectories);

        PrintResult(baseline);
        PrintResult(sdf);
        WriteCsv(cfg.csv_path, cfg.num_steps, cfg.step_size, baseline, sdf);
        if (!cfg.trajectory_csv_path.empty()) {
            WriteTrajectoryCsv(cfg.trajectory_csv_path, trajectories);
            std::cout << "Trajectory CSV written to: " << cfg.trajectory_csv_path << std::endl;
        }

        std::cout << "CSV written to: " << cfg.csv_path << std::endl;

        if (baseline.max_contacts == 0 || sdf.max_contacts == 0) {
            std::cerr << "Validation failed: one of the runs generated zero contacts." << std::endl;
            return 3;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    return 0;
}
