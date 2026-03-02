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
#include <vector>

#include "chrono/core/ChTypes.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/physics/ChContactMaterialSMC.h"
#include "chrono/physics/ChSystemSMC.h"

#include "SdfClustering.h"
#include "SdfCustomCollisionCallback.h"
#include "SdfNarrowphase.h"
#include "SdfNarrowphaseCallback.h"
#include "SdfOctreeVolume.h"
#include "SdfPairCache.h"
#include "SdfRegistry.h"
#include "SdfVolume.h"

namespace {

class BoxSdfVolume : public chrono_sdf_contact::SdfVolume {
  public:
    explicit BoxSdfVolume(const chrono::ChVector3d& half_lengths) : h_(half_lengths) {}

    Sample SampleLocal(const chrono::ChVector3d& p_local) const override {
        const chrono::ChVector3d q(std::abs(p_local.x()) - h_.x(), std::abs(p_local.y()) - h_.y(),
                                   std::abs(p_local.z()) - h_.z());

        const chrono::ChVector3d q_pos(std::max(q.x(), 0.0), std::max(q.y(), 0.0), std::max(q.z(), 0.0));
        const double outside_len = q_pos.Length();
        const double inside_term = std::min(std::max(q.x(), std::max(q.y(), q.z())), 0.0);
        const double phi = outside_len + inside_term;

        chrono::ChVector3d grad_local(1.0, 0.0, 0.0);
        if (outside_len > 1e-12) {
            grad_local = chrono::ChVector3d(Sign(p_local.x()) * q_pos.x(), Sign(p_local.y()) * q_pos.y(),
                                            Sign(p_local.z()) * q_pos.z()) /
                         outside_len;
        } else {
            // Inside the box: choose dominant axis of signed distance to faces.
            if (q.y() >= q.x() && q.y() >= q.z()) {
                grad_local = chrono::ChVector3d(0.0, Sign(p_local.y()), 0.0);
            } else if (q.z() >= q.x() && q.z() >= q.y()) {
                grad_local = chrono::ChVector3d(0.0, 0.0, Sign(p_local.z()));
            } else {
                grad_local = chrono::ChVector3d(Sign(p_local.x()), 0.0, 0.0);
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
    double step_size = 0.0;
    int steps = 0;
    double wall_time_s = 0.0;
    double avg_contacts = 0.0;
    unsigned int max_contacts = 0;
    double max_penetration = 0.0;
    double min_height = 0.0;
    double final_height = 0.0;
    double max_relative_energy_drift = 0.0;
    bool clustering_enabled = true;
    bool octree_enabled = false;
    bool exploded = false;
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

struct StabilityConfig {
    std::string csv_path = "sdf_sdf_stability_compare.csv";
    double duration = 1.5;
    std::vector<double> step_sizes = {1e-3, 2e-3, 4e-3};
    std::size_t kmax = 64;
    double cell_size = 0.05;
    double octree_eps_force = 500.0;
    int octree_max_depth = 8;
    double octree_c_error = 0.5;
    double octree_narrow_band = 0.25;
    double octree_half_extent_margin = 0.15;
    std::uint64_t seed = 42;
};

bool ParseDouble(const std::string& text, double& out) {
    try {
        out = std::stod(text);
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseInt(const std::string& text, int& out) {
    try {
        out = std::stoi(text);
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

void BuildTopInitialState(std::uint64_t seed, chrono::ChVector3d& pos, chrono::ChQuaternion<>& rot) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> xy_jitter(-0.02, 0.02);
    std::uniform_real_distribution<double> z_jitter(-0.015, 0.015);
    std::uniform_real_distribution<double> tilt_jitter(-0.05, 0.05);

    pos = chrono::ChVector3d(0.08 + xy_jitter(rng), 0.55 + 0.5 * xy_jitter(rng), z_jitter(rng));
    rot = chrono::QuatFromAngleZ(0.35 + tilt_jitter(rng)) * chrono::QuatFromAngleX(0.1 + tilt_jitter(rng)) *
          chrono::QuatFromAngleY(tilt_jitter(rng));
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
              << " [--csv <path>] [--duration <seconds>] [--dts <h1,h2,...>] [--kmax <int>] [--cell-size <meters>]"
                 " [--octree-eps-force <double>] [--octree-max-depth <int>] [--octree-c-error <double>]"
                 " [--octree-band <double>] [--octree-margin <double>] [--seed <u64>]\n";
}

bool ParseArgs(int argc, char* argv[], StabilityConfig& cfg) {
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
                std::cerr << "Invalid --dts value. Example: --dts 0.001,0.002,0.004" << std::endl;
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

        if (arg == "--octree-eps-force") {
            const auto* val = need_value("--octree-eps-force");
            if (!val || !ParseDouble(val, cfg.octree_eps_force)) {
                std::cerr << "Invalid --octree-eps-force value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--octree-max-depth") {
            const auto* val = need_value("--octree-max-depth");
            if (!val || !ParseInt(val, cfg.octree_max_depth)) {
                std::cerr << "Invalid --octree-max-depth value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--octree-c-error") {
            const auto* val = need_value("--octree-c-error");
            if (!val || !ParseDouble(val, cfg.octree_c_error)) {
                std::cerr << "Invalid --octree-c-error value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--octree-band") {
            const auto* val = need_value("--octree-band");
            if (!val || !ParseDouble(val, cfg.octree_narrow_band)) {
                std::cerr << "Invalid --octree-band value" << std::endl;
                return false;
            }
            continue;
        }

        if (arg == "--octree-margin") {
            const auto* val = need_value("--octree-margin");
            if (!val || !ParseDouble(val, cfg.octree_half_extent_margin)) {
                std::cerr << "Invalid --octree-margin value" << std::endl;
                return false;
            }
            continue;
        }

        std::cerr << "Unknown option: " << arg << std::endl;
        return false;
    }

    if (cfg.duration <= 0.0 || cfg.kmax == 0 || cfg.cell_size <= 0.0 || cfg.step_sizes.empty() ||
        cfg.octree_eps_force <= 0.0 || cfg.octree_max_depth < 0 || cfg.octree_c_error <= 0.0 ||
        cfg.octree_narrow_band <= 0.0 || cfg.octree_half_extent_margin < 0.0) {
        std::cerr << "Invalid configuration: check duration/kmax/cell-size/octree parameters." << std::endl;
        return false;
    }

    return true;
}

std::shared_ptr<chrono_sdf_contact::SdfVolume> BuildBoxProxyVolume(const chrono::ChVector3d& half_lengths,
                                                                    const StabilityConfig& cfg,
                                                                    bool use_octree) {
    using namespace chrono_sdf_contact;

    auto analytic = std::make_shared<BoxSdfVolume>(half_lengths);
    if (!use_octree) {
        return analytic;
    }

    SdfOctreeVolume::BuildConfig oct_cfg;
    oct_cfg.half_extent = std::max({half_lengths.x(), half_lengths.y(), half_lengths.z()}) + cfg.octree_half_extent_margin;
    oct_cfg.narrow_band = cfg.octree_narrow_band;
    oct_cfg.max_depth = cfg.octree_max_depth;
    oct_cfg.kn = 6.0e4;
    oct_cfg.eps_force = cfg.octree_eps_force;
    oct_cfg.c_error = cfg.octree_c_error;
    return std::make_shared<SdfOctreeVolume>(analytic, oct_cfg);
}

RunResult RunCase(chrono_sdf_contact::SdfSdfSamplingMode sampling_mode,
                  const std::string& mode_name,
                  double step_size,
                  const StabilityConfig& cfg,
                  bool clustering_enabled,
                  bool octree_enabled,
                  std::uint64_t seed) {
    using namespace chrono;
    using namespace chrono_sdf_contact;

    ChSystemSMC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    sys.SetGravitationalAcceleration(ChVector3d(0.0, -9.81, 0.0));

    auto material = chrono_types::make_shared<ChContactMaterialSMC>();
    material->SetFriction(0.4f);
    material->SetRestitution(0.0f);
    material->SetYoungModulus(3e6f);
    material->SetKn(6e4f);
    material->SetGn(250.0f);

    const ChVector3d kBaseHalf(0.6, 0.1, 0.6);
    const ChVector3d kTopHalf(0.25, 0.15, 0.2);

    auto base = chrono_types::make_shared<ChBodyEasyBox>(2.0 * kBaseHalf.x(), 2.0 * kBaseHalf.y(), 2.0 * kBaseHalf.z(),
                                                          1000.0, true, true, material);
    base->SetFixed(true);
    base->SetPos(ChVector3d(0.0, -0.1, 0.0));
    sys.AddBody(base);

    ChVector3d top_init_pos;
    ChQuaternion<> top_init_rot;
    const auto dt_code = static_cast<std::uint64_t>(std::llround(step_size * 1e9));
    BuildTopInitialState(seed ^ (dt_code * 0x9e3779b97f4a7c15ULL), top_init_pos, top_init_rot);

    auto top = chrono_types::make_shared<ChBodyEasyBox>(2.0 * kTopHalf.x(), 2.0 * kTopHalf.y(), 2.0 * kTopHalf.z(),
                                                         1000.0, true, true, material);
    top->SetPos(top_init_pos);
    top->SetRot(top_init_rot);
    sys.AddBody(top);

    auto registry = std::make_shared<SdfRegistry>();
    auto pair_cache = std::make_shared<SdfPairCache>();
    auto narrowphase = std::make_shared<SdfNarrowphase>();
    auto clustering = std::make_shared<SdfClustering>(cfg.kmax, cfg.cell_size);
    narrowphase->SetSdfSdfSamplingMode(sampling_mode);

    auto narrow_cb = std::make_shared<SdfNarrowphaseCallback>(registry, pair_cache);
    auto custom_cb = std::make_shared<SdfCustomCollisionCallback>(registry, pair_cache, narrowphase, clustering);
    custom_cb->SetClusteringEnabled(clustering_enabled);

    SdfProxy base_proxy;
    base_proxy.sdf = BuildBoxProxyVolume(kBaseHalf, cfg, octree_enabled);
    base_proxy.material = material;
    registry->Register(base->GetCollisionModel().get(), base_proxy);

    SdfProxy top_proxy;
    top_proxy.sdf = BuildBoxProxyVolume(kTopHalf, cfg, octree_enabled);
    top_proxy.material = material;
    registry->Register(top->GetCollisionModel().get(), top_proxy);

    sys.GetCollisionSystem()->RegisterNarrowphaseCallback(narrow_cb);
    sys.RegisterCustomCollisionCallback(custom_cb);

    const int num_steps = std::max(1, static_cast<int>(cfg.duration / step_size));
    auto reporter = std::make_shared<MinDistanceReporter>();
    double min_dist = std::numeric_limits<double>::infinity();
    double min_height = top->GetPos().y();
    const auto gravity = sys.GetGravitationalAcceleration();
    const double energy0 = ComputeBodyMechanicalEnergy(top.get(), gravity);
    double max_relative_energy_drift = 0.0;
    unsigned int max_contacts = 0;
    unsigned long long total_contacts = 0;
    bool exploded = false;

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; ++step) {
        pair_cache->BeginStep();
        sys.DoStepDynamics(step_size);

        const auto pos = top->GetPos();
        if (!std::isfinite(pos.x()) || !std::isfinite(pos.y()) || !std::isfinite(pos.z()) || std::abs(pos.y()) > 2.0) {
            exploded = true;
            break;
        }

        min_height = std::min(min_height, pos.y());
        const double energy = ComputeBodyMechanicalEnergy(top.get(), gravity);
        const double rel_drift = std::abs((energy - energy0) / (std::abs(energy0) + 1e-12));
        max_relative_energy_drift = std::max(max_relative_energy_drift, rel_drift);
        const auto ncontacts = sys.GetNumContacts();
        total_contacts += ncontacts;
        max_contacts = std::max(max_contacts, ncontacts);

        if (ncontacts > 0) {
            reporter->Reset();
            sys.GetContactContainer()->ReportAllContacts(reporter);
            if (reporter->HasContact()) {
                min_dist = std::min(min_dist, reporter->GetMinDistance());
            }
        }
    }
    const auto t1 = std::chrono::high_resolution_clock::now();

    RunResult out;
    out.mode = mode_name;
    out.step_size = step_size;
    out.steps = num_steps;
    out.wall_time_s = std::chrono::duration<double>(t1 - t0).count();
    out.avg_contacts = (num_steps > 0) ? static_cast<double>(total_contacts) / static_cast<double>(num_steps) : 0.0;
    out.max_contacts = max_contacts;
    out.max_penetration = (min_dist < 0.0) ? -min_dist : 0.0;
    out.min_height = min_height;
    out.final_height = top->GetPos().y();
    out.max_relative_energy_drift = max_relative_energy_drift;
    out.clustering_enabled = clustering_enabled;
    out.octree_enabled = octree_enabled;
    out.exploded = exploded;
    return out;
}

void WriteCsv(const std::string& path, const std::vector<RunResult>& results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open CSV file: " + path);
    }

    out << "mode,step_size,steps,wall_time_s,avg_contacts,max_contacts,max_penetration,min_height,final_height,max_relative_energy_drift,clustering_enabled,octree_enabled,exploded\n";
    for (const auto& r : results) {
        out << r.mode << "," << r.step_size << "," << r.steps << "," << r.wall_time_s << "," << r.avg_contacts << ","
            << r.max_contacts << "," << r.max_penetration << "," << r.min_height << "," << r.final_height << ","
            << r.max_relative_energy_drift << "," << (r.clustering_enabled ? 1 : 0) << ","
            << (r.octree_enabled ? 1 : 0) << "," << (r.exploded ? 1 : 0) << "\n";
    }
}

void PrintResult(const RunResult& r) {
    std::cout << "[" << r.mode << ", dt=" << r.step_size << "] "
              << "wall_time_s=" << std::fixed << std::setprecision(6) << r.wall_time_s
              << ", avg_contacts=" << r.avg_contacts << ", max_contacts=" << r.max_contacts
              << ", max_penetration=" << r.max_penetration << ", min_height=" << r.min_height
              << ", final_height=" << r.final_height
              << ", max_relative_energy_drift=" << r.max_relative_energy_drift
              << ", clustering=" << (r.clustering_enabled ? 1 : 0)
              << ", octree=" << (r.octree_enabled ? 1 : 0)
              << ", exploded=" << (r.exploded ? 1 : 0) << std::endl;
}

struct ModeSpec {
    std::string name;
    chrono_sdf_contact::SdfSdfSamplingMode sampling_mode = chrono_sdf_contact::SdfSdfSamplingMode::Bidirectional;
    bool clustering_enabled = true;
    bool octree_enabled = false;
};

std::vector<ModeSpec> BuildDefaultModeSpecs() {
    using chrono_sdf_contact::SdfSdfSamplingMode;
    return {
        {"a_to_b_only", SdfSdfSamplingMode::AtoBOnly, true, false},
        {"b_to_a_only", SdfSdfSamplingMode::BtoAOnly, true, false},
        {"bidirectional", SdfSdfSamplingMode::Bidirectional, true, false},
        {"a_to_b_only_no_cluster", SdfSdfSamplingMode::AtoBOnly, false, false},
        {"b_to_a_only_no_cluster", SdfSdfSamplingMode::BtoAOnly, false, false},
        {"bidirectional_no_cluster", SdfSdfSamplingMode::Bidirectional, false, false},
        {"bidirectional_octree", SdfSdfSamplingMode::Bidirectional, true, true},
    };
}

}  // namespace

int main(int argc, char* argv[]) {
    StabilityConfig cfg;
    if (!ParseArgs(argc, argv, cfg)) {
        return 1;
    }

    try {
        std::vector<RunResult> results;
        const auto mode_specs = BuildDefaultModeSpecs();
        results.reserve(cfg.step_sizes.size() * mode_specs.size());

        for (double dt : cfg.step_sizes) {
            for (const auto& mode : mode_specs) {
                results.push_back(
                    RunCase(mode.sampling_mode, mode.name, dt, cfg, mode.clustering_enabled, mode.octree_enabled, cfg.seed));
            }
        }

        for (const auto& r : results) {
            PrintResult(r);
        }
        WriteCsv(cfg.csv_path, results);
        std::cout << "CSV written to: " << cfg.csv_path << std::endl;

        bool any_bidir_contact = false;
        for (const auto& r : results) {
            if (r.mode.find("bidirectional") == 0 && r.max_contacts > 0 && !r.exploded) {
                any_bidir_contact = true;
                break;
            }
        }
        if (!any_bidir_contact) {
            std::cerr << "Validation failed: bidirectional mode did not produce stable contacts." << std::endl;
            return 3;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    return 0;
}
