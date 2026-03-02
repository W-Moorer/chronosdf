#include <algorithm>
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

struct ScaleupResult {
    std::string mode;
    int body_count = 0;
    int steps = 0;
    double step_size = 0.0;
    double wall_time_s = 0.0;
    double avg_contacts = 0.0;
    unsigned int max_contacts = 0;
    double max_penetration = 0.0;
    double max_relative_energy_drift = 0.0;
    bool exploded = false;
};

struct ScaleupConfig {
    std::string csv_path = "sdf_p2_scaleup.csv";
    double duration = 0.8;
    double step_size = 1.5e-3;
    std::vector<int> body_counts = {8, 16, 32, 64};
    std::size_t kmax = 64;
    double cell_size = 0.05;
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

bool ParseIntList(const std::string& text, std::vector<int>& out_counts) {
    out_counts.clear();
    std::size_t start = 0;
    while (start < text.size()) {
        std::size_t end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        auto token = text.substr(start, end - start);
        try {
            const int value = std::stoi(token);
            if (value <= 0) {
                return false;
            }
            out_counts.push_back(value);
        } catch (...) {
            return false;
        }
        start = end + 1;
    }
    return !out_counts.empty();
}

void PrintUsage(const char* exe_name) {
    std::cout << "Usage: " << exe_name
              << " [--csv <path>] [--duration <seconds>] [--dt <seconds>] [--counts <c1,c2,...>] [--kmax <int>]"
                 " [--cell-size <meters>] [--seed <u64>]\n";
}

bool ParseArgs(int argc, char* argv[], ScaleupConfig& cfg) {
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
        if (arg == "--dt") {
            const auto* val = need_value("--dt");
            if (!val || !ParseDouble(val, cfg.step_size)) {
                std::cerr << "Invalid --dt value" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--counts") {
            const auto* val = need_value("--counts");
            if (!val || !ParseIntList(val, cfg.body_counts)) {
                std::cerr << "Invalid --counts value. Example: --counts 8,16,32,64" << std::endl;
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

    if (cfg.duration <= 0.0 || cfg.step_size <= 0.0 || cfg.body_counts.empty() || cfg.kmax == 0 || cfg.cell_size <= 0.0) {
        std::cerr << "Invalid configuration: duration/dt/counts/kmax/cell-size must be positive." << std::endl;
        return false;
    }
    return true;
}

ScaleupResult RunCase(int body_count, bool use_sdf_contact, const ScaleupConfig& cfg) {
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

    auto ground = chrono_types::make_shared<ChBodyEasyBox>(8.0, 0.2, 8.0, 1000.0, true, true, material);
    ground->SetFixed(true);
    ground->SetPos(ChVector3d(0.0, -0.1, 0.0));
    sys.AddBody(ground);

    const ChVector3d half(0.12, 0.12, 0.12);
    const int per_row = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(body_count)))));
    const int layer_size = per_row * per_row;
    const double spacing = 0.30;
    const double layer_height = 2.0 * half.y() + 0.04;

    std::mt19937_64 rng(cfg.seed ^ static_cast<std::uint64_t>(body_count) * 0x9e3779b97f4a7c15ULL);
    std::uniform_real_distribution<double> jitter(-0.012, 0.012);

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.reserve(static_cast<std::size_t>(body_count));

    std::shared_ptr<SdfRegistry> registry;
    std::shared_ptr<SdfPairCache> pair_cache;
    if (use_sdf_contact) {
        registry = std::make_shared<SdfRegistry>();
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

    for (int i = 0; i < body_count; ++i) {
        auto box = chrono_types::make_shared<ChBodyEasyBox>(2.0 * half.x(), 2.0 * half.y(), 2.0 * half.z(), 1000.0, true,
                                                            true, material);

        const int layer = i / layer_size;
        const int id = i % layer_size;
        const int ix = id % per_row;
        const int iz = id / per_row;
        const double x = (static_cast<double>(ix) - 0.5 * static_cast<double>(per_row - 1)) * spacing + jitter(rng);
        const double z = (static_cast<double>(iz) - 0.5 * static_cast<double>(per_row - 1)) * spacing + jitter(rng);
        const double y = half.y() + 0.01 + static_cast<double>(layer) * layer_height;

        box->SetPos(ChVector3d(x, y, z));
        box->SetRot(QuatFromAngleZ(0.12 * jitter(rng)));
        if ((i % 7) == 0) {
            box->SetPosDt(ChVector3d(0.45, 0.0, 0.0));
        }
        sys.AddBody(box);
        bodies.push_back(box);

        if (use_sdf_contact) {
            SdfProxy proxy;
            proxy.sdf = std::make_shared<BoxSdfVolume>(half);
            proxy.material = material;
            registry->Register(box->GetCollisionModel().get(), proxy);
        }
    }

    const int num_steps = std::max(1, static_cast<int>(cfg.duration / cfg.step_size));
    const auto gravity = sys.GetGravitationalAcceleration();
    const double energy0 = ComputeBodiesMechanicalEnergy(bodies, gravity);
    double max_relative_energy_drift = 0.0;
    bool exploded = false;

    auto reporter = std::make_shared<MinDistanceReporter>();
    double min_dist = std::numeric_limits<double>::infinity();
    unsigned int max_contacts = 0;
    unsigned long long total_contacts = 0;

    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int step = 0; step < num_steps; ++step) {
        if (pair_cache) {
            pair_cache->BeginStep();
        }
        sys.DoStepDynamics(cfg.step_size);

        for (const auto& body : bodies) {
            const auto p = body->GetPos();
            if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()) || std::abs(p.y()) > 8.0) {
                exploded = true;
                break;
            }
        }
        if (exploded) {
            break;
        }

        const auto ncontacts = sys.GetNumContacts();
        total_contacts += ncontacts;
        max_contacts = std::max(max_contacts, ncontacts);

        const double energy = ComputeBodiesMechanicalEnergy(bodies, gravity);
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

    ScaleupResult out;
    out.mode = use_sdf_contact ? "sdf_contact" : "chrono_baseline";
    out.body_count = body_count;
    out.steps = num_steps;
    out.step_size = cfg.step_size;
    out.wall_time_s = std::chrono::duration<double>(t1 - t0).count();
    out.avg_contacts = (num_steps > 0) ? static_cast<double>(total_contacts) / static_cast<double>(num_steps) : 0.0;
    out.max_contacts = max_contacts;
    out.max_penetration = (min_dist < 0.0) ? -min_dist : 0.0;
    out.max_relative_energy_drift = max_relative_energy_drift;
    out.exploded = exploded;
    return out;
}

void PrintResult(const ScaleupResult& r) {
    std::cout << "[" << r.mode << ", body_count=" << r.body_count << "] "
              << "wall_time_s=" << std::fixed << std::setprecision(6) << r.wall_time_s
              << ", avg_contacts=" << r.avg_contacts << ", max_contacts=" << r.max_contacts
              << ", max_penetration=" << r.max_penetration
              << ", max_relative_energy_drift=" << r.max_relative_energy_drift
              << ", exploded=" << (r.exploded ? 1 : 0) << std::endl;
}

void WriteCsv(const std::string& path, const std::vector<ScaleupResult>& rows) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open CSV file: " + path);
    }
    out << "mode,body_count,steps,step_size,wall_time_s,avg_contacts,max_contacts,max_penetration,max_relative_energy_drift,exploded\n";
    for (const auto& r : rows) {
        out << r.mode << "," << r.body_count << "," << r.steps << "," << r.step_size << "," << r.wall_time_s << ","
            << r.avg_contacts << "," << r.max_contacts << "," << r.max_penetration << "," << r.max_relative_energy_drift
            << "," << (r.exploded ? 1 : 0) << "\n";
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    ScaleupConfig cfg;
    if (!ParseArgs(argc, argv, cfg)) {
        return 1;
    }

    try {
        std::vector<ScaleupResult> results;
        results.reserve(cfg.body_counts.size() * 2);
        for (int body_count : cfg.body_counts) {
            results.push_back(RunCase(body_count, false, cfg));
            results.push_back(RunCase(body_count, true, cfg));
        }

        for (const auto& r : results) {
            PrintResult(r);
        }
        WriteCsv(cfg.csv_path, results);
        std::cout << "CSV written to: " << cfg.csv_path << std::endl;

        for (int count : cfg.body_counts) {
            bool ok = false;
            for (const auto& r : results) {
                if (r.mode == "sdf_contact" && r.body_count == count && !r.exploded && r.max_contacts > 0) {
                    ok = true;
                    break;
                }
            }
            if (!ok) {
                std::cerr << "Validation failed: no stable SDF contacts for body_count=" << count << std::endl;
                return 3;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    return 0;
}
