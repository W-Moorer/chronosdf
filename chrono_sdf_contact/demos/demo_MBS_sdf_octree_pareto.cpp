#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "chrono/core/ChTypes.h"

#include "SdfOctreeVolume.h"
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

struct M3Config {
    std::string csv_path = "sdf_octree_pareto.csv";
    std::string pareto_csv_path = "sdf_octree_pareto_front.csv";
    std::vector<double> eps_force_list = {4000.0, 2000.0, 1000.0, 500.0, 250.0, 125.0};
    int query_count = 40000;
    int max_depth = 8;
    double kn = 6.0e4;
    double c_error = 0.5;
    double narrow_band = 0.25;
    double half_extent = 1.0;
    double sphere_radius = 0.45;
    std::uint64_t seed = 42;
};

struct M3CaseResult {
    double eps_force = 0.0;
    int max_depth = 0;
    double kn = 0.0;
    double c_error = 0.0;
    double narrow_band = 0.0;
    double half_extent = 0.0;
    std::size_t node_count = 0;
    std::size_t leaf_count = 0;
    double min_leaf_size = 0.0;
    double max_leaf_size = 0.0;
    double est_delta_d = 0.0;
    double est_delta_f = 0.0;
    double build_time_s = 0.0;
    double query_time_s = 0.0;
    double rmse_phi = 0.0;
    double max_abs_phi_err = 0.0;
    double rmse_force_err = 0.0;
    double max_force_err = 0.0;
    double rmse_force_err_linear = 0.0;
    double max_force_err_linear = 0.0;
    double force_map_mae = 0.0;
    double force_map_rel_mae = 0.0;
    double force_map_r2_fixed = 0.0;
    double force_map_slope = 0.0;
    double force_map_bound_violation_rate = 0.0;
};

bool ParseDouble(const std::string& s, double& out) {
    try {
        out = std::stod(s);
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseInt(const std::string& s, int& out) {
    try {
        out = std::stoi(s);
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseUInt64(const std::string& s, std::uint64_t& out) {
    try {
        out = static_cast<std::uint64_t>(std::stoull(s));
        return true;
    } catch (...) {
        return false;
    }
}

bool ParseDoubleList(const std::string& s, std::vector<double>& out) {
    out.clear();
    std::size_t start = 0;
    while (start < s.size()) {
        std::size_t end = s.find(',', start);
        if (end == std::string::npos) {
            end = s.size();
        }
        const auto token = s.substr(start, end - start);
        double value = 0.0;
        if (!ParseDouble(token, value) || value <= 0.0) {
            return false;
        }
        out.push_back(value);
        start = end + 1;
    }
    return !out.empty();
}

void PrintUsage(const char* exe) {
    std::cout << "Usage: " << exe
              << " [--csv <path>] [--pareto-csv <path>] [--eps-list <e1,e2,...>] [--queries <int>] [--max-depth <int>]"
                 " [--kn <double>] [--c-error <double>] [--band <double>] [--half-extent <double>] [--radius <double>] [--seed <u64>]\n";
}

bool ParseArgs(int argc, char* argv[], M3Config& cfg) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return false;
        }

        auto need_value = [&](const char* opt) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << opt << std::endl;
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "--csv") {
            const auto* v = need_value("--csv");
            if (!v) {
                return false;
            }
            cfg.csv_path = v;
            continue;
        }
        if (arg == "--pareto-csv") {
            const auto* v = need_value("--pareto-csv");
            if (!v) {
                return false;
            }
            cfg.pareto_csv_path = v;
            continue;
        }
        if (arg == "--eps-list") {
            const auto* v = need_value("--eps-list");
            if (!v || !ParseDoubleList(v, cfg.eps_force_list)) {
                std::cerr << "Invalid --eps-list (example: 4000,2000,1000,500,250)" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--queries") {
            const auto* v = need_value("--queries");
            if (!v || !ParseInt(v, cfg.query_count)) {
                std::cerr << "Invalid --queries" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--max-depth") {
            const auto* v = need_value("--max-depth");
            if (!v || !ParseInt(v, cfg.max_depth)) {
                std::cerr << "Invalid --max-depth" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--kn") {
            const auto* v = need_value("--kn");
            if (!v || !ParseDouble(v, cfg.kn)) {
                std::cerr << "Invalid --kn" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--c-error") {
            const auto* v = need_value("--c-error");
            if (!v || !ParseDouble(v, cfg.c_error)) {
                std::cerr << "Invalid --c-error" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--band") {
            const auto* v = need_value("--band");
            if (!v || !ParseDouble(v, cfg.narrow_band)) {
                std::cerr << "Invalid --band" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--half-extent") {
            const auto* v = need_value("--half-extent");
            if (!v || !ParseDouble(v, cfg.half_extent)) {
                std::cerr << "Invalid --half-extent" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--radius") {
            const auto* v = need_value("--radius");
            if (!v || !ParseDouble(v, cfg.sphere_radius)) {
                std::cerr << "Invalid --radius" << std::endl;
                return false;
            }
            continue;
        }
        if (arg == "--seed") {
            const auto* v = need_value("--seed");
            if (!v || !ParseUInt64(v, cfg.seed)) {
                std::cerr << "Invalid --seed" << std::endl;
                return false;
            }
            continue;
        }

        std::cerr << "Unknown option: " << arg << std::endl;
        return false;
    }

    if (cfg.query_count <= 0 || cfg.max_depth < 0 || cfg.kn <= 0.0 || cfg.c_error <= 0.0 || cfg.narrow_band <= 0.0 ||
        cfg.half_extent <= 0.0 || cfg.sphere_radius <= 0.0 || cfg.eps_force_list.empty()) {
        std::cerr << "Invalid configuration values." << std::endl;
        return false;
    }

    return true;
}

chrono::ChVector3d RandomUnit(std::mt19937_64& rng) {
    std::normal_distribution<double> norm(0.0, 1.0);
    chrono::ChVector3d v(norm(rng), norm(rng), norm(rng));
    double len = v.Length();
    if (len <= 1e-12) {
        return chrono::ChVector3d(1.0, 0.0, 0.0);
    }
    return v / len;
}

M3CaseResult RunCase(const M3Config& cfg, double eps_force) {
    using namespace chrono_sdf_contact;

    auto analytic = std::make_shared<SphereSdfVolume>(cfg.sphere_radius);

    SdfOctreeVolume::BuildConfig oct_cfg;
    oct_cfg.half_extent = cfg.half_extent;
    oct_cfg.narrow_band = cfg.narrow_band;
    oct_cfg.max_depth = cfg.max_depth;
    oct_cfg.kn = cfg.kn;
    oct_cfg.eps_force = eps_force;
    oct_cfg.c_error = cfg.c_error;

    const auto t0 = std::chrono::high_resolution_clock::now();
    auto octree = std::make_shared<SdfOctreeVolume>(analytic, oct_cfg);
    const auto t1 = std::chrono::high_resolution_clock::now();

    std::mt19937_64 rng(cfg.seed + static_cast<std::uint64_t>(eps_force));
    std::uniform_real_distribution<double> offset_dist(-cfg.narrow_band, cfg.narrow_band);

    double sum_sq_phi = 0.0;
    double max_abs_phi_err = 0.0;
    double sum_sq_force = 0.0;
    double max_force_err = 0.0;
    double sum_sq_force_linear = 0.0;
    double max_force_err_linear = 0.0;

    double sum_abs_map_diff = 0.0;
    double sum_rel_map_diff = 0.0;
    double sse_map_fixed = 0.0;
    double sum_force_err = 0.0;
    double sum_force_err_sq = 0.0;
    double sum_map_xx = 0.0;
    double sum_map_xy = 0.0;
    int bound_violation_count = 0;

    const auto q0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < cfg.query_count; ++i) {
        const auto dir = RandomUnit(rng);
        const double offset = offset_dist(rng);
        const auto p = (cfg.sphere_radius + offset) * dir;

        const auto truth = analytic->SampleLocal(p);
        const auto approx = octree->SampleLocal(p);

        const double dphi = approx.phi - truth.phi;
        const double abs_phi = std::abs(dphi);
        sum_sq_phi += dphi * dphi;
        max_abs_phi_err = std::max(max_abs_phi_err, abs_phi);

        const double pen_truth = std::max(0.0, -truth.phi);
        const double pen_approx = std::max(0.0, -approx.phi);
        const double force_truth = cfg.kn * pen_truth;
        const double force_approx = cfg.kn * pen_approx;
        const double force_err_actual = std::abs(force_approx - force_truth);

        const double force_err_linear = cfg.kn * abs_phi;

        sum_sq_force += force_err_actual * force_err_actual;
        max_force_err = std::max(max_force_err, force_err_actual);

        sum_sq_force_linear += force_err_linear * force_err_linear;
        max_force_err_linear = std::max(max_force_err_linear, force_err_linear);

        const double map_diff = force_err_actual - force_err_linear;
        sum_abs_map_diff += std::abs(map_diff);
        sum_rel_map_diff += std::abs(map_diff) / (force_err_linear + 1e-12);
        sse_map_fixed += map_diff * map_diff;

        sum_force_err += force_err_actual;
        sum_force_err_sq += force_err_actual * force_err_actual;

        sum_map_xx += force_err_linear * force_err_linear;
        sum_map_xy += force_err_linear * force_err_actual;

        if (force_err_actual > force_err_linear + 1e-10) {
            bound_violation_count++;
        }
    }
    const auto q1 = std::chrono::high_resolution_clock::now();

    const auto stats = octree->GetStats();

    M3CaseResult out;
    out.eps_force = eps_force;
    out.max_depth = cfg.max_depth;
    out.kn = cfg.kn;
    out.c_error = cfg.c_error;
    out.narrow_band = cfg.narrow_band;
    out.half_extent = cfg.half_extent;
    out.node_count = stats.node_count;
    out.leaf_count = stats.leaf_count;
    out.min_leaf_size = stats.min_leaf_size;
    out.max_leaf_size = stats.max_leaf_size;
    out.est_delta_d = stats.est_max_delta_d_band;
    out.est_delta_f = cfg.kn * stats.est_max_delta_d_band;
    out.build_time_s = std::chrono::duration<double>(t1 - t0).count();
    out.query_time_s = std::chrono::duration<double>(q1 - q0).count();
    out.rmse_phi = std::sqrt(sum_sq_phi / static_cast<double>(cfg.query_count));
    out.max_abs_phi_err = max_abs_phi_err;
    out.rmse_force_err = std::sqrt(sum_sq_force / static_cast<double>(cfg.query_count));
    out.max_force_err = max_force_err;
    out.rmse_force_err_linear = std::sqrt(sum_sq_force_linear / static_cast<double>(cfg.query_count));
    out.max_force_err_linear = max_force_err_linear;
    out.force_map_mae = sum_abs_map_diff / static_cast<double>(cfg.query_count);
    out.force_map_rel_mae = sum_rel_map_diff / static_cast<double>(cfg.query_count);
    const double mean_force_err = sum_force_err / static_cast<double>(cfg.query_count);
    const double sst_force = std::max(0.0, sum_force_err_sq - static_cast<double>(cfg.query_count) * mean_force_err * mean_force_err);
    if (sst_force > 1e-16) {
        out.force_map_r2_fixed = 1.0 - sse_map_fixed / sst_force;
    } else {
        out.force_map_r2_fixed = (sse_map_fixed <= 1e-16) ? 1.0 : 0.0;
    }
    out.force_map_slope = (sum_map_xx > 1e-16) ? (sum_map_xy / sum_map_xx) : 0.0;
    out.force_map_bound_violation_rate = static_cast<double>(bound_violation_count) / static_cast<double>(cfg.query_count);
    return out;
}

bool Dominates(const M3CaseResult& a, const M3CaseResult& b) {
    const bool better_or_equal_time = a.query_time_s <= b.query_time_s;
    const bool better_or_equal_error = a.rmse_phi <= b.rmse_phi;
    const bool strictly_better = (a.query_time_s < b.query_time_s) || (a.rmse_phi < b.rmse_phi);
    return better_or_equal_time && better_or_equal_error && strictly_better;
}

std::vector<M3CaseResult> ExtractParetoFront(const std::vector<M3CaseResult>& all) {
    std::vector<M3CaseResult> front;
    for (std::size_t i = 0; i < all.size(); ++i) {
        bool dominated = false;
        for (std::size_t j = 0; j < all.size(); ++j) {
            if (i == j) {
                continue;
            }
            if (Dominates(all[j], all[i])) {
                dominated = true;
                break;
            }
        }
        if (!dominated) {
            front.push_back(all[i]);
        }
    }
    std::sort(front.begin(), front.end(), [](const M3CaseResult& lhs, const M3CaseResult& rhs) {
        return lhs.query_time_s < rhs.query_time_s;
    });
    return front;
}

void WriteCsv(const std::string& path, const std::vector<M3CaseResult>& rows) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open CSV file: " + path);
    }
    out << "eps_force,kn,c_error,narrow_band,half_extent,max_depth,node_count,leaf_count,min_leaf_size,max_leaf_size,"
           "est_delta_d,est_delta_f,build_time_s,query_time_s,rmse_phi,max_abs_phi_err,rmse_force_err,max_force_err,"
           "rmse_force_err_linear,max_force_err_linear,force_map_mae,force_map_rel_mae,force_map_r2_fixed,force_map_slope,force_map_bound_violation_rate\n";
    for (const auto& r : rows) {
        out << r.eps_force << "," << r.kn << "," << r.c_error << "," << r.narrow_band << "," << r.half_extent << ","
            << r.max_depth << "," << r.node_count << "," << r.leaf_count << "," << r.min_leaf_size << ","
            << r.max_leaf_size << "," << r.est_delta_d << "," << r.est_delta_f << "," << r.build_time_s << ","
            << r.query_time_s << "," << r.rmse_phi << "," << r.max_abs_phi_err << "," << r.rmse_force_err << ","
            << r.max_force_err << "," << r.rmse_force_err_linear << "," << r.max_force_err_linear << ","
            << r.force_map_mae << "," << r.force_map_rel_mae << "," << r.force_map_r2_fixed << "," << r.force_map_slope
            << "," << r.force_map_bound_violation_rate << "\n";
    }
}

void PrintCase(const M3CaseResult& r) {
    std::cout << "eps_F=" << r.eps_force << ", nodes=" << r.node_count << ", leaves=" << r.leaf_count
              << ", est_dF=" << r.est_delta_f << ", query_time_s=" << std::fixed << std::setprecision(6)
              << r.query_time_s << ", rmse_phi=" << r.rmse_phi << ", max_phi_err=" << r.max_abs_phi_err
              << ", map_r2=" << r.force_map_r2_fixed << ", map_slope=" << r.force_map_slope << std::endl;
}

}  // namespace

int main(int argc, char* argv[]) {
    M3Config cfg;
    if (!ParseArgs(argc, argv, cfg)) {
        return 1;
    }

    try {
        std::vector<M3CaseResult> all_results;
        all_results.reserve(cfg.eps_force_list.size());

        for (double eps_f : cfg.eps_force_list) {
            auto result = RunCase(cfg, eps_f);
            all_results.push_back(result);
            PrintCase(result);
        }

        std::sort(all_results.begin(), all_results.end(),
                  [](const M3CaseResult& lhs, const M3CaseResult& rhs) { return lhs.eps_force > rhs.eps_force; });

        const auto pareto = ExtractParetoFront(all_results);
        WriteCsv(cfg.csv_path, all_results);
        WriteCsv(cfg.pareto_csv_path, pareto);

        std::cout << "CSV written to: " << cfg.csv_path << std::endl;
        std::cout << "Pareto CSV written to: " << cfg.pareto_csv_path << std::endl;

        if (pareto.empty()) {
            std::cerr << "Validation failed: empty Pareto front." << std::endl;
            return 3;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    return 0;
}
