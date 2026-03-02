#pragma once

#include <array>
#include <memory>
#include <vector>

#include "SdfVolume.h"

namespace chrono_sdf_contact {

class SdfOctreeVolume : public SdfVolume {
  public:
    struct BuildConfig {
        double half_extent = 1.0;
        double narrow_band = 0.25;
        int max_depth = 8;
        double kn = 6.0e4;
        double eps_force = 500.0;
        double c_error = 0.5;
    };

    struct Stats {
        std::size_t node_count = 0;
        std::size_t leaf_count = 0;
        double min_leaf_size = 0.0;
        double max_leaf_size = 0.0;
        double est_max_delta_d_all = 0.0;
        double est_max_delta_d_band = 0.0;
    };

    SdfOctreeVolume(std::shared_ptr<SdfVolume> source, const BuildConfig& cfg);

    Sample SampleLocal(const chrono::ChVector3d& p_local) const override;

    const BuildConfig& GetConfig() const { return cfg_; }
    const Stats& GetStats() const { return stats_; }

  private:
    struct Node {
        chrono::ChVector3d center = chrono::VNULL;
        double half_size = 0.0;
        std::array<double, 8> phi_corner{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
        std::array<int, 8> children{{-1, -1, -1, -1, -1, -1, -1, -1}};
        bool is_leaf = true;
    };

    static int CornerIndex(int ix, int iy, int iz);
    static chrono::ChVector3d CornerOffset(int idx);
    static bool InsideNodeBounds(const Node& node, const chrono::ChVector3d& p);

    int BuildNode(const chrono::ChVector3d& center, double half_size, int depth);
    bool ShouldRefine(const Node& node, int depth) const;
    bool IntersectsNarrowBand(const std::array<double, 8>& phi_corner, double half_size) const;
    Sample SampleLeaf(const Node& leaf, const chrono::ChVector3d& p_local) const;
    int LocateLeafIndex(const chrono::ChVector3d& p_local) const;
    double EvaluatePhiAt(const chrono::ChVector3d& p_local) const;

    std::shared_ptr<SdfVolume> source_;
    BuildConfig cfg_;
    std::vector<Node> nodes_;
    Stats stats_;
};

}  // namespace chrono_sdf_contact
