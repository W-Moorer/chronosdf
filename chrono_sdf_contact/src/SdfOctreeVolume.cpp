#include "SdfOctreeVolume.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace chrono_sdf_contact {

namespace {

double Clamp01(double x) {
    if (x < 0.0) {
        return 0.0;
    }
    if (x > 1.0) {
        return 1.0;
    }
    return x;
}

double Lerp(double a, double b, double t) {
    return a + t * (b - a);
}

}  // namespace

SdfOctreeVolume::SdfOctreeVolume(std::shared_ptr<SdfVolume> source, const BuildConfig& cfg)
    : source_(std::move(source)), cfg_(cfg) {
    if (!source_) {
        source_ = std::make_shared<SdfVolume>();
    }

    if (cfg_.half_extent <= 0.0) {
        cfg_.half_extent = 1.0;
    }
    if (cfg_.narrow_band <= 0.0) {
        cfg_.narrow_band = 0.25;
    }
    if (cfg_.max_depth < 0) {
        cfg_.max_depth = 0;
    }
    if (cfg_.kn <= 0.0) {
        cfg_.kn = 6.0e4;
    }
    if (cfg_.eps_force <= 0.0) {
        cfg_.eps_force = 500.0;
    }
    if (cfg_.c_error <= 0.0) {
        cfg_.c_error = 0.5;
    }

    nodes_.reserve(1024);
    stats_.min_leaf_size = std::numeric_limits<double>::infinity();
    stats_.max_leaf_size = 0.0;
    stats_.est_max_delta_d_all = 0.0;
    stats_.est_max_delta_d_band = 0.0;

    BuildNode(chrono::VNULL, cfg_.half_extent, 0);

    if (stats_.leaf_count == 0) {
        stats_.min_leaf_size = 0.0;
    }
    if (stats_.est_max_delta_d_band == 0.0) {
        stats_.est_max_delta_d_band = stats_.est_max_delta_d_all;
    }
    stats_.node_count = nodes_.size();
}

SdfVolume::Sample SdfOctreeVolume::SampleLocal(const chrono::ChVector3d& p_local) const {
    if (!source_ || nodes_.empty()) {
        return source_ ? source_->SampleLocal(p_local) : SdfVolume::Sample{1.0, chrono::ChVector3d(1.0, 0.0, 0.0)};
    }

    const int leaf_idx = LocateLeafIndex(p_local);
    if (leaf_idx < 0) {
        return source_->SampleLocal(p_local);
    }
    return SampleLeaf(nodes_[static_cast<std::size_t>(leaf_idx)], p_local);
}

int SdfOctreeVolume::CornerIndex(int ix, int iy, int iz) {
    return ix + (iy << 1) + (iz << 2);
}

chrono::ChVector3d SdfOctreeVolume::CornerOffset(int idx) {
    return chrono::ChVector3d((idx & 1) ? 1.0 : -1.0, (idx & 2) ? 1.0 : -1.0, (idx & 4) ? 1.0 : -1.0);
}

bool SdfOctreeVolume::InsideNodeBounds(const Node& node, const chrono::ChVector3d& p) {
    return std::abs(p.x() - node.center.x()) <= node.half_size && std::abs(p.y() - node.center.y()) <= node.half_size &&
           std::abs(p.z() - node.center.z()) <= node.half_size;
}

int SdfOctreeVolume::BuildNode(const chrono::ChVector3d& center, double half_size, int depth) {
    Node node;
    node.center = center;
    node.half_size = half_size;
    node.is_leaf = true;

    for (int i = 0; i < 8; ++i) {
        const auto corner = center + half_size * CornerOffset(i);
        node.phi_corner[static_cast<std::size_t>(i)] = source_->SampleLocal(corner).phi;
    }

    const int this_idx = static_cast<int>(nodes_.size());
    nodes_.push_back(node);

    if (ShouldRefine(nodes_[static_cast<std::size_t>(this_idx)], depth)) {
        nodes_[static_cast<std::size_t>(this_idx)].is_leaf = false;

        const double child_half = half_size * 0.5;
        for (int i = 0; i < 8; ++i) {
            const auto child_center = center + child_half * CornerOffset(i);
            const int child_idx = BuildNode(child_center, child_half, depth + 1);
            nodes_[static_cast<std::size_t>(this_idx)].children[static_cast<std::size_t>(i)] = child_idx;
        }
    } else {
        const double leaf_size = 2.0 * half_size;
        stats_.leaf_count++;
        stats_.min_leaf_size = std::min(stats_.min_leaf_size, leaf_size);
        stats_.max_leaf_size = std::max(stats_.max_leaf_size, leaf_size);
        const double delta_d_est = cfg_.c_error * leaf_size;
        stats_.est_max_delta_d_all = std::max(stats_.est_max_delta_d_all, delta_d_est);
        if (IntersectsNarrowBand(nodes_[static_cast<std::size_t>(this_idx)].phi_corner, half_size)) {
            stats_.est_max_delta_d_band = std::max(stats_.est_max_delta_d_band, delta_d_est);
        }
    }

    return this_idx;
}

bool SdfOctreeVolume::ShouldRefine(const Node& node, int depth) const {
    if (depth >= cfg_.max_depth) {
        return false;
    }

    if (!IntersectsNarrowBand(node.phi_corner, node.half_size)) {
        return false;
    }

    const double h = 2.0 * node.half_size;
    const double delta_d_est = cfg_.c_error * h;
    const double delta_f_est = cfg_.kn * delta_d_est;
    return delta_f_est > cfg_.eps_force;
}

bool SdfOctreeVolume::IntersectsNarrowBand(const std::array<double, 8>& phi_corner, double half_size) const {
    double min_phi = phi_corner[0];
    double max_phi = phi_corner[0];
    double min_abs_phi = std::abs(phi_corner[0]);
    for (double phi : phi_corner) {
        min_phi = std::min(min_phi, phi);
        max_phi = std::max(max_phi, phi);
        min_abs_phi = std::min(min_abs_phi, std::abs(phi));
    }

    if (min_phi <= cfg_.narrow_band && max_phi >= -cfg_.narrow_band) {
        return true;
    }
    if (min_phi <= 0.0 && max_phi >= 0.0) {
        return true;
    }
    // Conservative Lipschitz-based test (|grad phi| <= 1 for signed distance fields):
    // if the best corner value can still reach the narrow band inside this cell, refine.
    const double max_variation = 2.0 * half_size * std::sqrt(3.0);
    if (min_abs_phi - max_variation <= cfg_.narrow_band) {
        return true;
    }

    return false;
}

SdfVolume::Sample SdfOctreeVolume::SampleLeaf(const Node& leaf, const chrono::ChVector3d& p_local) const {
    const double h = 2.0 * leaf.half_size;
    if (h <= 1e-12) {
        return source_->SampleLocal(p_local);
    }

    const auto min_corner = leaf.center - chrono::ChVector3d(leaf.half_size, leaf.half_size, leaf.half_size);
    const double u = Clamp01((p_local.x() - min_corner.x()) / h);
    const double v = Clamp01((p_local.y() - min_corner.y()) / h);
    const double w = Clamp01((p_local.z() - min_corner.z()) / h);

    const auto& c = leaf.phi_corner;
    const double c000 = c[CornerIndex(0, 0, 0)];
    const double c100 = c[CornerIndex(1, 0, 0)];
    const double c010 = c[CornerIndex(0, 1, 0)];
    const double c110 = c[CornerIndex(1, 1, 0)];
    const double c001 = c[CornerIndex(0, 0, 1)];
    const double c101 = c[CornerIndex(1, 0, 1)];
    const double c011 = c[CornerIndex(0, 1, 1)];
    const double c111 = c[CornerIndex(1, 1, 1)];

    const double c00 = Lerp(c000, c100, u);
    const double c10 = Lerp(c010, c110, u);
    const double c01 = Lerp(c001, c101, u);
    const double c11 = Lerp(c011, c111, u);

    const double c0 = Lerp(c00, c10, v);
    const double c1 = Lerp(c01, c11, v);
    const double phi = Lerp(c0, c1, w);

    const double dphi_du = (1.0 - v) * (1.0 - w) * (c100 - c000) + v * (1.0 - w) * (c110 - c010) +
                           (1.0 - v) * w * (c101 - c001) + v * w * (c111 - c011);
    const double dphi_dv = (1.0 - u) * (1.0 - w) * (c010 - c000) + u * (1.0 - w) * (c110 - c100) +
                           (1.0 - u) * w * (c011 - c001) + u * w * (c111 - c101);
    const double dphi_dw = (1.0 - u) * (1.0 - v) * (c001 - c000) + u * (1.0 - v) * (c101 - c100) +
                           (1.0 - u) * v * (c011 - c010) + u * v * (c111 - c110);

    chrono::ChVector3d grad(dphi_du / h, dphi_dv / h, dphi_dw / h);
    if (grad.Length2() > 1e-16) {
        grad.Normalize();
    } else {
        grad = source_->SampleLocal(p_local).grad;
        if (grad.Length2() > 1e-16) {
            grad.Normalize();
        }
    }

    return {phi, grad};
}

int SdfOctreeVolume::LocateLeafIndex(const chrono::ChVector3d& p_local) const {
    if (nodes_.empty()) {
        return -1;
    }

    int idx = 0;
    if (!InsideNodeBounds(nodes_[0], p_local)) {
        return -1;
    }

    while (true) {
        const auto& node = nodes_[static_cast<std::size_t>(idx)];
        if (node.is_leaf) {
            return idx;
        }

        const int ix = (p_local.x() >= node.center.x()) ? 1 : 0;
        const int iy = (p_local.y() >= node.center.y()) ? 1 : 0;
        const int iz = (p_local.z() >= node.center.z()) ? 1 : 0;
        const int child_slot = CornerIndex(ix, iy, iz);
        const int child_idx = node.children[static_cast<std::size_t>(child_slot)];
        if (child_idx < 0) {
            return idx;
        }
        idx = child_idx;
    }
}

double SdfOctreeVolume::EvaluatePhiAt(const chrono::ChVector3d& p_local) const {
    const int leaf_idx = LocateLeafIndex(p_local);
    if (leaf_idx < 0) {
        return source_->SampleLocal(p_local).phi;
    }
    return SampleLeaf(nodes_[static_cast<std::size_t>(leaf_idx)], p_local).phi;
}

}  // namespace chrono_sdf_contact
