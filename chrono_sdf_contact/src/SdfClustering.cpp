#include "SdfClustering.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>

namespace chrono_sdf_contact {

namespace {

struct Key3 {
    std::int64_t ix;
    std::int64_t iy;
    std::int64_t iz;

    bool operator==(const Key3& other) const noexcept {
        return ix == other.ix && iy == other.iy && iz == other.iz;
    }
};

struct Key3Hash {
    std::size_t operator()(const Key3& key) const noexcept {
        const auto h1 = std::hash<std::int64_t>{}(key.ix);
        const auto h2 = std::hash<std::int64_t>{}(key.iy);
        const auto h3 = std::hash<std::int64_t>{}(key.iz);
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        seed ^= h3 + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};

struct ClusterAgg {
    chrono::ChVector3d sum_vpA = chrono::VNULL;
    chrono::ChVector3d sum_vpB = chrono::VNULL;
    chrono::ChVector3d sum_vN = chrono::VNULL;
    double min_distance = std::numeric_limits<double>::infinity();
    std::size_t count = 0;

    chrono::ChCollisionModel* modelA = nullptr;
    chrono::ChCollisionModel* modelB = nullptr;
    std::shared_ptr<chrono::ChContactMaterial> matA;
    std::shared_ptr<chrono::ChContactMaterial> matB;
};

Key3 Quantize(const chrono::ChVector3d& p, double h) {
    return {static_cast<std::int64_t>(std::floor(p.x() / h)), static_cast<std::int64_t>(std::floor(p.y() / h)),
            static_cast<std::int64_t>(std::floor(p.z() / h))};
}

void KeepDeepest(std::vector<GeneratedContact>& contacts, std::size_t max_contacts) {
    if (contacts.size() <= max_contacts) {
        return;
    }

    std::sort(contacts.begin(), contacts.end(), [](const GeneratedContact& lhs, const GeneratedContact& rhs) {
        return lhs.cinfo.distance < rhs.cinfo.distance;
    });

    contacts.erase(contacts.begin() + static_cast<std::ptrdiff_t>(max_contacts), contacts.end());
}

}  // namespace

SdfClustering::SdfClustering(std::size_t max_contacts, double cell_size)
    : max_contacts_(max_contacts), cell_size_(cell_size) {}

void SdfClustering::SetMaxContacts(std::size_t max_contacts) {
    max_contacts_ = max_contacts;
}

void SdfClustering::SetCellSize(double cell_size) {
    cell_size_ = cell_size;
}

void SdfClustering::Apply(std::vector<GeneratedContact>& contacts) const {
    if (contacts.empty()) {
        return;
    }

    if (cell_size_ <= 0.0) {
        KeepDeepest(contacts, max_contacts_);
        return;
    }

    std::unordered_map<Key3, ClusterAgg, Key3Hash> clusters;
    clusters.reserve(contacts.size());

    for (const auto& contact : contacts) {
        const auto p = 0.5 * (contact.cinfo.vpA + contact.cinfo.vpB);
        const auto key = Quantize(p, cell_size_);
        auto& agg = clusters[key];

        if (agg.count == 0) {
            agg.modelA = contact.cinfo.modelA;
            agg.modelB = contact.cinfo.modelB;
            agg.matA = contact.matA;
            agg.matB = contact.matB;
            agg.min_distance = contact.cinfo.distance;
        }

        auto n = contact.cinfo.vN;
        if (agg.sum_vN.Length2() > 1e-16 && n.Dot(agg.sum_vN) < 0.0) {
            n = -n;
        }

        agg.sum_vpA += contact.cinfo.vpA;
        agg.sum_vpB += contact.cinfo.vpB;
        agg.sum_vN += n;
        agg.count++;

        if (contact.cinfo.distance < agg.min_distance) {
            agg.min_distance = contact.cinfo.distance;
            agg.matA = contact.matA;
            agg.matB = contact.matB;
        }
    }

    std::vector<GeneratedContact> reduced;
    reduced.reserve(clusters.size());

    for (const auto& kv : clusters) {
        const auto& agg = kv.second;
        if (agg.count == 0 || !agg.modelA || !agg.modelB || !agg.matA || !agg.matB) {
            continue;
        }

        chrono::ChCollisionInfo cinfo;
        cinfo.modelA = agg.modelA;
        cinfo.modelB = agg.modelB;
        cinfo.shapeA = nullptr;
        cinfo.shapeB = nullptr;
        cinfo.vpA = agg.sum_vpA / static_cast<double>(agg.count);
        cinfo.vpB = agg.sum_vpB / static_cast<double>(agg.count);
        cinfo.vN = agg.sum_vN;
        if (cinfo.vN.Length2() > 1e-16) {
            cinfo.vN.Normalize();
        } else {
            cinfo.vN = chrono::ChVector3d(0.0, 1.0, 0.0);
        }
        cinfo.distance = agg.min_distance;
        cinfo.eff_radius = chrono::ChCollisionInfo::GetDefaultEffectiveCurvatureRadius();
        cinfo.reaction_cache = nullptr;

        GeneratedContact out;
        out.cinfo = cinfo;
        out.matA = agg.matA;
        out.matB = agg.matB;
        reduced.push_back(std::move(out));
    }

    contacts.swap(reduced);
    KeepDeepest(contacts, max_contacts_);
}

}  // namespace chrono_sdf_contact
