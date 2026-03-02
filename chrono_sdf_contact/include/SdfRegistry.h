#pragma once

#include <atomic>
#include <cstdint>
#include <memory>
#include <shared_mutex>
#include <unordered_map>
#include <vector>

#include "chrono/collision/ChCollisionModel.h"
#include "chrono/core/ChVector3.h"
#include "chrono/physics/ChContactMaterial.h"

#include "SdfVolume.h"

namespace chrono_sdf_contact {

struct SdfProxy {
    std::shared_ptr<SdfVolume> sdf;
    std::shared_ptr<chrono::ChContactMaterial> material;
    double envelope = 0.0;
    std::vector<chrono::ChVector3d> boundary_samples_local;
};

class SdfRegistry {
  public:
    void Register(chrono::ChCollisionModel* model, const SdfProxy& proxy);
    void QueryPresence(chrono::ChCollisionModel* modelA, chrono::ChCollisionModel* modelB, bool& hasA, bool& hasB) const;
    void GetPairPtrs(chrono::ChCollisionModel* modelA,
                     chrono::ChCollisionModel* modelB,
                     const SdfProxy*& out_proxyA,
                     const SdfProxy*& out_proxyB) const;
    void GetPair(chrono::ChCollisionModel* modelA,
                 chrono::ChCollisionModel* modelB,
                 bool& hasA,
                 SdfProxy& out_proxyA,
                 bool& hasB,
                 SdfProxy& out_proxyB) const;
    bool Has(chrono::ChCollisionModel* model) const;
    bool Get(chrono::ChCollisionModel* model, SdfProxy& out_proxy) const;
    std::uint64_t GetRevision() const noexcept;
    void Clear();

  private:
    mutable std::shared_mutex mutex_;
    std::unordered_map<chrono::ChCollisionModel*, SdfProxy> entries_;
    std::atomic<std::uint64_t> revision_{0};
};

}  // namespace chrono_sdf_contact
