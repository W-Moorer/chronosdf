#pragma once

#include <memory>

#include "chrono/physics/ChSystem.h"

#include "SdfClustering.h"
#include "SdfNarrowphase.h"
#include "SdfPairCache.h"
#include "SdfRegistry.h"

namespace chrono_sdf_contact {

class SdfCustomCollisionCallback : public chrono::ChSystem::CustomCollisionCallback {
  public:
    SdfCustomCollisionCallback(std::shared_ptr<SdfRegistry> registry,
                               std::shared_ptr<SdfPairCache> pair_cache,
                               std::shared_ptr<SdfNarrowphase> narrowphase = nullptr,
                               std::shared_ptr<SdfClustering> clustering = nullptr);

    void SetClusteringEnabled(bool enabled) { clustering_enabled_ = enabled; }
    bool IsClusteringEnabled() const { return clustering_enabled_; }

    void OnCustomCollision(chrono::ChSystem* sys) override;

  private:
    std::shared_ptr<SdfRegistry> registry_;
    std::shared_ptr<SdfPairCache> pair_cache_;
    std::shared_ptr<SdfNarrowphase> narrowphase_;
    std::shared_ptr<SdfClustering> clustering_;
    bool clustering_enabled_ = true;
};

}  // namespace chrono_sdf_contact
