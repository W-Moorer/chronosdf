#pragma once

#include <cstdint>
#include <memory>
#include <unordered_map>

#include "chrono/collision/ChCollisionSystem.h"

#include "SdfPairCache.h"
#include "SdfRegistry.h"

namespace chrono_sdf_contact {

class SdfNarrowphaseCallback : public chrono::ChCollisionSystem::NarrowphaseCallback {
  public:
    SdfNarrowphaseCallback(std::shared_ptr<SdfRegistry> registry, std::shared_ptr<SdfPairCache> pair_cache);

    bool OnNarrowphase(chrono::ChCollisionInfo& cinfo) override;

  private:
    void RefreshCachesIfRegistryChanged();
    bool IsSdfModelCached(chrono::ChCollisionModel* model);
    bool HasSupportedShapeCached(chrono::ChCollisionModel* model);

    std::shared_ptr<SdfRegistry> registry_;
    std::shared_ptr<SdfPairCache> pair_cache_;
    std::unordered_map<chrono::ChCollisionModel*, bool> sdf_presence_cache_;
    std::unordered_map<chrono::ChCollisionModel*, bool> supported_shape_cache_;
    std::uint64_t cached_registry_revision_ = 0;
};

}  // namespace chrono_sdf_contact
