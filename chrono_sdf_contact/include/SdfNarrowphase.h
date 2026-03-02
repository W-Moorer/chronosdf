#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

#include "chrono/collision/ChCollisionInfo.h"
#include "chrono/collision/ChCollisionModel.h"
#include "chrono/physics/ChContactMaterial.h"

namespace chrono_sdf_contact {

class SdfRegistry;
struct SdfProxy;

enum class SdfSdfSamplingMode {
    Bidirectional,
    AtoBOnly,
    BtoAOnly
};

struct GeneratedContact {
    chrono::ChCollisionInfo cinfo;
    std::shared_ptr<chrono::ChContactMaterial> matA;
    std::shared_ptr<chrono::ChContactMaterial> matB;
};

class SdfNarrowphase {
  public:
    void SetSdfSdfSamplingMode(SdfSdfSamplingMode mode);
    SdfSdfSamplingMode GetSdfSdfSamplingMode() const;

    void BuildContactsForPair(const SdfRegistry& registry,
                              chrono::ChCollisionModel* modelA,
                              chrono::ChCollisionModel* modelB,
                              std::vector<GeneratedContact>& out_contacts) const;

  private:
    struct OtherModelSupportQuery {
        chrono::ChVector3d p_model = chrono::VNULL;
        double radius = 0.0;
        bool sphere = false;
        std::shared_ptr<chrono::ChContactMaterial> material;
    };

    struct OtherModelSupportCache {
        std::size_t shape_count = 0;
        double bound_radius = 0.0;
        std::vector<OtherModelSupportQuery> queries;
    };

    const std::vector<chrono::ChVector3d>& GetBoundarySamplesCached(chrono::ChCollisionModel* model,
                                                                     const SdfProxy& proxy) const;
    double GetModelBoundRadiusCached(chrono::ChCollisionModel* model) const;
    const OtherModelSupportCache& GetOtherModelSupportCached(chrono::ChCollisionModel* model) const;

    SdfSdfSamplingMode sdf_sdf_sampling_mode_ = SdfSdfSamplingMode::Bidirectional;
    mutable std::unordered_map<chrono::ChCollisionModel*, std::vector<chrono::ChVector3d>> boundary_samples_cache_;
    mutable std::unordered_map<chrono::ChCollisionModel*, double> model_bound_radius_cache_;
    mutable std::unordered_map<chrono::ChCollisionModel*, OtherModelSupportCache> other_support_cache_;
};

}  // namespace chrono_sdf_contact
