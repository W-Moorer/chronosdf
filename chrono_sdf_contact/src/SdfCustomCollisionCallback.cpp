#include "SdfCustomCollisionCallback.h"

#include <utility>
#include <vector>

#include "chrono/physics/ChContactContainer.h"

namespace chrono_sdf_contact {

SdfCustomCollisionCallback::SdfCustomCollisionCallback(std::shared_ptr<SdfRegistry> registry,
                                                       std::shared_ptr<SdfPairCache> pair_cache,
                                                       std::shared_ptr<SdfNarrowphase> narrowphase,
                                                       std::shared_ptr<SdfClustering> clustering)
    : registry_(std::move(registry)),
      pair_cache_(std::move(pair_cache)),
      narrowphase_(std::move(narrowphase)),
      clustering_(std::move(clustering)) {
    if (!narrowphase_) {
        narrowphase_ = std::make_shared<SdfNarrowphase>();
    }
    if (!clustering_) {
        clustering_ = std::make_shared<SdfClustering>();
    }
}

void SdfCustomCollisionCallback::OnCustomCollision(chrono::ChSystem* sys) {
    if (!sys || !registry_ || !pair_cache_ || !narrowphase_) {
        return;
    }

    auto contact_container = sys->GetContactContainer();
    if (!contact_container) {
        return;
    }

    auto pairs = pair_cache_->ConsumePairs();
    if (pairs.empty()) {
        return;
    }

    std::vector<GeneratedContact> contacts;
    contacts.reserve(128);
    for (const auto& pair : pairs) {
        contacts.clear();
        narrowphase_->BuildContactsForPair(*registry_, pair.first, pair.second, contacts);

        if (contacts.empty()) {
            continue;
        }

        // Apply clustering/K_max per body pair, not globally across all pairs.
        if (clustering_enabled_ && clustering_) {
            clustering_->Apply(contacts);
        }

        for (const auto& contact : contacts) {
            if (!contact.cinfo.modelA || !contact.cinfo.modelB) {
                continue;
            }
            if (!contact.matA || !contact.matB) {
                continue;
            }
            contact_container->AddContact(contact.cinfo, contact.matA, contact.matB);
        }
    }
}

}  // namespace chrono_sdf_contact
