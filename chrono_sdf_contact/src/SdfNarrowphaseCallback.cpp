#include "SdfNarrowphaseCallback.h"

namespace chrono_sdf_contact {

namespace {

bool HasSupportedShape(const chrono::ChCollisionModel* model) {
    if (!model) {
        return false;
    }

    const auto& shape_instances = model->GetShapeInstances();
    for (const auto& instance : shape_instances) {
        if (!instance.shape) {
            continue;
        }

        const auto type = instance.shape->GetType();
        if (type == chrono::ChCollisionShape::Type::SPHERE || type == chrono::ChCollisionShape::Type::BOX ||
            type == chrono::ChCollisionShape::Type::TRIANGLEMESH) {
            return true;
        }
    }

    return false;
}

}  // namespace

SdfNarrowphaseCallback::SdfNarrowphaseCallback(std::shared_ptr<SdfRegistry> registry,
                                               std::shared_ptr<SdfPairCache> pair_cache)
    : registry_(std::move(registry)), pair_cache_(std::move(pair_cache)) {
    if (registry_) {
        cached_registry_revision_ = registry_->GetRevision();
    }
}

void SdfNarrowphaseCallback::RefreshCachesIfRegistryChanged() {
    if (!registry_) {
        return;
    }

    const auto revision = registry_->GetRevision();
    if (revision == cached_registry_revision_) {
        return;
    }

    cached_registry_revision_ = revision;
    sdf_presence_cache_.clear();
    supported_shape_cache_.clear();
}

bool SdfNarrowphaseCallback::IsSdfModelCached(chrono::ChCollisionModel* model) {
    if (!model || !registry_) {
        return false;
    }

    auto it = sdf_presence_cache_.find(model);
    if (it != sdf_presence_cache_.end()) {
        return it->second;
    }

    const bool has_sdf = registry_->Has(model);
    sdf_presence_cache_[model] = has_sdf;
    return has_sdf;
}

bool SdfNarrowphaseCallback::HasSupportedShapeCached(chrono::ChCollisionModel* model) {
    if (!model) {
        return false;
    }

    auto it = supported_shape_cache_.find(model);
    if (it != supported_shape_cache_.end()) {
        return it->second;
    }

    const bool supported = HasSupportedShape(model);
    supported_shape_cache_[model] = supported;
    return supported;
}

bool SdfNarrowphaseCallback::OnNarrowphase(chrono::ChCollisionInfo& cinfo) {
    if (!registry_ || !pair_cache_) {
        return true;
    }

    if (!cinfo.modelA || !cinfo.modelB) {
        return true;
    }

    RefreshCachesIfRegistryChanged();

    const bool hasA = IsSdfModelCached(cinfo.modelA);
    const bool hasB = IsSdfModelCached(cinfo.modelB);

    if (hasA && hasB) {
        pair_cache_->PushPair(cinfo.modelA, cinfo.modelB);
        return false;
    }

    if (!(hasA || hasB)) {
        return true;
    }

    chrono::ChCollisionModel* model_other = hasA ? cinfo.modelB : cinfo.modelA;
    if (!HasSupportedShapeCached(model_other)) {
        return true;
    }

    pair_cache_->PushPair(cinfo.modelA, cinfo.modelB);
    return false;
}

}  // namespace chrono_sdf_contact
