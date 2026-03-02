#include "SdfPairCache.h"

namespace chrono_sdf_contact {

void SdfPairCache::BeginStep() {
    std::lock_guard<std::mutex> lock(mutex_);
    pairs_.clear();
    seen_.clear();
}

void SdfPairCache::PushPair(chrono::ChCollisionModel* modelA, chrono::ChCollisionModel* modelB) {
    if (!modelA || !modelB) {
        return;
    }

    auto pair = Canonicalize(modelA, modelB);

    std::lock_guard<std::mutex> lock(mutex_);
    if (seen_.insert(pair).second) {
        pairs_.push_back(pair);
    }
}

std::vector<CollisionPair> SdfPairCache::ConsumePairs() {
    std::lock_guard<std::mutex> lock(mutex_);
    auto out = std::move(pairs_);
    pairs_.clear();
    seen_.clear();
    return out;
}

std::size_t SdfPairCache::PairHash::operator()(const CollisionPair& pair) const noexcept {
    const auto h1 = std::hash<const void*>{}(pair.first);
    const auto h2 = std::hash<const void*>{}(pair.second);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U));
}

CollisionPair SdfPairCache::Canonicalize(chrono::ChCollisionModel* modelA, chrono::ChCollisionModel* modelB) {
    if (modelB < modelA) {
        return {modelB, modelA};
    }
    return {modelA, modelB};
}

}  // namespace chrono_sdf_contact

