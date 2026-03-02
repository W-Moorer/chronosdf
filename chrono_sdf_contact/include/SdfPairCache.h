#pragma once

#include <mutex>
#include <unordered_set>
#include <utility>
#include <vector>

#include "chrono/collision/ChCollisionModel.h"

namespace chrono_sdf_contact {

using CollisionPair = std::pair<chrono::ChCollisionModel*, chrono::ChCollisionModel*>;

class SdfPairCache {
  public:
    void BeginStep();
    void PushPair(chrono::ChCollisionModel* modelA, chrono::ChCollisionModel* modelB);
    std::vector<CollisionPair> ConsumePairs();

  private:
    struct PairHash {
        std::size_t operator()(const CollisionPair& pair) const noexcept;
    };

    static CollisionPair Canonicalize(chrono::ChCollisionModel* modelA, chrono::ChCollisionModel* modelB);

    mutable std::mutex mutex_;
    std::vector<CollisionPair> pairs_;
    std::unordered_set<CollisionPair, PairHash> seen_;
};

}  // namespace chrono_sdf_contact

