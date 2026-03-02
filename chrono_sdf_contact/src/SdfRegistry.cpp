#include "SdfRegistry.h"

namespace chrono_sdf_contact {

void SdfRegistry::Register(chrono::ChCollisionModel* model, const SdfProxy& proxy) {
    if (!model) {
        return;
    }

    std::unique_lock<std::shared_mutex> lock(mutex_);
    entries_[model] = proxy;
    revision_.fetch_add(1, std::memory_order_relaxed);
}

void SdfRegistry::QueryPresence(chrono::ChCollisionModel* modelA,
                                chrono::ChCollisionModel* modelB,
                                bool& hasA,
                                bool& hasB) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    hasA = (entries_.find(modelA) != entries_.end());
    hasB = (entries_.find(modelB) != entries_.end());
}

void SdfRegistry::GetPairPtrs(chrono::ChCollisionModel* modelA,
                              chrono::ChCollisionModel* modelB,
                              const SdfProxy*& out_proxyA,
                              const SdfProxy*& out_proxyB) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);

    out_proxyA = nullptr;
    out_proxyB = nullptr;

    auto itA = entries_.find(modelA);
    if (itA != entries_.end()) {
        out_proxyA = &itA->second;
    }

    auto itB = entries_.find(modelB);
    if (itB != entries_.end()) {
        out_proxyB = &itB->second;
    }
}

void SdfRegistry::GetPair(chrono::ChCollisionModel* modelA,
                          chrono::ChCollisionModel* modelB,
                          bool& hasA,
                          SdfProxy& out_proxyA,
                          bool& hasB,
                          SdfProxy& out_proxyB) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);

    auto itA = entries_.find(modelA);
    hasA = (itA != entries_.end());
    if (hasA) {
        out_proxyA = itA->second;
    }

    auto itB = entries_.find(modelB);
    hasB = (itB != entries_.end());
    if (hasB) {
        out_proxyB = itB->second;
    }
}

bool SdfRegistry::Has(chrono::ChCollisionModel* model) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return entries_.find(model) != entries_.end();
}

bool SdfRegistry::Get(chrono::ChCollisionModel* model, SdfProxy& out_proxy) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    auto it = entries_.find(model);
    if (it == entries_.end()) {
        return false;
    }
    out_proxy = it->second;
    return true;
}

std::uint64_t SdfRegistry::GetRevision() const noexcept {
    return revision_.load(std::memory_order_relaxed);
}

void SdfRegistry::Clear() {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    entries_.clear();
    revision_.fetch_add(1, std::memory_order_relaxed);
}

}  // namespace chrono_sdf_contact
