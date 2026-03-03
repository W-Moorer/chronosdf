#include "SdfNarrowphase.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <unordered_set>
#include <vector>

#include "chrono/collision/ChCollisionShapeBox.h"
#include "chrono/collision/ChCollisionShapeSphere.h"
#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChContactable.h"

#include "SdfRegistry.h"

namespace chrono_sdf_contact {

namespace {

constexpr std::size_t kMaxMeshTrianglesPerShape = 256;
constexpr std::size_t kMaxBoundarySamples = 256;
constexpr std::size_t kBoxDirectionalSamples = 8;
constexpr std::size_t kBoxFaceSamples = 6;
constexpr std::size_t kMaxMeshQueriesPerShape = 128;
constexpr std::size_t kMaxQueriesPerModel = 512;
constexpr double kMeshQueryCellSize = 0.04;
// Clamp extreme one-shot penetration to avoid unstable impulses from outlier samples.
constexpr double kMaxPenetrationDepth = 0.05;

struct BodyPose {
    chrono::ChVector3d pos;
    chrono::ChQuaternion<> rot;
};

struct QuantKey3 {
    std::int64_t ix;
    std::int64_t iy;
    std::int64_t iz;

    bool operator==(const QuantKey3& other) const noexcept {
        return ix == other.ix && iy == other.iy && iz == other.iz;
    }
};

struct QuantKey3Hash {
    std::size_t operator()(const QuantKey3& key) const noexcept {
        const auto h1 = std::hash<std::int64_t>{}(key.ix);
        const auto h2 = std::hash<std::int64_t>{}(key.iy);
        const auto h3 = std::hash<std::int64_t>{}(key.iz);
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        seed ^= h3 + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        return seed;
    }
};

QuantKey3 QuantizePoint(const chrono::ChVector3d& p, double h) {
    return {static_cast<std::int64_t>(std::floor(p.x() / h)), static_cast<std::int64_t>(std::floor(p.y() / h)),
            static_cast<std::int64_t>(std::floor(p.z() / h))};
}

bool ExtractBodyPose(chrono::ChCollisionModel* model, BodyPose& pose) {
    if (!model || !model->GetContactable()) {
        return false;
    }

    auto* body = dynamic_cast<chrono::ChBody*>(model->GetContactable());
    if (!body) {
        return false;
    }

    pose.pos = body->GetPos();
    pose.rot = body->GetRot();
    return true;
}

chrono::ChVector3d ToModelPoint(const chrono::ChCollisionShapeInstance& shape_instance,
                                const chrono::ChVector3d& p_shape_local) {
    return shape_instance.frame.GetPos() + shape_instance.frame.GetRot().Rotate(p_shape_local);
}

void AppendBoundarySample(std::vector<chrono::ChVector3d>& samples, const chrono::ChVector3d& p_model) {
    if (samples.size() >= kMaxBoundarySamples) {
        return;
    }
    samples.push_back(p_model);
}

void AppendSphereBoundarySamples(const chrono::ChCollisionShapeInstance& instance,
                                 double radius,
                                 std::vector<chrono::ChVector3d>& out_samples) {
    static const std::array<chrono::ChVector3d, 7> kDirs = {
        chrono::VNULL, chrono::ChVector3d(1.0, 0.0, 0.0),  chrono::ChVector3d(-1.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 1.0, 0.0), chrono::ChVector3d(0.0, -1.0, 0.0), chrono::ChVector3d(0.0, 0.0, 1.0),
        chrono::ChVector3d(0.0, 0.0, -1.0)};

    for (const auto& d : kDirs) {
        AppendBoundarySample(out_samples, ToModelPoint(instance, radius * d));
    }
}

void AppendBoxBoundarySamples(const chrono::ChCollisionShapeInstance& instance,
                              const chrono::ChVector3d& h,
                              std::vector<chrono::ChVector3d>& out_samples) {
    static const std::array<chrono::ChVector3d, kBoxDirectionalSamples> kDirs = {
        chrono::ChVector3d(-1.0, -1.0, -1.0), chrono::ChVector3d(1.0, -1.0, -1.0),
        chrono::ChVector3d(-1.0, 1.0, -1.0),  chrono::ChVector3d(1.0, 1.0, -1.0),
        chrono::ChVector3d(-1.0, -1.0, 1.0),  chrono::ChVector3d(1.0, -1.0, 1.0),
        chrono::ChVector3d(-1.0, 1.0, 1.0),   chrono::ChVector3d(1.0, 1.0, 1.0)};
    static const std::array<chrono::ChVector3d, kBoxFaceSamples> kFaceDirs = {
        chrono::ChVector3d(1.0, 0.0, 0.0),  chrono::ChVector3d(-1.0, 0.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0),
        chrono::ChVector3d(0.0, -1.0, 0.0), chrono::ChVector3d(0.0, 0.0, 1.0),  chrono::ChVector3d(0.0, 0.0, -1.0)};

    AppendBoundarySample(out_samples, ToModelPoint(instance, chrono::VNULL));

    for (const auto& d : kDirs) {
        const chrono::ChVector3d p_shape_local(d.x() * h.x(), d.y() * h.y(), d.z() * h.z());
        AppendBoundarySample(out_samples, ToModelPoint(instance, p_shape_local));
    }
    for (const auto& d : kFaceDirs) {
        const chrono::ChVector3d p_shape_local(d.x() * h.x(), d.y() * h.y(), d.z() * h.z());
        AppendBoundarySample(out_samples, ToModelPoint(instance, p_shape_local));
    }
}

void AppendMeshBoundarySamples(const chrono::ChCollisionShapeInstance& instance,
                               const std::shared_ptr<chrono::ChTriangleMesh>& mesh,
                               std::vector<chrono::ChVector3d>& out_samples) {
    if (!mesh) {
        return;
    }

    const std::size_t max_triangles = std::min<std::size_t>(mesh->GetNumTriangles(), kMaxMeshTrianglesPerShape);
    if (max_triangles == 0) {
        return;
    }

    const std::size_t target_samples = std::max<std::size_t>(1, kMaxBoundarySamples / 2);
    const std::size_t stride = std::max<std::size_t>(1, max_triangles / target_samples);

    for (std::size_t i = 0; i < max_triangles && out_samples.size() < kMaxBoundarySamples; i += stride) {
        auto tri = mesh->GetTriangle(static_cast<unsigned int>(i));
        const auto centroid_shape = (tri.p1 + tri.p2 + tri.p3) / 3.0;
        AppendBoundarySample(out_samples, ToModelPoint(instance, centroid_shape));
    }
}

void BuildBoundarySamplesFromModel(chrono::ChCollisionModel* model, std::vector<chrono::ChVector3d>& out_samples) {
    out_samples.clear();
    if (!model) {
        out_samples.push_back(chrono::VNULL);
        return;
    }

    const auto& shape_instances = model->GetShapeInstances();
    for (const auto& instance : shape_instances) {
        if (!instance.shape || out_samples.size() >= kMaxBoundarySamples) {
            continue;
        }

        if (instance.shape->GetType() == chrono::ChCollisionShape::Type::SPHERE) {
            auto sphere_shape = std::static_pointer_cast<chrono::ChCollisionShapeSphere>(instance.shape);
            AppendSphereBoundarySamples(instance, sphere_shape->GetRadius(), out_samples);
            continue;
        }

        if (instance.shape->GetType() == chrono::ChCollisionShape::Type::BOX) {
            auto box_shape = std::static_pointer_cast<chrono::ChCollisionShapeBox>(instance.shape);
            AppendBoxBoundarySamples(instance, box_shape->GetHalflengths(), out_samples);
            continue;
        }

        if (instance.shape->GetType() == chrono::ChCollisionShape::Type::TRIANGLEMESH) {
            auto mesh_shape = std::static_pointer_cast<chrono::ChCollisionShapeTriangleMesh>(instance.shape);
            AppendMeshBoundarySamples(instance, mesh_shape->GetMesh(), out_samples);
            continue;
        }
    }

    if (out_samples.empty()) {
        out_samples.push_back(chrono::VNULL);
    }
}

void AppendPointContact(const chrono::ChVector3d& pA_world,
                        chrono::ChCollisionModel* modelA,
                        chrono::ChCollisionModel* modelB,
                        const SdfProxy& proxy_sdf,
                        const BodyPose& sdf_pose,
                        const std::shared_ptr<chrono::ChContactMaterial>& matA,
                        std::vector<GeneratedContact>& out_contacts) {
    if (!proxy_sdf.sdf || !proxy_sdf.material || !matA) {
        return;
    }

    auto sample = proxy_sdf.sdf->SampleBodyPose(sdf_pose.pos, sdf_pose.rot, pA_world);
    if (sample.grad.Length2() < 1e-16) {
        return;
    }

    const double distance = std::max(sample.phi, -kMaxPenetrationDepth);
    if (distance >= 0.0) {
        return;
    }

    chrono::ChCollisionInfo cinfo;
    cinfo.modelA = modelA;
    cinfo.modelB = modelB;
    cinfo.shapeA = nullptr;
    cinfo.shapeB = nullptr;

    const auto nB_world = sample.grad;
    cinfo.vN = -nB_world;  // respect to A, A->B
    cinfo.vpA = pA_world;
    cinfo.vpB = pA_world - distance * nB_world;
    cinfo.distance = distance;
    cinfo.eff_radius = chrono::ChCollisionInfo::GetDefaultEffectiveCurvatureRadius();
    cinfo.reaction_cache = nullptr;

    GeneratedContact generated;
    generated.cinfo = cinfo;
    generated.matA = matA;
    generated.matB = proxy_sdf.material;
    out_contacts.push_back(std::move(generated));
}

void AppendSphereContact(const chrono::ChVector3d& center_world,
                         double radius,
                         chrono::ChCollisionModel* modelA,
                         chrono::ChCollisionModel* modelB,
                         const SdfProxy& proxy_sdf,
                         const BodyPose& sdf_pose,
                         const std::shared_ptr<chrono::ChContactMaterial>& matA,
                         std::vector<GeneratedContact>& out_contacts) {
    if (!proxy_sdf.sdf || !proxy_sdf.material || !matA) {
        return;
    }

    auto sample = proxy_sdf.sdf->SampleBodyPose(sdf_pose.pos, sdf_pose.rot, center_world);
    if (sample.grad.Length2() < 1e-16) {
        return;
    }

    const double distance = std::max(sample.phi - radius, -kMaxPenetrationDepth);
    if (distance >= 0.0) {
        return;
    }

    chrono::ChCollisionInfo cinfo;
    cinfo.modelA = modelA;
    cinfo.modelB = modelB;
    cinfo.shapeA = nullptr;
    cinfo.shapeB = nullptr;

    const auto nB_world = sample.grad;
    cinfo.vN = -nB_world;  // respect to A, A->B
    cinfo.vpA = center_world + radius * cinfo.vN;
    cinfo.vpB = center_world - (distance + radius) * nB_world;
    cinfo.distance = distance;
    cinfo.eff_radius = chrono::ChCollisionInfo::GetDefaultEffectiveCurvatureRadius();
    cinfo.reaction_cache = nullptr;

    GeneratedContact generated;
    generated.cinfo = cinfo;
    generated.matA = matA;
    generated.matB = proxy_sdf.material;
    out_contacts.push_back(std::move(generated));
}

void AppendSdfSdfContactAtoB(const chrono::ChVector3d& pA_world,
                             chrono::ChCollisionModel* modelA,
                             chrono::ChCollisionModel* modelB,
                             const SdfProxy& proxyA,
                             const SdfProxy& proxyB,
                             const BodyPose& poseB,
                             std::vector<GeneratedContact>& out_contacts) {
    auto sampleB = proxyB.sdf->SampleBodyPose(poseB.pos, poseB.rot, pA_world);
    if (sampleB.grad.Length2() < 1e-16) {
        return;
    }

    if (sampleB.phi >= 0.0) {
        return;
    }
    const double distance = std::max(sampleB.phi, -kMaxPenetrationDepth);

    chrono::ChCollisionInfo cinfo;
    cinfo.modelA = modelA;
    cinfo.modelB = modelB;
    cinfo.shapeA = nullptr;
    cinfo.shapeB = nullptr;
    cinfo.vN = -sampleB.grad;
    cinfo.vpA = pA_world;
    cinfo.vpB = pA_world - distance * sampleB.grad;
    cinfo.distance = distance;
    cinfo.eff_radius = chrono::ChCollisionInfo::GetDefaultEffectiveCurvatureRadius();
    cinfo.reaction_cache = nullptr;

    GeneratedContact out;
    out.cinfo = cinfo;
    out.matA = proxyA.material;
    out.matB = proxyB.material;
    out_contacts.push_back(std::move(out));
}

void AppendSdfSdfContactBtoA(const chrono::ChVector3d& pB_world,
                             chrono::ChCollisionModel* modelA,
                             chrono::ChCollisionModel* modelB,
                             const SdfProxy& proxyA,
                             const SdfProxy& proxyB,
                             const BodyPose& poseA,
                             std::vector<GeneratedContact>& out_contacts) {
    auto sampleA = proxyA.sdf->SampleBodyPose(poseA.pos, poseA.rot, pB_world);
    if (sampleA.grad.Length2() < 1e-16) {
        return;
    }

    if (sampleA.phi >= 0.0) {
        return;
    }
    const double distance = std::max(sampleA.phi, -kMaxPenetrationDepth);

    chrono::ChCollisionInfo cinfo;
    cinfo.modelA = modelA;
    cinfo.modelB = modelB;
    cinfo.shapeA = nullptr;
    cinfo.shapeB = nullptr;
    cinfo.vN = sampleA.grad;
    cinfo.vpA = pB_world - distance * sampleA.grad;
    cinfo.vpB = pB_world;
    cinfo.distance = distance;
    cinfo.eff_radius = chrono::ChCollisionInfo::GetDefaultEffectiveCurvatureRadius();
    cinfo.reaction_cache = nullptr;

    GeneratedContact out;
    out.cinfo = cinfo;
    out.matA = proxyA.material;
    out.matB = proxyB.material;
    out_contacts.push_back(std::move(out));
}

void BuildSdfSdfContacts(chrono::ChCollisionModel* modelA,
                         chrono::ChCollisionModel* modelB,
                         const SdfProxy& proxyA,
                         const SdfProxy& proxyB,
                         const std::vector<chrono::ChVector3d>& samplesA,
                         const std::vector<chrono::ChVector3d>& samplesB,
                         SdfSdfSamplingMode mode,
                         std::vector<GeneratedContact>& out_contacts) {
    if (!modelA || !modelB || !modelA->GetContactable() || !modelB->GetContactable()) {
        return;
    }
    if (!proxyA.sdf || !proxyB.sdf || !proxyA.material || !proxyB.material) {
        return;
    }

    BodyPose poseA;
    BodyPose poseB;
    if (!ExtractBodyPose(modelA, poseA) || !ExtractBodyPose(modelB, poseB)) {
        return;
    }

    if (mode == SdfSdfSamplingMode::Bidirectional || mode == SdfSdfSamplingMode::AtoBOnly) {
        for (const auto& pA_local : samplesA) {
            const auto pA_world = poseA.pos + poseA.rot.Rotate(pA_local);
            AppendSdfSdfContactAtoB(pA_world, modelA, modelB, proxyA, proxyB, poseB, out_contacts);
        }
    }

    if (mode == SdfSdfSamplingMode::Bidirectional || mode == SdfSdfSamplingMode::BtoAOnly) {
        for (const auto& pB_local : samplesB) {
            const auto pB_world = poseB.pos + poseB.rot.Rotate(pB_local);
            AppendSdfSdfContactBtoA(pB_world, modelA, modelB, proxyA, proxyB, poseA, out_contacts);
        }
    }
}

}  // namespace

const std::vector<chrono::ChVector3d>& SdfNarrowphase::GetBoundarySamplesCached(chrono::ChCollisionModel* model,
                                                                                  const SdfProxy& proxy) const {
    if (!proxy.boundary_samples_local.empty()) {
        return proxy.boundary_samples_local;
    }

    if (!cache_enabled_) {
        BuildBoundarySamplesFromModel(model, boundary_samples_scratch_);
        return boundary_samples_scratch_;
    }

    auto it = boundary_samples_cache_.find(model);
    if (it != boundary_samples_cache_.end()) {
        return it->second;
    }

    std::vector<chrono::ChVector3d> samples;
    BuildBoundarySamplesFromModel(model, samples);
    auto inserted = boundary_samples_cache_.emplace(model, std::move(samples));
    return inserted.first->second;
}

double SdfNarrowphase::GetModelBoundRadiusCached(chrono::ChCollisionModel* model) const {
    if (!cache_enabled_) {
        return ComputeModelBoundRadius(model);
    }

    auto it = model_bound_radius_cache_.find(model);
    if (it != model_bound_radius_cache_.end()) {
        return it->second;
    }

    const double radius = ComputeModelBoundRadius(model);
    model_bound_radius_cache_[model] = radius;
    return radius;
}

double SdfNarrowphase::ComputeModelBoundRadius(chrono::ChCollisionModel* model) const {
    double radius = 0.0;
    if (model) {
        const auto& shape_instances = model->GetShapeInstances();
        for (const auto& instance : shape_instances) {
            if (!instance.shape) {
                continue;
            }

            const auto center_model = instance.frame.GetPos();
            const auto type = instance.shape->GetType();

            if (type == chrono::ChCollisionShape::Type::SPHERE) {
                auto sphere_shape = std::static_pointer_cast<chrono::ChCollisionShapeSphere>(instance.shape);
                radius = std::max(radius, center_model.Length() + sphere_shape->GetRadius());
                continue;
            }

            if (type == chrono::ChCollisionShape::Type::BOX) {
                auto box_shape = std::static_pointer_cast<chrono::ChCollisionShapeBox>(instance.shape);
                radius = std::max(radius, center_model.Length() + box_shape->GetHalflengths().Length());
                continue;
            }

            if (type == chrono::ChCollisionShape::Type::TRIANGLEMESH) {
                auto mesh_shape = std::static_pointer_cast<chrono::ChCollisionShapeTriangleMesh>(instance.shape);
                auto mesh = mesh_shape->GetMesh();
                if (!mesh) {
                    continue;
                }

                const std::size_t num_triangles = mesh->GetNumTriangles();
                for (std::size_t i = 0; i < num_triangles; ++i) {
                    auto tri = mesh->GetTriangle(static_cast<unsigned int>(i));
                    radius = std::max(radius, ToModelPoint(instance, tri.p1).Length());
                    radius = std::max(radius, ToModelPoint(instance, tri.p2).Length());
                    radius = std::max(radius, ToModelPoint(instance, tri.p3).Length());
                }
            }
        }
    }
    return radius;
}

const SdfNarrowphase::OtherModelSupportCache& SdfNarrowphase::GetOtherModelSupportCached(
    chrono::ChCollisionModel* model) const {
    static const OtherModelSupportCache kEmptyCache{};
    if (!model) {
        return kEmptyCache;
    }

    if (!cache_enabled_) {
        other_support_scratch_ = BuildOtherModelSupport(model);
        return other_support_scratch_;
    }

    const auto& shape_instances = model->GetShapeInstances();
    auto it = other_support_cache_.find(model);
    if (it != other_support_cache_.end() && it->second.shape_count == shape_instances.size()) {
        return it->second;
    }

    auto inserted = other_support_cache_.insert_or_assign(model, BuildOtherModelSupport(model));
    return inserted.first->second;
}

SdfNarrowphase::OtherModelSupportCache SdfNarrowphase::BuildOtherModelSupport(chrono::ChCollisionModel* model) const {
    OtherModelSupportCache cache;
    if (!model) {
        return cache;
    }

    const auto& shape_instances = model->GetShapeInstances();
    cache.shape_count = shape_instances.size();
    cache.queries.reserve(std::min<std::size_t>(kMaxQueriesPerModel, shape_instances.size() * 8));

    auto append_point_query = [&](const chrono::ChVector3d& p_model,
                                  const std::shared_ptr<chrono::ChContactMaterial>& material) {
        if (!material || cache.queries.size() >= kMaxQueriesPerModel) {
            return;
        }
        OtherModelSupportQuery q;
        q.p_model = p_model;
        q.radius = 0.0;
        q.sphere = false;
        q.material = material;
        cache.queries.push_back(std::move(q));
    };

    auto append_sphere_query = [&](const chrono::ChVector3d& center_model,
                                   double radius,
                                   const std::shared_ptr<chrono::ChContactMaterial>& material) {
        if (!material || radius <= 0.0 || cache.queries.size() >= kMaxQueriesPerModel) {
            return;
        }
        OtherModelSupportQuery q;
        q.p_model = center_model;
        q.radius = radius;
        q.sphere = true;
        q.material = material;
        cache.queries.push_back(std::move(q));
    };

    auto append_mesh_queries = [&](const chrono::ChCollisionShapeInstance& instance,
                                   const std::shared_ptr<chrono::ChTriangleMesh>& mesh,
                                   const std::shared_ptr<chrono::ChContactMaterial>& material) {
        if (!mesh || !material || cache.queries.size() >= kMaxQueriesPerModel) {
            return;
        }

        const std::size_t max_triangles = std::min<std::size_t>(mesh->GetNumTriangles(), kMaxMeshTrianglesPerShape);
        if (max_triangles == 0) {
            return;
        }

        const std::size_t max_remaining = kMaxQueriesPerModel - cache.queries.size();
        if (max_remaining == 0) {
            return;
        }

        const std::size_t target_samples = std::min<std::size_t>(kMaxMeshQueriesPerShape, max_remaining);
        const std::size_t stride = std::max<std::size_t>(1, max_triangles / std::max<std::size_t>(1, target_samples));

        std::unordered_set<QuantKey3, QuantKey3Hash> seen_cells;
        seen_cells.reserve(target_samples * 2);

        for (std::size_t i = 0; i < max_triangles; i += stride) {
            if (cache.queries.size() >= kMaxQueriesPerModel) {
                break;
            }
            auto tri = mesh->GetTriangle(static_cast<unsigned int>(i));
            const auto centroid_shape = (tri.p1 + tri.p2 + tri.p3) / 3.0;
            const auto centroid_model = ToModelPoint(instance, centroid_shape);
            if (!seen_cells.insert(QuantizePoint(centroid_model, kMeshQueryCellSize)).second) {
                continue;
            }

            append_point_query(centroid_model, material);
        }
    };

    static const std::array<chrono::ChVector3d, kBoxDirectionalSamples> kBoxDirs = {
        chrono::ChVector3d(-1.0, -1.0, -1.0), chrono::ChVector3d(1.0, -1.0, -1.0),
        chrono::ChVector3d(-1.0, 1.0, -1.0),  chrono::ChVector3d(1.0, 1.0, -1.0),
        chrono::ChVector3d(-1.0, -1.0, 1.0),  chrono::ChVector3d(1.0, -1.0, 1.0),
        chrono::ChVector3d(-1.0, 1.0, 1.0),   chrono::ChVector3d(1.0, 1.0, 1.0)};
    static const std::array<chrono::ChVector3d, kBoxFaceSamples> kBoxFaceDirs = {
        chrono::ChVector3d(1.0, 0.0, 0.0),  chrono::ChVector3d(-1.0, 0.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0),
        chrono::ChVector3d(0.0, -1.0, 0.0), chrono::ChVector3d(0.0, 0.0, 1.0),  chrono::ChVector3d(0.0, 0.0, -1.0)};

    for (const auto& instance : shape_instances) {
        if (!instance.shape || cache.queries.size() >= kMaxQueriesPerModel) {
            continue;
        }

        const auto material = instance.shape->GetMaterial();
        if (!material) {
            continue;
        }

        if (instance.shape->GetType() == chrono::ChCollisionShape::Type::SPHERE) {
            auto sphere_shape = std::static_pointer_cast<chrono::ChCollisionShapeSphere>(instance.shape);
            append_sphere_query(ToModelPoint(instance, chrono::VNULL), sphere_shape->GetRadius(), material);
            continue;
        }

        if (instance.shape->GetType() == chrono::ChCollisionShape::Type::BOX) {
            auto box_shape = std::static_pointer_cast<chrono::ChCollisionShapeBox>(instance.shape);
            const auto h = box_shape->GetHalflengths();
            append_point_query(ToModelPoint(instance, chrono::VNULL), material);
            for (const auto& d : kBoxDirs) {
                const chrono::ChVector3d p_shape_local(d.x() * h.x(), d.y() * h.y(), d.z() * h.z());
                append_point_query(ToModelPoint(instance, p_shape_local), material);
            }
            for (const auto& d : kBoxFaceDirs) {
                const chrono::ChVector3d p_shape_local(d.x() * h.x(), d.y() * h.y(), d.z() * h.z());
                append_point_query(ToModelPoint(instance, p_shape_local), material);
            }
            continue;
        }

        if (instance.shape->GetType() == chrono::ChCollisionShape::Type::TRIANGLEMESH) {
            auto mesh_shape = std::static_pointer_cast<chrono::ChCollisionShapeTriangleMesh>(instance.shape);
            append_mesh_queries(instance, mesh_shape->GetMesh(), material);
            continue;
        }
    }

    cache.bound_radius = GetModelBoundRadiusCached(model);
    return cache;
}

void SdfNarrowphase::SetSdfSdfSamplingMode(SdfSdfSamplingMode mode) {
    sdf_sdf_sampling_mode_ = mode;
}

SdfSdfSamplingMode SdfNarrowphase::GetSdfSdfSamplingMode() const {
    return sdf_sdf_sampling_mode_;
}

void SdfNarrowphase::SetCacheEnabled(bool enabled) {
    if (cache_enabled_ == enabled) {
        return;
    }
    cache_enabled_ = enabled;
    boundary_samples_cache_.clear();
    model_bound_radius_cache_.clear();
    other_support_cache_.clear();
}

bool SdfNarrowphase::IsCacheEnabled() const {
    return cache_enabled_;
}

void SdfNarrowphase::SetBoundCullEnabled(bool enabled) {
    bound_cull_enabled_ = enabled;
}

bool SdfNarrowphase::IsBoundCullEnabled() const {
    return bound_cull_enabled_;
}

void SdfNarrowphase::BuildContactsForPair(const SdfRegistry& registry,
                                          chrono::ChCollisionModel* modelA,
                                          chrono::ChCollisionModel* modelB,
                                          std::vector<GeneratedContact>& out_contacts) const {
    if (!modelA || !modelB) {
        return;
    }

    const SdfProxy* proxyA_ptr = nullptr;
    const SdfProxy* proxyB_ptr = nullptr;
    registry.GetPairPtrs(modelA, modelB, proxyA_ptr, proxyB_ptr);
    const bool hasA = (proxyA_ptr != nullptr);
    const bool hasB = (proxyB_ptr != nullptr);

    if (hasA && hasB) {
        const auto& proxyA = *proxyA_ptr;
        const auto& proxyB = *proxyB_ptr;
        const auto& samplesA = GetBoundarySamplesCached(modelA, proxyA);
        const auto& samplesB = GetBoundarySamplesCached(modelB, proxyB);
        BuildSdfSdfContacts(modelA, modelB, proxyA, proxyB, samplesA, samplesB, sdf_sdf_sampling_mode_, out_contacts);
        return;
    }

    // Current scope: SDF-vs-(sphere/box/trianglemesh).
    if (!(hasA || hasB)) {
        return;
    }

    chrono::ChCollisionModel* model_other = hasA ? modelB : modelA;
    chrono::ChCollisionModel* model_sdf = hasA ? modelA : modelB;
    const SdfProxy& proxy_sdf = hasA ? *proxyA_ptr : *proxyB_ptr;

    if (!proxy_sdf.sdf || !proxy_sdf.material || !model_sdf->GetContactable()) {
        return;
    }
    if (!model_other || !model_other->GetContactable()) {
        return;
    }

    BodyPose sdf_pose;
    BodyPose other_pose;
    if (!ExtractBodyPose(model_sdf, sdf_pose) || !ExtractBodyPose(model_other, other_pose)) {
        return;
    }

    const auto& support = GetOtherModelSupportCached(model_other);
    if (support.queries.empty()) {
        return;
    }

    if (bound_cull_enabled_ && support.bound_radius > 0.0) {
        const auto center_sample = proxy_sdf.sdf->SampleBodyPose(sdf_pose.pos, sdf_pose.rot, other_pose.pos);
        if (center_sample.phi > support.bound_radius) {
            return;
        }
    }

    for (const auto& q : support.queries) {
        if (!q.material) {
            continue;
        }

        const auto p_world = other_pose.pos + other_pose.rot.Rotate(q.p_model);
        if (q.sphere) {
            AppendSphereContact(p_world, q.radius, model_other, model_sdf, proxy_sdf, sdf_pose, q.material,
                                out_contacts);
        } else {
            AppendPointContact(p_world, model_other, model_sdf, proxy_sdf, sdf_pose, q.material, out_contacts);
        }
    }
}

}  // namespace chrono_sdf_contact
