#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        Color3f li(0.0f);
        auto lightSources = scene->getAllEmmiterMeshes();

        if (!scene->rayIntersect(ray, its)) {
            return Color3f(0.0f);
        }

        if (its.mesh->getBSDF() && its.mesh->getBSDF()->isDiffuse()) {
            Normal3f n = its.shFrame.n;

            auto discretePDF = DiscretePDF(lightSources.size());
            for (size_t i = 0; i < lightSources.size(); i++)
            {
                discretePDF.append(1);
            }
            discretePDF.normalize();

            size_t index = discretePDF.sample(sampler->next1D());
            auto lightSource = lightSources[index];
            auto lightSourceEmitter = lightSource->getEmitter();

            auto itsMeshEmitter = its.mesh->getEmitter();
            if (itsMeshEmitter != nullptr) {
                li += itsMeshEmitter->getRadiance(its.mesh);
            }

            EmitterQueryRecord q = EmitterQueryRecord(its.p, lightSource);
            auto sample = lightSourceEmitter->sample(sampler->next2D(), q);
            int isVisible = 0;
            Intersection its2;
            if (scene->rayIntersect(Ray3f(its.p, q.wo.normalized()), its2)) {
                if (its2.mesh == lightSource) {
                    isVisible = 1;
                }
            }

            float geometricTerm = isVisible * fabs(n.dot(q.wo.normalized()) * q.n.normalized().dot(-q.wo.normalized())) / q.wo.squaredNorm();
            Color3f le = lightSourceEmitter->eval(q);

            Vector3f d = its.shFrame.toLocal(Vector3f(ray.d));
            BSDFQueryRecord bsdfQ(-d.normalized());
            bsdfQ.wo = its.shFrame.toLocal(q.wo.normalized());
            bsdfQ.measure = ESolidAngle;
            const BSDF* bsdf = its.mesh->getBSDF();
            Color3f diffuse(0.0f);
            if (bsdf != nullptr) {
                diffuse = bsdf->eval(bsdfQ);
            }

            li += ((geometricTerm * le.cwiseProduct(diffuse)) / (lightSourceEmitter->pdf(q) * discretePDF[index]));
            //li += ((geometricTerm * le.cwiseProduct(diffuse)) / (q.pdf * discretePDF[index]));
        }
        else if (its.mesh->getBSDF()) {
            auto bsdf = its.mesh->getBSDF();
            Vector3f d = its.shFrame.toLocal(Vector3f(ray.d));
            BSDFQueryRecord bsdfQ(-d.normalized());
            bsdfQ.measure = EDiscrete;
            auto c = bsdf->sample(bsdfQ, sampler->next2D());
            auto end = sampler->next1D();
            if (end < 0.95) {
                li = c * Li(scene, sampler, Ray3f(its.p, its.shFrame.toWorld(bsdfQ.wo))) / 0.95f;
                //li = bsdfQ.eta * Li(scene, sampler, Ray3f(its.p, its.shFrame.toWorld(bsdfQ.wo))) / 0.95f;
            }
        }

        return li;
        //return LiRecursive(scene, sampler, ray, nullptr);
    }

    Color3f LiRecursive(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Mesh* currMesh) const {
    
    }



    std::string toString() const {
        return "SimpleIntegrator[]";
    }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END