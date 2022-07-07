#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class MaterialSamplingIntegrator : public Integrator {
public:
    MaterialSamplingIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        Color3f li(0.0f);
        Color3f L(1.0f);


        Ray3f rayToCheck = ray;
        float acceta = 1.0f;

        while (true) {

            if (!scene->rayIntersect(rayToCheck, its)) {
                break;
            }

            Normal3f n = its.shFrame.n;

            auto itsMeshEmitter = its.mesh->getEmitter();
            if (itsMeshEmitter != nullptr) {
                if (n.dot(-rayToCheck.d.normalized()) > 0) {
                    li += (L * itsMeshEmitter->getRadiance());
                }
            }

            float probability = fmin(L.maxCoeff(), 0.99);

            auto end = sampler->next1D();
            if (end < probability) {
                Vector3f d = its.shFrame.toLocal(Vector3f(rayToCheck.d));
                BSDFQueryRecord bsdfQ(-d.normalized());
                const BSDF* bsdf = its.mesh->getBSDF();
                L /= probability;
                if (bsdf != nullptr) {
                    L *= bsdf->sample(bsdfQ, sampler->next2D());
                    rayToCheck = Ray3f(its.p, its.shFrame.toWorld(bsdfQ.wo));
                }
                else {
                    break;
                }
            }
            else {
                break;
            }
        }

        return li;
    }



    std::string toString() const {
        return "MaterialSamplingIntegrator[]";
    }
};

NORI_REGISTER_CLASS(MaterialSamplingIntegrator, "path_mats");
NORI_NAMESPACE_END