#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class MultipleImportanceSamplingIntegrator : public Integrator {
public:
    MultipleImportanceSamplingIntegrator(const PropertyList& props) {
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        Color3f li(0.0f);
        Color3f L(1.0f);


        Ray3f rayToCheck = ray;
        float acceta = 1.0f;

        if (!scene->rayIntersect(rayToCheck, its)) {
            return li;
        }
        size_t k = 1;
        bool wasRefractive = true;
        bool firstTime = true;
        while (true) {

            Normal3f n = its.shFrame.n;

            auto itsMeshEmitter = its.mesh->getEmitter();
            if (itsMeshEmitter != nullptr && wasRefractive) {
                if (n.dot(-rayToCheck.d.normalized()) > 0) {
                    li += (L * itsMeshEmitter->getRadiance());
                }
            }

            const BSDF* bsdf = its.mesh->getBSDF();
            auto lightSources = scene->getAllEmmiterMeshes();
            auto discretePDF = DiscretePDF(lightSources.size());
            for (size_t i = 0; i < lightSources.size(); i++)
            {
                discretePDF.append(1);
            }
            discretePDF.normalize();

            size_t index = discretePDF.sample(sampler->next1D());
            auto lightSource = lightSources[index];
            auto lightSourceEmitter = lightSource->getEmitter();

            if (bsdf) {
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
                Color3f diffuse(0.0f);
                if (bsdf != nullptr) {
                    diffuse = bsdf->eval(bsdfQ);
                }

                Color3f directIllumination = ((geometricTerm * le.cwiseProduct(diffuse)) / (lightSourceEmitter->pdf(q) * discretePDF[index]));
                //Color3f directIllumination = ((geometricTerm * le.cwiseProduct(diffuse)) / (q.pdf * discretePDF[index]));
                float weight = 1.0f;
                if (bsdf) {
                    auto plight = (lightSourceEmitter->pdf(q));
                    auto pbsdf = bsdf->pdf(bsdfQ);
                    if (fabs(q.n.normalized().dot(-q.wo.normalized())) != 0) {
                        plight = plight / fabs(q.n.normalized().dot(-q.wo.normalized())) * q.wo.squaredNorm();
                        weight = plight / (plight + pbsdf);
                    }
                }
                li += (L * directIllumination * weight);
            }

            float probability = fmin(L.maxCoeff() * acceta * acceta, 0.99);
            if (firstTime) {
                probability = 1.1f;
                firstTime = false;
            }


            auto end = sampler->next1D();
            if (end < probability) {
                Vector3f d = its.shFrame.toLocal(Vector3f(rayToCheck.d));
                BSDFQueryRecord bsdfQ(-d.normalized());
                if (bsdf != nullptr) {
                    auto brdfTemp = bsdf->sample(bsdfQ, sampler->next2D());
                    rayToCheck = Ray3f(its.p, its.shFrame.toWorld(bsdfQ.wo));
                    if (bsdf) {
                        wasRefractive = !bsdf->isDiffuse();
                    }
                    if (!scene->rayIntersect(rayToCheck, its)) {
                        break;
                    }
                    L /= probability;
                    L *= brdfTemp;
                    acceta *= bsdfQ.eta;
                    if (its.mesh->getEmitter()) {
                        EmitterQueryRecord q = EmitterQueryRecord(its.p, its.mesh);
                        q.wo = bsdfQ.wo;
                        Color3f directIllumination = ((its.mesh->getEmitter()->eval(q).cwiseProduct(brdfTemp) / (discretePDF[index])));
                        float weight = 1.0f;
                        if (bsdf) {
                            auto plight = (its.mesh->getEmitter()->pdf(q));
                            auto pbsdf = bsdf->pdf(bsdfQ);
                            if (fabs(q.n.normalized().dot(-q.wo.normalized())) != 0) {
                                plight = plight / fabs(q.n.normalized().dot(-q.wo.normalized())) * q.wo.squaredNorm();
                                weight = pbsdf / (plight + pbsdf);
                            }
                        }
                        li += (L * directIllumination * weight);
                    }
                }
                else {
                    break;
                }
            }
            else {
                break;
            }
        }

        return li / k;
    }

    std::string toString() const {
        return "MultipleImportanceSamplingIntegrator[]";
    }
};

NORI_REGISTER_CLASS(MultipleImportanceSamplingIntegrator, "path_mis");
NORI_NAMESPACE_END