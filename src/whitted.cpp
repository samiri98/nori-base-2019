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


        EmitterQueryRecord q = EmitterQueryRecord(its.p, lightSource);
        auto sample = lightSourceEmitter->sample(Point3f(sampler->next1D(), sampler->next1D(), sampler->next1D()), q);
        int isVisible = 0;
        Intersection its2;
        
        if (scene->rayIntersect(Ray3f(its.p, q.wo), its2)) {
            if (its2.mesh == lightSource) {
                isVisible = 1;
                count++;
            }
        }

        auto geometricTerm = isVisible * fabs(n.dot(q.wo.normalized())) * fmax(0.0f, q.n.normalized().dot(-q.wo.normalized())) / q.wo.squaredNorm();
        auto le = lightSourceEmitter->eval(q);

        Vector3f d = its.shFrame.toLocal(Vector3f(ray.d));
        BSDFQueryRecord bsdfQ(-d);
        bsdfQ.wo = its.shFrame.toLocal(q.wo.normalized());
        bsdfQ.measure = ESolidAngle;
        const BSDF* bsdf = its.mesh->getBSDF();
        Color3f diffuse(0.0f);
        if (bsdf != nullptr) {
            diffuse = bsdf->eval(bsdfQ);
        }
        
        li = geometricTerm * le * diffuse / (lightSourceEmitter->pdf(q) * discretePDF[index]);

        auto itsMeshEmitter = its.mesh->getEmitter();
        if (itsMeshEmitter != nullptr) {
            li = li + itsMeshEmitter->getRadiance(its.mesh);
        }

        return li;
    }

    std::string toString() const {
        return "SimpleIntegrator[]";
    }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END