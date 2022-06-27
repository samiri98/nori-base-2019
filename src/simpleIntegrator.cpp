#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList& props) {
        position = props.getPoint("position");
        energy = props.getColor("energy");
    }

    Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        /* Return the component-wise absolute
           value of the shading normal as a color */
        Normal3f n = its.shFrame.n;
        auto diff = (position - its.p);
        float cosTh = n.normalized().dot(diff.normalized());
        cosTh = cosTh < 0 ? 0 : cosTh;
        float squaredNorm = diff.squaredNorm();
        int isVisible = 0;
        if (!scene->rayIntersect(Ray3f(its.p, diff))) {
            isVisible = 1;
        }

        return ((energy * cosTh * isVisible)) / (squaredNorm*4* M_PI * M_PI);
    }

    std::string toString() const {
        return "SimpleIntegrator[]";
    }
protected:
    Point3f position;
    Color3f energy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END