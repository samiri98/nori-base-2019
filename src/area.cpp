#include <nori/emitter.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList& props) {
        radiance = props.getColor("radiance");
    }

    Color3f sample(const Point2f& sample, EmitterQueryRecord& query) const {
        auto mesh = query.mesh;
        Point3f p;
        Vector3f n;
        mesh->squareToArea(sample, p, n, query.pdf);
        query.wo = p - query.p;
        query.n = n;
        return radiance;
    }

    Color3f eval(EmitterQueryRecord& query) const {
        if ((-query.wo).normalized().dot(query.n.normalized()) <= 0) return Color3f(0.f);

        return radiance;
    }

    float pdf(EmitterQueryRecord& query) const {
        auto mesh = query.mesh;
        auto pdf = mesh->squareToAreaPDF();
        return pdf;
    }

    virtual Color3f getRadiance() const {
        return radiance;
    }

    std::string toString() const {
        return "AreaLight[]";
    }
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END