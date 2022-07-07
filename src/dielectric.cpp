/*
    This file is part of Nori, a simple educational ray tracer
    Copyright (c) 2015 by Wenzel Jakob
    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.
    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList& propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord&) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord&) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord& bRec, const Point2f& sample) const {
        // if incident light is parallel to the surface, then then ignore it
        if (Frame::cosTheta(bRec.wi) == 0.0f)
            return Color3f(0.0f);

        float ind1, ind2;
        Vector3f normal;
        bRec.measure = EDiscrete; // type of measurement used 

        // in the case that the ray comes from behing the mesh
        // this first case.
        // intIOR is the index of refraction inside of the object
        // extIOR is the index of refraction outside the object
        if (Frame::cosTheta(bRec.wi) <= 0.0f) {
            ind1 = m_intIOR;
            ind2 = m_extIOR;
            normal = Vector3f(0, 0, -1.0f);
        }
        else {
            ind1 = m_extIOR;
            ind2 = m_intIOR;
            normal = Vector3f(0, 0, 1.0f);
        }

        // This will calculate the fraction of instnsity of the light that is reflected and 
        // refracted based on the interior and exterior index values
        float fresPos = fresnel(fabs(Frame::cosTheta(bRec.wi)), ind1, ind2);

        if (sample.x() > fresPos) {
            float eta = ind1 / ind2; // snells law

            // From slide 26 of lecture 4
            float cost = bRec.wi.dot(normal);
            cost = sqrt(1.0f - eta * eta * (1 - cost * cost));

            Vector3f reflection = -cost * normal; // This is the vector form to obtain refracted (transmitted) direction
            reflection -= eta * (bRec.wi - bRec.wi.dot(normal) * normal); // From slide 31 of lecture 4

            bRec.wo = reflection;
            bRec.eta = eta;
        }
        else {
            bRec.wo = 2 * normal * bRec.wi.dot(normal) - bRec.wi; // From slide 31 of lecture 4
            bRec.eta = 1.0f;
        }

        // if at right angle, return 0
        if (Frame::cosTheta(bRec.wo) == 0.0f) return Color3f(0.0f);
        return Color3f(bRec.eta * bRec.eta);
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END