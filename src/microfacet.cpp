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
#include <nori/warp.h>
#include <nori/dpdf.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    float getG1(Vector3f w1, Vector3f w2) const {
        float chi = w1.dot(w2) / Frame::cosTheta(w1);
        chi = chi > 0 ? 1 : 0;
        auto b = m_alpha * Frame::tanTheta(w1);
        b = 1.0f / b;
        auto coeff = b < 1.6 ? ((3.535f*b+2.181f*b*b)/ (1+2.276f * b + 2.577f * b * b)) : 1.0f;
        return chi * coeff;
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        auto wi = bRec.wi;
        auto wo = bRec.wo;
        auto wh = (wi + wo).normalized();
        auto F = fresnel(wh.dot(wi), m_intIOR, m_extIOR);
        auto denum = 4.0f * Frame::cosTheta(wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh);
        auto D = Warp::squareToBeckmannPdf(wh, m_alpha);
        auto G = getG1(wi, wh) * getG1(wo, wh);
        return (m_kd / M_PI) + (Color3f(1.0f) * m_ks * D * F * G)/(denum);
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
    	//throw NoriException("MicrofacetBRDF::pdf(): not implemented!");
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        auto wi = bRec.wi;
        auto wo = bRec.wo;
        auto wh = (wi + wo).normalized();
        auto Jh = 1.0f / (4 * (wh.dot(wo)));
        return (m_ks * Warp::squareToBeckmannPdf(wh, m_alpha) * Jh) + ((1-m_ks) * Frame::cosTheta(wo)/M_PI);
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
    	//throw NoriException("MicrofacetBRDF::sample(): not implemented!");
        if (Frame::cosTheta(bRec.wi) <= 0.0f) return Color3f(0.0f);

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;
        bRec.measure = ESolidAngle;

        auto discretePDF = DiscretePDF(2);
        discretePDF.append(1 - m_ks);
        discretePDF.append(m_ks);
        discretePDF.normalize();
        float reusedSample = _sample[0];
        size_t isSpecular = discretePDF.sampleReuse(reusedSample);
        if (isSpecular) {
            auto normal = Warp::squareToBeckmann(Point2f(reusedSample,_sample[1]), m_alpha).normalized();
            bRec.wo = (- bRec.wi) + 2 * normal * bRec.wi.dot(normal);
        }
        else {

            /* Warp a uniformly distributed sample on [0,1]^2
               to a direction on a cosine-weighted hemisphere */
            bRec.wo = Warp::squareToCosineHemisphere(Point2f(reusedSample, _sample[1]));
        }


        float pdfVal = pdf(bRec);

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        //return eval(bRec) * Frame::cosTheta(bRec.wo) / pdfVal;
        if (pdfVal > 0) return eval(bRec) * Frame::cosTheta(bRec.wo) / pdfVal;
        else return Color3f(0.0f);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
