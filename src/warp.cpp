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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    throw NoriException("Warp::squareToTent() is not yet implemented!");
}

float Warp::squareToTentPdf(const Point2f &p) {
    throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = sqrt(sample[0]);
    float th = 2 * M_PI * sample[1];

    return Point2f(r * cosf(th), r * sinf(th));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return p.squaredNorm() < 1.0f ? INV_PI : 0;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    throw NoriException("Warp::squareToUniformSphere() is not yet implemented!");
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    throw NoriException("Warp::squareToUniformSpherePdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    throw NoriException("Warp::squareToUniformHemisphere() is not yet implemented!");
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    throw NoriException("Warp::squareToUniformHemispherePdf() is not yet implemented!");
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f diskPoint = squareToUniformDisk(sample);
    return Vector3f(diskPoint[0], diskPoint[1], sqrt(1 - sample[0]));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    float cosTh = v.dot(Vector3f(0, 0, 1));
    cosTh = cosTh < 0 ? 0 : cosTh;
    return cosTh / M_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float phi = 2.0f * M_PI * sample[0];
    float logSample = log(1.0f-sample[1]);
    if (std::isinf(logSample)) {
        logSample = 0;
    }

    float tan2Theta = -alpha * alpha * logSample;
    float cosTheta = 1 / sqrt(1+tan2Theta);
    float sinTheta = std::sqrt(std::max(float(0), 1 - cosTheta * cosTheta));
    //float sin2Theta = tan2Theta / (1.0f + tan2Theta);
    //float cos2Theta = 1.0f / (1.0f + tan2Theta);
    //return Vector3f(sqrt(sin2Theta) * cosf(phi), sqrt(sin2Theta) * sinf(phi), sqrt(cos2Theta));
    return Vector3f(sinTheta * cosf(phi), sinTheta * sinf(phi), cosTheta);
}

float Warp::squareToBeckmannPdf(const Vector3f& m, float alpha) {
    float cos2Th = m[2] * m[2];
    float tan2Th = (1- cos2Th)/cos2Th;
    float cos3Th = cos2Th * m[2];
    float num = -tan2Th / (alpha * alpha);
    float denum = alpha * alpha * cos3Th;
    float result = (1.0f / (2.0f * M_PI)) * (2 * exp(num)) / (denum);
    result = result < 0 ? 0 : result;
    result = isfinite(result) ? result : 0;
    return result;
}

NORI_NAMESPACE_END
