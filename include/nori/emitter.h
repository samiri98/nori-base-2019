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

#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

struct EmitterQueryRecord {
    /// Incident position
    Point3f p;

    /// Outgoing direction (in the local frame)
    Vector3f wo;

    Vector3f n;

    float pdf;

    const Mesh* mesh;

    EmitterQueryRecord(Point3f p, const Mesh *mesh) : p(p), mesh(mesh) {}
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:

    virtual Color3f sample(const Point2f& sample, EmitterQueryRecord& query) const = 0;

    virtual Color3f eval(EmitterQueryRecord& query) const = 0;

    virtual float pdf(EmitterQueryRecord& query) const = 0;

    virtual Color3f getRadiance(const Mesh* mesh) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

protected:
    Color3f radiance;
};

NORI_NAMESPACE_END
