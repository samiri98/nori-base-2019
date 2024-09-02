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

#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    m_mesh.push_back(mesh);
    for (auto mesh : m_mesh) {
        m_bbox.expandBy(mesh->getBoundingBox());
    }
}

// todo: parallelize this code
Accel::Node* Accel::build(BoundingBox3f box, std::vector<std::vector<uint32_t>> triangles, uint32_t depth) {
    size_t sizeOfTriangles = 0;
    for (size_t i = 0; i < triangles.size(); i++)
    {
        sizeOfTriangles += triangles[i].size();
    }
    
    if (sizeOfTriangles == 0)
        return nullptr;

    if (sizeOfTriangles <= 10 || depth > 10) {
        Accel::ChildNode* leaf = new Accel::ChildNode();
        leaf->type = 1;
        leaf->box = box;
        leaf->triangles = triangles;
        return leaf;
    }

    std::vector<std::vector<uint32_t>> list[8];
    std::vector<Point3f> dividors = { {0,0,0},
                                      {0,1,0},
                                      {1,0,0},
                                      {0,0,1},
                                      {1,1,0},
                                      {0,1,1},
                                      {1,0,1},
                                      {1,1,1} };

    std::vector<BoundingBox3f> sub_bbox;
    auto extent = box.getExtents();
    

    for (int i = 0; i < 8; i++) {
        auto newMax = box.max - (Point3f(extent[0] / 2, extent[1] / 2, extent[2] / 2).cwiseProduct(dividors[i]));
        auto newMin = newMax - (extent / 2);
        BoundingBox3f i_bbox(newMin, newMax);
        sub_bbox.push_back(i_bbox);
    }

    for (int j = 0; j < 8; j++) {
        std::vector<std::vector<uint32_t>> trianglesPerMesh;
        for (int i = 0; i < m_mesh.size(); i++) {
            std::vector<uint32_t> trianglesTemp;
            for (auto triangle : triangles[i]) {
                BoundingBox3f t_bbox = m_mesh[i]->getBoundingBox(triangle);
                if (sub_bbox[j].overlaps(t_bbox, false)) {
                    trianglesTemp.push_back(triangle);
                }
            }
            trianglesPerMesh.push_back(trianglesTemp);
        }
        list[j] = trianglesPerMesh;
    }
    
    Accel::ParentNode* node = new ParentNode();
    node->box = box;
    node->type = 0;
    for (int i = 0; i < 8; ++i) {
        node->children[i] = build(sub_bbox[i], list[i], depth+1);
    }

    return node;
}


void Accel::build() {
    /* Nothing to do here for now */
    std::vector<std::vector<uint32_t>> v;
    for (size_t i = 0; i < m_mesh.size(); i++) {
        std::vector<uint32_t> temp;
        for (size_t j = 0; j < m_mesh[i]->getTriangleCount(); j++) {
            temp.push_back(j);
        }
        v.push_back(temp);
    }
    m_root = build(m_bbox, v, 0);
}

// todo: early stopping by ordered traversal
bool Accel::octree_traversal(bool shadowRay, Ray3f &ray, Intersection& its, uint32_t &f, Node *node) const {
    bool foundIntersection = false;
    if (node == nullptr) {
        return false;
    }

    if (node->box.rayIntersect(ray)) {
        if (node->type == 0) {
            ParentNode* p = (ParentNode *)node;
            for (auto c : p->children) {
                foundIntersection = octree_traversal(shadowRay, ray, its, f, c) || foundIntersection;
                if (foundIntersection && shadowRay) {
                    return foundIntersection;
                }
            }
        }
        else {
            ChildNode* c = (ChildNode*)node;
            for (size_t i = 0; i < m_mesh.size(); i++)
            {
                for (auto triangle : c->triangles[i]) {
                    float u, v, t;
                    if (m_mesh[i]->rayIntersect(triangle, ray, u, v, t)) {
                        if (shadowRay) {
                            return true;
                        }
                        ray.maxt = its.t = t;
                        its.uv = Point2f(u, v);
                        its.mesh = m_mesh[i];
                        f = triangle;
                        foundIntersection = true;
                    }
                }
            }
        }
    }
    
    return foundIntersection;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    foundIntersection = octree_traversal(shadowRay, ray, its, f, m_root);
    
    if (foundIntersection) {
        count++;
    }

    if (foundIntersection && shadowRay) {
        return foundIntersection;
    }

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

