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
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

Accel::Node* Accel::build(BoundingBox3f box, std::vector<uint32_t> triangles, uint32_t depth) {
    if (triangles.size() == 0)
        return nullptr;

    if (triangles.size() <= 10 || depth > 10) {
        /*if (depth > 500) {
            cout << "depth problem" << endl;
            cout << "number of triangles:" << endl;
            cout << triangles.size() << endl;
        }*/
        Accel::ChildNode* leaf = new Accel::ChildNode();
        leaf->type = 1;
        leaf->box = box;
        leaf->triangles = triangles;
        for (size_t i = 0; i < triangles.size(); i++)
        {
            count.insert(triangles[i]);
        }
        return leaf;
    }

    std::vector<uint32_t> list[8];
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
        /*cout << (Point3f(extent[0] / 2, extent[1] / 2, extent[2] / 2).cwiseProduct(dividors[i])) << endl;
        cout << "+++++++++++++" << endl;*/
        auto newMax = box.max - (Point3f(extent[0] / 2, extent[1] / 2, extent[2] / 2).cwiseProduct(dividors[i]));
        auto newMin = newMax - (extent / 2);
        BoundingBox3f i_bbox(newMin, newMax);
        sub_bbox.push_back(i_bbox);
    }

    for (auto triangle : triangles) {
        BoundingBox3f t_bbox = m_mesh->getBoundingBox(triangle);
        for (int i = 0; i < 8; i++) {
            //cout << "enter" << endl;
            if (sub_bbox[i].overlaps(t_bbox, false)) {
                list[i].push_back(triangle);
            }
        }
        //cout << "exit" << endl;
    }

    Accel::ParentNode* node = new ParentNode();
    node->box = box;
    node->type = 0;
    //node->children = new ParentNode;
    for (int i = 0; i < 8; ++i) {
        //ParentNode* p = (ParentNode*)node->children;
        //p->children[i] = build(sub_bbox[i], list[i]);
        node->children[i] = build(sub_bbox[i], list[i], depth+1);
    }

    return node;
}


void Accel::build() {
    /* Nothing to do here for now */
    uint32_t size = m_mesh->getTriangleCount();
    //count = 0;
    std::vector<uint32_t> v(size);
    for (int i = 0; i < size; i++) {
        v[i] = i;
    }

    m_root = build(m_bbox, v, 0);
    std::cout << "hello" << std::endl;
}

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
            }
        }
        else {
            ChildNode* c = (ChildNode*)node;
            for (auto triangle : c->triangles) {
                float u, v, t;
                if (m_mesh->rayIntersect(triangle, ray, u, v, t)) {
                    if (shadowRay) {
                        return true;
                    }
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = triangle;
                    foundIntersection = true;
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

    /* Brute force search through all triangles */
    //for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //    float u, v, t;
    //    if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //        /* An intersection was found! Can terminate
    //           immediately if this is a shadow ray query */
    //        if (shadowRay)
    //            return true;
    //        ray.maxt = its.t = t;
    //        its.uv = Point2f(u, v);
    //        its.mesh = m_mesh;
    //        f = idx;
    //        foundIntersection = true;
    //    }
    //}
    foundIntersection = octree_traversal(shadowRay, ray, its, f, m_root);


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

