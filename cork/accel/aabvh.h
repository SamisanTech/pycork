// +-------------------------------------------------------------------------
// | aabvh.h  (optimized: templated queries, no std::function, flat stack)
// +-------------------------------------------------------------------------
#ifndef CORK_AABVH_H_HEADER_HAS_BEEN_INCLUDED
#define CORK_AABVH_H_HEADER_HAS_BEEN_INCLUDED

#include <vector>
#include <algorithm>

#include <cork/math/bbox.h>
#include <cork/util/shortVec.h>
#include <cork/util/iterPool.h>

// larger leaves = fewer nodes, better for dense meshes / SIMD-ish traversal
static const uint LEAF_SIZE = 16;

template<class GeomIdx>
struct GeomBlob
{
    BBox3d  bbox;
    Vec3d   point; // representative point, usually the box midpoint
    GeomIdx id;
};

template<class GeomIdx>
struct AABVHNode
{
    BBox3d                          bbox;
    AABVHNode                       *left;
    AABVHNode                       *right;
    ShortVec<uint, LEAF_SIZE>       blobids;
    inline bool isLeaf() const { return left == nullptr; }
};

template<class GeomIdx>
class AABVH
{
public:
    AABVH(const std::vector< GeomBlob<GeomIdx> > &geoms) : root(nullptr),
                                                           blobs(geoms),
                                                           tmpids(geoms.size())
    {
        ENSURE(blobs.size() > 0);

        for(uint k=0; k<tmpids.size(); k++)
            tmpids[k] = k;

        root = constructTree(0, tmpids.size(), 2);
    }

    ~AABVH() {}

    // Hot path: templated action — no std::function / heap alloc
    template<class Action>
    inline void for_each_in_box(const BBox3d &bbox, Action &&action) const {
        AABVHNode<GeomIdx>* stack[128];
        int sp = 0;
        stack[sp++] = root;

        while(sp > 0) {
            AABVHNode<GeomIdx> *node = stack[--sp];

            if(!hasIsct(node->bbox, bbox))
                continue;

            if(node->isLeaf()) {
                for(uint bid : node->blobids) {
                    if(hasIsct(bbox, blobs[bid].bbox))
                        action(blobs[bid].id);
                }
            } else {
                // push farther child first (optional); both needed
                if(sp + 2 <= 128) {
                    stack[sp++] = node->right;
                    stack[sp++] = node->left;
                }
            }
        }
    }

private:
    AABVHNode<GeomIdx>* constructTree(uint begin, uint end, uint last_dim)
    {
        ENSURE(end - begin > 0);
        if(end-begin <= LEAF_SIZE) {
            AABVHNode<GeomIdx> *node = node_pool.alloc();
            node->left = nullptr;
            node->right = nullptr;
            node->blobids.resize(end-begin);
            for(uint k=0; k<end-begin; k++) {
                uint blobid = node->blobids[k] = tmpids[begin + k];
                node->bbox = convex(node->bbox, blobs[blobid].bbox);
            }
            return node;
        }

        // split along longest axis of aggregate bbox of points
        BBox3d extent;
        for(uint i=begin; i<end; i++)
            extent = convex(extent, blobs[tmpids[i]].bbox);
        Vec3d diag = extent.maxp - extent.minp;
        uint dim = 0;
        if(diag[1] > diag[dim]) dim = 1;
        if(diag[2] > diag[dim]) dim = 2;
        (void)last_dim;

        uint mid = (begin + end) / 2;
        quickSelect(mid, begin, end, dim);

        AABVHNode<GeomIdx> *node = node_pool.alloc();
        node->left = constructTree(begin, mid, dim);
        node->right = constructTree(mid, end, dim);
        node->bbox = convex(node->left->bbox, node->right->bbox);
        return node;
    }

    void quickSelect(uint select, uint begin, uint end, uint dim)
    {
        while(true) {
            if(end-1 == select)
                return;

            uint pi = randMod(end-begin) + begin;
            double pv = blobs[tmpids[pi]].point[dim];

            int front = (int)begin;
            int back  = (int)end-1;
            while(front < back) {
                if(blobs[tmpids[front]].point[dim] < pv) {
                    front++;
                }
                else if(blobs[tmpids[back]].point[dim] > pv) {
                    back--;
                }
                else {
                    std::swap(tmpids[front], tmpids[back]);
                    front++;
                    back--;
                }
            }
            if(front == back && blobs[tmpids[front]].point[dim] <= pv) {
                front++;
            }

            if(select < uint(front)) {
                end = (uint)front;
            } else {
                begin = (uint)front;
            }
        }
    }

private:
    AABVHNode<GeomIdx>                  *root;

    IterPool< AABVHNode<GeomIdx> >      node_pool;
    std::vector< GeomBlob<GeomIdx> >    blobs;
    std::vector<uint>                   tmpids;
};

#endif
