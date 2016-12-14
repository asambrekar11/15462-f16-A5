#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CMU462 { namespace StaticScene {
  

  void BVHAccel::Swap(Primitive* &a, Primitive* &b)
  {
    Primitive* temp = a;
    a = b;
    b = temp;
  }

  void BVHAccel::PrintAxis(Axis &currAxis) const
  {
    switch (currAxis)
    {
      case XAXIS:
          printf("X Axis\n");
          break;
      case YAXIS:
          printf("Y Axis\n");
          break;
      case ZAXIS:
          printf("Z Axis\n");
          break;
    }
  }

  bool BVHAccel::IsLessThanEquals(Axis &currAxis, Primitive * &left, Primitive * &right) const
  {
    bool retValue;
    switch (currAxis)
    {
      case XAXIS:
          retValue = left->get_bbox().centroid().x <= right->get_bbox().centroid().x;
          break;
      case YAXIS:
          retValue = left->get_bbox().centroid().y <= right->get_bbox().centroid().y;
          break;
      case ZAXIS:
          retValue = left->get_bbox().centroid().z <= right->get_bbox().centroid().z;
          break;
    }
    return retValue;
  }

  size_t BVHAccel::Partition(Axis &currAxis, size_t start, size_t end)
  {
    // Primitive * pivot = primitives[end];
    // size_t i = start-1;
    // for (size_t j  = start; j < end; j++)
    // {
    //   if (IsLessThan(currAxis,primitives[j],pivot))
    //   {
    //     i++;
    //     Swap(primitives[i], primitives[j]);
    //   }
    // }
    // Swap(primitives[i+1],primitives[end]);
    // return (i+1);

    size_t mid = start + (end-start)/2;
    Primitive * pivot = primitives[mid];
    Swap(primitives[mid],primitives[start]);
    size_t i = start + 1;
    size_t j = end;
    while (i <= j)
    {
      while (i <= j && IsLessThanEquals(currAxis,primitives[i],pivot))
      {
        i++;
      }
      while (i <= j && !IsLessThanEquals(currAxis,primitives[j],pivot))
      {
        j--;
      }
      if (i<j)
      {
        Swap(primitives[i],primitives[j]);
      }
    }
    Swap(primitives[i-1],primitives[start]);
    return i-1;

  }

  void BVHAccel::QuickSort(Axis &currAxis, size_t start, size_t end)
  {
    
    // printf("start = %d, end = %d\n",(int)start,(int)end);
    if ((long long int)start < (long long int)end)
    {
      // printf("Begin Partition with start = %d, end = %d\n",(int)start,(int)end);
      size_t p = Partition(currAxis,start,end);
      QuickSort(currAxis,start,p-1);
      QuickSort(currAxis,p+1,end);
    }
  }

  void BVHAccel::SplitNode(BVHNode* currNode, size_t max_leaf_size)
  {
    if (currNode->isLeaf())
    {
      // Number of buckets
      const int nBuckets = 8;
      Bucket globalBest(XAXIS);

      // For each axis X ->0, Y->1, Z->2 
      for (int axis = 0; axis <3; axis++)
      {
        // Initialize buckets
        Bucket buckets[nBuckets];
        Axis currentAxis = (Axis)axis;
        
        // PrintAxis(currentAxis);
        Bucket axialBest(currentAxis);
        for (int i = 0; i < nBuckets; i++)
        {
          buckets[i].myAxis = currentAxis;
        }
        // printf("Buckets initialized\n");

        // Put primitives in buckets
        for (int i = currNode->start; i<(currNode->start+currNode->range); i++)
        {
          Vector3D centroid = primitives[i]->get_bbox().centroid();
          int currBucket;
          double factor = double(nBuckets);
          switch (currentAxis)
          {
            case XAXIS:
                currBucket = floor(factor * (centroid.x - currNode->bb.min.x) / currNode->bb.extent.x);
                break;
            case YAXIS:
                currBucket = floor(factor * (centroid.y - currNode->bb.min.y) / currNode->bb.extent.y);
                break;
            case ZAXIS:
                currBucket = floor(factor * (centroid.z - currNode->bb.min.z) / currNode->bb.extent.z);
                break;
          }
          if (currBucket >= nBuckets)
          {
            currBucket = nBuckets-1;
          }
          // printf("Current Bucket = %d\n",currBucket);
          for (int j = 0; j < nBuckets; j++)
          {
            if (j < currBucket)
            {
              buckets[j].r_bb.expand(primitives[i]->get_bbox());
              buckets[j].r_count++;
            }
            else
            {
              buckets[j].l_bb.expand(primitives[i]->get_bbox());
              buckets[j].l_count++;
            }
          }
        }
        // printf("Primitives are in buckets!\n");

        // Calculate cost of each bucket
        for (int i = 0; i<nBuckets; i++)
        {
          double currentCost = buckets[i].UpdateCost(currNode->bb.surface_area());
          // printf("Cost = %lf\n",currentCost);
          if (buckets[i].l_count != 0 && buckets[i].r_count != 0 && (axialBest.cost == INF_D || axialBest.cost >= currentCost))
          {
            axialBest = buckets[i];
          }
        }
        // printf("Cost calculated\n");

        if (globalBest.cost == INF_D || globalBest.cost >= axialBest.cost)
        {
          globalBest = axialBest;
          // printf("New global Best! l_count = %d, r_count = %d\n", (int)globalBest.l_count,(int)globalBest.r_count);
        }
      }

      // Sort along global best
      // printf("Start = %d, Range = %d\n",(int)currNode->start, (int)currNode->range);
      QuickSort(globalBest.myAxis, currNode->start, (currNode->start + currNode->range-1));
      // printf("Sorted!\n");

      // Initialize l node
      currNode->l = new BVHNode(globalBest.l_bb,currNode->start,globalBest.l_count);

      // Split l node
      if (currNode->l->range > max_leaf_size)
      {
        SplitNode(currNode->l, max_leaf_size);
      }

      // Initialize r node
      currNode->r = new BVHNode(globalBest.r_bb,currNode->start + globalBest.l_count, globalBest.r_count);

      // Split r node
      if (currNode->r->range > max_leaf_size)
      {
        SplitNode(currNode->r, max_leaf_size);
      }
    }
  }

  void BVHAccel::DestroyNode(BVHNode* currNode)
  {
    if (currNode != NULL)
    {
      if (currNode->isLeaf())
      {
        return;
      }
      if (currNode->l != NULL)
      {
        DestroyNode(currNode->l);
      }
      if (currNode->r != NULL)
      {
        DestroyNode(currNode->r);
      }
      delete currNode;
      currNode = NULL;
    }
  }

  BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
      size_t max_leaf_size) {

    this->primitives = _primitives;

    // TODO:
    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration. The starter code build a BVH aggregate with a
    // single leaf node (which is also the root) that encloses all the
    // primitives.

    BBox bb;
    for (size_t i = 0; i < primitives.size(); ++i) {
      bb.expand(primitives[i]->get_bbox());
    }    
    root = new BVHNode(bb, 0, primitives.size());
    // cout<<bb.min<<endl<<bb.max<<endl;
    if (primitives.size() > max_leaf_size)
      SplitNode(root,max_leaf_size);

  }

  BVHAccel::~BVHAccel() {

    // TODO:
    // Implement a proper destructor for your BVH accelerator aggregate
    DestroyNode(root);

  }

  BBox BVHAccel::get_bbox() const {
    return root->bb;
  }

  bool BVHAccel::intersectWithNodePrimitives(const Ray &ray, BVHNode* currNode) const
  {
    bool hit = false;
    if (currNode->isLeaf())
    {
      for (int j = currNode->start; j < (currNode->start + currNode->range); j++)
      {
        if (primitives[j]->intersect(ray))
        {
          hit = true;
        }
      }
    }
    return hit;
  }

  bool BVHAccel::intersectWithNodePrimitives(const Ray &ray, BVHNode* currNode, Intersection* i) const
  {
    bool hit = false;
    double min_t = INF_D;
    Intersection temp;
    // if (i == NULL)
    // {
    //   i = new Intersection();
    // }
    
    if (currNode->isLeaf())
    {
      for (int j = currNode->start; j < (currNode->start + currNode->range); j++)
      {
        if (primitives[j]->intersect(ray,&temp))
        {
          // printf("Yay! Found an intersecting primitive!\n");
          // printf("start = %d, range = %d, j = %d, t = %lf\n",(int)currNode->start, (int)currNode->range, j, temp.t);
          hit = true;
          
          // printf("Checking min_t\n");
          if (min_t == INF_D || min_t >= temp.t)
          {
            min_t = temp.t;
            *i = temp;
          }
        }
      }
    }
    // printf("Getting out of primitives\n");
    return hit;
  }

  bool BVHAccel::intersectWithNode(const Ray &ray, BVHNode *currNode) const
  {
    double ray_mint = 0, ray_maxt = INF_D;
    // if (!currNode->bb.intersect(ray,ray_mint,ray_maxt))
    // {
    //   // printf("Not in bbox\n");
    //   return false;
    // }
    // ray.min_t = ray_mint;
    // ray.max_t = ray_maxt;
    if (currNode->isLeaf())
    {
      return intersectWithNodePrimitives(ray,currNode);
    }
    bool hit = false;
    // double min_t1 = ray.min_t, max_t1 = ray.max_t, min_t2 = ray.min_t, max_t2 = ray.max_t;
    double min_t1 = 0.0, max_t1 = INF_D, min_t2 = 0.0, max_t2 = INF_D;

    bool hit1 = currNode->l->bb.intersect(ray,min_t1,max_t1);
    bool hit2 = currNode->r->bb.intersect(ray,min_t2,max_t2);
    // printf("min_t1 = %lf; max_t1 = %lf; min_t2 = %lf; max_t2 = %lf\n",min_t1, max_t1, min_t2, max_t2);

    if (hit1 == true && hit2 == true)
    {
      BVHNode* first = min_t1 <= min_t2 ? currNode->l : currNode->r;
      BVHNode* second = min_t1 <= min_t2 ? currNode->r : currNode->l;
      // printf("Checking for first\n");
      // ray.min_t = min(min_t1,min_t2);
      // ray.max_t = ray.min_t == min_t1?max_t1:max_t2;
      if (intersectWithNode(ray,first))
      {
        hit = true;
        // printf("First child intersection!\n");
      }
      else
      {
        // printf("Checking for second\n");
        // ray.min_t = max(min_t1,min_t2);
        // ray.max_t = ray.min_t == min_t1?max_t1:max_t2;
        if (intersectWithNode(ray,second))
        {
          hit  = true;
          // printf("Second child intersection!\n");
        }
      }
    }
    else if (hit1)
    {
      // ray.min_t = min_t1;
      // ray.max_t = max_t1;
      if (intersectWithNode(ray,currNode->l))
      {
        hit = true;
      }
    }
    else if (hit2)
    {
      // ray.min_t = min_t2;
      // ray.max_t = max_t2;
      if (intersectWithNode(ray,currNode->r))
      {
        hit = true;
      }
    }
    return hit;
  }

  bool BVHAccel::intersectWithNode(const Ray &ray, BVHNode *currNode, Intersection *i) const
  {
    // double ray_mint = 0.0, ray_maxt = INF_D;
    // if (!currNode->bb.intersect(ray,ray_mint,ray_maxt))
    // {
    //   return false;
    // }
    // ray.min_t = ray_mint;
    // ray.max_t = ray_maxt;
    
    if (currNode->isLeaf())
    {
      return intersectWithNodePrimitives(ray,currNode,i);
    }
    // printf("Not a leaf\n");
    bool hit = false;
    double min_t1 = ray.min_t, max_t1 = ray.max_t, min_t2 = ray.min_t, max_t2 = ray.max_t;
    // double min_t1 = 0.0, max_t1 = INF_D, min_t2 = 0.0, max_t2 = INF_D;
    bool hit1 = currNode->l->bb.intersect(ray,min_t1,max_t1);
    bool hit2 = currNode->r->bb.intersect(ray,min_t2,max_t2);
    // printf("min_t1 = %lf; max_t1 = %lf; min_t2 = %lf; max_t2 = %lf\n",min_t1, max_t1, min_t2, max_t2);
      
    if (hit1 == true && hit2 == true)
    {
      BVHNode* first = min_t1 <= min_t2 ? currNode->l : currNode->r;
      BVHNode* second = min_t1 <= min_t2 ? currNode->r : currNode->l;
      // ray.min_t = min(min_t1,min_t2);
      // ray.max_t = ray.min_t == min_t1?max_t1:max_t2;
      // printf("Checking for first\n");
      if (intersectWithNode(ray,first,i))
      {
        hit = true;
        // printf("intersection t = %lf\n",i->t);
        if (i->t >= max(min_t1,min_t2))
        {
            // printf("Checking for second\n");
            // printf("old t = %lf\n",i->t);
            // ray.min_t = max(min_t1,min_t2);
            // ray.max_t = ray.min_t == min_t1?max_t1:max_t2;
            Intersection temp;
            if (intersectWithNode(ray,second,&temp))
            {
              // printf("new t = %lf\n",temp.t);
              if (temp.t < i->t)
              {
                *i = temp;
                hit  = true;
              }
            }
        }
      }
      else
      {
        // ray.min_t = max(min_t1,min_t2);
        // ray.max_t = ray.min_t == min_t1?max_t1:max_t2;
        if (intersectWithNode(ray,second,i))
        {
          hit  = true;
        }
      }
    }
    else if (hit1)
    {
      // ray.min_t = min_t1;
      // ray.max_t = max_t1;
      if (intersectWithNode(ray,currNode->l,i))
      {
        hit = true;
        // printf("First child intersection\n");
      }
    }
    else if (hit2)
    {
      // ray.min_t = min_t2;
      // ray.max_t = max_t2;
      if (intersectWithNode(ray,currNode->r,i))
      {
        hit = true;
        // printf("Second child intersection\n");
      }
    }
    return hit;
  }

  bool BVHAccel::intersect(const Ray &ray) const {

    // TODO:
    // Implement ray - bvh aggregate intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.
    // printf("Without intersection\n");
    // printf("Check shadow ray intersection\n");
    if (root->bb.intersect(ray,ray.min_t,ray.max_t))
      return intersectWithNode(ray, root);
    return false;

  }

  bool BVHAccel::intersect(const Ray &ray, Intersection *i) const {

    // TODO:
    // Implement ray - bvh aggregate intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate. When an intersection does happen.
    // You should store the non-aggregate primitive in the intersection data
    // and not the BVH aggregate itself.
    // printf("With intersection\n");
    ray.min_t = 0.0; ray.max_t = INF_D;
    if (root->bb.intersect(ray,ray.min_t,ray.max_t))
      return intersectWithNode(ray, root,i);
    return false;
  }

}  // namespace StaticScene
}  // namespace CMU462
