#include <Eigen/Core>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#pragma once

/// Null node flag.
constexpr uint32_t NULL_NODE = std::numeric_limits<uint32_t>::max();


// Axis-aligned bounding box
class AABB
{
  public:
   AABB();
   AABB(const Eigen::Vector3d& lower, const Eigen::Vector3d& upper);

   // Compute the surface area of the box.
   double computeSurfaceArea() const;

   // Get the surface area of the box.
   const double surfaceArea() const { return m_surfaceArea; }
   double& surfaceArea() { return m_surfaceArea; }

   // Merge two AABBs into this one.
   void merge(const AABB&, const AABB&);

   // Test whether the AABB is contained within this one.
   bool contains(const AABB&) const;

   // Test whether the AABB overlaps this one.
   bool overlaps(const AABB&, bool touchIsOverlap) const;

   // Compute the center of the AABB.
   const Eigen::Vector3d computeCenter() const;
   const Eigen::Vector3d center() const { return m_center; }
   Eigen::Vector3d& center() { return m_center; }

   // Lower bound of AABB
   Eigen::Vector3d lowerBound;

   // Upper bound of AABB
   Eigen::Vector3d upperBound;

  private:
   // The position of the AABB center.
   Eigen::Vector3d m_center;

   /// The AABB's surface area.
   double m_surfaceArea;
};


// AABB Tree Node
//
// The AABB objects are "fattened" before they are stored to avoid having to
// continually update and rebalance the tree when displacements are small.
//
// Nodes are aware of their position within in the tree. The isLeaf member
// function allows the tree to query whether the node is a leaf, i.e. to
// determine whether it holds a single object.
//
class AABBNode
{
  public:
   // Constructor.
   AABBNode();

   // The fattened axis-aligned bounding box.
   AABB aabb;

   // Index of the parent node.
   uint32_t parent;

   // Index of the next node.
   uint32_t next;

   // Index of the left-hand child.
   uint32_t left;

   // Index of the right-hand child.
   uint32_t right;

   // Height of the node. This is 0 for a leaf and -1 for a free node.
   int height;

   // The index of the original object that the node contains (leaf nodes only).
   // Used for indexing into some external array; not used internall
   uint32_t origin;

   // Test whether the node is a leaf.
   bool isLeaf() const { return (left == NULL_NODE); }
};


// Dynamic AABB tree.
//
// The dynamic AABB tree is a hierarchical data structure that can be used to
// efficiently query overlaps between objects of arbitrary shape and size that
// lie inside of a simulation box.
//
class AABBTree
{
  public:
   // Constructor
   // skinThickness - The skin thickness for fattened AABBs, as a fraction
   //                 of the AABB base length.
   // numObjects - The number of objects this AABB tree will contain
   // touchIsOverlap - Does touching count as overlapping in query operations?
   AABBTree(double skinThickness = 0.05, uint32_t numObjects = 16, bool touchIsOverlap=true);

   // Insert a sphere into the tree
   void addSphere(uint32_t origin, const Eigen::Vector3d& position, double radius);

   // Insert a bounding box into the tree
   void insertAABB(uint32_t origin, const Eigen::Vector3d& lowerBound,
                   const Eigen::Vector3d& upperBound);
   void insertAABB(uint32_t, const AABB& aabb);

   // Return the number of objects stored in the tree.
   const uint32_t numObjects() const;

   //! Remove an object from the tree.
   void removeObject(uint32_t index);

   /// Remove all objects from the tree.
   void removeAll();

   // Update the tree if an object moves
   // alwaysReinsert - Always reinsert the object, even if it's within its
   //   old AABB (default:false)
   // Returns true if the object was reinserted.
   bool updateSphere(uint32_t index, const Eigen::Vector3d& position, double radius,
                     bool alwaysReinsert=false);

   // Update the tree if an object moves
   // alwaysReinsert - Always reinsert the object, even if it's within its old
   //    AABB (default: false)

   bool updateAABB(uint32_t index,  const Eigen::Vector3d& lowerBound,
                   const Eigen::Vector3d& upperBound, bool alwaysReinsert=false);

   // Query the tree to find candidate intersections for an object (i.e. does
   // the AABB of this object intersect any of its neighbors?
   std::vector<uint32_t> query(uint32_t index) const;

   // Query the tree to find candidate interactions for a given AABB, ignoring
   // intersections with object 'index'
   std::vector<uint32_t> query(uint32_t index, const AABB& aabb) const;

   // Query the tree to find candidate interactions for an AABB.
   std::vector<uint32_t> query(const AABB& aabb) const;

   // Intersect this AABB with a second, returning a list of pairs of
   // box intersections
   std::vector< std::pair<uint32_t, uint32_t> > intersect(const AABBTree& tree);

   // Get the object's AABB.
   const AABB& getAABB(uint32_t index) const;

   // Get the height of the tree.
   uint32_t getHeight() const;

   // Get the number of nodes in the tree.
   uint32_t getNodeCount() const;

   // Compute the maximum balancance (the maximum difference between the height
   // of two children of a node) of the tree
   uint32_t computeMaximumBalance() const;

   // Compute the surface area ratio (the ratio of the sum of the node surface
   // area to the surface area of the root node) of the tree.
   double computeSurfaceAreaRatio() const;

   // Validate the tree.
   void validate() const;

   // Rebuild an optimal tree.
   void rebuild();

  private:
   // The index of the root node.
   uint32_t root;

   // The dynamic tree of AABBNodes
   std::vector<AABBNode> nodes;

   // The current number of nodes in the tree.
   uint32_t nodeCount;

   // The current node capacity.
   uint32_t nodeCapacity;

   // The position of node at the top of the free list.
   uint32_t freeList;

   // The skin thickness of the fattened AABBs, as a fraction of the AABB base length.
   double skinThickness;

   // A map between an external object list and node indices.
   std::unordered_map<uint32_t, uint32_t> externalMap;

   // Does touching count as overlapping in tree queries?
   bool touchIsOverlap;

   // Allocate a new node. Returns the index of the allocated node.
   uint32_t allocateNode();

   // Free an existing node.
   void freeNode(uint32_t index);

   // Insert a leaf into the tree.
   void insertLeaf(uint32_t index);

   // Remove a leaf from the tree.
   void removeLeaf(uint32_t index);

   // Balance the tree.
   uint32_t balance(uint32_t index);

   // Compute the height of the tree.
   uint32_t computeHeight() const;

   // Compute the height of a sub-tree, using 'index' as the sub-tree root node
   uint32_t computeHeight(uint32_t index) const;

   // Assert that a sub-tree (with root 'root_index') has a valid structure.
   void validateStructure(uint32_t root_index) const;

   // Assert that the sub-tree (with root 'root_index') has valid metrics.
   void validateMetrics(uint32_t root_index) const;
};

std::vector< std::pair<uint32_t, uint32_t> > intersectAABBTrees(const AABBTree& treeA,
                                                                const AABBTree& treeB);


