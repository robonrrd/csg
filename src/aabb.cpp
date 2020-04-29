// GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999

#include "aabb.h"


AABB::AABB()
{
}

AABB::AABB(const Eigen::Vector3d& lower, const Eigen::Vector3d& upper)
    : lowerBound(lower)
    , upperBound(upper)
{
    // Validate that the upper bounds exceed the lower bounds.
    for (uint32_t ii=0; ii<3; ++ii)
        if (lowerBound[ii] > upperBound[ii])
            throw std::invalid_argument("AABB lower bound is greater than the upper bound");

    m_surfaceArea = computeSurfaceArea();
    m_center = computeCenter();
}

double
AABB::computeSurfaceArea() const
{
    // Sum of area of all sides.
    double sum = 0.0;

    // General formula for one side: hold one dimension constant and multiply by
    // all the other ones.
    for (uint32_t d1=0; d1<3; ++d1)
    {
        // Area of current side.
        double product = 1.0;
        for (uint32_t d2=0; d2<3; ++d2)
        {
            if (d1 == d2)
                continue;

            double dx = upperBound[d2] - lowerBound[d2];
            product *= dx;
        }
        sum += product;
    }

    return 2.0*sum;
}

void AABB::merge(const AABB& aabb1, const AABB& aabb2)
{
    // TODO: This can be done with component-wise operations in Eigen
    for (uint32_t ii=0; ii<3; ++ii)
    {
        lowerBound[ii] = std::min(aabb1.lowerBound[ii], aabb2.lowerBound[ii]);
        upperBound[ii] = std::max(aabb1.upperBound[ii], aabb2.upperBound[ii]);
    }

    m_surfaceArea = computeSurfaceArea();
    m_center = computeCenter();
}

bool AABB::contains(const AABB& aabb) const
{
    for (uint32_t ii=0; ii<3; ++ii)
    {
        if (aabb.lowerBound[ii] < lowerBound[ii])
            return false;
        if (aabb.upperBound[ii] > upperBound[ii])
            return false;
    }

    return true;
}

bool AABB::overlaps(const AABB& aabb, bool touchIsOverlap=true) const
{
    if (touchIsOverlap)
    {
        for (uint32_t ii=0; ii<3; ++ii)
            if (aabb.upperBound[ii] < lowerBound[ii] || aabb.lowerBound[ii] > upperBound[ii])
                return false;
    }
    else
    {
        for (uint32_t ii=0; ii<3; ++ii)
            if (aabb.upperBound[ii] <= lowerBound[ii] || aabb.lowerBound[ii] >= upperBound[ii])
                return false;
    }

    return true;
}

const Eigen::Vector3d AABB::computeCenter() const
{
    return 0.5*(lowerBound + upperBound);
}



AABBNode::AABBNode()
    : parent(NULL_NODE)
    , next(NULL_NODE)
    , left(NULL_NODE)
    , right(NULL_NODE)
    , height(0)
{
    // empty
}



AABBTree::AABBTree(double skin_thickness, uint32_t numObjects, bool touch_is_overlap)
    : skinThickness(skin_thickness)
    , touchIsOverlap(touch_is_overlap)
{
    // Initialise the tree.
    root = NULL_NODE;
    nodeCount = 0;
    nodeCapacity = numObjects;
    nodes.resize(nodeCapacity);

    // Build a linked list for the list of free nodes.
    for (uint32_t ii=0; ii<nodeCapacity-1; ++ii)
    {
        nodes[ii].next = ii + 1;
        nodes[ii].height = -1;
    }
    nodes[nodeCapacity-1].next = NULL_NODE;
    nodes[nodeCapacity-1].height = -1;

    // Assign the index of the first free node.
    freeList = 0;
}

uint32_t AABBTree::allocateNode()
{
    // Exand the node pool as needed.
    if (freeList == NULL_NODE)
    {
        assert(nodeCount == nodeCapacity);

        // The free list is empty. Rebuild a bigger pool.
        nodeCapacity *= 2;
        nodes.resize(nodeCapacity);

        // Build a linked list for the list of free nodes.
        for (uint32_t i=nodeCount;i<nodeCapacity-1;i++)
        {
            nodes[i].next = i + 1;
            nodes[i].height = -1;
        }
        nodes[nodeCapacity-1].next = NULL_NODE;
        nodes[nodeCapacity-1].height = -1;

        // Assign the index of the first free node.
        freeList = nodeCount;
    }

    // Peel a node off the free list.
    uint32_t node = freeList;
    freeList = nodes[node].next;
    nodes[node].parent = NULL_NODE;
    nodes[node].left = NULL_NODE;
    nodes[node].right = NULL_NODE;
    nodes[node].height = 0;
    nodeCount++;

    return node;
}

void AABBTree::freeNode(uint32_t node)
{
    assert(node < nodeCapacity);
    assert(0 < nodeCount);

    nodes[node].next = freeList;
    nodes[node].height = -1;
    freeList = node;
    nodeCount--;
}

void AABBTree::addSphere(uint32_t origin, const Eigen::Vector3d& position, double radius)
{
    // Make sure this object doesn't already exist.
    if (externalMap.count(origin) != 0)
    {
        throw std::invalid_argument("[ERROR]: Object already exists in tree!");
    }

    // Allocate a new node for the object.
    uint32_t node = allocateNode();

    // AABB size
    Eigen::Vector3d size;

    // Compute the AABB limits.
    for (uint32_t ii=0; ii<3; ii++)
    {
        nodes[node].aabb.lowerBound[ii] = position[ii] - radius;
        nodes[node].aabb.upperBound[ii] = position[ii] + radius;
        size[ii] = nodes[node].aabb.upperBound[ii] - nodes[node].aabb.lowerBound[ii];
    }

    // Fatten the AABB.
    for (uint32_t ii=0; ii<3; ii++)
    {
        nodes[node].aabb.lowerBound[ii] -= skinThickness * size[ii];
        nodes[node].aabb.upperBound[ii] += skinThickness * size[ii];
    }
    nodes[node].aabb.surfaceArea(nodes[node].aabb.computeSurfaceArea());
    nodes[node].aabb.center(nodes[node].aabb.computeCenter());

    // Zero the height.
    nodes[node].height = 0;

    // Insert a new leaf into the tree.
    insertLeaf(node);

    // Add the new object to the map.
    externalMap.insert(std::unordered_map<uint32_t, uint32_t>::value_type(origin, node));

    // Store the object index.
    nodes[node].origin = origin;
}

void AABBTree::insertAABB(uint32_t index, const AABB& aabb)
{
    insertAABB(index, aabb.lowerBound, aabb.upperBound);
}


void AABBTree::insertAABB(uint32_t index, const Eigen::Vector3d& lowerBound,
                          const  Eigen::Vector3d& upperBound)
{
    // Make sure the object doesn't already exist.
    if (externalMap.count(index) != 0)
    {
        throw std::invalid_argument("[ERROR]: Object already exists in tree!");
    }

    // Allocate a new node for the object.
    uint32_t node = allocateNode();

    // AABB size
    Eigen::Vector3d size;

    // Compute the AABB limits.
    for (uint32_t ii=0; ii<3; ii++)
    {
        // Validate the bound.
        if (lowerBound[ii] > upperBound[ii])
        {
            throw std::invalid_argument("AABB lower bound is greater than the upper bound!");
        }

        nodes[node].aabb.lowerBound[ii] = lowerBound[ii];
        nodes[node].aabb.upperBound[ii] = upperBound[ii];
        size[ii] = upperBound[ii] - lowerBound[ii];
    }

    // Fatten the AABB.
    for (uint32_t ii=0; ii<3; ii++)
    {
        nodes[node].aabb.lowerBound[ii] -= skinThickness * size[ii];
        nodes[node].aabb.upperBound[ii] += skinThickness * size[ii];
    }
    nodes[node].aabb.surfaceArea(nodes[node].aabb.computeSurfaceArea());
    nodes[node].aabb.center(nodes[node].aabb.computeCenter());

    // Zero the height.
    nodes[node].height = 0;

    // Insert a new leaf into the tree.
    insertLeaf(node);

    // Add the new object to the map.
    externalMap.insert(std::unordered_map<uint32_t, uint32_t>::value_type(index, node));

    // Store the object index.
    nodes[node].origin = index;
}

const uint32_t AABBTree::numObjects() const
{
    return externalMap.size();
}

void AABBTree::removeObject(uint32_t index)
{
    // Map iterator.
    std::unordered_map<uint32_t, uint32_t>::iterator it;

    // Find the original object
    it = externalMap.find(index);
    if (it == externalMap.end())
    {
        throw std::invalid_argument("[ERROR]: Invalid particle index!");
    }

    // Extract the node index.
    uint32_t node = it->second;

    // Erase the node from the map.
    externalMap.erase(it);

    assert(node < nodeCapacity);
    assert(nodes[node].isLeaf());

    removeLeaf(node);
    freeNode(node);
}

void AABBTree::removeAll()
{
    // Iterator pointing to the start of the map.
    std::unordered_map<uint32_t, uint32_t>::iterator it = externalMap.begin();

    // Iterate over the map.
    while (it != externalMap.end())
    {
        // Extract the node index.
        uint32_t node = it->second;

        assert(node < nodeCapacity);
        assert(nodes[node].isLeaf());

        removeLeaf(node);
        freeNode(node);

        ++it;
    }

    // Clear the map.
    externalMap.clear();
}

bool AABBTree::updateSphere(uint32_t index, const Eigen::Vector3d& position, double radius,
                            bool alwaysReinsert)
{
    // AABB bounds vectors.
    Eigen::Vector3d lowerBound;
    Eigen::Vector3d upperBound;

    // Compute the AABB limits.
    for (uint32_t ii=0; ii<3; ii++)
    {
        lowerBound[ii] = position[ii] - radius;
        upperBound[ii] = position[ii] + radius;
    }

    // Update the object
    return updateAABB(index, lowerBound, upperBound, alwaysReinsert);
}

bool AABBTree::updateAABB(uint32_t index, const Eigen::Vector3d& lowerBound,
                          const Eigen::Vector3d& upperBound, bool alwaysReinsert)
{
    // Map iterator.
    std::unordered_map<uint32_t, uint32_t>::iterator it;

    // Find the AABBNode of the object
    it = externalMap.find(index);

    // The object doesn't exist.
    if (it == externalMap.end())
    {
        throw std::invalid_argument("[ERROR]: Invalid external object index!");
    }

    // Extract the node index.
    uint32_t node = it->second;

    assert(node < nodeCapacity);
    assert(nodes[node].isLeaf());

    // AABB size
    Eigen::Vector3d size;

    // Compute the AABB limits.
    for (uint32_t ii=0; ii<3; ii++)
    {
        // Validate the bound.
        if (lowerBound[ii] > upperBound[ii])
        {
            throw std::invalid_argument("AABB lower bound is greater than the upper bound!");
        }
    }
    size = upperBound - lowerBound;

    // Create the new AABB.
    AABB aabb(lowerBound, upperBound);

    // No need to update if the object is still within its fattened AABB.
    if (!alwaysReinsert && nodes[node].aabb.contains(aabb))
        return false;

    // Remove the current leaf.
    removeLeaf(node);

    // Fatten the new AABB.
    aabb.lowerBound -= skinThickness * size;
    aabb.upperBound += skinThickness * size;

    // Assign the new AABB.
    nodes[node].aabb = aabb;

    // Update the surface area and centroid.
    nodes[node].aabb.surfaceArea(nodes[node].aabb.computeSurfaceArea());
    nodes[node].aabb.center(nodes[node].aabb.computeCenter());

    // Insert a new leaf node.
    insertLeaf(node);

    return true;
}

std::vector<uint32_t> AABBTree::query(uint32_t index) const
{
    // Make sure that this is a valid object.
    if (externalMap.count(index) == 0)
    {
        throw std::invalid_argument("[ERROR]: Invalid external index index!");
    }

    // Test overlap of object's AABB against all other objects
    return query(index, nodes[externalMap.find(index)->second].aabb);
}

std::vector<uint32_t> AABBTree::query(uint32_t index, const AABB& aabb) const
{
    std::vector<uint32_t> stack;
    stack.reserve(256);
    stack.push_back(root);

    std::vector<uint32_t> indices;

    while (stack.size() > 0)
    {
        uint32_t node = stack.back();
        stack.pop_back();

        // Copy the AABB.
        AABB nodeAABB = nodes[node].aabb;

        if (node == NULL_NODE) continue;


        // Test for overlap between the AABBs.
        if (aabb.overlaps(nodeAABB, touchIsOverlap))
        {
            // Check that we're at a leaf node.
            if (nodes[node].isLeaf())
            {
                // Can't interact with itself.
                if (nodes[node].origin != index)
                {
                    indices.push_back(nodes[node].origin);
                }
            }
            else
            {
                stack.push_back(nodes[node].left);
                stack.push_back(nodes[node].right);
            }
        }
    }

    return indices;
}

std::vector<uint32_t> AABBTree::query(const AABB& aabb) const
{
    // Make sure the tree isn't empty.
    if (externalMap.size() == 0)
    {
        return std::vector<uint32_t>();
    }

    // Test overlap of AABB against all objects
    return query(std::numeric_limits<uint32_t>::max(), aabb);
}

const AABB& AABBTree::getAABB(uint32_t index) const
{
    auto it = externalMap.find(index);
    if (it == externalMap.end())
    {
        throw std::invalid_argument("[ERROR]: Invalid particle index!");
    }

    // Extract the node index.
    return nodes[it->second].aabb;
    // return nodes[externalMap[index]].aabb;
}

void AABBTree::insertLeaf(uint32_t leaf)
{
    if (root == NULL_NODE)
    {
        root = leaf;
        nodes[root].parent = NULL_NODE;
        return;
    }

    // Find the best sibling for the node.

    AABB leafAABB = nodes[leaf].aabb;
    uint32_t index = root;

    while (!nodes[index].isLeaf())
    {
        // Extract the children of the node.
        uint32_t left  = nodes[index].left;
        uint32_t right = nodes[index].right;

        double surfaceArea = nodes[index].aabb.surfaceArea();

        AABB combinedAABB;
        combinedAABB.merge(nodes[index].aabb, leafAABB);
        double combinedSurfaceArea = combinedAABB.surfaceArea();

        // Cost of creating a new parent for this node and the new leaf.
        double cost = 2.0 * combinedSurfaceArea;

        // Minimum cost of pushing the leaf further down the tree.
        double inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

        // Cost of descending to the left.
        double costLeft;
        if (nodes[left].isLeaf())
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[left].aabb);
            costLeft = aabb.surfaceArea() + inheritanceCost;
        }
        else
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[left].aabb);
            double oldArea = nodes[left].aabb.surfaceArea();
            double newArea = aabb.surfaceArea();
            costLeft = (newArea - oldArea) + inheritanceCost;
        }

        // Cost of descending to the right.
        double costRight;
        if (nodes[right].isLeaf())
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[right].aabb);
            costRight = aabb.surfaceArea() + inheritanceCost;
        }
        else
        {
            AABB aabb;
            aabb.merge(leafAABB, nodes[right].aabb);
            double oldArea = nodes[right].aabb.surfaceArea();
            double newArea = aabb.surfaceArea();
            costRight = (newArea - oldArea) + inheritanceCost;
        }

        // Descend according to the minimum cost.
        if ((cost < costLeft) && (cost < costRight))
            break;

        // Descend.
        if (costLeft < costRight)
            index = left;
        else
            index = right;
    }

    uint32_t sibling = index;

    // Create a new parent.
    uint32_t oldParent = nodes[sibling].parent;
    uint32_t newParent = allocateNode();
    nodes[newParent].parent = oldParent;
    nodes[newParent].aabb.merge(leafAABB, nodes[sibling].aabb);
    nodes[newParent].height = nodes[sibling].height + 1;

    // The sibling was not the root.
    if (oldParent != NULL_NODE)
    {
        if (nodes[oldParent].left == sibling) nodes[oldParent].left = newParent;
        else                                  nodes[oldParent].right = newParent;

        nodes[newParent].left = sibling;
        nodes[newParent].right = leaf;
        nodes[sibling].parent = newParent;
        nodes[leaf].parent = newParent;
    }
    // The sibling was the root.
    else
    {
        nodes[newParent].left = sibling;
        nodes[newParent].right = leaf;
        nodes[sibling].parent = newParent;
        nodes[leaf].parent = newParent;
        root = newParent;
    }

    // Walk back up the tree fixing heights and AABBs.
    index = nodes[leaf].parent;
    while (index != NULL_NODE)
    {
        index = balance(index);

        uint32_t left = nodes[index].left;
        uint32_t right = nodes[index].right;

        assert(left != NULL_NODE);
        assert(right != NULL_NODE);

        nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);
        nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);

        index = nodes[index].parent;
    }
}

void AABBTree::removeLeaf(uint32_t leaf)
{
    if (leaf == root)
    {
        root = NULL_NODE;
        return;
    }

    uint32_t parent = nodes[leaf].parent;
    uint32_t grandParent = nodes[parent].parent;
    uint32_t sibling;

    if (nodes[parent].left == leaf) sibling = nodes[parent].right;
    else                            sibling = nodes[parent].left;

    // Destroy the parent and connect the sibling to the grandparent.
    if (grandParent != NULL_NODE)
    {
        if (nodes[grandParent].left == parent) nodes[grandParent].left = sibling;
        else                                   nodes[grandParent].right = sibling;

        nodes[sibling].parent = grandParent;
        freeNode(parent);

        // Adjust ancestor bounds.
        uint32_t index = grandParent;
        while (index != NULL_NODE)
        {
            index = balance(index);

            uint32_t left = nodes[index].left;
            uint32_t right = nodes[index].right;

            nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);
            nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);

            index = nodes[index].parent;
        }
    }
    else
    {
        root = sibling;
        nodes[sibling].parent = NULL_NODE;
        freeNode(parent);
    }
}

uint32_t AABBTree::balance(uint32_t node)
{
    assert(node != NULL_NODE);

    if (nodes[node].isLeaf() || (nodes[node].height < 2))
        return node;

    uint32_t left = nodes[node].left;
    uint32_t right = nodes[node].right;

    assert(left < nodeCapacity);
    assert(right < nodeCapacity);

    int currentBalance = nodes[right].height - nodes[left].height;

    // Rotate right branch up.
    if (currentBalance > 1)
    {
        uint32_t rightLeft = nodes[right].left;
        uint32_t rightRight = nodes[right].right;

        assert(rightLeft < nodeCapacity);
        assert(rightRight < nodeCapacity);

        // Swap node and its right-hand child.
        nodes[right].left = node;
        nodes[right].parent = nodes[node].parent;
        nodes[node].parent = right;

        // The node's old parent should now point to its right-hand child.
        if (nodes[right].parent != NULL_NODE)
        {
            if (nodes[nodes[right].parent].left == node) nodes[nodes[right].parent].left = right;
            else
            {
                assert(nodes[nodes[right].parent].right == node);
                nodes[nodes[right].parent].right = right;
            }
        }
        else root = right;

        // Rotate.
        if (nodes[rightLeft].height > nodes[rightRight].height)
        {
            nodes[right].right = rightLeft;
            nodes[node].right = rightRight;
            nodes[rightRight].parent = node;
            nodes[node].aabb.merge(nodes[left].aabb, nodes[rightRight].aabb);
            nodes[right].aabb.merge(nodes[node].aabb, nodes[rightLeft].aabb);

            nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightRight].height);
            nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightLeft].height);
        }
        else
        {
            nodes[right].right = rightRight;
            nodes[node].right = rightLeft;
            nodes[rightLeft].parent = node;
            nodes[node].aabb.merge(nodes[left].aabb, nodes[rightLeft].aabb);
            nodes[right].aabb.merge(nodes[node].aabb, nodes[rightRight].aabb);

            nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightLeft].height);
            nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightRight].height);
        }

        return right;
    }

    // Rotate left branch up.
    if (currentBalance < -1)
    {
        uint32_t leftLeft = nodes[left].left;
        uint32_t leftRight = nodes[left].right;

        assert(leftLeft < nodeCapacity);
        assert(leftRight < nodeCapacity);

        // Swap node and its left-hand child.
        nodes[left].left = node;
        nodes[left].parent = nodes[node].parent;
        nodes[node].parent = left;

        // The node's old parent should now point to its left-hand child.
        if (nodes[left].parent != NULL_NODE)
        {
            if (nodes[nodes[left].parent].left == node) nodes[nodes[left].parent].left = left;
            else
            {
                assert(nodes[nodes[left].parent].right == node);
                nodes[nodes[left].parent].right = left;
            }
        }
        else root = left;

        // Rotate.
        if (nodes[leftLeft].height > nodes[leftRight].height)
        {
            nodes[left].right = leftLeft;
            nodes[node].left = leftRight;
            nodes[leftRight].parent = node;
            nodes[node].aabb.merge(nodes[right].aabb, nodes[leftRight].aabb);
            nodes[left].aabb.merge(nodes[node].aabb, nodes[leftLeft].aabb);

            nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftRight].height);
            nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftLeft].height);
        }
        else
        {
            nodes[left].right = leftRight;
            nodes[node].left = leftLeft;
            nodes[leftLeft].parent = node;
            nodes[node].aabb.merge(nodes[right].aabb, nodes[leftLeft].aabb);
            nodes[left].aabb.merge(nodes[node].aabb, nodes[leftRight].aabb);

            nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftLeft].height);
            nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftRight].height);
        }

        return left;
    }

    return node;
}

uint32_t AABBTree::computeHeight() const
{
    return computeHeight(root);
}

uint32_t AABBTree::computeHeight(uint32_t node) const
{
    assert(node < nodeCapacity);

    if (nodes[node].isLeaf()) return 0;

    uint32_t height1 = computeHeight(nodes[node].left);
    uint32_t height2 = computeHeight(nodes[node].right);

    return 1 + std::max(height1, height2);
}

uint32_t AABBTree::getHeight() const
{
    if (root == NULL_NODE) return 0;
    return nodes[root].height;
}

uint32_t AABBTree::getNodeCount() const
{
    return nodeCount;
}

uint32_t AABBTree::computeMaximumBalance() const
{
    uint32_t maxBalance = 0;
    for (uint32_t i=0; i<nodeCapacity; i++)
    {
        if (nodes[i].height <= 1)
            continue;

        assert(nodes[i].isLeaf() == false);

        uint32_t balance = std::abs(nodes[nodes[i].left].height - nodes[nodes[i].right].height);
        maxBalance = std::max(maxBalance, balance);
    }

    return maxBalance;
}

double AABBTree::computeSurfaceAreaRatio() const
{
    if (root == NULL_NODE) return 0.0;

    double rootArea = nodes[root].aabb.computeSurfaceArea();
    double totalArea = 0.0;

    for (uint32_t i=0; i<nodeCapacity;i++)
    {
        if (nodes[i].height < 0) continue;

        totalArea += nodes[i].aabb.computeSurfaceArea();
    }

    return totalArea / rootArea;
}

void AABBTree::validate() const
{
    validateStructure(root);
    validateMetrics(root);

    uint32_t freeCount = 0;
    uint32_t freeIndex = freeList;

    while (freeIndex != NULL_NODE)
    {
        assert(freeIndex < nodeCapacity);
        freeIndex = nodes[freeIndex].next;
        freeCount++;
    }

    assert(getHeight() == computeHeight());
    assert((nodeCount + freeCount) == nodeCapacity);
}

void AABBTree::rebuild()
{
    std::vector<uint32_t> nodeIndices(nodeCount);
    uint32_t count = 0;

    for (uint32_t i=0;i<nodeCapacity;i++)
    {
        // Free node.
        if (nodes[i].height < 0) continue;

        if (nodes[i].isLeaf())
        {
            nodes[i].parent = NULL_NODE;
            nodeIndices[count] = i;
            count++;
        }
        else freeNode(i);
    }

    while (count > 1)
    {
        double minCost = std::numeric_limits<double>::max();
        int iMin = -1, jMin = -1;

        for (uint32_t i=0;i<count;i++)
        {
            AABB aabbi = nodes[nodeIndices[i]].aabb;

            for (uint32_t j=i+1;j<count;j++)
            {
                AABB aabbj = nodes[nodeIndices[j]].aabb;
                AABB aabb;
                aabb.merge(aabbi, aabbj);
                double cost = aabb.surfaceArea();

                if (cost < minCost)
                {
                    iMin = i;
                    jMin = j;
                    minCost = cost;
                }
            }
        }

        uint32_t index1 = nodeIndices[iMin];
        uint32_t index2 = nodeIndices[jMin];

        uint32_t parent = allocateNode();
        nodes[parent].left = index1;
        nodes[parent].right = index2;
        nodes[parent].height = 1 + std::max(nodes[index1].height, nodes[index2].height);
        nodes[parent].aabb.merge(nodes[index1].aabb, nodes[index2].aabb);
        nodes[parent].parent = NULL_NODE;

        nodes[index1].parent = parent;
        nodes[index2].parent = parent;

        nodeIndices[jMin] = nodeIndices[count-1];
        nodeIndices[iMin] = parent;
        count--;
    }

    root = nodeIndices[0];

    validate();
}

void AABBTree::validateStructure(uint32_t node) const
{
    if (node == NULL_NODE) return;

    if (node == root) assert(nodes[node].parent == NULL_NODE);

    uint32_t left = nodes[node].left;
    uint32_t right = nodes[node].right;

    if (nodes[node].isLeaf())
    {
        assert(left == NULL_NODE);
        assert(right == NULL_NODE);
        assert(nodes[node].height == 0);
        return;
    }

    assert(left < nodeCapacity);
    assert(right < nodeCapacity);

    assert(nodes[left].parent == node);
    assert(nodes[right].parent == node);

    validateStructure(left);
    validateStructure(right);
}

void AABBTree::validateMetrics(uint32_t node) const
{
    if (node == NULL_NODE) return;

    uint32_t left = nodes[node].left;
    uint32_t right = nodes[node].right;

    if (nodes[node].isLeaf())
    {
        assert(left == NULL_NODE);
        assert(right == NULL_NODE);
        assert(nodes[node].height == 0);
        return;
    }

    assert(left < nodeCapacity);
    assert(right < nodeCapacity);

    int height1 = nodes[left].height;
    int height2 = nodes[right].height;
    int height = 1 + std::max(height1, height2);
    (void)height; // Unused variable in Release build
    assert(nodes[node].height == height);

    AABB aabb;
    aabb.merge(nodes[left].aabb, nodes[right].aabb);

    for (uint32_t ii=0; ii<3; ii++)
    {
        assert(aabb.lowerBound[ii] == nodes[node].aabb.lowerBound[ii]);
        assert(aabb.upperBound[ii] == nodes[node].aabb.upperBound[ii]);
    }

    validateMetrics(left);
    validateMetrics(right);
}

std::unordered_set<std::pair<uint32_t, uint32_t>> AABBTree::intersect(const AABBTree& tree)
{
    std::unordered_set<std::pair<uint32_t, uint32_t>> intersections;
    for (uint32_t ii=0; ii<numObjects(); ++ii)
    {
        std::vector<uint32_t> interx = tree.query(getAABB(ii));
        for (const auto& ix : interx)
            intersections.emplace(std::make_pair(ii, ix));
    }

    std::cout << intersections.size() << " potential intersections found"
              << std::endl;
    return intersections;
}
