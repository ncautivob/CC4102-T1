#include <bits/stdc++.h>
using namespace std;
typedef long double ld;

#define B 4096

using Point = pair<ld, ld>;

typedef struct {
  Point point;
  ld coveringRadius;
  set<Entry>* children;  // If not set<Entry>, then it is a single point.
} Entry;

using Node = set<Entry>;  // Conjunto de puntos
using Nodes = set<Node>;  // set<set<Entry>>

/**
 * Computes the euclidean distance between firstPoint and secondPoint.
 */
ld distance(Point& firstPoint, Point& secondPoint) {
  ld delta_x = firstPoint.first - secondPoint.first;
  ld delta_y = firstPoint.second - secondPoint.second;

  return (ld)sqrt(delta_x * delta_x + delta_y * delta_y);
}

Nodes nearestNodes(Nodes input) {
  ld minDistance = 1e18;
  Nodes candidateNearest = {NULL, NULL};

  /**
   * It is known that euclidean distance is symmetric, id. est., d(x, y) = d(y,
   * x). Even more, d(x, x) = 0. So we discard the diagonal and one half of the
   * total pairs. Here, we are doing binom(n, 2) comparisons ~ O(n^2), but it's
   * way better that use the naive approach to do n^2 comparisons.
   */
  for (auto entry = input.begin(); entry != input.end(); entry++) {
    Node node = *entry;
    Entry medioid = computeMedioid(node);
    Point point = medioid.point;
    for (auto otherEntry = next(entry); otherEntry != input.end();
         otherEntry++) {
      Node otherNode = *otherEntry;
      Entry otherMedioid = computeMedioid(otherNode);
      Point otherPoint = otherMedioid.point;

      // As we are using <=, we return the very last pair of clusters, if there
      // are more than one.
      ld dist = distance(point, otherPoint);
      if (dist <= minDistance) {
        minDistance = dist;
        candidateNearest = {node, otherNode};
      }
    }
  }

  return candidateNearest;
};

Node nearestNeighbour(Node objective, Nodes nodes) {
  Entry objectiveMedioid = computeMedioid(objective);
  Point objectivePoint = objectiveMedioid.point;

  Node candidateNearest = {};
  ld minDistance = 1e18;

  for (auto entry = nodes.begin(); entry != nodes.end(); entry++) {
    Node node = *entry;
    Entry nodeMedioid = computeMedioid(node);
    Point nodePoint = nodeMedioid.point;

    ld dist = distance(objectivePoint, nodePoint);
    /**
     * <= is a relajation to choose the very last node that
     * verifies this condition. Dist must not be 0, in other
     * case we are talking about the same node.
     */
    if (dist <= minDistance && dist != 0) {
      minDistance = dist;
      candidateNearest = node;
    }
  }

  return candidateNearest;
}

/**
 * Computes a medioid given a set of points (cluster).
 * If multiple medioids, it computes the last in the set order.
 */
Entry computeMedioid(Node input) {
  if (input.empty()) return {{0.0, 0.0}, NULL, NULL};

  Entry initial = *(input.begin());
  Entry candidateMedioid = initial;

  // A very large number, because we want to minimize.
  ld minTotalDistance = (ld)1e18;

  for (auto entry = input.begin(); entry != input.end(); entry++) {
    // Initialize the total distance for this point.
    ld totalDistance = 0.0;
    // The (p, c_r, a) tuple.
    Entry tuple = *entry;
    Point point = tuple.point;

    for (auto otherEntry = input.begin(); otherEntry != input.end();
         otherEntry++) {
      Entry otherTuple = *otherEntry;
      Point otherPoint = otherTuple.point;

      totalDistance += distance(point, otherPoint);
    }

    /**
     * <= is a relajation to choose the very last point that
     * verifies this condition.
     */
    if (totalDistance <= minTotalDistance) {
      minTotalDistance = totalDistance;
      candidateMedioid = tuple;
    }
  }

  return candidateMedioid;
}

Entry outputHoja(Node input) {
  Entry primaryMedioid = computeMedioid(input);
  Point g = primaryMedioid.point;
  ld r = 0.0;
  Node cluster;

  for (Entry entry : input) {
    Point p = entry.point;
    Entry newEntry = {p, NULL, NULL};

    cluster.insert(newEntry);
    r = max(r, distance(g, p));
  }

  Node* a = &cluster;
  return {g, r, a};
}

Entry outputInterno(Node input) {
  Entry primaryMedioid = computeMedioid(input);
  Point G = primaryMedioid.point;
  ld R = 0.0;
  Node cluster;

  for (Entry entry : input) {
    Point g = entry.point;
    ld r = entry.coveringRadius;

    cluster.insert(entry);
    R = max(R, distance(G, g) + r);
  }

  Node* A = &cluster;
  return {G, R, A};
}

Entry SSAlgorithm(Node input) {
  if (sizeof(input) <= B) {
    Entry leaf = outputHoja(input);
    return leaf;
  }

  Nodes out = cluster(input);
  Node C;

  for (Node Node : out) C.insert(outputHoja(Node));

  // Falta del paso 4 en abajo
  while (sizeof(C) > B) {
  }
}

Nodes cluster(Node input) {
  Nodes cluster_output, aux_cluster;

  for (Entry p : input) {
    Node newNode;
    newNode.insert(p);  // Insertamos el singleton {p}
    aux_cluster.insert(newNode);
  }

  while (aux_cluster.size() > 1) {
    Nodes nearest = nearestNodes(aux_cluster);

    // Nearest must only be two, as nearestNodes returns two nodes.
    // We retrieve this data to erase the iterators from the set.
    auto first = nearest.begin();
    auto second = prev(nearest.end());

    // It must hold the inequality |firstNode| >= |secondNode|, by instructions.
    if (sizeof(*first) <= sizeof(*second)) {
      swap(first, second);
    }

    Node firstNode = *first;
    Node secondNode = *second;

    Node unionNode;
    unionNode.insert(firstNode.begin(), firstNode.end());
    unionNode.insert(secondNode.begin(), secondNode.end());

    if (sizeof(unionNode) <= B) {
      aux_cluster.erase(first);
      aux_cluster.erase(second);

      aux_cluster.insert(unionNode);
    } else {
      aux_cluster.erase(first);
      cluster_output.insert(firstNode);
    }
  }

  Node lastNode = *(prev(aux_cluster.end()));
  Node nearestToLast = {};
  if (sizeof(cluster_output) > 0) {
    nearestToLast = nearestNeighbour(lastNode, cluster_output);
    cluster_output.erase(cluster_output.find(nearestToLast));
  }

  // Falta del paso 6 en adelante
};

int main() {
  int a = 1;
  return 0;
}
