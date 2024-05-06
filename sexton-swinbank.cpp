#include <bits/stdc++.h>
using namespace std;

using Point = pair<double, double>;

typedef struct {
  Point point;
  double coveringRadius;
  set<Entry>* children;  // If not set<Entry>, then it is a single point.
} Entry;

#define ENTRY_SIZE sizeof(Entry)
#define B 4096 / ENTRY_SIZE
#define b 0.5 * B / ENTRY_SIZE

using Node = set<Entry>;            // Entries set.
using Nodes = set<Node>;            // Set that have entries set as elements.
using Points = set<Point>;          // Set of points.
using PointClusters = set<Points>;  // Set that have points set as elements.

/**
 * Computes the squared euclidean distance between firstPoint and secondPoint.
 * As the square function strictly increases in the positive X-axis, then the
 * result doesn't change.
 */
double distance(Point& firstPoint, Point& secondPoint) {
  double delta_x = firstPoint.first - secondPoint.first;
  double delta_y = firstPoint.second - secondPoint.second;

  return (double)delta_x * delta_x + delta_y * delta_y;
}

/**
 * Uses the MinMax split policy to separate a cluster into two.
 */
PointClusters splittingPolicy(Points input) {
  vector<Point> pointIndexing;
  // distances[i] = (d, j), where d is distance(i, j).
  vector<set<pair<double, int>>> distances;
  // Marks if the point was already used or not.
  vector<bool> used(input.size(), false);

  // Point mapping to index
  for (auto entry = input.begin(); entry != input.end(); entry++) {
    Point point = *entry;
    pointIndexing.push_back(point);
  }

  int numberOfPoints = pointIndexing.size();
  // binom(n, 2) because distance is symmetric.
  for (int i = 0; i < numberOfPoints; i++) {
    Point first = pointIndexing[i];
    for (int j = i + 1; j < numberOfPoints; j++) {
      Point second = pointIndexing[j];

      double dist = distance(first, second);
      distances[i].insert({dist, j});
      distances[j].insert({dist, i});
    }
  }

  double minMaxRadius = DBL_MAX;
  pair<int, int> candidateClusters = {NULL, NULL};

  for (int i = 0; i < numberOfPoints; i++) {
    for (int j = i + 1; j < numberOfPoints; j++) {
      used[i] = true;
      used[j] = true;

      auto firstItr = distances[i].begin();
      auto secondItr = distances[j].begin();

      double firstRadius = 0.0;
      double secondRadius = 0.0;

      int usedPoints = 2;
      while (usedPoints < numberOfPoints) {
        // If used points are even, then it's i's turn. Otherwise, it's j's.
        int turn = usedPoints % 2 == 0 ? i : j;

        pair<double, int> selected;
        if (turn == i) {
          selected = *firstItr;

          // Move the iterator while the selected point was already used...
          while (used[selected.second]) {
            firstItr++;
            selected = *firstItr;
          }

          // selected.first is the distance between i and the selected point.
          firstRadius = max(firstRadius, selected.first);
          /**
           * Mark as used and move the iterator as we don't care about this
           * element anymore.
           */
          used[selected.second] = true;
          firstItr++;
        } else {
          selected = *secondItr;

          while (used[selected.second]) {
            secondItr++;
            selected = *secondItr;
          }

          secondRadius = max(secondRadius, selected.first);
          used[selected.second] = true;
          secondItr++;
        }

        usedPoints++;
      }

      double maxRadius = max(firstRadius, secondRadius);

      if (maxRadius < minMaxRadius) {
        minMaxRadius = maxRadius;
        candidateClusters = {i, j};
      }

      // Reset used points for the next iteration.
      used.assign(numberOfPoints, false);
    }
  }

  /**
   * As candidateClusters must have the best indices, we use them to
   * generate the two set of points.
   */
  Points p, q;

  int pIndex = candidateClusters.first;
  int qIndex = candidateClusters.second;

  /**
   * These (pIndex, qIndex) are the points that generate this cluster
   * configuration (p, q).
   */
  p.insert(pointIndexing[pIndex]);
  q.insert(pointIndexing[qIndex]);

  used[pIndex] = true;
  used[qIndex] = true;

  auto pItr = distances[pIndex].begin();
  auto qItr = distances[qIndex].begin();

  int usedPoints = 2;
  while (usedPoints < numberOfPoints) {
    int turn = usedPoints % 2 == 0 ? pIndex : qIndex;

    pair<double, int> selected;
    if (turn == pIndex) {
      selected = *pItr;

      while (used[selected.second]) {
        pItr++;
        selected = *pItr;
      }

      // selected.second is the point marked by this iteration.
      p.insert(pointIndexing[selected.second]);
      used[selected.second] = true;
      pItr++;
    } else {
      selected = *qItr;

      while (used[selected.second]) {
        qItr++;
        selected = *qItr;
      }

      q.insert(pointIndexing[selected.second]);
      used[selected.second] = true;
      qItr++;
    }

    usedPoints++;
  }

  return {p, q};
}

/**
 * Joins two points set first and second.
 */
Points pointsUnion(Points first, Points second) {
  Points result;
  result.insert(first.begin(), first.end());
  result.insert(second.begin(), second.end());

  return result;
}

PointClusters nearestPointClusters(PointClusters input) {
  double minDistance = DBL_MAX;
  PointClusters candidateNearest = {NULL, NULL};

  /**
   * It is known that euclidean distance is symmetric, id. est., d(x, y) = d(y,
   * x). Even more, d(x, x) = 0. So we discard the diagonal and one half of the
   * total pairs. Here, we are doing binom(n, 2) comparisons ~ O(n^2), but it's
   * way better that use the naive approach to do n^2 comparisons.
   */
  for (auto entry = input.begin(); entry != input.end(); entry++) {
    Points points = *entry;
    Point medioid = computeMedioid(points);

    for (auto otherEntry = next(entry); otherEntry != input.end();
         otherEntry++) {
      Points otherPoints = *otherEntry;
      Point otherMedioid = computeMedioid(otherPoints);

      // As we are using <=, we return the very last pair of point clusters, if
      // there are more than one that satisfies the inequality.
      double dist = distance(medioid, otherMedioid);
      if (dist <= minDistance) {
        minDistance = dist;
        candidateNearest = {points, otherPoints};
      }
    }
  }

  return candidateNearest;
};

Points nearestNeighbour(Points objective, PointClusters clusters) {
  Point objectiveMedioid = computeMedioid(objective);

  Points candidateNearest = {};
  double minDistance = DBL_MAX;

  for (auto entry = clusters.begin(); entry != clusters.end(); entry++) {
    Points points = *entry;
    Point pointsMedioid = computeMedioid(points);

    double dist = distance(objectiveMedioid, pointsMedioid);
    /**
     * <= is a relajation to choose the very last node that
     * verifies this condition. Dist must not be 0, in other
     * case we are talking about the same node.
     */
    if (dist <= minDistance && dist != 0) {
      minDistance = dist;
      candidateNearest = points;
    }
  }

  return candidateNearest;
}

/**
 * Computes a medioid given a set of points.
 * If multiple medioids, it computes the last in the set order.
 */
Point computeMedioid(Points input) {
  if (input.empty()) return {0.0, 0.0};

  Point initial = *(input.begin());
  Point candidateMedioid = initial;

  // A very large number, because we want to minimize.
  double minTotalDistance = DBL_MAX;

  for (auto entry = input.begin(); entry != input.end(); entry++) {
    // Initialize the total distance for this point.
    double totalDistance = 0.0;
    // The current point to be compared.
    Point point = *entry;

    for (auto otherEntry = input.begin(); otherEntry != input.end();
         otherEntry++) {
      Point otherPoint = *otherEntry;

      totalDistance += distance(point, otherPoint);
    }

    /**
     * <= is a relajation to choose the very last point that
     * verifies this condition.
     */
    if (totalDistance <= minTotalDistance) {
      minTotalDistance = totalDistance;
      candidateMedioid = point;
    }
  }

  return candidateMedioid;
}

/**
 * Returns a tuple (g, r, a) where g is the primary medioid of the input,
 * r is the covering radius, and a the node direction.
 */
Entry OutputLeafPage(Points input) {
  Point g = computeMedioid(input);
  double r = 0.0;
  Node cluster;

  for (Point p : input) {
    Entry newEntry = {p, NULL, NULL};

    cluster.insert(newEntry);
    r = max(r, distance(g, p));
  }

  Node* a = &cluster;
  return {g, r, a};
}

/**
 * Returns a tuple (G, R, A) where G is the primary medioid of the points
 * that belongs to the input node, R is the covering radius and A the node
 * direction.
 */
Entry OutputInternalPage(Node input) {
  Points C_in;
  // We add the points in the input node to C_in.
  for (Entry entry : input) C_in.insert(entry.point);

  Point G = computeMedioid(C_in);
  double R = 0.0;
  Node cluster;

  for (Entry entry : input) {
    Point g = entry.point;
    double r = entry.coveringRadius;
    Entry* a = entry.children;

    cluster.insert(entry);
    R = max(R, distance(G, g) + r);
  }

  Node* A = &cluster;
  return {G, R, A};
}

Entry* SSAlgorithm(Points input) {
  // There are at most B points.
  if (input.size() <= B) {
    Entry rootRepresentation = OutputLeafPage(input);
    return rootRepresentation.children;
  }

  // We have > B points. Then, we must do clustering to reduce size.
  PointClusters C_out = cluster(input);
  Node C = {};

  for (Points points : C_out) {
    C.insert(OutputLeafPage(points));
  }

  /**
   * We do this in order to retrieve an entry given a point, to optimize
   * search. Warning: Saving is global, and it doesn't reset in every while
   * iteration, this may be a problem.
   */
  map<Point, Entry> pointToEntry;

  while (C.size() > B) {
    Points C_in;

    // We insert into C_in the points contained in C.
    for (Entry entry : C) {
      Point point = entry.point;
      C_in.insert(point);

      // Only add to the mapping if is not present.
      if (pointToEntry.find(point) == pointToEntry.end()) {
        pointToEntry[point] = entry;
      }
    }

    // Clustering again to reduce set size of C.
    C_out = cluster(C_in);
    Nodes C_mra = {};

    for (Points points : C_out) {
      Node s;
      for (Point point : points) {
        s.insert(pointToEntry[point]);
      }

      C_mra.insert(s);
    }

    C = {};
    for (Node node : C_mra) {
      C.insert(OutputInternalPage(node));
    }
  }

  Entry rootRepresentation = OutputInternalPage(C);
  return rootRepresentation.children;
}

/**
 * Split a set of points into clusters of size at least b and size at most B.
 */
PointClusters cluster(Points input) {
  PointClusters C_out = {}, C = {};

  for (Point p : input) {
    // We insert the singleton {p}.
    C.insert({p});
  }

  while (C.size() > 1) {
    PointClusters nearest = nearestPointClusters(C);

    // Nearest must only be two, as nearestNodes returns two nodes.
    // We retrieve this data to erase the iterators from the set.
    auto first = nearest.begin();
    auto second = prev(nearest.end());

    // It must hold the inequality |firstNode| >= |secondNode|, by instructions.
    if (sizeof(*first) <= sizeof(*second)) {
      swap(first, second);
    }

    Points firstPointSet = *first;
    Points secondPointSet = *second;

    Points firstUnion = pointsUnion(firstPointSet, secondPointSet);

    if (firstUnion.size() <= B) {
      C.erase(first);
      C.erase(second);

      C.insert(firstUnion);
    } else {
      C.erase(first);
      C_out.insert(firstPointSet);
    }
  }

  Points lastPointSet = *(prev(C.end()));
  Points nearestToLast = {};
  if (C_out.size() > 0) {
    nearestToLast = nearestNeighbour(lastPointSet, C_out);
    C_out.erase(C_out.find(nearestToLast));
  }

  Points secondUnion = pointsUnion(lastPointSet, nearestToLast);
  if (secondUnion.size() <= B) {
    C_out.insert(secondUnion);
  } else {
    PointClusters separated = splittingPolicy(secondUnion);

    Points newFirst = *(separated.begin());
    Points newSecond = *(prev(separated.end()));

    C_out.insert(newFirst);
    C_out.insert(newSecond);
  }

  return C_out;
};

int main() {
  // The test set goes here...
  return 0;
}
