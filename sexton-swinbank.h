/**
 * Sexton-Swinbank algorithm to construct M-trees.
 *
 * Valentina Alarcón Y. [https://github.com/valentina-alarcon-y]
 * Naomí Cautivo B. [https://github.com/ncautivob]
 * Máximo Flores Valenzuela [https://github.com/maxfloresv]
 */

#include <cfloat>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <utility>
#include <vector>

using namespace std;
using namespace std::chrono;

using Point = pair<double, double>;

struct Node;

typedef struct {
  Point point;
  double coveringRadius;
  Node* children;
} Entry;

#define ENTRY_SIZE sizeof(Entry)
#define B 4096 / ENTRY_SIZE
#define b 0.5 * B / ENTRY_SIZE
#define N_ITER 100
#define query_radius 0.02

using Entries = vector<Entry>;
typedef struct Node {
  Entries entries;
} Node;
using Nodes = vector<Node>;
using Points = set<Point>;
using PointClusters = set<Points>;

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
  int numberOfPoints = (int)input.size();

  vector<Point> pointIndexing;
  // distances[i] = (d, j), where d is distance(i, j).
  vector<set<pair<double, int>>> distances(numberOfPoints);
  // Marks if the point was already used or not.
  vector<bool> used(numberOfPoints, false);

  // Point mapping to index
  for (Point point : input) {
    pointIndexing.push_back(point);
  }

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
  pair<int, int> candidateClusters = {-1, -1};

  for (int i = 0; i < numberOfPoints; i++) {
    for (int j = i + 1; j < numberOfPoints; j++) {
      used[i] = true;
      used[j] = true;

      auto firstItr = distances[i].begin();
      auto secondItr = distances[j].begin();

      double firstRadius = 0.0;
      double secondRadius = 0.0;

      int usedPoints = 2;
      while (usedPoints < numberOfPoints && firstItr != distances[i].end() &&
             secondItr != distances[j].end()) {
        // If used points are even, then it's i's turn. Otherwise, it's j's.
        int turn = usedPoints % 2 == 0 ? i : j;

        pair<double, int> selected;
        if (turn == i) {
          selected = *firstItr;

          // Move the iterator while the selected point was already used...
          while (used[selected.second] && firstItr != distances[i].end()) {
            firstItr++;
            selected = *firstItr;
          }

          if (firstItr != distances[i].end()) {
            // selected.first is the distance between i and the selected point.
            firstRadius = max(firstRadius, selected.first);
            /**
             * Mark as used and move the iterator as we don't care about this
             * element anymore.
             */
            used[selected.second] = true;
            firstItr++;
          }
        } else {
          selected = *secondItr;

          while (used[selected.second] && secondItr != distances[j].end()) {
            secondItr++;
            selected = *secondItr;
          }

          if (secondItr != distances[j].end()) {
            secondRadius = max(secondRadius, selected.first);
            used[selected.second] = true;
            secondItr++;
          }
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
  while (usedPoints < numberOfPoints && pItr != distances[pIndex].end() &&
         qItr != distances[qIndex].end()) {
    int turn = usedPoints % 2 == 0 ? pIndex : qIndex;

    pair<double, int> selected;
    if (turn == pIndex) {
      selected = *pItr;

      while (used[selected.second] && pItr != distances[pIndex].end()) {
        pItr++;
        selected = *pItr;
      }

      if (pItr != distances[pIndex].end()) {
        // selected.second is the point marked by this iteration.
        p.insert(pointIndexing[selected.second]);
        used[selected.second] = true;
        pItr++;
      }
    } else {
      selected = *qItr;

      while (used[selected.second] && qItr != distances[qIndex].end()) {
        qItr++;
        selected = *qItr;
      }

      if (qItr != distances[qIndex].end()) {
        q.insert(pointIndexing[selected.second]);
        used[selected.second] = true;
        qItr++;
      }
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

  for (Point point : input) {
    // Initialize the total distance for this point.
    double totalDistance = 0.0;

    for (Point otherPoint : input) {
      totalDistance += distance(point, otherPoint);
    }

    /**
     * As it's strictly less, it chooses the very first point
     * that satisfies this condition.
     */
    if (totalDistance < minTotalDistance) {
      minTotalDistance = totalDistance;
      candidateMedioid = point;
    }
  }

  return candidateMedioid;
}

PointClusters nearestPointClusters(PointClusters input) {
  double minDistance = DBL_MAX;
  PointClusters candidateNearest = {};

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

      double dist = distance(medioid, otherMedioid);
      /**
       * As it's strictly less, it chooses the very first PointClusters
       * that satisfies this condition.
       */
      if (dist < minDistance) {
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
     * As it's strictly less, it chooses the very first point set
     * that satisfies this condition.
     */
    if (dist < minDistance && dist != 0) {
      minDistance = dist;
      candidateNearest = points;
    }
  }

  return candidateNearest;
}

/**
 * Returns a tuple (g, r, a) where g is the primary medioid of the input,
 * r is the covering radius, and a the node direction.
 */
Entry OutputLeafPage(Points input) {
  Point g = computeMedioid(input);
  double r = 0.0;
  Node* cluster = new Node;

  for (Point p : input) {
    // 0.0 as NULL in double entry gives a warning.
    Entry newEntry = {p, 0.0, nullptr};

    cluster->entries.push_back(newEntry);
    // Remember that distance(., .) returns the squared euclidean distance.
    r = max(r, sqrt(distance(g, p)));
  }

  return {g, r, cluster};
}

/**
 * Returns a tuple (G, R, A) where G is the primary medioid of the points
 * that belongs to the input node, R is the covering radius and A the node
 * direction.
 */
Entry OutputInternalPage(Node input) {
  Points C_in;
  // We add the points in the input node to C_in.
  for (Entry entry : input.entries) {
    C_in.insert(entry.point);
  }

  Point G = computeMedioid(C_in);
  double R = 0.0;
  Node* cluster = new Node;

  for (Entry entry : input.entries) {
    Point g = entry.point;
    double r = entry.coveringRadius;

    (*cluster).entries.push_back(entry);
    // Remember that distance(., .) returns the squared euclidean distance.
    R = max(R, sqrt(distance(G, g)) + r);
  }

  return {G, R, cluster};
}

/**
 * Split a set of points into clusters of size at least b and size at most B.
 */
PointClusters cluster(Points input) {
  PointClusters C_out, C;

  for (Point p : input) {
    // We insert the singleton {p}.
    C.insert({p});
  }

  while (C.size() > 1) {
    PointClusters nearest = nearestPointClusters(C);

    // Nearest must only be two, as nearestNodes returns two nodes.
    // We retrieve this data to erase the iterators from the set.
    auto first = C.find(*nearest.begin());
    auto second = C.find(*prev(nearest.end()));

    // It must hold the inequality |firstNode| >= |secondNode|, by instructions.
    if ((*first).size() <= (*second).size()) {
      swap(first, second);
    }

    Points firstUnion = pointsUnion(*first, *second);

    if (firstUnion.size() <= B) {
      C.insert(firstUnion);

      C.erase(first);
      C.erase(second);
    } else {
      C_out.insert(*first);
      C.erase(first);
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

Node* SSAlgorithm(Points input) {
  // There are at most B points.
  if (input.size() <= B) {
    Entry rootRepresentation = OutputLeafPage(input);
    return rootRepresentation.children;
  }

  // We have > B points. Then, we must do clustering to reduce size.
  PointClusters C_out = cluster(input);
  Node C;

  for (Points points : C_out) {
    C.entries.push_back(OutputLeafPage(points));
  }

  /**
   * We do this in order to retrieve an entry given a point, to optimize
   * search. Warning: Saving is global, and it doesn't reset in every while
   * iteration, this may be a problem.
   */
  map<Point, Entry> pointToEntry;

  while (C.entries.size() > B) {
    Points C_in;

    // We insert into C_in the points contained in C.
    for (Entry entry : C.entries) {
      Point point = entry.point;
      C_in.insert(point);

      // Only add to the mapping if is not present.
      if (pointToEntry.find(point) == pointToEntry.end()) {
        pointToEntry[point] = entry;
      }
    }

    // Clustering again to reduce set size of C.
    C_out = cluster(C_in);
    Nodes C_mra;

    for (Points points : C_out) {
      Node s;
      for (Point point : points) {
        s.entries.push_back(pointToEntry[point]);
      }

      C_mra.push_back(s);
    }

    C = {};
    for (Node node : C_mra) {
      C.entries.push_back(OutputInternalPage(node));
    }
  }

  Entry rootRepresentation = OutputInternalPage(C);
  return rootRepresentation.children;
};

Points createSet(int n) {
  Points result;

  while (n--) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0);

    double first = dis(gen);
    double second = dis(gen);

    result.insert({first, second});
  }

  return result;
}

/**
 * Does a 1-by-1 deletion for each entry that propagates to its correspondent
 * recursive subtree roots.
 */
void freeMem(Node* root, Entry parent) {
  if (root == nullptr) return;
  parent.children = nullptr;

  for (Entry entry : root->entries) {
    // If it is a leaf, we delete it directly.
    if (entry.children != nullptr) {
      freeMem(entry.children, entry);
    }

    delete root;
  }
}

/** The vector v has the accesses. */
void search(Node* root, Point q, double r, vector<int>& v, int index) {
  Entries entries = root->entries;
  if (entries.size() == 0) return;

  // Given that distance returns the squared euclidean distance.
  double r_square = pow(r, 2);
  // This is a external node, and else condition is internal.
  if (entries[0].children == nullptr) {
    v[index]++;
    /**
     * Verify for each entry if it satisfies (d(p, q))^2 <= r^2.
     * This is possible because square function increases strictly in >= 0.
     */
    for (Entry entry : entries) {
      Point p = entry.point;
      if (distance(p, q) <= r_square) {
        // ans.insert(p);
      }
    }
  } else {
    v[index]++;
    for (Entry entry : entries) {
      Point p = entry.point;
      double cr = entry.coveringRadius;
      double cr_square = pow(cr, 2);
      double binom_square = r_square + cr_square + 2 * r * cr;

      if (distance(p, q) <= binom_square) {
        // Search recursively in childs
        search(entry.children, q, r, v, index);
      }
    }
  }
}

Points parsePoints(int n) {
  string line;

  string fileName =
      n == N_ITER ? "queries.txt" : "pow2" + to_string(n) + ".txt";
  ifstream File("./input/" + fileName);

  Points result;
  if (File.is_open()) {
    while (getline(File, line)) {
      // Resets the index. The format is (x, y).
      int i = 0;
      string x, y;

      while (line[i] != '(') i++;
      i++;
      while (line[i] != ',') {
        x += line[i];
        i++;
      }

      i += 2;
      while (line[i] != ')') {
        y += line[i];
        i++;
      }

      result.insert({stod(x), stod(y)});
    }
    File.close();
  }

  return result;
}

int main() {
  // The argument is the power of 2, or N_ITER if query.
  Points testSet = parsePoints(10);
  Points testQueries = parsePoints(N_ITER);

  vector<Point> queryPoints;
  vector<pair<int, double>> response;
  vector<int> accesses(N_ITER);

  for (Point point : testQueries) {
    queryPoints.push_back(point);
  }

  auto start = high_resolution_clock::now();
  Node* root = SSAlgorithm(testSet);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);

  cout << "Construcción M-Tree (SS): " << duration.count() << " s\n";

  //  Frees the memory used in the M-Tree.
  // freeMem(root, {{NULL, NULL}, NULL, nullptr});

  for (int i = 0; i < N_ITER; i++) {
    Point q = queryPoints[i];

    auto start = high_resolution_clock::now();
    search(root, q, query_radius, accesses, i);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<nanoseconds>(stop - start);
    response.push_back({accesses[i], duration.count()});
  }

  for (int i = 0; i < N_ITER; i++) {
    cout << i + 1 << " (accesos, tiempo): (" << response[i].first << ", "
         << response[i].second << " ns)\n";
  }

  return 0;
}