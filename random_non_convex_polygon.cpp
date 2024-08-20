// Copyright 2024 Tier IV, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "autoware/universe_utils/geometry/boost_geometry.hpp"
#include "autoware/universe_utils/geometry/random_non_convex_polygon.hpp"
#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/algorithms/is_valid.hpp>
#include <boost/geometry/algorithms/is_simple.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <algorithm>
#include <random>
#include <vector>
#include <limits>
#include <cmath>
#include <thread>
#include <mutex>
#include <iostream>


namespace autoware::universe_utils
{
namespace
{
    // define Edge as a pair of Points
    struct Edge {
        Point2d first;
        Point2d second;
        bool valid = true;
    };

    struct PolygonWithEdges {
    Polygon2d polygon;
    std::vector<Edge> edges;
    };

    struct VectorsWithMin
    {
        std::vector<double> vectors;
        double min;
    };

    VectorsWithMin prepare_coordinate_vectors(
        const size_t nb_vertices, std::uniform_real_distribution<double>& random_double,
        std::uniform_int_distribution<int>& random_bool, std::default_random_engine& random_engine)
    {
        std::vector<double> v;
        for (auto i = 0UL; i < nb_vertices; ++i) {
            v.push_back(random_double(random_engine));
        }
        std::sort(v.begin(), v.end());
        const auto min_v = v.front();
        const auto max_v = v.back();
        std::vector<double> v1;
        v1.push_back(min_v);
        std::vector<double> v2;
        v2.push_back(min_v);
        for (auto i = 1UL; i + 1 < v.size(); ++i) {
            if (random_bool(random_engine) == 0) {
                v1.push_back((v[i]));
            } else {
                v2.push_back((v[i]));
            }
        }
        v1.push_back(max_v);
        v2.push_back(max_v);
        std::vector<double> diffs;
        for (auto i = 0UL; i + 1 < v1.size(); ++i) {
            diffs.push_back(v1[i + 1] - v1[i]);
        }
        for (auto i = 0UL; i + 1 < v2.size(); ++i) {
            diffs.push_back(v2[i] - v2[i + 1]);
        }
        VectorsWithMin vectors;
        vectors.vectors = diffs;
        vectors.min = min_v;
        return vectors;
    }

    double dist(const Edge& e, const Point2d& p)
    {
        double x = p.x();
        double y = p.y();
        double x1 = e.first.x();
        double y1 = e.first.y();
        double x2 = e.second.x();
        double y2 = e.second.y();

        double dx = x2 - x1;
        double dy = y2 - y1;

        if (dx == 0.0 && dy == 0.0) {
            dx = x - x1;
            dy = y - y1;
            return std::sqrt(dx * dx + dy * dy);
        }

        double t = ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy);
        t = std::max(0.0, std::min(1.0, t));
        dx = x - (x1 + t * dx);
        dy = y - (y1 + t * dy);
        return std::sqrt(dx * dx + dy * dy);
    }

    bool intersecting(const Edge& e, const Polygon2d& polygon)
    {
        Segment2d edge_segment{e.first, e.second};

        for (size_t i = 0; i < polygon.outer().size(); ++i)
        {
            const Point2d& p1 = polygon.outer()[i];
            const Point2d& p2 = polygon.outer()[(i + 1) % polygon.outer().size()];
            Segment2d poly_segment{p1, p2};

            if (boost::geometry::intersects(edge_segment, poly_segment))
            {
                return true;
            }
        }
        return false;
    }

bool isValid(const Edge& e, const Polygon2d& P, const std::vector<Point2d>& Q) {
    bool valid = false;
    size_t i = 0;
    
    while (!valid && i < Q.size()) {
        Edge e1 = {e.first, Q[i]};
        Edge e2 = {Q[i], e.second};

        // Debug output
        std::cout << "Checking edge from (" << e1.first.x() << ", " << e1.first.y() << ") to (" << e1.second.x() << ", " << e1.second.y() << ")" << std::endl;
        std::cout << "  vs (" << e2.first.x() << ", " << e2.first.y() << ") to (" << e2.second.x() << ", " << e2.second.y() << ")" << std::endl;

        bool intersects_e1 = intersecting(e1, P);
        bool intersects_e2 = intersecting(e2, P);

        std::cout << "  Intersection status: e1 = " << intersects_e1 << ", e2 = " << intersects_e2 << std::endl;

        if (!intersects_e1 && !intersects_e2) {
            valid = true;  // Edge is valid
            std::cout << "  Edge is valid" << std::endl;
        } else {
            std::cout << "  Edge is not valid" << std::endl;
        }
        
        ++i;
    }
    
    return valid;  // Edge is not valid if no valid configuration is found
}

Point2d getNearestNode(const std::vector<Point2d>& Q, const Edge& e) {
    double min_distance = std::numeric_limits<double>::max();  // Use a large initial value
    Point2d nearest_node(0, 0);  // Initialize to a default value

    std::cout << "Points in Q:" << std::endl;
    for (const auto& point : Q) {
        std::cout << "(" << point.x() << ", " << point.y() << ")" << std::endl;
    }

    for (const auto& node : Q) {
        double distance = dist(e, node);  // Calculate the distance from the node to the edge
        std::cout << "Distance from (" << node.x() << ", " << node.y() << ") to edge: " << distance << std::endl;

        if (distance < min_distance) {
            min_distance = distance;
            nearest_node = node;
        }
    }

    std::cout << "Nearest node found: (" << nearest_node.x() << ", " << nearest_node.y() << ")" << std::endl;
    return nearest_node;
}


    Edge getBreakingEdge(PolygonWithEdges& P_with_edges, const std::vector<Point2d>& Q) {
        double min_distance = std::numeric_limits<double>::max();
        Edge eBreaking;
        eBreaking.valid = false;

        for (const auto& edge : P_with_edges.edges) {
            if (edge.valid) {
                Point2d nearest_node = getNearestNode(Q, edge);
                double distance = dist(edge, nearest_node);
                if (distance < min_distance) {
                    min_distance = distance;
                    eBreaking = edge;
                }
            }
        }

        return eBreaking;
    }


    bool isValidInsertion(const LinearRing2d& new_ring) {
        for (size_t j = 0; j < new_ring.size(); ++j) {
            const Point2d& a = new_ring[j];
            const Point2d& b = new_ring[(j + 1) % new_ring.size()];
            const Point2d& c = new_ring[(j + 2) % new_ring.size()];

            double cross_product = (b.x() - a.x()) * (c.y() - b.y()) - (b.y() - a.y()) * (c.x() - b.x());
            if (cross_product <= 0) {
                return true;
            }
        }
        return false;
    }

    void updateEdges(PolygonWithEdges& P_with_edges) {
        auto& ring = P_with_edges.polygon.outer();
        P_with_edges.edges.clear();  // Clear the old edges
        for (size_t i = 0; i < ring.size(); ++i) {
            Edge e = {ring[i], ring[(i + 1) % ring.size()]};
            P_with_edges.edges.push_back(e);  // Add the new edges
        }
    }

    void insertNode(PolygonWithEdges& P_with_edges, const Point2d& w) {
        auto& ring = P_with_edges.polygon.outer();  // Get the outer ring of the polygon
        
        // Check if the point `w` already exists in the ring
        auto it = std::find_if(ring.begin(), ring.end(), [&](const Point2d& p) {
            return p.x() == w.x() && p.y() == w.y();
        });

        if (it != ring.end()) {
            // Point already exists, so skip insertion
            return;
        }

        // If not found, insert the new node `w`
        for (size_t i = 0; i < ring.size(); ++i) {
            LinearRing2d new_ring = ring;
            new_ring.insert(new_ring.begin() + (i + 1), w);  // Insert the new node `w`
            ring = std::move(new_ring);  // Update the polygon's ring with the new ring
            updateEdges(P_with_edges);  // Update the edges in `P_with_edges`
            return;
        }
    }

    

    void removeNode(std::vector<Point2d>& Q, const Point2d& w) {
        const double epsilon = 1e-9;  // Tolerance for comparing points

        Q.erase(std::remove_if(Q.begin(), Q.end(), [&](const Point2d& p) {
            return std::abs(p.x() - w.x()) < epsilon && std::abs(p.y() - w.y()) < epsilon;
        }), Q.end());
    }

    void markValidEdges(PolygonWithEdges& P_with_edges, const std::vector<Point2d>& Q) {
        for (auto& edge : P_with_edges.edges) {
            edge.valid = isValid(edge, P_with_edges.polygon, Q);
        }
    }

    Polygon2d inwardDenting(Polygon2d& polygon)
    {
        Polygon2d P = polygon;

        // Compute convex hull and initialize the point sets
        const auto& outer_ring = P.outer();
        std::vector<Point2d> points(outer_ring.begin(), outer_ring.end());
        std::vector<Point2d> Q(points.begin(), points.end());

        PolygonWithEdges P_with_edges;
        P_with_edges.polygon = P;
        P_with_edges.edges.resize(P_with_edges.polygon.outer().size());

        // Populate initial edges
        for (size_t i = 0; i < P_with_edges.edges.size(); ++i) {
            P_with_edges.edges[i] = {P_with_edges.polygon.outer()[i], 
                                    P_with_edges.polygon.outer()[(i + 1) % P_with_edges.polygon.outer().size()]};
        }
        // Iterative inward denting
        while (!Q.empty()) {
            Edge e = getBreakingEdge(P_with_edges, Q);
            Point2d w = getNearestNode(Q, e);
            insertNode(P_with_edges, w);
            std::cout << "(" << w.x() << ", " << w.y() << ")" << std::endl;
            removeNode(Q, w);
            markValidEdges(P_with_edges, Q);
        }
        return P_with_edges.polygon;
}

}

bool is_convex(const autoware::universe_utils::Polygon2d & polygon)
{
  const auto & outer_ring = polygon.outer();
  size_t num_points = outer_ring.size();

  if (num_points < 4) {
    return true; 
  }

  bool is_positive = false;
  bool is_negative = false;

  for (size_t i = 0; i < num_points; ++i) {
    auto p1 = outer_ring[i];
    auto p2 = outer_ring[(i + 1) % num_points];
    auto p3 = outer_ring[(i + 2) % num_points];

    double cross_product =
      (p2.x() - p1.x()) * (p3.y() - p2.y()) - (p2.y() - p1.y()) * (p3.x() - p2.x());

    if (cross_product > 0) {
      is_positive = true;
    } else if (cross_product < 0) {
      is_negative = true;
    }

    if (is_positive && is_negative) {
      return false;  
    }
  }

  return true; 
}

const double TOLERANCE = 1e-6;
// Function to check if two points are approximately equal
bool arePointsEqual(const autoware::universe_utils::Point2d& p1, const autoware::universe_utils::Point2d& p2) {
    return std::fabs(p1.x() - p2.x()) < TOLERANCE && std::fabs(p1.y() - p2.y()) < TOLERANCE;
}

// Function to remove duplicate points from a polygon
void removeDuplicatePoints(autoware::universe_utils::Polygon2d& polygon) {
    std::vector<autoware::universe_utils::Point2d> uniquePoints;
    
    for (const auto& point : polygon.outer()) {
        // Check if the point is already in the uniquePoints vector
        auto it = std::find_if(uniquePoints.begin(), uniquePoints.end(),
            [&](const autoware::universe_utils::Point2d& p) { return arePointsEqual(p, point); });
        
        if (it == uniquePoints.end()) {
            uniquePoints.push_back(point);
        }
    }

    // Replace the polygon's points with the unique points
    polygon.outer().clear();
    polygon.outer().assign(uniquePoints.begin(), uniquePoints.end());
}



Polygon2d random_non_convex_polygon(const size_t vertices, const double max)
{
    if (vertices < 4) {
        return Polygon2d();
    }

    std::random_device r;
    std::default_random_engine random_engine(r());
    std::uniform_real_distribution<double> uniform_dist(-max, max);
    std::uniform_int_distribution<int> random_bool(0, 1);

    Polygon2d poly;
    bool is_non_convex = false;
    int max_attempts = 100;
    int attempt = 0;

    while (!is_non_convex && attempt < max_attempts) {
        auto xs = prepare_coordinate_vectors(vertices, uniform_dist, random_bool, random_engine);
        auto ys = prepare_coordinate_vectors(vertices, uniform_dist, random_bool, random_engine);
        std::shuffle(ys.vectors.begin(), ys.vectors.end(), random_engine);

        LinearRing2d vectors;
        for (size_t i = 0; i < xs.vectors.size(); ++i) {
            vectors.emplace_back(xs.vectors[i], ys.vectors[i]);
        }

        std::sort(vectors.begin(), vectors.end(), [](const Point2d& p1, const Point2d& p2) {
            return std::atan2(p1.y(), p1.x()) < std::atan2(p2.y(), p2.x());
        });
        auto min_x = max;
        auto min_y = max;
        auto x = 0.0;
        auto y = 0.0;
        LinearRing2d points;
        for (const auto & p : vectors) {
            points.emplace_back(x, y);
            x += p.x();
            y += p.y();
            min_x = std::min(p.x(), min_x);
            min_y = std::min(p.y(), min_y);
        }
        const auto shift_x = min_x - xs.min;
        const auto shift_y = min_y - ys.min;
        for (auto & p : points) {
            p.x() += shift_x;
            p.y() += shift_y;
        }

        LinearRing2d poly_points;
        for (const auto& p : points) {
            poly_points.emplace_back(p.x(), p.y());
        }

        for (const auto& triangle : poly_points) {
        std::cout << "before: " << triangle << std::endl;
        }

        poly.outer() = poly_points;
        // apply inward denting algorithm
        poly = inwardDenting(poly);
        // check for convexity
        if (!is_convex(poly)) {
            is_non_convex = true;
        }

        // Ensure all points are positive after inward denting
        double min_x_final = std::numeric_limits<double>::max();
        double min_y_final = std::numeric_limits<double>::max();
        for (const auto& p : poly.outer()) {
            min_x_final = std::min(p.x(), min_x_final);
            min_y_final = std::min(p.y(), min_y_final);
        }

        if (min_x_final < 0 || min_y_final < 0) {
            double shift_x_final = std::max(0.0, -min_x_final);
            double shift_y_final = std::max(0.0, -min_y_final);
            for (auto& p : poly.outer()) {
                p.x() += shift_x_final;
                p.y() += shift_y_final;
            }
        }

        for (const auto& triangle : poly.outer()) {
        std::cout << "after: " << triangle << std::endl;
        }

        boost::geometry::correct(poly);


        ++attempt;
    }
    return poly;
}
}  // namespace autoware::universe_utils

