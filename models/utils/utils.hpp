////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILS_H
#define UTILS_H

#include "utilLibs.hpp"
#include "utilIO.hpp"

#define X_POSITIVE_BOUNDARY 0
#define X_NEGATIVE_BOUNDARY 1
#define Y_POSITIVE_BOUNDARY 2
#define Y_NEGATIVE_BOUNDARY 3
#define Z_POSITIVE_BOUNDARY 4
#define Z_NEGATIVE_BOUNDARY 5

#define OUTER_BOUNDARY_NODE 0
#define INNER_BOUNDARY_NODE 1
#define INNER_NETWORK_NODE 2

namespace util {

/*!
 * @brief Return the Square of numbers and ignores values lower than 0 and
 * higher than 1
 *
 * @param x  Number
 * @return value Value between 0 and 1
 */
double square(double x);

/*!
 * @brief Ignores values lower than 0 and higher than 1
 *
 * @param x  Number
 * @return value Value between 0 and 1
 */
double linear(double x);

/*!
 * @brief Regularized Heaviside function
 *
 * @param x  Number
 * @return value Heaveside function
 */
double heaviside(double x);

/*!
 * @brief Modulus of a vector (norm)
 *
 * @param xg X coordinate
 * @param yg Y coordinate
 * @return value Norm
 */
double modulus(double xg, double yg);

/*!
 * @brief Determinant of matrix
 *
 * @param M Matrix
 * @return value Determinant of M
 */
double determinant(const std::vector<std::vector<double>> &M);

/*!
 * @brief Inverse of matrix
 *
 * @param M Matrix
 * @return Minv Inverse of matrix M
 */
std::vector<std::vector<double>> inverse(const
std::vector<std::vector<double>> &M);

/*!
 * @brief Transpose of matrix
 *
 * @param M Matrix
 * @return Mt Transpose of matrix M
 */
std::vector<std::vector<double>> transpose(const
std::vector<std::vector<double>> &M);

/*!
 * @brief Add element to list
 *
 * @param i Element to be added
 * @param data List
 * @return true True if i is added to list
 */
bool addToList(unsigned int i, std::vector< unsigned int> *data);

/*!
 * @brief compute the mass of the qoi
 *
 * @param es Equation system
 * @param system_name Name of system
 * @param var_name Name of variable
 * @param value_mass Value of mass
 */
void computeMass(EquationSystems &es, const std::string &system_name,
                 const std::string &var_name, double &value_mass);

/*! @brief Position index of given point in 2d
   * @param x X coordinate of point
   * @param y Y coordinate of point
   * @param Nelx Number of elements in along x-axis (and y-axis)
   * @return Index Unique integer corresponding to the point
 */
int positionIndex(double x, double y, int Nelx);

/*!
 * @brief Projects concentration to physical range i.e. [0,1]
 *
 * @param x Number
 * @return value Value between 0 and 1
 */
double project_concentration(double x);

/*!
 * @brief Prints message to std::cout and also writes to the file
 */
class Logger {

public:

  Logger() : d_comm_p(nullptr), d_screen_out(true) {};

  Logger(const std::string &log_file, Parallel::Communicator *comm, bool
                                                                        screen_out = true)
      : d_comm_p(comm), d_screen_out(screen_out) {

    d_dbg_file.open(log_file, std::ios_base::out);
  };

  // overload function call
  void operator()(std::ostringstream &oss, bool screen_out) {
    log(oss, screen_out);
  }
  void operator()(std::ostringstream &oss) {
    log(oss, d_screen_out);
  }

  void operator()(const std::string &str, bool screen_out) {
    log(str, screen_out);
  }
  void operator()(const std::string &str) {
    log(str, d_screen_out);
  }

  void log(std::ostringstream &oss, bool screen_out) {

    if (screen_out)
      out << oss.str();
    if (d_comm_p->rank() == 0)
      d_dbg_file << oss.str();

    oss.str("");
    oss.clear();
  };

  void log(std::ostringstream &oss) {
    log(oss, d_screen_out);
  };

  void log(const std::string &str, bool screen_out) {

    if (screen_out)
      out << str;
    if (d_comm_p->rank() == 0)
      d_dbg_file << str;
  };

  void log(const std::string &str) {

    log(str, d_screen_out);
  };

  ~Logger() { d_dbg_file.close(); }

  void init(const std::string &log_file, Parallel::Communicator *comm) {
    d_comm_p = comm;
    d_dbg_file.open(log_file, std::ios_base::out);
  }

private:
  Parallel::Communicator *d_comm_p;
  std::ofstream d_dbg_file;
  bool d_screen_out{};
};

/*!
 * @brief Find node closer to given point
 *
 * @param x Point
 * @param mesh Mesh
 * @return id Id of node
 */
unsigned int locate_node(const Point &x, const MeshBase &mesh);

/*!
 * @brief Returns point in line formed by points p1 and p2
 *
 * @param p1 Point 1
 * @param p2 Point 2
 * @param s Parametric coordinate
 * @return p Point
 */
Point point_on_line(const Point &p1, const Point &p2, const double &s);

/*!
 * @brief Returns direction from point 1 to 2
 *
 * @param p1 Point 1
 * @param p2 Point 2
 * @return p Unit vector pointing at point 2 from point 1
 */
Point get_direction(const Point &p1, const Point &p2);

/*!
 * @brief Returns true if point is inside the box
 *
 * @param p Point
 * @param box Box
 * @return True If point inside box otherwise false
 */
bool is_inside_box(const Point &p, const std::pair<Point, Point> &box, double
tol = 0.);
bool is_inside_box(const Point &p, const double &box_size, double
tol = 0.);

/*!
 * @brief Returns true if point is inside the cylinder
 *
 * @param p Point
 * @param length Length of cylinder
 * @param radius Radius of cylinder
 * @param axis Axis of cylinder
 * @return True If point inside cylinder otherwise false
 */
bool is_inside_cylinder(const Point &p, const double &length, const double
&radius, const Point &axis);

/*!
 * @brief Returns true if point is inside the cylinder
 *
 * @param p Point
 * @param R Radius of cylinder
 * @param x1 Point at the center of cross-section at s=0
 * @param x2 Point at the center of cross-section at s=L
 * @param
 * @return True If point inside cylinder otherwise false
 */
bool is_inside_cylinder(const Point &p, const double &radius, const Point &x1,
                        const Point &x2);

/*!
 * @brief Returns true if point is inside the ellipsoid
 *
 * @param p Point
 * @param center Center of ellipse
 * @param radius_vec Vector of radius describing ellipse
 * @param dim Dimension
 * @param d
 * @return True If point inside otherwise false
 */
bool is_inside_ellipse(const Point &p, const Point &center, const
std::vector<double> &radius_vec, unsigned int dim);

/*!
 * @brief Returns true if point is inside the ellipsoid
 *
 * Also computes
 * d = x^2 / r1^2 + y^2 / r2^2 + z^2 / r3^2
 *
 * @param p Point
 * @param center Center of ellipse
 * @param radius_vec Vector of radius describing ellipse
 * @param dim Dimension
 * @param d
 * @return True If point inside otherwise false
 */
bool is_inside_ellipse(const Point &p, const Point &center, const
std::vector<double> &radius_vec, unsigned int dim, double &d);

/*!
 * @brief Transforms point in ellipse to point in ball of given radius
 *
 * @param p Point
 * @param center Center of ellipse
 * @param radius_vec Vector of radius describing ellipse
 * @param dim Dimension
 * @param ball_r Radius of ball
 * @return Point Point in ball
 */
Point ellipse_to_ball(const Point &p, const Point &center, const
std::vector<double> &radius_vec, unsigned int dim, const double &ball_r);

/*!
 * @brief Returns the vector after rotating by desired angle
 *
 * @param p Vector
 * @param theta Angle of rotation
 * @param axis Axis of rotation
 */
Point rotate(const Point &p, const double &theta, const Point &axis);

/*!
 * @brief Returns the unit vector corresponding to cross product
 *
 * Note: Better to use inbuilt function in Point - p1.cross(p2)
 *
 * @param p1 Vector 1
 * @param p2 Vector 2
 * @return vector Unit vector along the cross product
 */
Point cross_product(const Point &p1, const Point &p2);

/*!
 * @brief Computes angle between two vectors
 * @param a Vector 1
 * @param b Vector 2
 */
double angle(Point a, Point b);

/*!
 * @brief Computes angle between two vectors
 * @param a Vector 1
 * @param b Vector 2
 * @param axis Axis of rotation
 * @param is_axis If true then axis is the axis of orientation, otherwise
 * axis specifies the +ve side of the plane in which a and b are
 */
double angle(Point a, Point b, Point axis, bool is_axis = true);

/*!
 * @brief Computes following function
 *
 * exp[1 - 1 / (1 - r^a)]
 *
 * @param r Argument of function
 * @param exponent Value of power a in above function
 */
double exp_decay_function(double r, double exponent = 4.);

std::vector<Point> discretize_ball_surface(const unsigned int
                                                 &disc_num, const double
                                                 &ball_r, unsigned int dim);

std::vector<Point> discretize_cube(const unsigned int
                                           &disc_num, const double
                                           &cube_size, unsigned int dim);


/*!
 * @brief Do lines intersect
 *
 * @param line_1 Line 1
 * @param line_2 Line 2
 * @return True If lines intersect, else false
 */
bool lines_intersect(const std::pair<Point, Point> &line_1, const std::pair<Point,
    Point> &line_2);

/*!
 * @brief Compute distance between lines
 *
 * @param line_1 Line 1
 * @param line_2 Line 2
 * @return Value Distance
 */
double distance_between_lines(const std::pair<Point, Point> &line_1,
                              const std::pair<Point, Point> &line_2);
double distance_between_segments(const std::pair<Point, Point> &line_1,
                              const std::pair<Point, Point> &line_2);

/*!
 * @brief Compute distance between planes
 *
 * @param plane_1 Plane 1 given by pair of normal and one point which it
 * contains
 * @param plane_2 Plane 2 given by pair of normal and one point which it
 * contains
 * @return Value Distance
 */
double distance_between_planes(const std::pair<Point, Point> &plane_1,
                              const std::pair<Point, Point> &plane_2);

/*!
 * @brief Compute distance between point and line
 *
 * @param p Point
 * @param line Line
 * @return Value Distance
 */
double point_distance_line(const Point &p,
                              const std::pair<Point, Point> &line);
double point_distance_segment(const Point &p,
                           const std::pair<Point, Point> &line);

/*!
 * @brief Compute distance between point and plane
 *
 * @param p Point
 * @param plane Plane given by pair of normal and one point which it
 * contains
 * @return Value Distance
 */
double point_distance_plane(const Point &p,
                           const std::pair<Point, Point> &plane);

/*!
 * @brief Get list of elements which have this node as vertex
 *
 * @param node_id Node id
 * @param mesh Mesh
 * @return List
 */
std::vector<unsigned int> find_elems(const unsigned int node_id,
                                     const MeshBase &mesh);

/*!
 * @brief Get list of elements which have this node as vertex
 *
 * @param node_id Node id
 * @param mesh Mesh
 * @return List
 */
void add_unique(unsigned int dof, Real val,
                std::vector<unsigned int> &list, std::vector<Real>
                &list_val);

/*!
 * @brief Get the major axis which is parallel to given vector within some
 * tolerance
 *
 * @param vec Given vector
 * @return axis Parallel to vec
 */
std::string get_vec_major_axis(Point vec);

/*!
 * @brief Get the length between two points
 *
 * @param p1
 * @param p2
 * @return distance
 */
double dist_between_points(const std::vector<double> &p1, const
std::vector<double> &p2);

void get_unit_vector(std::vector<double> &p1);

unsigned int get_elem_id(const std::vector<double> &p, const double &mesh_size,
                         const unsigned int &num_elems, const unsigned int &dim);
unsigned int get_elem_id(const Point &p, const double &mesh_size,
                         const unsigned int &num_elems, const unsigned int &dim);

Point to_point(const std::vector<double> &p);

Point determineRotator( const Point & dir );
std::vector<double> determineRotator( const std::vector<double> &dir );
} // namespace util

#endif // UTILS_H
