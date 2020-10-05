////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "utils.hpp"
#include <math.h>

double util::square(double x) {
  if (x < 0) {
    return 0.;
  } else if (x > 1) {
    return 1.;
  } else
    return x * x;
}

double util::linear(double x) {
  if (x < 0) {
    return 0.;
  } else if (x > 1) {
    return 1.;
  } else
    return x;
}

//double util::heaviside(double x) {
//  if (fabs(x) >= 0.064) {
//    return (x < 0) ? 0. : 1.;
//  } else {
//    return (1. / (2. * 0.064)) * (x + 0.064);
//  }
//}

double util::heaviside(double x) {
  if (x < 1.0E-10)
    return 0.;
  else {
    if (x < 0.0001 - 1.0E-10)
      return 1. - std::exp(1. - 1. / (1. - std::pow(x / 0.0001, 4)));
    else
      return 1.;
  }
}

double util::modulus(double xg, double yg) {
  return sqrt(pow(xg, 2) + pow(yg, 2));
}

double util::determinant(const std::vector<std::vector<double>> &M) {

  if (M.size() == 2)
    return M[0][0] * M[1][1] - M[1][0] * M[0][1];
  else if (M.size() == 3)
    return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
           M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
           M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  else {
    std::cerr << "Error: Determinant function only supports square matrix of "
                 "dimension 2 and 3.\n";
    exit(1);
  }
}

std::vector<std::vector<double>>
util::inverse(const std::vector<std::vector<double>> &M) {

  std::vector<std::vector<double>> Minv;
  for (size_t i = 0; i < M.size(); i++)
    Minv.emplace_back(M.size(), 0.);

  double detM = util::determinant(M);

  if (M.size() == 2) {

    Minv[0][0] = M[1][1] / detM;
    Minv[0][1] = -M[0][1] / detM;
    Minv[1][0] = -M[1][0] / detM;
    Minv[1][1] = M[0][0] / detM;

    return Minv;
  } else if (M.size() == 3) {

    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++)
        Minv[j][i] =
          (M[(i + 1) % 3][(j + 1) % 3] * M[(i + 2) % 3][(j + 2) % 3] -
           M[(i + 1) % 3][(j + 2) % 3] * M[(i + 2) % 3][(j + 1) % 3]) /
          detM;
    }

    return Minv;
  } else {
    std::cerr << "Error: Inverse function only supports square matrix of "
                 "dimension 2 and 3.\n";
    exit(1);
  }
}

std::vector<std::vector<double>>
util::transpose(const std::vector<std::vector<double>> &M) {

  std::vector<std::vector<double>> Mt;
  for (size_t i = 0; i < M.size(); i++)
    Mt.emplace_back(M.size(), 0.);

  for (size_t i = 0; i < M.size(); i++) {
    for (size_t j = 0; j < M.size(); j++)
      Mt[i][j] = M[j][i];
  }

  return Mt;
}

bool util::addToList(unsigned int i, std::vector<unsigned int> *data) {
  for (auto d : *data)
    if (i == d)
      return false;

  data->push_back(i);
  return true;
}

void util::computeMass(EquationSystems &es, const std::string &system_name,
                       const std::string &var_name, double &value_mass) {

  value_mass = 0.;

  TransientLinearImplicitSystem &sys =
    es.get_system<TransientLinearImplicitSystem>(system_name);

  const unsigned int var = sys.variable_number(var_name);
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_type = sys.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const DofMap &sys_map = sys.get_dof_map();
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_var;
  Number qoi_mass = 0.;
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();
  for (; el != end_el; ++el) {
    const Elem *elem = *el;
    sys_map.dof_indices(elem, dof_indices);
    sys_map.dof_indices(elem, dof_indices_var, var);
    fe->reinit(elem);
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++) {
      Number partial_var_mass = 0.;
      for (unsigned int l = 0; l < phi.size(); l++) {
        double var_val = sys.current_solution(dof_indices_var[l]);
        if (var_val > 1.0)
          var_val = 1.0;
        else if (var_val < 0.0)
          var_val = 0.0;
        partial_var_mass += JxW[qp] * phi[l][qp] * var_val;
      }
      qoi_mass += partial_var_mass;
    }
  }
  value_mass = qoi_mass;
}

int util::positionIndex(double x, double y, int Nelx) {

  return std::floor(y + 0.5) * (Nelx + 1) + std::floor(x + 0.5);
}

double util::project_concentration(double x) {

  auto y = 1. < x ? 1. : x;
  return 0. > y ? 0. : y;
}

unsigned int util::locate_node(const Point &x, const MeshBase &mesh) {

  if (mesh.n_nodes() == 0)
    return 0;

  long int locate = -1;
  double dist = (x - mesh.point(0)).norm();
  for (const auto &node : mesh.node_ptr_range()) {

    auto id = node->id();
    Point dx = x - mesh.point(id);
    if (dx.norm() < dist) {
      locate = id;
      dist = dx.norm();
    }
  }

  if (locate == -1) {
    out << "x = (" << x(0) << ", " << x(1) << ", " << x(2) << ")\n";
    libmesh_error_msg("Can not find node closer to point");
    exit(1);
  }

  return locate;
}

Point util::point_on_line(const Point &p1, const Point &p2, const double &s) {

  return Point(p1(0) + s * (p2(0) - p1(0)), p1(1) + s * (p2(1) - p1(1)),
               p1(2) + s * (p2(2) - p1(2)));
}

Point util::get_direction(const Point &p1, const Point &p2) {

  Point dp = p2 - p1;
  return dp / dp.norm();
}

bool util::is_inside_box(const Point &p, const std::pair<Point, Point> &box,
                         double tol) {

  const double eps = 1.0E-10;
  return !(
    p(0) < box.first(0) + tol - eps || p(0) > box.second(0) - tol + eps ||
    p(1) < box.first(1) + tol - eps || p(1) > box.second(1) - tol + eps ||
    p(2) < box.first(2) + tol - eps || p(2) > box.second(2) - tol + eps);
}

bool util::is_inside_box(const Point &p, const double &box_size, double tol) {

  const double eps = 1.0E-10;
  return !(p(0) < tol - eps || p(0) > box_size - tol + eps ||
           p(1) < tol - eps || p(1) > box_size - tol + eps ||
           p(2) < tol - eps || p(2) > box_size - tol + eps);
}

bool util::is_inside_cylinder(const Point &p, const double &length,
                              const double &radius, const Point &axis) {

  double p_dot_a = p * axis;
  if (p_dot_a > length or p_dot_a < 0.)
    return false;
  else {

    Point p_parallel = p - p_dot_a * axis;

    return p_parallel.norm_sq() < radius * radius;
  }
}

bool util::is_inside_cylinder(const Point &p, const double &radius,
                              const Point &x1, const Point &x2) {

  Point p_new = p - x1;
  Point a = x2 - x1;
  double p_dot_a = p_new * a;

  // note here we should 1 if a is not normalized
  if (p_dot_a > 1. or p_dot_a < 0.)
    return false;
  else {

    Point p_parallel = p_new - p_dot_a * a;

    return p_parallel.norm_sq() < radius * radius;
  }
}

bool util::is_inside_ellipse(const Point &p, const Point &center,
                             const std::vector<double> &radius_vec,
                             unsigned int dim) {

  double d = 0.;
  Point x = p - center;
  for (unsigned int i = 0; i < dim; i++)
    d += x(i) * x(i) / (radius_vec[i] * radius_vec[i]);

  return d < 1.;
}

bool util::is_inside_ellipse(const Point &p, const Point &center,
                             const std::vector<double> &radius_vec,
                             unsigned int dim, double &d) {

  d = 0.;
  Point x = p - center;
  for (unsigned int i = 0; i < dim; i++)
    d += x(i) * x(i) / (radius_vec[i] * radius_vec[i]);

  return d < 1.;
}

Point util::ellipse_to_ball(const Point &p, const Point &center,
                            const std::vector<double> &radius_vec,
                            unsigned int dim, const double &ball_r) {

  Point x = p - center;
  Point xb = Point();
  for (unsigned int i = 0; i < dim; i++)
    xb(i) = ball_r * x(0) / radius_vec[i];

  return xb;
}

Point util::rotate(const Point &p, const double &theta, const Point &axis) {

  auto ct = std::cos(theta);
  auto st = std::sin(theta);

  // dot
  double p_dot_n = p * axis;

  // cross
  Point n_cross_p = axis.cross(p);

  return (1. - ct) * p_dot_n * axis + ct * p + st * n_cross_p;
}

std::vector<double> util::rotate(std::vector<double> &p, double theta, std::vector<double> &axis) {

  std::vector<double> rotated_point = std::vector<double>(3, 0.0);

  auto ct = std::cos(theta);
  auto st = std::sin(theta);

  double p_dot_n = 0.0;

  for (int i = 0; i < 3; i++) {

    p_dot_n = p_dot_n + p[i] * axis[i];
  }

  std::vector<double> n_cross_p = cross_prod(axis, p);

  for (int i = 0; i < 3; i++) {

    rotated_point[i] = (1. - ct) * p_dot_n * axis[i] + ct * p[i] + st * n_cross_p[i];
  }

  return rotated_point;
}

Point util::cross_product(const Point &p1, const Point &p2) {

  //  auto p = Point(p1(1)*p2(2)-p1(2)*p2(1),
  //               p1(2)*p2(0)-p1(0)*p2(2),
  //               p1(0)*p2(1)-p1(1)*p2(0));

  Point p = p1.cross(p2);

  if (p.norm() > 1.0E-8)
    return p / p.norm();
  else
    return Point();
}

std::vector<double> util::cross_prod(std::vector<double> &p1, std::vector<double> &p2) {

  std::vector<double> c_product = std::vector<double>(3, 0.0);

  c_product[0] = p1[1] * p2[2] - p1[2] * p2[1];
  c_product[1] = p1[2] * p2[0] - p1[0] * p2[2];
  c_product[2] = p1[0] * p2[1] - p1[1] * p2[0];

  double norm_c_product = std::sqrt(c_product[0] * c_product[0] + c_product[1] * c_product[1] + c_product[2] * c_product[2]);

  if (norm_c_product > 0.0) {

    for (int i = 0; i < 3; i++) {

      c_product[i] = c_product[i] / norm_c_product;
    }

  } else {

    for (int i = 0; i < 3; i++) {

      c_product[i] = 0.0;
    }
  }

  return c_product;
}

double util::angle(Point a, Point b) {

  if ((a - b).norm_sq() < 1.0E-12)
    return 0.;

  // since we do not know which side of plane given by normal
  // a x b / |a x b| is +ve, we compute the angle using cosine
  return std::acos(b * a / (b.norm() * a.norm()));
}

double util::angle(Point a, Point b, Point axis, bool is_axis) {

  if ((a - b).norm_sq() < 1.0E-12)
    return 0.;

  if (is_axis) {

    // normal to plane of rotation
    Point n = axis / axis.norm();

    Point na = n.cross(a);

    double theta = std::atan(b * na / (a * b - (b * n) * (a * n)));
    if (theta < 0.)
      theta += M_PI;

    if (b * na < 0.)
      theta = M_PI + theta;

    return theta;
  } else {

    auto theta = angle(a, b);

    // TODO below only works in specific cases such as when vectors in xy
    //  plane and vector x gives the positive plane direction, i.e. whether
    //  (0, 0, 1) is +ve or (0, 0, -1) is +ve. Similar is true for yz, zx
    //  planes.

    // normal to a and b
    Point n_ab = a.cross(b);

    double orient = axis * n_ab;
    if (orient < 0.)
      return 2. * M_PI - theta;
    else
      return theta;
  }
}

double util::exp_decay_function(double r, double exponent) {

  return std::exp(1. - 1. / (1. - std::pow(std::abs(r), exponent)));
}

std::vector<Point> util::discretize_ball_surface(const unsigned int &disc_num,
                                                 const double &ball_r,
                                                 unsigned int dim) {

  std::vector<Point> ball_points;
  // discretization along each axis
  int N = disc_num;
  unsigned int counter = 0;
  for (int i = -N; i <= N; i++)
    for (int j = -N; j <= N; j++)
      for (int k = -N; k <= N; k++) {

        if (dim == 2 && k != 0)
          continue;

        // we want points on the face of cube
        // +x
        bool skip_x = true;
        if (i == N or i == -N)
          skip_x = false;

        if ((j == N or j == -N) or !skip_x)
          skip_x = false;

        if (dim == 3 and ((k == N or k == -N) or !skip_x))
          skip_x = false;

        if (skip_x)
          continue;

        Point x = Point(double(i), double(j), double(k)) * 1. / double(N);

        // transform x to point in unit ball
        Point y = ball_r * x / x.norm();

        ball_points.push_back(y);
      }

  return ball_points;
}

std::vector<Point> util::discretize_cube(const unsigned int &disc_num,
                                         const double &cube_size,
                                         unsigned int dim) {

  std::vector<Point> points;

  // discretization along each axis
  int N = disc_num;
  unsigned int counter = 0;
  for (int i = -N; i <= N; i++)
    for (int j = -N; j <= N; j++)
      for (int k = -N; k <= N; k++) {

        if (dim == 2 && k != 0)
          continue;

        Point x =
          Point(double(i), double(j), double(k)) * cube_size / double(N);

        points.push_back(x);
      }

  return points;
}

bool util::lines_intersect(const std::pair<Point, Point> &line_1,
                           const std::pair<Point, Point> &line_2) {

  // change of variable so that first point of line_1 is at origin
  // After change of variables:
  // a is the second point of line_1
  // b is the first point of line_2
  // c is the difference of second and first point of line_2
  Point a = line_1.second - line_1.first;
  Point b = line_2.first - line_1.first;
  Point c = line_2.second - line_2.first;

  // check if the two lines are parallel
  if (util::angle(a / a.norm(), c / c.norm()) < 1.0E-8)
    return false;

  double a_dot_a = a.norm_sq();
  double a_dot_b = a * b;
  double a_dot_c = a * c;
  double b_dot_c = b * c;
  double c_dot_c = c.norm_sq();

  double r = (a_dot_a * b_dot_c - a_dot_b * a_dot_c) /
             (a_dot_c * a_dot_c - c_dot_c * a_dot_a);

  // if r is in (0,1) then this gives the intersection point
  // otherwise b + r c gives the point where two vectors originating from
  // a and originating from b would intersect
  return r > 0. and r < 1.;
}

double util::distance_between_lines(const std::pair<Point, Point> &line_1,
                                    const std::pair<Point, Point> &line_2) {

  // let line 1 is l1(s) = p + s u
  // and line 2 is l2(r) = q + r v
  // and let normal to the plane containing line 1 and 2 is
  // n = u x v / |u x v|

  Point u = line_1.second - line_1.first;
  Point v = line_2.second - line_2.first;
  Point w0 = line_1.first - line_2.first;

  double a = u * u;
  double b = u * v;
  double c = v * v;
  double d = u * w0;
  double e = v * w0;

  Point dp = w0 + ((b * e - c * d) * u + (a * e - b * d) * v) / (a * c - b * b);

  return dp.norm();
}

double util::distance_between_segments(const std::pair<Point, Point> &line_1,
                                       const std::pair<Point, Point> &line_2) {

  // let line 1 is l1(s) = p + s u
  // and line 2 is l2(r) = q + r v
  // and let normal to the plane containing line 1 and 2 is
  // n = u x v / |u x v|

  Point u = line_1.second - line_1.first;
  Point v = line_2.second - line_2.first;
  Point w0 = line_1.first - line_2.first;

  double a = u * u;
  double b = u * v;
  double c = v * v;
  double d = u * w0;
  double e = v * w0;
  double D = a * c - b * b;
  double sc, sN, sD = D;
  double tc, tN, tD = D;

  // compute line parameters of two closest points
  if (D < 1.0E-12) {

    sN = 0.;
    sD = 1.;
    tN = e;
    tD = c;
  } else {

    sN = b * e - c * d;
    tN = a * e - b * d;

    if (sN < 0.) {
      sN = 0.;
      tN = e;
      tD = c;
    } else if (sN > sD) {
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }

  if (tN < 0.) {

    tN = 0.;

    if (-d < 0.)
      sN = 0.;
    else if (-d > a)
      sN = sD;
    else {
      sN = -d;
      sD = a;
    }
  } else if (tN > tD) {

    tN = tD;

    if (-d + b < 0.)
      sN = 0.;
    else if (-d + b > a)
      sN = sD;
    else {
      sN = -d + b;
      sD = a;
    }
  }

  sc = abs(sN) < 1.0E-12 ? 0. : sN / sD;
  tc = abs(tN) < 1.0E-12 ? 0. : tN / tD;

  Point dp = w0 + sc * u - tc * v;

  return dp.norm();
}

double util::distance_between_planes(const std::pair<Point, Point> &plane_1,
                                     const std::pair<Point, Point> &plane_2) {

  // check if planes are parallel
  if (util::angle(plane_1.first, plane_2.first) < 1.0E-8)
    return 0.;

  return std::abs(plane_1.first * (plane_1.second - plane_2.second)) /
         plane_1.first.norm();
}

double util::point_distance_line(const Point &p,
                                 const std::pair<Point, Point> &line) {

  // line vector
  Point v = line.second - line.first;

  // vector from 1st point of line to p
  Point w = p - line.first;

  // project w onto v and add 1st point to get projected point's location
  Point w_on_line = line.first + (w * v) * v / v.norm_sq();

  return (p - w_on_line).norm();
}

double util::point_distance_segment(const Point &p,
                                    const std::pair<Point, Point> &line) {

  // line vector
  Point v = line.second - line.first;

  // vector from 1st point of line to p
  Point w = p - line.first;

  // determine if w is on left side or right side of line
  double w_dot_v = w * v;
  if (w_dot_v < 1.0E-12)
    return (p - line.first).norm();

  if (w_dot_v > v.norm_sq() - 1.0E-12)
    return (p - line.second).norm();

  // project w onto v and add 1st point to get projected point's location
  Point w_on_line = line.first + w_dot_v * v / v.norm_sq();

  return (p - w_on_line).norm();
}

double util::point_distance_plane(const Point &p,
                                  const std::pair<Point, Point> &plane) {

  // if plane is given by unit normal n and a point a which is contained in it
  // then the distance of point p from plane is
  // |(p - a) dot n| / |n|

  Point pa = p - plane.second;
  return std::abs(pa * plane.first) / plane.first.norm();
}

std::vector<unsigned int> util::find_elems(const unsigned int node_id,
                                           const MeshBase &mesh) {

  std::vector<unsigned int> list;
  for (const auto *elem : mesh.element_ptr_range())
    for (unsigned int i = 0; i < elem->n_nodes(); i++)
      if (elem->node_id(i) == node_id)
        list.emplace_back(elem->id());

  return list;
}

void util::add_unique(unsigned int dof, Real val,
                      std::vector<unsigned int> &list,
                      std::vector<Real> &list_val) {

  for (unsigned int i = 0; i < list.size(); i++) {

    if (list[i] == dof) {
      list_val[i] += val;

      return;
    }
  }

  list.push_back(dof);
  list_val.push_back(val);
}

std::string util::get_vec_major_axis(Point vec) {

  // get dot product and compare with the norm of vec
  auto norm = vec.norm();
  if (std::abs(norm - std::abs(vec(0))) < 1.0e-3 * norm)
    return "x";
  else if (std::abs(norm - std::abs(vec(1))) < 1.0e-3 * norm)
    return "y";
  else if (std::abs(norm - std::abs(vec(2))) < 1.0e-3 * norm)
    return "z";
  else
    return "";
}

double util::dist_between_points(const std::vector<double> &p1,
                                 const std::vector<double> &p2) {
  double r = 0;
  for (unsigned int i = 0; i < p1.size(); i++)
    r += (p1[i] - p2[i]) * (p1[i] - p2[i]);

  return std::sqrt(r);
}

void util::get_unit_vector(std::vector<double> &v) {
  double l = 0.;
  for (double i : v)
    l += i;

  l = std::sqrt(l);
  for (double &i : v)
    i = i / l;
}

unsigned int util::get_elem_id(const std::vector<double> &p,
                               const double &mesh_size,
                               const unsigned int &num_elems,
                               const unsigned int &dim) {

  if (dim == 2) {
    unsigned int i = p[0] / mesh_size;
    unsigned int j = p[1] / mesh_size;

    return j * num_elems + i;
  } else {
    unsigned int i = p[0] / mesh_size;
    unsigned int j = p[1] / mesh_size;
    unsigned int k = p[2] / mesh_size;

    return k * num_elems * num_elems + j * num_elems + i;
  }
}

unsigned int util::get_elem_id(const Point &p,
                               const double &mesh_size,
                               const unsigned int &num_elems,
                               const unsigned int &dim) {

  if (dim == 2) {
    unsigned int i = p(0) / mesh_size;
    unsigned int j = p(1) / mesh_size;

    return j * num_elems + i;
  } else {
    unsigned int i = p(0) / mesh_size;
    unsigned int j = p(1) / mesh_size;
    unsigned int k = p(2) / mesh_size;

    return k * num_elems * num_elems + j * num_elems + i;
  }
}

Point util::to_point(const std::vector<double> &p) {

  auto q = Point();
  for (unsigned int i = 0; i < p.size(); i++)
    q(i) = p[i];

  return q;
}

std::vector<double> util::determineRotator(const std::vector<double> &dir) {

  std::vector<double> rotator;

  if (std::abs(dir[0]) < 1.0e-11 && std::abs(dir[1]) < 1.0e-11) {

    rotator.push_back(1.0);

    rotator.push_back(0.0);

    rotator.push_back(0.0);

  } else {

    rotator.push_back(-dir[1]);

    rotator.push_back(dir[0]);

    rotator.push_back(0.0);
  }

  return rotator;
}

Point util::determineRotator(const Point &dir) {

  if (std::abs(dir(0)) < 1.0e-11 && std::abs(dir(1)) < 1.0e-11) {

    return Point(1., 0., 0.);
  } else {

    return Point(-dir(1), dir(0), 0.);
  }
}
