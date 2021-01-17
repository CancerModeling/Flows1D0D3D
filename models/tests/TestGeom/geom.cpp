////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "libmesh/point_locator_tree.h"
#include "tests.hpp"
#include <utils.hpp>

namespace {
void random_init() { srand(time(nullptr)); }

// random number between 0 and 1
double get_random_number() { return double(std::rand()) / double(RAND_MAX); }

std::pair<Point, Point> bounding_box(const std::vector<Point> &nodes) {

  auto box = std::pair<Point, Point>(Point(), Point());

  for (const auto &x : nodes) {

    for (unsigned int i = 0; i < 3; i++) {
      if (x(i) < box.first(i))
        box.first(i) = x(i);

      if (x(i) > box.second(i))
        box.second(i) = x(i);
    }
  }

  return box;
}
} // namespace

void test::geom::rotation(int argc, char **argv, Parallel::Communicator *comm) {

  if (comm->rank() == 0) {
    printf("********** TestGeom: Rotation **************\n");
  }

  //
  // Test description:
  // We consider a arbitrary vector in 3d, and set this as axis of rotation.
  // We then apply the rotation of vector (1,0,0) about this axis and write
  // the points to file
  //
  std::ofstream file;
  file.open("rotate.csv");
  // header
  //  file << "x, y, z, r, theta\n";

  {
    Point ax = Point(1., 1., 1.);
    ax = ax / ax.norm();
    Point x0 = Point(3., 3., 0.);

    file << 0. << ", " << 0. << ", " << 0. << ", " << 10. << ", " << -1.
         << std::endl;
    file << ax(0) << ", " << ax(1) << ", " << ax(2) << ", " << 5. << ", " << -1.
         << std::endl;

    unsigned int N = 100;
    for (unsigned int i = 0; i < N; i++) {
      double theta = i * 2. * M_PI / N;

      auto x = util::rotate(x0, theta, ax);
      file << x(0) << ", " << x(1) << ", " << x(2) << ", " << 1. << ", "
           << theta << std::endl;
    }
  }

  {
    Point ax = Point(1., 1., 0.);
    ax = ax / ax.norm();
    Point x0 = Point(6., 0., 0.);

    file << ax(0) << ", " << ax(1) << ", " << ax(2) << ", " << 5. << ", "
         << -10. << std::endl;
    // file << x0(0) << ", " << x0(1) << ", " << x0(2) << ", " << -1. <<
    // std::endl;

    unsigned int N = 100;
    for (unsigned int i = 0; i < N; i++) {
      double theta = i * 2. * M_PI / N;

      auto x = util::rotate(x0, -theta, ax);
      file << x(0) << ", " << x(1) << ", " << x(2) << ", " << 1. << ", "
           << 10. * theta << std::endl;
    }
  }

  file.close();
}

void test::geom::cylinder(int argc, char **argv, Parallel::Communicator *comm, Point x0,
                          Point a, double R, double L) {

  if (comm->rank() == 0) {
    printf("********** TestGeom: Cylinder **************\n");
  }

  //
  // Test description:
  // Render cross-section of cylinder with given axis and point x0
  //
  std::ofstream file;
  file.open("cylinder.csv");
  // header
  //  file << "x, y, z, r, theta\n";

  // must normalize axis of cylinder
  a = a / a.norm();
  {
    file << 0. << ", " << 0. << ", " << 0. << ", " << 10. << ", " << -1.
         << std::endl;
    file << x0(0) << ", " << x0(1) << ", " << x0(2) << ", " << 10. << ", "
         << -1. << std::endl;
    file << a(0) << ", " << a(1) << ", " << a(2) << ", " << 10. << ", "
         << -1. << std::endl;
  }

  // render cross-section at fixed intervals
  unsigned int N = 10;
  for (unsigned int j = 0; j <= N; j++) {

    // parametric coordinate
    double s = double(j) * L / (double(N));

    // coordinate
    Point xs = x0 + s * a;
    file << xs(0) << ", " << xs(1) << ", " << xs(2) << ", " << 10. << ", "
         << -1. << std::endl;

    // get basis vector for plane parallel to cross-section and passing
    // through origin
    Point e1 = xs - (xs * a) * a;
    e1 = e1 / e1.norm();

    // rotate R*e1 from 0 to 2 pi to generate disk in plane passing through
    // origin
    // then translate it to xs
    unsigned int Nr = 100;
    for (unsigned int i = 0; i < Nr; i++) {
      double theta = i * 2. * M_PI / Nr;

      auto x = xs + util::rotate(R * e1, theta, a);
      file << x(0) << ", " << x(1) << ", " << x(2) << ", "
           << 10. * double(j + 1) << ", " << theta << std::endl;
    }
  }

  file.close();
}

void test::geom::point_in_cylinder(int argc, char **argv, Parallel::Communicator *comm, Point x0,
                                   Point a, double R, double L) {

  if (comm->rank() == 0) {
    printf("********** TestGeom: Point In Cylinder **************\n");
  }

  //
  // Test description:
  // Render crosso-section of cylinder with given axis and point x0
  //
  std::ofstream file;
  file.open("point_in_cylinder.csv");
  // header
  //  file << "x, y, z, r, theta\n";

  a = a / a.norm();
  // get basis vector for plane parallel to cross-section and passing
  // through origin
  Point e1 = x0 - (x0 * a) * a;
  e1 = e1 / e1.norm();

  // get random point between (0,L) and (0,R) and check if random point lies
  // in cylinder
  unsigned int N = 100;
  for (unsigned int j = 0; j <= N; j++) {

    // get random r
    double r = get_random_number() * 2. * R;
    double l = get_random_number() * 2. * L;

    Point x = x0 + l * a + r * e1;

    int in_cylinder = -1.;

    if (util::is_inside_cylinder(x - x0, L, R, a))
      in_cylinder = 1.;

    // write
    file << x(0) << ", " << x(1) << ", " << x(2) << ", " << in_cylinder << ", "
         << 0. << std::endl;
  }

  file.close();
}

void test::geom::angle_test(int argc, char **argv, Parallel::Communicator *comm) {

  if (comm->rank() == 0) {
    printf("********** TestGeom: Angle Test **************\n");
  }

  //
  // Test description:
  // We consider a arbitrary vector in 3d, and set this as axis of rotation.
  // We then apply the rotation of vector x0 about this axis, and use the
  // resulting point to determine the angle using angle() function util
  //
  std::ofstream file;
  file.open("angle_test.csv");

  {
    std::pair<Point, Point> test_setup_1;
    std::vector<std::pair<Point, double>> test_points_1;

    Point ax = Point(0., 0., 1.);
    ax = ax / ax.norm();
    Point x0 = Point(1., 0., 0.);
    test_setup_1 = std::make_pair(ax, x0);

    //    file << 0. << ", " << 0. << ", " << 0. << ", " << 10. << ", " << -1.
    //         << std::endl;
    //    file << ax(0) << ", " << ax(1) << ", " << ax(2) << ", " << 5. << ", " << -1.
    //         << std::endl;

    unsigned int N = 20;
    for (unsigned int i = 0; i < N; i++) {
      double theta = i * 2. * M_PI / N;

      auto x = util::rotate(x0, theta, ax);
      test_points_1.emplace_back(x, theta);

      // normal
      Point normal_1 = ax;
      Point normal_2 = normal_1.cross(x0);
      double p_upper = x * normal_2;
      if (p_upper < 0.)
        p_upper = -1.;
      else
        p_upper = 1.;

      file << x(0) << ", " << x(1) << ", " << x(2) << ", " << p_upper << ", "
           << theta << std::endl;
    }

    // find angle without orientation
    for (auto p : test_points_1) {

      Point x = p.first;
      double theta_find = util::angle(x0, x);
      x += Point(2., 2., 2.);

      file << x(0) << ", " << x(1) << ", " << x(2) << ", " << 1. << ", "
           << theta_find << std::endl;
    }

    // find angle with axis of rotation
    for (auto p : test_points_1) {

      Point x = p.first;
      double theta_find = util::angle(x0, x, ax);
      x += Point(4., 4., 4.);

      file << x(0) << ", " << x(1) << ", " << x(2) << ", " << 1. << ", "
           << theta_find << std::endl;
    }

    // find angle with orientation
    for (auto p : test_points_1) {

      Point x = p.first;
      double theta_find = util::angle(x0, x, Point(0., 0., 1.), false);
      x += Point(6., 6., 6.);

      file << x(0) << ", " << x(1) << ", " << x(2) << ", " << 1. << ", "
           << theta_find << std::endl;
    }
  }

  file.close();
}

void test::geom::cube_to_ball(int argc, char **argv, unsigned int dim, double a,
                              Parallel::Communicator *comm) {

  if (comm->rank() == 0) {
    printf("********** TestGeom: Cube to Ball **************\n");
  }

  //
  // Test description:
  // We consider a uniform discretization of the cube of given size and map
  // all the points to the unit ball
  //
  std::ofstream file;
  file.open("cube_to_ball_" + std::to_string(dim) + ".csv");

  {
    // discretization along each axis
    int N = 2;
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

          Point x = Point(double(i), double(j), double(k)) * a / double(N);

          {
            auto y = x;
            // y += Point(2. * a, 2. * a, 2. * a);
            // file << y(0) << ", " << y(1) << ", " << y(2) << ", " << 0 << ", "
            //      << counter << std::endl;
          }

          // transform x to point in unit ball
          Point y = 1. * x / x.norm();

          // translate for visualization
          // y += Point(2. * a, 2. * a, 2. * a);

          file << y(0) << ", " << y(1) << ", " << y(2) << ", " << 1 << ", "
               << counter << std::endl;

          counter++;
        }
  }

  file.close();
}

void test::geom::cylinder_intersection_mesh(int argc, char **argv, Parallel::Communicator *comm, Point x0,
                                            Point a, double R, double L) {

  if (comm->rank() == 0) {
    printf("********** TestGeom: Cylinder intersecting Mesh **************\n");
  }

  //
  // Test description:
  // Discretize cylinder's surface and find the elements in mesh containing
  // discrete points on cylinder
  //
  std::ofstream file;
  file.open("cylinder_intersecting_mesh.csv");

  // mesh
  ReplicatedMesh mesh(*comm);

  // must normalize axis of cylinder
  a = a / a.norm();

  // render cross-section at fixed intervals
  std::vector<Point> cyl_points;
  std::vector<double> cyl_weights;
  std::vector<unsigned int> cyl_elems;

  unsigned int N = 10; // 10 elements along length
  for (unsigned int j = 0; j < N; j++) {

    // parametric coordinate
    double s = (double(j) + 0.5) * L / (double(N));

    // coordinate at center of cross-section
    Point xs = x0 + s * a;

    // get basis vector for plane parallel to cross-section and passing
    // through origin
    Point e1 = xs - (xs * a) * a;
    e1 = e1 / e1.norm();

    // now discretize the theta direction
    // rotate R*e1 from 0 to 2 pi to generate disk in plane passing through
    // origin then translate it to xs
    unsigned int Nr = 10; // 10 intervals in [0, 2pi]
    for (unsigned int i = 0; i < Nr; i++) {

      double theta = (double(i) + 0.5) * 2. * M_PI / Nr;

      auto x = xs + util::rotate(R * e1, theta, a);

      // compute weight of the point x
      double w = (L / double(N)) * 2. * M_PI / Nr * R;

      // save
      cyl_points.push_back(x);
      cyl_weights.push_back(w);
    }
  }

  // find the bounding box
  auto box = bounding_box(cyl_points);
  box.first -= Point(L / 5, L / 5, L / 5);
  box.second += Point(L / 5, L / 5, L / 5);

  // create mesh
  MeshTools::Generation::build_cube(mesh, 10, 10, 10,
                                    box.first(0), box.second(0),
                                    box.first(1), box.second(1),
                                    box.first(2), box.second(2),
                                    HEX8);

  // get point locator
  const auto &mesh_locator = mesh.point_locator();

  // loop over points on cylinder and find element in the mesh

  for (const auto &x : cyl_points) {

    const auto *e = mesh_locator(x);
    cyl_elems.push_back(e->id());
  }

  // write the data to file
  for (unsigned int i = 0; i < cyl_points.size(); i++)
    file << cyl_points[i](0) << ", " << cyl_points[i](1)
         << ", " << cyl_points[i](2) << ", " << cyl_weights[i]
         << ", " << cyl_elems[i] << "\n";

  file.close();

  // compute sum of weights and see if it is equal to area of cylinder
  double sum = 0.;
  for (auto w : cyl_weights)
    sum += w;
  out << "Sum of weights = " << sum
      << ", surface area of cylinder = " << 2. * M_PI * R * L << "\n";

  mesh.write("mesh.e");
}