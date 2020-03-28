////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2019 Prashant K. Jha
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include "tests.hpp"

// local anonymous namespace
namespace {

struct Tip {
  unsigned int d_n;
  unsigned int d_id;
  unsigned int d_pid;
  std::vector<unsigned int> d_nodes;

  Tip() : d_n(0), d_id(0), d_pid(0), d_nodes(std::vector<unsigned int>()){};

  Tip(const unsigned int &n, const unsigned int &id, const unsigned int &pid)
      : d_n(n), d_id(id), d_pid(pid), d_nodes(std::vector<unsigned int>(n)){};

  void print(const int &info_level) {
    out << "Tip id = " << d_id << ", node = " << d_n << ", pid = " << d_pid
        << ", old nodes = ";
    for (const auto &n : d_nodes)
      out << n << "; ";
    out << std::endl;
  }
};

bool findTip(unsigned int l, unsigned int &loc, const std::vector<Tip> &list) {
  for (size_t i = 0; i < list.size(); i++)
    if (list[i].d_id == l) {
      loc = i;
      return true;
    }

  return false;
}

unsigned int
get_unique_id(const std::vector<Tip> &list,
              const std::pair<unsigned int, unsigned int> &id_range) {

  unsigned int dummy = 0;
  for (size_t i = id_range.first; i < id_range.second; i++) {

    if (findTip(i, dummy, list))
      continue;

    return i;
  }

  //  libmesh_error();
}

unsigned int getPid(unsigned int l, const unsigned int &num_local_nodes) {

  return l / num_local_nodes;
}

std::string get_print_string(const std::vector<std::string> &message) {

  std::ostringstream oss;

  unsigned int i = 0;
  for (const auto &s : message) {

    oss << "Processor = " << i << std::endl;
    oss << "  " << s << std::endl;

    i++;
  }

  return oss.str();
}

std::string get_string(std::ostringstream &oss) {
  auto s = oss.str();

  oss.str("");
  oss.clear();

  return s;
}

} // namespace

void test::comm::run(int argc, char **argv, Parallel::Communicator *comm) {

  if (comm->rank() == 0) {
    printf("********** TestComm **************\n");
  }

  // dummy mesh
  ReplicatedMesh mesh(*comm);
  MeshTools::Generation::build_square(mesh, 10, 10, 0., 1., 0., 1., QUAD4);

  // for debug output
  std::ostringstream oss;

  // assign each processor a range from which they can give global ids to
  // their local tips
  unsigned int num_local_nodes = 1000;
  std::pair<unsigned int, unsigned int> id_range = {
      comm->rank() * num_local_nodes, (comm->rank() + 1) * num_local_nodes};

  // generate dummy local nodes
  std::vector<unsigned int> nodes;
  for (size_t i = 0; i < num_local_nodes; i++)
    nodes.emplace_back(comm->rank() * num_local_nodes + i);

  // create some tips
  std::vector<Tip> local_tips;
  for (size_t i = 0; i < 2; i++)
    local_tips.emplace_back(
        Tip(nodes[i], get_unique_id(local_tips, id_range), comm->rank()));

  oss << "  Local tips = ";
  for (const auto &t : local_tips)
    oss << "(" << t.d_id << ", " << t.d_pid << ", " << t.d_n << "); ";

  //
  // Test 1: Pass the string which has information of local tips it owns to
  // all other processor and print it on processor 0
  //
  //
  // To make below work, note the size of global_message
  // Size is 1 and not comm->size()
  //
  auto global_message = std::vector<std::string>(1, "");
  global_message[0] = get_string(oss);
  comm->allgather(global_message);

  // out << "----- received strings -------\n";
  // for (size_t i=0; i<global_message.size(); i++)
  //   out << "string id = " << i << ", string = " << global_message[i] << "\n";

  // print message
  {
    std::ostringstream oss1;
    oss1 << "Allgather message strings test\n";
    oss1 << "In processor " << comm->rank() << ", received messages are :\n";
    oss1 << get_print_string(global_message) << std::endl;

    // print message in 0 processor
    out << oss1.str();

    // also print message in 1 processor (use printf as std::cout does not work)
    if (comm->rank() == 1) {
      printf("***\nDebug processor 1 received data \n");
      printf("%s", oss1.str().c_str());
      printf("\n ***\n");
    }
  }

  // communicate with all the processor about new nodes
  std::vector<Tip> global_tips;
  std::vector<unsigned int> global_tip_ids;
  std::vector<unsigned int> global_tip_nodes;
  for (const auto &l : local_tips) {
    global_tip_ids.push_back(l.d_id);
    global_tip_nodes.push_back(l.d_n);
  }
  comm->allgather(global_tip_ids);
  comm->allgather(global_tip_nodes);

  // create global tip from above received information
  for (size_t i = 0; i < global_tip_ids.size(); i++)
    global_tips.emplace_back(global_tip_nodes[i], global_tip_ids[i],
                             getPid(global_tip_ids[i], num_local_nodes));

  // output to screen
  oss << "  Global tips = ";
  for (const auto &t : global_tips)
    oss << "(" << t.d_id << ", " << t.d_pid << ", " << t.d_n << "); ";
  oss << std::endl;
  global_message = {get_string(oss)};
  comm->allgather(global_message);

  // print message
  {
    std::ostringstream oss1;
    oss1 << "Allgather global tips test\n";
    oss1 << "In processor " << comm->rank()
         << ", updated tips in all processors are :\n";
    oss1 << get_print_string(global_message) << std::endl;

    // print message in 0 processor
    out << oss1.str();

    // also print message in 1 processor (use printf as std::cout does not work)
    if (comm->rank() == 1) {
      printf("***\nDebug processor 1 received data \n");
      printf("%s", oss1.str().c_str());
      printf("\n ***\n");
    }
  }

  // end
}
