//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/Test_Mesh_Interface.cc
 * \author Ryan Wollaeger <wollaeger@lanl.gov>
 * \date   Thursday, Jun 07, 2018, 15:43 pm
 * \brief  Helper class for generating test meshes.
 * \note   Copyright (C) 2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_mesh_test_Test_Mesh_Interface_hh
#define rtt_mesh_test_Test_Mesh_Interface_hh

#include "mesh/Draco_Mesh.hh"
#include <algorithm>
#include <bitset>

namespace rtt_mesh_test {

//===========================================================================//
/*!
 * \class Test_Mesh_Interface
 *
 * \brief Class to help generate and contain serialized data needed to construct
 * meshes.
 *
 * This class is currently restricted to 2D meshes with quadrilateral cells.
 */
//===========================================================================//

class Test_Mesh_Interface {
public:
  typedef rtt_mesh::Draco_Mesh::Layout Layout;

  // >>> DATA

  const unsigned dim;
  const size_t num_cells;
  const size_t num_nodes;
  const size_t num_sides;
  const unsigned num_nodes_per_cell;
  std::vector<unsigned> cell_type;
  std::vector<unsigned> cell_to_node_linkage;
  std::vector<unsigned> side_set_flag;
  std::vector<unsigned> side_node_count;
  std::vector<unsigned> side_to_node_linkage;
  std::vector<double> coordinates;
  std::vector<unsigned> global_node_number;

public:
  //! Constructor.
  Test_Mesh_Interface(const size_t num_xdir_, const size_t num_ydir_,
                      const std::vector<unsigned> &global_node_number_ = {},
                      const double xdir_offset_ = 0.0,
                      const double ydir_offset_ = 0.0);

  // >>> SERVICES

  std::vector<unsigned> flatten_cn_linkage(const Layout &layout,
                                           const Layout &bd_layout,
                                           const Layout &go_layout) const;

  std::vector<unsigned> flatten_sn_linkage(const Layout &bd_layout) const;
};

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

Test_Mesh_Interface::Test_Mesh_Interface(
    const size_t num_xdir_, const size_t num_ydir_,
    const std::vector<unsigned> &global_node_number_, const double xdir_offset_,
    const double ydir_offset_)
    : dim(2), num_cells(num_xdir_ * num_ydir_),
      num_nodes((num_xdir_ + 1) * (num_ydir_ + 1)),
      num_sides(2 * (num_xdir_ + num_ydir_)), num_nodes_per_cell(4),
      global_node_number(global_node_number_) {

  // set the number of cells and nodes
  const size_t num_xdir = num_xdir_;
  const size_t num_ydir = num_ydir_;

  // set the number of sides and side flags
  const size_t poff = num_xdir + num_ydir; // parallel side offset

  // use two dimensions and Cartesian geometry
  cell_type.resize(num_cells, num_nodes_per_cell);

  // set the cell-to-node linkage counterclockwise about each cell
  cell_to_node_linkage.resize(num_cells * num_nodes_per_cell);
  for (size_t j = 0; j < num_ydir; ++j) {
    for (size_t i = 0; i < num_xdir; ++i) {

      // get a serial cell index
      size_t cell = i + num_xdir * j;

      // set each node entry per cell
      cell_to_node_linkage[4 * cell] = cell + j;
      cell_to_node_linkage[4 * cell + 1] = cell + j + 1;
      cell_to_node_linkage[4 * cell + 2] = cell + num_xdir + 1 + j + 1;
      cell_to_node_linkage[4 * cell + 3] = cell + num_xdir + 1 + j;
    }
  }

  // set two nodes per each side cell (only possibility in 2D)
  side_node_count.resize(num_sides, 2);

  // calculate arrays storing data about the sides of the mesh
  side_set_flag.resize(num_sides);
  side_to_node_linkage.resize(2 * num_sides);
  // ... over top and bottom faces
  for (size_t i = 0; i < num_xdir; ++i) {
    // bottom face
    side_set_flag[i] = 1;
    side_to_node_linkage[2 * i] = i;
    side_to_node_linkage[2 * i + 1] = i + 1;
    // top face
    side_set_flag[i + poff] = 3;
    side_to_node_linkage[2 * (i + poff)] = num_nodes - 1 - i;
    side_to_node_linkage[2 * (i + poff) + 1] = num_nodes - 1 - i - 1;
  }
  // ... over left and right faces
  for (size_t j = 0; j < num_ydir; ++j) {
    // right face
    side_set_flag[j + num_xdir] = 2;
    side_to_node_linkage[2 * (j + num_xdir)] = num_xdir * (j + 1) + j;
    side_to_node_linkage[2 * (j + num_xdir) + 1] = num_xdir * (j + 2) + j + 1;
    // left face
    side_set_flag[j + poff + num_xdir] = 4;
    side_to_node_linkage[2 * (j + poff + num_xdir)] =
        num_nodes - 1 - num_xdir * (j + 1) - j;
    side_to_node_linkage[2 * (j + poff + num_xdir) + 1] =
        num_nodes - 1 - num_xdir * (j + 2) - (j + 1);
  }

  // set some coordinates and global node indices
  coordinates.resize(dim * num_nodes);
  for (unsigned i = 0; i < num_nodes; ++i) {

    // TODO: generalize coordinate generation
    std::bitset<8> b2(i);
    coordinates[dim * i] = xdir_offset_ + static_cast<double>(b2[0]);
    coordinates[1 + dim * i] = ydir_offset_ + static_cast<double>(b2[1]);
  }

  if (global_node_number.size() == 0) {

    global_node_number.resize(num_nodes);

    // set global node index to local index
    for (unsigned i = 0; i < num_nodes; ++i)
      global_node_number[i] = i;
  }
}

//---------------------------------------------------------------------------//
// SERVICES
//---------------------------------------------------------------------------//

std::vector<unsigned>
Test_Mesh_Interface::flatten_cn_linkage(const Layout &layout,
                                        const Layout &bd_layout,
                                        const Layout &go_layout) const {
  // this function assumes a fixed number of nodes per cell
  std::vector<unsigned> cn_linkage;
  cn_linkage.reserve(num_cells * num_nodes_per_cell);
  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // node vector
    std::vector<unsigned> node_vec;
    // insert node indices from interior cell-to-cell linkage
    if (layout.count(cell) > 0) {
      for (auto lpair : layout.at(cell)) {
        node_vec.insert(node_vec.end(), lpair.second.begin(),
                        lpair.second.end());
      }
    }

    // insert node indices from boundary cell-to-side linkage
    if (bd_layout.count(cell) > 0) {
      for (auto lpair : bd_layout.at(cell)) {
        node_vec.insert(node_vec.end(), lpair.second.begin(),
                        lpair.second.end());
      }
    }

    // insert node indices from cell-to-ghost-cell linkage
    if (go_layout.count(cell) > 0) {
      for (auto lpair : go_layout.at(cell)) {
        node_vec.insert(node_vec.end(), lpair.second.begin(),
                        lpair.second.end());
      }
    }

    // sort and unique
    std::sort(node_vec.begin(), node_vec.end());
    auto last = std::unique(node_vec.begin(), node_vec.end());
    node_vec.erase(last, node_vec.end());

    // insert unique nodes into linkage array
    cn_linkage.insert(cn_linkage.end(), node_vec.begin(), node_vec.end());
  }

  // return the flattened cell-to-node linkage
  return cn_linkage;
}

//---------------------------------------------------------------------------//

std::vector<unsigned>
Test_Mesh_Interface::flatten_sn_linkage(const Layout &bd_layout) const {

  // this function assumes 2 nodes per side
  std::vector<unsigned> sn_linkage(num_sides * 2);

  for (unsigned cell = 0; cell < num_cells; ++cell) {

    // insert node indices from boundary cell-to-side linkage
    if (bd_layout.count(cell) > 0) {
      for (auto lpair : bd_layout.at(cell)) {

        // get the side
        unsigned side = lpair.first;

        // set the vertices
        sn_linkage[2 * side] = lpair.second[0];
        sn_linkage[2 * side + 1] = lpair.second[1];
      }
    }
  }

  // return the flattened side-to-node linage
  return sn_linkage;
}

} // end namespace rtt_mesh_test

#endif // rtt_mesh_Test_Mesh_Interface_hh
