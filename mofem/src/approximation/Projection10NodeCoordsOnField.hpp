/** \file Projection10NodeCoordsOnField.hpp

FIXME: Move code to cpp file.

Project displacements/coordinates from 10 node tetrahedra on hierarchical
approximation base.

This is example how to use MoFEM::DofMethod when some operator for each node
need to be applied.

*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __PROJECTION10NODECOORDSONFIELD_HPP__
#define __PROJECTION10NODECOORDSONFIELD_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief Projection of edge entities with one mid-node on hierarchical basis
 */
struct Projection10NodeCoordsOnField : public DofMethod {

  Projection10NodeCoordsOnField(Interface &m_field, std::string field_name,
                                int verb = 0);

  MoFEMErrorCode preProcess();

  MoFEMErrorCode operator()();

  MoFEMErrorCode postProcess();

protected:

  Interface &mField;
  std::string fieldName;
  int vErbose;

  VectorDouble coords;
  VectorDouble3 aveMidCoord;
  VectorDouble3 midNodeCoord;
  VectorDouble3 diffNodeCoord;
  VectorDouble3 dOf;

};

struct ProjectionFieldOn10NodeTet : public Projection10NodeCoordsOnField {


  ProjectionFieldOn10NodeTet(Interface &m_field, std::string _fieldName,
                             bool set_nodes, bool on_coords,
                             std::string on_tag = "NoNE");

  MoFEMErrorCode preProcess();

  MoFEMErrorCode operator()(); 


  bool setNodes;
  bool onCoords;
  std::string onTag;

  const int maxApproximationOrder;

  Tag th;

protected:

  Field_multiIndex::index<FieldName_mi_tag>::type::iterator field_it;
  VectorDouble L;
  VectorDouble K;

};

} // namespace MoFEM

#endif // __PROJECTION10NODECOORDSONFIELD_HPP__
