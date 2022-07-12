/** \file Projection10NodeCoordsOnField.hpp

FIXME: Move code to cpp file.

Project displacements/coordinates from 10 node tetrahedra on hierarchical
approximation base.

This is example how to use MoFEM::DofMethod when some operator for each node
need to be applied.

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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
