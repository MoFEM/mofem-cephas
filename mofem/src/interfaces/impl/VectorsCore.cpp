/** \file Vectors.cpp
 * \brief Managing Vec, IS and Scatter
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <ISManager.hpp>
#include <VecManager.hpp>

namespace MoFEM {

// const static int debug = 1;

PetscErrorCode Core::VecCreateSeq(const std::string &name,RowColData rc,Vec *V) const {
  return VecManager(*this).vecCreateSeq(name,rc,V);
}
PetscErrorCode Core::VecCreateGhost(const std::string &name,RowColData rc,Vec *V) const {
  return VecManager(*this).vecCreateGhost(name,rc,V);
}
PetscErrorCode Core::ISCreateProblemOrder(
  const std::string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb
) const {
  return ISManager(*this).isCreateProblemOrder(problem,rc,min_order,max_order,is);
}
PetscErrorCode Core::ISCreateProblemFieldAndRank(
  const std::string &problem,
  RowColData rc,
  const std::string &field,
  int min_coeff_idx,
  int max_coeff_idx,
  IS *is,
  int verb
) const {
  return ISManager(*this).isCreateProblemFieldAndRank(problem,rc,field,min_coeff_idx,max_coeff_idx,is);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::ISCreateFromProblemFieldToOtherProblemField(
  const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
  const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
  std::vector<int> &idx,std::vector<int> &idy,int verb
) const {
  return ISManager(*this).isCreateFromProblemFieldToOtherProblemField(
    x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,idx,idy
  );
}

PetscErrorCode Core::ISCreateFromProblemFieldToOtherProblemField(
  const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
  const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
  IS *ix,IS *iy,int verb
) const {
  return ISManager(*this).isCreateFromProblemFieldToOtherProblemField(
    x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,ix,iy
  );
}
PetscErrorCode Core::VecScatterCreate(
  Vec xin,const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
  Vec yin,const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
  VecScatter *newctx,int verb
) const {
  return VecManager(*this).vecScatterCreate(
    xin,x_problem,x_field_name,x_rc,
    yin,y_problem,y_field_name,y_rc,
    newctx
  );
}
PetscErrorCode Core::ISCreateFromProblemToOtherProblem(
  const std::string &x_problem,RowColData x_rc,
  const std::string &y_problem,RowColData y_rc,
  std::vector<int> &idx,std::vector<int> &idy,
  int verb
) const {
  return ISManager(*this).isCreateFromProblemToOtherProblem(
    x_problem,x_rc,y_problem,y_rc,idx,idy
  );
}
PetscErrorCode Core::ISCreateFromProblemToOtherProblem(
  const std::string &x_problem,RowColData x_rc,
  const std::string &y_problem,RowColData y_rc,
  IS *ix,IS *iy,
  int verb
) const {
  return ISManager(*this).isCreateFromProblemToOtherProblem(
    x_problem,x_rc,y_problem,y_rc,ix,iy
  );
}
PetscErrorCode Core::VecScatterCreate(
  Vec xin,
  const std::string &x_problem,
  RowColData x_rc,
  Vec yin,
  const std::string &y_problem,
  RowColData y_rc,
  VecScatter *newctx,
  int verb
) const {
  return VecManager(*this).vecScatterCreate(
    xin,x_problem,x_rc,yin,y_problem,y_rc,newctx
  );
}
PetscErrorCode Core::set_local_ghost_vector(
  const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  return VecManager(*this).setLocalGhostVector(problem_ptr,rc,V,mode,scatter_mode);
}
PetscErrorCode Core::set_local_ghost_vector(
  const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  return VecManager(*this).setLocalGhostVector(name,rc,V,mode,scatter_mode);
}
PetscErrorCode Core::set_global_ghost_vector(
  const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  return VecManager(*this).setGlobalGhostVector(problem_ptr,rc,V,mode,scatter_mode);
}
PetscErrorCode Core::set_global_ghost_vector(
  const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
) const {
  return VecManager(*this).setGlobalGhostVector(name,rc,V,mode,scatter_mode);
}
PetscErrorCode Core::set_other_local_ghost_vector(
  const Problem *problem_ptr,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,
  Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
) {
  return VecManager(*this).setOtherLocalGhostVector(
    problem_ptr,field_name,cpy_field_name,rc,V,mode,scatter_mode
  );
}
PetscErrorCode Core::set_other_local_ghost_vector(
  const std::string &name,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,
  Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
) {
  return VecManager(*this).setOtherLocalGhostVector(
    name,field_name,cpy_field_name,rc,V,mode,scatter_mode
  );
}
PetscErrorCode Core::set_other_global_ghost_vector(
  const Problem *problem_ptr,
  const std::string& field_name,
  const std::string& cpy_field_name,
  RowColData rc,
  Vec V,
  InsertMode mode,
  ScatterMode scatter_mode,
  int verb
) {
  return VecManager(*this).setOtherGlobalGhostVector(
    problem_ptr,field_name,cpy_field_name,rc,V,mode,scatter_mode
  );
}
PetscErrorCode Core::set_other_global_ghost_vector(
  const std::string &name,
  const std::string& field_name,
  const std::string& cpy_field_name,
  RowColData rc,
  Vec V,
  InsertMode mode,
  ScatterMode scatter_mode,
  int verb
) {
  return VecManager(*this).setOtherGlobalGhostVector(
    name,field_name,cpy_field_name,rc,V,mode,scatter_mode
  );
}


}
