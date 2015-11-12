/** \file GetAdjacancies.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

// const static int debug = 1;

PetscErrorCode Core::get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) {
  PetscFunctionBegin;
  RefMoFEMEntity from_ref_entiti(moab,from_entiti);
  //cerr << "from:\n";
  //cerr << from_ref_entiti << endl;
  rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERR_PETSC(rval);
  Range::iterator eit = adj_entities.begin();
  //cerr << "to:\n";
  for(;eit!=adj_entities.end();) {
    RefMoFEMEntity adj_entiti(moab,*eit);
    //cerr << "\t" << adj_entiti << endl;
    if(from_ref_entiti.get_BitRefLevel() != adj_entiti.get_BitRefLevel()) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) {
  PetscFunctionBegin;
  RefMoFEMEntity from_ref_entiti(moab,from_entiti);
  //cerr << "from:\n";
  //cerr << from_ref_entiti << endl;
  rval = moab.get_adjacencies(&from_entiti,1,to_dimension,false,adj_entities); CHKERR_PETSC(rval);
  Range::iterator eit = adj_entities.begin();
  //cerr << "to:\n";
  for(;eit!=adj_entities.end();) {
    RefMoFEMEntity adj_entiti(moab,*eit);
    //cerr << "\t" << adj_entiti << endl;
    if(!(from_ref_entiti.get_BitRefLevel()&adj_entiti.get_BitRefLevel()).any()) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_adjacencies(
    const MoFEMProblem *problem_ptr,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type,
    const int verb) {
  PetscFunctionBegin;
  BitRefLevel bit = problem_ptr->get_BitRefLevel();
  ierr = get_adjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type,const int verb) {
  PetscFunctionBegin;
  if(verb>0) {
    ostringstream ss;
    ss << "from: " << bit << endl << "to: " << endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  rval = moab.get_adjacencies(from_entities,num_netities,to_dimension,false,adj_entities,operation_type); CHKERR_PETSC(rval);
  Range::iterator eit = adj_entities.begin();
  //cerr << "to:\n";
  for(;eit!=adj_entities.end();) {
    RefMoFEMEntity adj_entiti(moab,*eit);
    if(verb>0) {
      ostringstream ss;
      ss << "\t" << adj_entiti << endl;
      PetscPrintf(comm,ss.str().c_str());
    }
    if(!(adj_entiti.get_BitRefLevel()&bit).any() ) {
      eit = adj_entities.erase(eit);
    } else {
      eit++;
    }
  }
  PetscFunctionReturn(0);
}

}
