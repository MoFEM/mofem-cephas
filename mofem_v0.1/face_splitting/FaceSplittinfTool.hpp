/* Copyright (C) 2014, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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


#ifndef __FACESPLITTINGTOOL_HPP__
#define __FACESPLITTINGTOOL_HPP__

#include "FieldInterface.hpp"
#include "CoreDataStructures.hpp"

#include <moab/AdaptiveKDTree.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;


namespace MoFEM {

struct FaceSplittingTools {

  FieldInterface& mField;
  AdaptiveKDTree kdTree;
  
  Interface& moab_distance_from_crack_surface;
  Core mb_instance_distance_from_crack_surface;
  AdaptiveKDTree kdTree_DistanceFromCrackSurface;

  FaceSplittingTools(FieldInterface& _mField): 
    mField(_mField),
    kdTree(&_mField.get_moab(),true),
    moab_distance_from_crack_surface(mb_instance_distance_from_crack_surface),
    kdTree_DistanceFromCrackSurface(&moab_distance_from_crack_surface,true) {


    kdTree_rootMeshset_DistanceFromCrackSurface = 0;
    opositeFrontEdges = 0;
    nodesOnCrackSurface = 0;
    crackSurfaceCrossingEdges = 0;
    crackFrontTests = 0;
    chopTetsFaces = 0;

  }

  ~FaceSplittingTools() {}
 
  PetscErrorCode cleanMeshsets() {
    PetscFunctionBegin;

    rval = mField.get_moab().delete_entities(&kdTree_rootMeshset_DistanceFromCrackSurface,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&opositeFrontEdges,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&nodesOnCrackSurface,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&crackSurfaceCrossingEdges,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&crackFrontTests,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&chopTetsFaces,1); CHKERR_PETSC(rval);

    PetscFunctionReturn(0);
  }

  //Cloasest face to crack surface, crack surface perspective 

  EntityHandle kdTree_rootMeshset_DistanceFromCrackSurface;
  map<EntityHandle,EntityHandle> map_nodes;
  PetscErrorCode buildKDTreeForCrackSurface(
    Range &entities,const BitRefLevel bit_mesh);
  PetscErrorCode buildKDTreeForCrackSurface(const BitRefLevel bit_mesh);

  //Init bit level data
  Range mesh_level_nodes;
  Range mesh_level_edges;
  Range mesh_level_tris;
  Range mesh_level_tets;
  PetscErrorCode initBitLevelData(const BitRefLevel bit_mesh);

  //Calulte distance on mesh
  PetscErrorCode calculateDistanceFromCrackSurface(Range &nodes);
  PetscErrorCode calculateDistanceFromCrackSurface();

  //Front edges

  EntityHandle opositeFrontEdges;
  EntityHandle nodesOnCrackSurface;
  EntityHandle crackSurfaceCrossingEdges;

  PetscErrorCode getOpositeForntEdges(bool createMeshset);
  PetscErrorCode getCrackSurfaceCorssingEdges(bool createMeshset);

  //Front tets

  EntityHandle crackFrontTests;
  EntityHandle chopTetsFaces;

  PetscErrorCode getCrackFrontTets(bool createMeshset);
  PetscErrorCode chopTetsUntilNonOneLeftOnlyCrackSurfaceFaces(bool createMeshset);

  //Split new crack front faces

  PetscErrorCode addNewSurfaceFaces_to_Cubit_msId200();
  PetscErrorCode addcrackFront_to_Cubit201();
  PetscErrorCode meshRefine(const BitRefLevel bit_mesh,const BitRefLevel new_bit_mesh);
  PetscErrorCode splitFaces(const BitRefLevel bit_mesh,const BitRefLevel new_bit_mesh);

  //Cat mesh
  
  PetscErrorCode catMesh(const BitRefLevel bit_mesh,const BitRefLevel new_bit_mesh);

  private:
  ErrorCode rval;
  PetscErrorCode ierr;

  double diffNTET[4*3];
  Tag th_b;
  Tag th_distance;

  PetscErrorCode calculate_qualityAfterProjectingNodes(EntityHandle meshset);
  PetscErrorCode calculate_qualityAfterProjectingNodes(Range &option_nodes,double &current_b);


};

}

#endif // __FACESPLITTINGTOOL_HPP__

