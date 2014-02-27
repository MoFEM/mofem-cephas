/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#include<FaceSplittinfTool.hpp>

namespace MoFEM {

PetscErroCode FaceSplittingTools::loadMesh(const BitRefLevel &bit_mesh,EntityHandle meshset) {
  PetscFunctionBegin;
  Range tets;
  rval = mField.get_moab().get_entities_by_type(meshset,MBTET,tets,true); CHKERR_PETSC(rval);
  Range nodes;
  ierr = mField.get_moab().get_connectivity(tets,nodes,true); CHKERRQ(ierr);
  for(Range::iterate nit = nodes.begin();nit!=nodes.end();nit++) {
    double coords[3];
    rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
    ierr = moab_work.create_vertex(coords,map_nodes[*nit]);  CHKERRQ(ierr);
  }
  for(Range::iterator tit = tets.begin();tit!=tets.end();tit++) {
    const EntityHandle* conn;
    int num_nodes; 
    rval = mField.get_moab().get_connectivity(&*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle work_conn[num_nodes];
    for(int nn = 0;nn<num_nodes;nn++) {
      work_conn[nn] = map_nodes[conn[nn]];
    }
    rval = moab_work.create_element(MBTET,work_conn,4,map_tets[*tit]); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}

PetscErroCode FaceSplittingTools::findFaces(
  ublas::vector<double,ublas::bounded_array<double,3> > &face_normal) {
  PetscFunctionBegin;


  PetscFunctionReturn(0);
}




}

