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

/* This file is calculated of Plane wave impinged on a sphere as Dirichlet BC */

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <moab/ParallelComm.hpp>

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>
#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#include <DirichletBC.hpp>

struct PlaneWaveProblem: public DisplacementBCFEMethodPreAndPostProc {

	PlaneWaveProblem(
		FieldInterface& _mField,const string &_field_name,Mat _Aij,Vec _X,Vec _F): 
		DisplacementBCFEMethodPreAndPostProc(_mField,_field_name,_Aij,_X,_F) {}


	double reVal;
	PetscErrorCode analyticalFunction() {
		PetscFunctionBegin;
		PetscFunctionReturn(0);
	}

	double R,Phi;
	PetscErrorCode sphericalCoordinates(double x,double y,double z) {
		PetscFunctionBegin;
		R = sqrt(pow(x,2.0)+pow(u,2.0));
		PetscFunctionReturn(0);
	}

	PetscErrorCode iNitalize() {
		PetscFunctionBegin;

		if(map_zero_rows.empty()) {
			ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
			for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
				if(it->get_Cubit_name().compare(0,13,"ANALYTICAL_SOL") == 0) {
					DisplacementCubitBcData mydata;
					ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
					for(int dim = 0;dim<3;dim++) {
						Range ents;
						ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
						if(dim>1) {
							Range _edges;
							ierr = mField.get_moab().get_adjacencies(ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
							ents.insert(_edges.begin(),_edges.end());
						}
						if(dim>0) {
							Range _nodes;
							rval = mField.get_moab().get_connectivity(ents,_nodes,true); CHKERR_PETSC(rval);
							ents.insert(_nodes.begin(),_nodes.end());
						}
						for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
							for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof)) {
								if(dof->get_ent_type() == MBVERTEX) {
									EntityHandle ent = dof->get_ent();
									double coords[3];
									rval = mField.get_moab().get_coords(&ent,1,coords); CHKERR_PETSC(rval);
									ierr = sphericalCoordinates(coords[0],coords[1],coords[2]); CHKERRQ(ierr); //get spherical coordinates from Cartesian Coordinates
									ierr = analyticalFunction(); CHKERRQ(ierr); //calculate the funictional scalar value for acoustic potential
									map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = reVal;
								} else {
									map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;  //the scalar value defined on vertex nodes of surface
								}
							}
						}
					}
				}
			}
			dofsIndices.resize(map_zero_rows.size());
			dofsValues.resize(map_zero_rows.size());
			int ii = 0;
			map<DofIdx,FieldData>::iterator mit = map_zero_rows.begin();
			for(;mit!=map_zero_rows.end();mit++,ii++) { 
				dofsIndices[ii] = mit->first;  //get the key of the interiator pointer for map container  
				dofsValues[ii] = mit->second;
			}
		}


		PetscFunctionReturn(0);
	}

   //preprocess() fill in the solution vector with scalar Dirichlet values and zero otherwise.
   //postprocess() fill in 1 on maindiagonal and zero on row and column in stiffness matrix.
   //postprocess() fill zeros to right hand side vector F with Dirichlet position.

}

