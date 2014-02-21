/* Copyright (C) 2013, Michael Cortis (mikecortis AT gmail.com)
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


#ifndef __FUNTIONSFORFIELDDATA_HPP__
#define __FUNTIONSFORFIELDDATA_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "FEMethod_UpLevelStudent.hpp"

namespace MoFEM {

struct ComputeFieldGradients: public FEMethod_UpLevelStudent {
		
		Tag th_field_grad;
		
    FieldInterface& mField;
		string field_name;
		bool normalize;
    ComputeFieldGradients(FieldInterface& _mField,string _field_name = "POTENTIAL_FIELD",bool _normalize = false):
		FEMethod_UpLevelStudent(_mField.get_moab()),mField(_mField), field_name(_field_name), normalize(_normalize) {
			const int def_val_len = 0;
			rval = moab.tag_get_handle("FIELD_GRAD_ON_GAUSS_POINTS",def_val_len,MB_TYPE_DOUBLE,th_field_grad,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);

    };
		
		PetscErrorCode ComputeFibreDirection(vector<ublas::matrix<double> > &normalized_phi) {
			PetscFunctionBegin;
			
			MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator field_it;
			field_it = moabfields->get<FieldName_mi_tag>().find(field_name);
			const int field_rank = field_it->get_max_rank();
			
			vector< ublas::matrix< FieldData > > phi;
			ierr = GetGaussDiffDataVector(field_name,phi); CHKERRQ(ierr);
			

				for (unsigned int gg=0; gg<phi.size(); gg++) {
					normalized_phi[gg].resize(field_rank,3);
					for (int rank = 0; rank < field_rank; rank ++) {
						for (int ii=0; ii<3; ii++) {
							if (normalize == true) {
								normalized_phi[gg](rank,ii) = -phi[gg](rank,ii)/sqrt(pow(phi[gg](rank,0),2)+pow(phi[gg](rank,1),2)+pow(phi[gg](rank,2),2));}
							else{
								normalized_phi[gg](rank,ii) = phi[gg](rank,ii);
							}}}}
			PetscFunctionReturn(0);
		}
		
		vector<double> g_NTET,g_NTRI;
		const double* G_W_TET;
	
		PetscErrorCode preProcess() {
			PetscFunctionBegin;
			g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_W_TET = G_TET_W45;
			
			PetscFunctionReturn(0);
		}
		
		PetscErrorCode operator()() {
			PetscFunctionBegin;
			ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

			MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator field_it;
			field_it = moabfields->get<FieldName_mi_tag>().find(field_name);
			const int rank_size = field_it->get_max_rank();
			
			vector< ublas::matrix< double > > normalized_phi;
			normalized_phi.resize(coords_at_Gauss_nodes.size());
			ierr = ComputeFibreDirection(normalized_phi); CHKERRQ(ierr);

			double fibreVectorArray[3*rank_size*coords_at_Gauss_nodes.size()];
			for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
				for (int rank = 0; rank<rank_size; rank++) {
					for (int ii=0; ii<3; ii++) {
						fibreVectorArray[3*rank_size*gg+3*rank+ii]=normalized_phi[gg](rank,ii);
					}
				}
			}

			EntityHandle fe_handle1 = fe_ptr->get_ent();
			EntityHandle fe_handle2 = fe_ptr->get_parent_ent();
			
			int fibreVectorArraySize = 3*rank_size*coords_at_Gauss_nodes.size();
			void const* tag_fibre_vector_data[] = {fibreVectorArray};
					
			rval = moab.tag_set_by_ptr(th_field_grad,&fe_handle1,1,tag_fibre_vector_data,&fibreVectorArraySize); CHKERR_PETSC(rval);
			rval = moab.tag_set_by_ptr(th_field_grad,&fe_handle2,1,tag_fibre_vector_data,&fibreVectorArraySize); CHKERR_PETSC(rval);
			
			//Access DATA
//			double *data;
//			int arraysize_fibre;
//			rval = moab.tag_get_by_ptr(th_field_grad,&fe_handle1,1,(const void**)&data,&arraysize_fibre); CHKERR_PETSC(rval);
			
			ierr = OpStudentEnd(); CHKERRQ(ierr);

			PetscFunctionReturn(0); }
		
		PetscErrorCode postProcess() {
			PetscFunctionBegin;
			PetscFunctionReturn(0);
		}
	};

};

#endif // __POTENTIALFLOWFEMETHOD_HPP__

