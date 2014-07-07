/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ElasticFE_RVELagrange_RigidBodyTranslation_HPP__
#define __ElasticFE_RVELagrange_RigidBodyTranslation__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"

namespace MoFEM {
  
  struct ElasticFE_RVELagrange_RigidBodyTranslation: public ElasticFE_RVELagrange_Disp {
    
    ElasticFE_RVELagrange_RigidBodyTranslation(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F, ublas::vector<FieldData> _applied_strain):
    ElasticFE_RVELagrange_Disp(_mField, _dirihlet_ptr,_Aij, _D, _F, _applied_strain){};
    
    virtual PetscErrorCode GetN_and_Indices() {
      PetscFunctionBegin;
      
      //Find out indices for row and column for nodes on the surface, i.e. triangles
      row_mat = 0;
      RowGlob[row_mat].resize(3);
      ColGlob[row_mat].resize(9);
      
      typedef FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator row_dofs_iterator;
      typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
      const EntityHandle* conn_face;
      int num_nodes;
      EntityHandle face_tri;  face_tri=fe_ptr->get_ent();
      rval = moab.get_connectivity(face_tri,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      //        cout<< "num_nodes ="<<num_nodes << endl;
      //        cout<< "conn_face ="<<conn_face << endl;
      
      //minimum and maximum rows indices for each node on the surface
      row_dofs_iterator niit,hi_niit;   //iterator for rows
      niit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("Lagrange_mul_disp_rigid_trans");
      hi_niit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("Lagrange_mul_disp_rigid_trans");
      int nn = 0;
      for(;niit!=hi_niit;niit++) {
        RowGlob[row_mat][nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
      }
      
      nn = 0;
      for(;nn<3;nn++) {
        dofs_iterator col_niit,hi_col_niit;  // iterator for columns
        
        //minimum and maximum row and column indices for each node on the surface
        col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
        hi_col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
        
        // two different loops, i.e. one for row and one for column (may be need it for multiphysics problems)
        for(;col_niit!=hi_col_niit;col_niit++) {
          ColGlob[row_mat][nn*col_niit->get_max_rank()+col_niit->get_dof_rank()] = col_niit->get_petsc_gloabl_dof_idx();
        }
      }
      //        cout<<"\nFor nodes "<<endl;
      //        cout<<"\n RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
      //        for(int ii=0; ii<RowGlob[row_mat].size(); ii++) cout<<RowGlob[row_mat][ii]<<" ";
      //        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
      //        for(int ii=0; ii<ColGlob[row_mat].size(); ii++) cout<<ColGlob[row_mat][ii]<<" ";
      //        cout<<"\n\n\n"<<endl;
      //        //Stop code
      //        std::string wait;
      //        std::cin >> wait;
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode Lhs() {
      PetscFunctionBegin;
      
      ublas::matrix<FieldData> Mat_face;          Mat_face.resize(3,9);           Mat_face.clear();
      ublas::matrix<FieldData> Mat_face_Tran;     Mat_face_Tran.resize(9,3);      Mat_face_Tran.clear();
      //cout<<"Mat_face "<<Mat_face<<endl;
      for(int nn=0; nn<3; nn++){
        Mat_face(0,3*nn+0)=1.0;  Mat_face(1,3*nn+1)=1.0;   Mat_face(2,3*nn+2)=1.0;
      }
      //        cout<<"Mat_face "<< Mat_face << endl<<endl;
      //Assembly C1 with size (3 x 9) for each node
      ierr = MatSetValues(Aij,RowGlob[0].size(),&(RowGlob[0])[0],ColGlob[0].size(),&(ColGlob[0])[0],&(Mat_face.data())[0],INSERT_VALUES ); CHKERRQ(ierr);
      
      //Assembly C1T with size (9 x 3) for each node
      Mat_face_Tran=trans(Mat_face);
      ierr = MatSetValues(Aij,ColGlob[0].size(),&(ColGlob[0])[0],RowGlob[0].size(),&(RowGlob[0])[0],&(Mat_face_Tran.data())[0],INSERT_VALUES ); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      //        cout<<"Hi from class RigidBodyMotion"<<endl;
      ierr = GetN_and_Indices(); CHKERRQ(ierr);
      ierr = Lhs(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    
  };
  
  
}

#endif //__ElasticFE_RVELagrange_Periodic_RigidBodyMotion
