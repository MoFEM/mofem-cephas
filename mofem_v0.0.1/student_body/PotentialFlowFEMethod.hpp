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


#ifndef __POTENTIALFLOWFEMETHOD_HPP__
#define __POTENTIALFLOWFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "FEMethod_UpLevelStudent.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

namespace MoFEM {

struct BernoullyEquations {

    //Look to Zienkiewicz Book for details, Volume 3, Section 4.2

    BernoullyEquations(): set_pressure(NULL) {
      rho = 1;
      a.resize(3);
      ublas::noalias(a) = ublas::zero_vector<FieldData>(3);
    }

    double rho;
    ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > a;
    PetscErrorCode set_params(double _rho,ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > _a) {
      PetscFunctionBegin;
      ublas::noalias(a) = _a;
      rho = _rho;
      PetscFunctionReturn(0);
    }

    typedef NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator dof_it;
    PetscErrorCode (*set_pressure)
      (dof_it dof,ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > &coords,
      FieldData &pressure);

    PetscErrorCode compute_pressure(
      FieldData potential,
      FieldData d_phi_dt,ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > &u,
      FieldData &pressure) {
	PetscFunctionBegin;
	pressure = rho*( d_phi_dt - 0.5*inner_prod(u,u) - potential );
	if(pressure != pressure) {
	  SETERRQ1(PETSC_COMM_SELF,1,"Huston we have a problem, pressure is %6.4e",pressure);
	}
	PetscFunctionReturn(0);
    }

};

struct LaplacianElem: public FEMethod_UpLevelStudent {

    FieldInterface& mField;
    Mat A;
    Vec F;
    LaplacianElem(FieldInterface& _mField,Mat _A,Vec _F): 
      FEMethod_UpLevelStudent(_mField.get_moab()),mField(_mField),A(_A),F(_F) {
      body_velocity = ublas::zero_vector<FieldData>(3);
    }; 

    const double *G_TET_W,*G_TRI_W;
    vector<double> g_NTET,g_NTRI;
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(45*4);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_TET_W = G_TET_W45;
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
      G_TRI_W = G_TRI_W13;
      ierr = VecSetOption(F,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    vector<vector<DofIdx> > RowGlobDofs;
    vector<vector< ublas::matrix<FieldData> > > diffRowNMatrix;

    PetscErrorCode get_ShapeFunctionsAndIndices() {
      PetscFunctionBegin;

      RowGlobDofs.resize(1+6+4+1);
      diffRowNMatrix.resize(1+6+4+1);
      ierr = GetRowGlobalIndices("POTENTIAL_FIELD",RowGlobDofs[0]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",diffRowNMatrix[0]); CHKERRQ(ierr);
      int ee = 0;
      for(;ee<6;ee++) {
	ierr = GetRowGlobalIndices("POTENTIAL_FIELD",MBEDGE,RowGlobDofs[1+ee],ee); CHKERRQ(ierr);
	if(RowGlobDofs[1+ee].size()>0) {
	  ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",MBEDGE,diffRowNMatrix[1+ee],ee); CHKERRQ(ierr);
	} 
      }
      int ff = 0;
      for(;ff<4;ff++) {
	ierr = GetRowGlobalIndices("POTENTIAL_FIELD",MBTRI,RowGlobDofs[1+6+ff],ff); CHKERRQ(ierr);
	if(RowGlobDofs[1+6+ff].size()>0) {
	  ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",MBTRI,diffRowNMatrix[1+6+ff],ff); CHKERRQ(ierr);
	}
      }
      ierr = GetRowGlobalIndices("POTENTIAL_FIELD",MBTET,RowGlobDofs[1+6+4]); CHKERRQ(ierr);
      if(RowGlobDofs[1+6+4].size()>0) {
	ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",MBTET,diffRowNMatrix[1+6+4]); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode compute_LHS(FieldData a = 1) {
      PetscFunctionBegin;

      ublas::matrix<FieldData> K;
      for(int rr = 0;rr<(1+6+4+1); rr++) {
	if(RowGlobDofs[rr].size()==0) continue;
	for(int cc = 0;cc<(1+6+4+1);cc++) {
	  if(RowGlobDofs[cc].size()==0) continue;
	  K.resize(RowGlobDofs[rr].size(),RowGlobDofs[cc].size());
	  for(unsigned int gg = 0;gg<g_NTET.size()/4;gg++) {
	    ublas::noalias(K) = prod( trans(diffRowNMatrix[rr][gg]),diffRowNMatrix[cc][gg] ); 
	    K *= a*V*G_TET_W[gg];
	    ierr = MatSetValues(A,
	      RowGlobDofs[rr].size(),&(RowGlobDofs[rr])[0],
	      RowGlobDofs[cc].size(),&(RowGlobDofs[cc])[0],
	      &(K.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }

      PetscFunctionReturn(0);
    }

    ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > body_velocity;

    PetscErrorCode set_bodyVelocity(FieldData ux,FieldData uy,FieldData uz) {
      PetscFunctionBegin;
      body_velocity.resize(3);
      body_velocity[0] = ux;
      body_velocity[1] = uy;
      body_velocity[2] = uz;
      ostringstream ss;
      ss << "Set body velocity: ";
      ss << body_velocity;
      ss << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      PetscFunctionReturn(0);
    }

    PetscErrorCode set_bodyVelocity(ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > _body_velocity) {
      PetscFunctionBegin;
      body_velocity = _body_velocity;
      ostringstream ss;
      ss << "Set body velocity: ";
      ss << body_velocity;
      ss << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
      PetscFunctionReturn(0);
    }

    PetscErrorCode compute_SurfaceRHS() {
      PetscFunctionBegin;

      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {

	string name;
	if( (it->get_CubitBCType()&Cubit_BC_bitset(PressureSet)).none() ) {
	  name = it->get_Cubit_name();
	  if( name!="NormalVelocity" ) {
	    continue;
	  } 
	}

	double flux = 0;
	if( (it->get_CubitBCType()&Cubit_BC_bitset(PressureSet)).any() ) {
	  pressure_cubit_bc_data mydata;
	  ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
	  flux = mydata.data.value1;
	}

	SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
	SideNumber_multiIndex::nth_index<1>::type::iterator sit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
	SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
	Range meshsets;
	rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
	meshsets.insert(it->meshset); // check if children of pressure bc triangles are not in main cubit meshset
	for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
	SideNumber_multiIndex::nth_index<1>::type::iterator siit = sit;
	for(;siit!=hi_siit;siit++) {

	  if(!moab.contains_entities(*mit,&siit->ent,1)) continue;

	  ierr = ShapeFunctions_TRI(siit->ent,g_NTRI);  CHKERRQ(ierr);
	  
	  const EntityHandle *conn_face;
	  int num_nodes_face;
	  rval = moab.get_connectivity(siit->ent,conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
	  double coords_face[9];
	  rval = moab.get_coords(conn_face,num_nodes_face,coords_face); CHKERR_PETSC(rval);
	  ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);
	  ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > normal(3);
	  ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,&*normal.data().begin()); CHKERRQ(ierr);
	  double area = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;

	  if( name == "NormalVelocity" ) {
	    flux = -inner_prod(body_velocity,normal)/(2*area);
	    //cerr << body_velocity << endl;
	    //cerr << flux << endl;
	  }

        
	  //nodes
	  vector<DofIdx>& RowGlob_nodes = RowGlobDofs[0];
	  vector< ublas::matrix<FieldData> > FaceNMatrix_nodes;
	  ierr = GetGaussRowFaceNMatrix(siit->ent,"POTENTIAL_FIELD",FaceNMatrix_nodes,MBVERTEX); CHKERRQ(ierr);

	  //edges
	  vector<vector<ublas::matrix<FieldData> > > FaceNMatrix_edges(6);
	  SideNumber_multiIndex::nth_index<1>::type::iterator siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siiit = side_table.get<1>().upper_bound(boost::make_tuple(MBEDGE,6));
	  for(;siiit!=hi_siiit;siiit++) {
	    if(RowGlobDofs[1+siiit->side_number].size()>0) {
	      ierr = GetGaussRowFaceNMatrix(siit->ent,"POTENTIAL_FIELD",FaceNMatrix_edges[siiit->side_number],MBEDGE,siiit->ent); CHKERRQ(ierr);
	    }
	  }

	  //faces
	  vector<DofIdx>& RowGlob_face = RowGlobDofs[1+6+siit->side_number];
	  vector< ublas::matrix<FieldData> > FaceNMatrix_face;
	  if(RowGlob_face.size()>0) {
	    ierr = GetGaussRowFaceNMatrix(siit->ent,"POTENTIAL_FIELD",FaceNMatrix_face,MBTRI); CHKERRQ(ierr);
	  }

	  //calulate & assemble

	  //nodes
	  ublas::vector<FieldData> f_ext_nodes = ublas::zero_vector<FieldData>(FaceNMatrix_nodes[0].size2());
	  int g_dim = get_dim_gNTRI();
	  for(int gg = 0;gg<g_dim;gg++) {
	    double w = flux*area*G_TRI_W[gg];
	    f_ext_nodes += w*ublas::matrix_row<ublas::matrix<FieldData> >(FaceNMatrix_nodes[gg],0);
	  }
	  ierr = VecSetValues(F,RowGlob_nodes.size(),&(RowGlob_nodes)[0],&(f_ext_nodes.data())[0],ADD_VALUES); CHKERRQ(ierr);

	  //edges
	  siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	  for(;siiit!=hi_siiit;siiit++) {
	    vector<DofIdx>& RowGlob_edge = RowGlobDofs[1+siiit->side_number];
	    if(RowGlob_edge.size()>0) {
	      vector<ublas::matrix<FieldData> >& FaceNMatrix_edge = FaceNMatrix_edges[siiit->side_number];
	      ublas::vector<FieldData> f_ext_edges = ublas::zero_vector<FieldData>(FaceNMatrix_edge[0].size2());
	      for(int gg = 0;gg<g_dim;gg++) {
		double w = flux*area*G_TRI_W[gg];
		f_ext_edges += w*ublas::matrix_row<ublas::matrix<FieldData> >(FaceNMatrix_edge[gg],0);
	      }
	      if(RowGlob_edge.size()!=f_ext_edges.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_edge.size()!=f_ext_edges.size()");
	      ierr = VecSetValues(F,RowGlob_edge.size(),&(RowGlob_edge[0]),&(f_ext_edges.data())[0],ADD_VALUES); CHKERRQ(ierr);
	    }
	  }

	  //face     
	  if(RowGlob_face.size()>0) {
	    ublas::vector<FieldData> f_ext_faces = ublas::zero_vector<FieldData>(FaceNMatrix_face[0].size2());
	    for(int gg = 0;gg<g_dim;gg++) {
	      double w = flux*area*G_TRI_W[gg];
	      f_ext_faces += w*ublas::matrix_row<ublas::matrix<FieldData> >(FaceNMatrix_face[gg],0);
	    }
	    if(RowGlob_face.size()!=f_ext_faces.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_face.size()!=f_ext_faces.size()");
	    ierr = VecSetValues(F,RowGlob_face.size(),&(RowGlob_face)[0],&(f_ext_faces.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	
	}

      }}

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      ierr = get_ShapeFunctionsAndIndices(); CHKERRQ(ierr);
      ierr = compute_LHS(); CHKERRQ(ierr);
      ierr = compute_SurfaceRHS(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      vector<DofIdx> zero_pressure_dofs;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|UnknownCubitName,it)) {
	if(it->get_Cubit_name() != "ZeroPressure") continue;
	Range nodes;
	rval = moab.get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
	Range edges;
	rval = moab.get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
	Range tris;
	rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
	Range adj;
	rval = moab.get_connectivity(tris,adj,true); CHKERR_PETSC(rval);
	nodes.insert(adj.begin(),adj.end());
	rval = moab.get_connectivity(edges,adj,true); CHKERR_PETSC(rval);
	nodes.insert(adj.begin(),adj.end());
	rval = moab.get_adjacencies(tris,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
	for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
	  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
	    ierr = VecSetValue(F,dof->get_petsc_gloabl_dof_idx(),0,INSERT_VALUES); CHKERRQ(ierr);
	    zero_pressure_dofs.push_back(dof->get_petsc_gloabl_dof_idx());
	}}
	for(Range::iterator eit = edges.begin();eit!=edges.end();eit++) {
	  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*eit,dof)) {
	    ierr = VecSetValue(F,dof->get_petsc_gloabl_dof_idx(),0,INSERT_VALUES); CHKERRQ(ierr);
	    zero_pressure_dofs.push_back(dof->get_petsc_gloabl_dof_idx());
	}}
	for(Range::iterator tit = tris.begin();tit!=tris.end();tit++) {
	  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*tit,dof)) {
	    ierr = VecSetValue(F,dof->get_petsc_gloabl_dof_idx(),0,INSERT_VALUES); CHKERRQ(ierr);
	    zero_pressure_dofs.push_back(dof->get_petsc_gloabl_dof_idx());
	}}
      }
      ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatZeroRowsColumns(A,zero_pressure_dofs.size(),&zero_pressure_dofs[0],1,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

struct SteadyBernoullyElem: public FEMethod_UpLevelStudent,BernoullyEquations {

    FieldInterface& mField;
    Mat A;
    Vec F;
    SteadyBernoullyElem(FieldInterface& _mField,Mat _A,Vec _F): 
      FEMethod_UpLevelStudent(_mField.get_moab()),mField(_mField),A(_A),F(_F) {}; 

    const double *G_TET_W,*G_TRI_W;
    vector<double> g_NTET,g_NTRI;
    FieldSpace pressure_space;
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(45*4);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_TET_W = G_TET_W45;
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
      G_TRI_W = G_TRI_W13;
      ierr = VecSetOption(F,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);


      MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find("PRESSURE_FIELD");
      if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
      pressure_space = fiit->get_space();
    
      PetscFunctionReturn(0);
    }

    vector<vector<DofIdx> > RowGlobDofs;
    vector<vector< ublas::matrix<FieldData> > > RowNMatrix;
    vector< ublas::matrix<FieldData> > negativeVelocities;
    
    PetscErrorCode get_ShapeFunctionsAndIndices() {
      PetscFunctionBegin;

      RowGlobDofs.resize(1+6+4+1);
      RowNMatrix.resize(1+6+4+1);
      for(int rr= 0;rr<1+6+4+1;rr++) {
	RowGlobDofs[rr].resize(0);
      }

      switch(pressure_space) {
	case L2:
	  ierr = GetRowGlobalIndices("PRESSURE_FIELD",MBTET,RowGlobDofs[1+6+4]); CHKERRQ(ierr);
	  if(RowGlobDofs[1+6+4].size()>0) {
	    ierr = GetGaussRowNMatrix("PRESSURE_FIELD",MBTET,RowNMatrix[1+6+4]); CHKERRQ(ierr);
	  }
	  break;
	case H1: {
	  ierr = GetRowGlobalIndices("PRESSURE_FIELD",RowGlobDofs[0]); CHKERRQ(ierr);
	  ierr = GetGaussRowNMatrix("PRESSURE_FIELD",RowNMatrix[0]); CHKERRQ(ierr);
	  int ee = 0;
	  for(;ee<6;ee++) {
	    ierr = GetRowGlobalIndices("PRESSURE_FIELD",MBEDGE,RowGlobDofs[1+ee],ee); CHKERRQ(ierr);
	    if(RowGlobDofs[1+ee].size()>0) {
	      ierr = GetGaussRowNMatrix("PRESSURE_FIELD",MBEDGE,RowNMatrix[1+ee],ee); CHKERRQ(ierr);
	    } 
	  }
	  int ff = 0;
	  for(;ff<4;ff++) {
	    ierr = GetRowGlobalIndices("PRESSURE_FIELD",MBTRI,RowGlobDofs[1+6+ff],ff); CHKERRQ(ierr);
	    if(RowGlobDofs[1+6+ff].size()>0) {
	      ierr = GetGaussRowNMatrix("PRESSURE_FIELD",MBTRI,RowNMatrix[1+6+ff],ff); CHKERRQ(ierr);
	    }
	  }
	  ierr = GetRowGlobalIndices("PRESSURE_FIELD",MBTET,RowGlobDofs[1+6+4]); CHKERRQ(ierr);
	  if(RowGlobDofs[1+6+4].size()>0) {
	    ierr = GetGaussRowNMatrix("PRESSURE_FIELD",MBTET,RowNMatrix[1+6+4]); CHKERRQ(ierr);
	  }
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
      }

      ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",negativeVelocities); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


    PetscErrorCode compute_LHS() {
      PetscFunctionBegin;

      ublas::matrix<FieldData> K;
      for(int rr = 0;rr<(1+6+4+1); rr++) {
	if(RowGlobDofs[rr].size()==0) continue;
	for(int cc = 0;cc<(1+6+4+1);cc++) {
	  if(RowGlobDofs[cc].size()==0) continue;
	  K.resize(RowGlobDofs[rr].size(),RowGlobDofs[cc].size());
	  for(unsigned int gg = 0;gg<g_NTET.size()/4;gg++) {
	    ublas::noalias(K) = prod( trans(RowNMatrix[rr][gg]),RowNMatrix[cc][gg] ); 
	    K *= V*G_TET_W[gg];
	    ierr = MatSetValues(A,
	      RowGlobDofs[rr].size(),&(RowGlobDofs[rr])[0],
	      RowGlobDofs[cc].size(),&(RowGlobDofs[cc])[0],
	      &(K.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode compute_RHS() {
      PetscFunctionBegin;

      for(int rr = 0;rr<(1+6+4+1); rr++) {
	if(RowGlobDofs[rr].size()==0) continue;
	ublas::vector<FieldData> f(RowGlobDofs[rr].size());
	ublas::noalias(f) = ublas::zero_vector<FieldData>(RowGlobDofs[rr].size());
	for(unsigned int gg = 0;gg<g_NTET.size()/4;gg++) {
	  ublas::matrix_row<ublas::matrix<FieldData> > mr(negativeVelocities[gg],0);
	  ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > velocities(3);
	  velocities = -mr;
	  //cerr << velocities << endl;
	  FieldData pressure ;
	  ierr = compute_pressure(0,0,velocities,pressure); CHKERRQ(ierr);
	  //cerr << pressure << endl;
	  for(unsigned int dd = 0;dd<RowGlobDofs[rr].size();dd++) {
	    f[dd] += V*G_TET_W[gg]*(RowNMatrix[rr][gg])(0,dd)*pressure;
	  }
	  if(f.size()!=RowGlobDofs[rr].size()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	}
	ierr = VecSetValues(F,RowGlobDofs[rr].size(),&(RowGlobDofs[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      ierr = get_ShapeFunctionsAndIndices(); CHKERRQ(ierr);
      ierr = compute_LHS(); CHKERRQ(ierr);
      ierr = compute_RHS(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

};

struct SurfaceForces: public FEMethod_UpLevelStudent,BernoullyEquations {

  FieldInterface& mField;
  Vec total_Force;

  SurfaceForces(FieldInterface& _mField): 
      FEMethod_UpLevelStudent(_mField.get_moab()),mField(_mField) {}; 

  PetscErrorCode create_totalForceVector(Vec *_total_Force) {
    PetscFunctionBegin;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      ierr = VecCreateGhost(PETSC_COMM_WORLD,3,3,0,PETSC_NULL,_total_Force); CHKERRQ(ierr);
    } else {
      int ghosts[3] = {0,1,2};
      ierr = VecCreateGhost(PETSC_COMM_WORLD,0,3,3,ghosts,_total_Force); CHKERRQ(ierr);
    }
    total_Force = *_total_Force;
    PetscFunctionReturn(0);
  }

  PetscErrorCode set_totalForceVector(Vec _total_Force) {
    PetscFunctionBegin;
    total_Force = _total_Force;
    PetscFunctionReturn(0);
  }

  Tag th_normal;
  EntityHandle out_meshset;
  const double *G_TET_W,*G_TRI_W;
  vector<double> g_NTET,g_NTRI;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    g_NTET.resize(45*4);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    G_TET_W = G_TET_W45;
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
    G_TRI_W = G_TRI_W13;
    ierr = VecZeroEntries(total_Force); CHKERRQ(ierr);

    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    double def[3] = {0,0,0};
    rval = moab.tag_get_handle("BODY_NORMAL",3,MB_TYPE_DOUBLE,th_normal,MB_TAG_CREAT|MB_TAG_SPARSE,def); CHKERR_PETSC(rval);

    PetscFunctionReturn(0);
  }

  PetscErrorCode computeResultantSurfaceForces(
    vector<ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > > &faces_force) {
    PetscFunctionBegin;
  
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,it)) {

      string name;
      name = it->get_Cubit_name();
      if( name!="BodySurface" ) {
	continue;
      } 

      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator sit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
      Range meshsets;
      rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
      meshsets.insert(it->meshset); // check if children of pressure bc triangles are not in main cubit meshset
      for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
	SideNumber_multiIndex::nth_index<1>::type::iterator siit = sit;
	for(;siit!=hi_siit;siit++) {

	  if(!moab.contains_entities(*mit,&siit->ent,1)) continue;

	  ierr = ShapeFunctions_TRI(siit->ent,g_NTRI);  CHKERRQ(ierr);
	  
	  const EntityHandle *conn_face;
	  int num_nodes_face;
	  rval = moab.get_connectivity(siit->ent,conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
	  double coords_face[9];
	  rval = moab.get_coords(conn_face,num_nodes_face,coords_face); CHKERR_PETSC(rval);
	  ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);
	  ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > normal(3);
	  ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,&*normal.data().begin()); CHKERRQ(ierr);
	  double area = ublas::norm_2(normal)*0.5;
	  normal /= ublas::norm_2(normal);
	  //normal *= siit->sense
	  rval = moab.tag_set_data(th_normal,&(siit->ent),1,&*normal.data().begin()); CHKERR_PETSC(rval);
	  rval = moab.add_entities(out_meshset,&(siit->ent),1); CHKERR_PETSC(rval);

	  vector<ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > > coords_at_gauss_points;
	  int g_dim = get_dim_gNTRI();
	  coords_at_gauss_points.resize(g_dim);
	  for(int gg = 0;gg<g_dim;gg++) {
	    coords_at_gauss_points[gg].resize(3);
	    coords_at_gauss_points[gg][0] = cblas_ddot(3,&coords_face[0],3,&g_NTRI[3*gg],1);
	    coords_at_gauss_points[gg][1] = cblas_ddot(3,&coords_face[1],3,&g_NTRI[3*gg],1);
	    coords_at_gauss_points[gg][2] = cblas_ddot(3,&coords_face[2],3,&g_NTRI[3*gg],1);
	  }

	  vector<ublas::vector<FieldData> > pressure_at_gauss_points;
	  Data_at_FaceGaussPoints(siit->ent,"PRESSURE_FIELD",pressure_at_gauss_points);

	  for(int gg = 0;gg<g_dim;gg++) {
	    double w = G_TRI_W[gg];
	    //cerr << pressure_at_gauss_points[gg] << endl;
	    faces_force[siit->side_number] += w*pressure_at_gauss_points[gg][0]*normal*area;
	  }

	}
      }
    }

    PetscFunctionReturn(0);

  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

    vector<ublas::vector<FieldData,ublas::bounded_array<FieldData,3> > > faces_force;
    faces_force.resize(4);
    for(int ff = 0;ff<4;ff++) {
      faces_force[ff] = ublas::zero_vector<FieldData>(3);
    }
    ierr = computeResultantSurfaceForces(faces_force); CHKERRQ(ierr);
    
    int forces_idx[] = { 0,1,2 };
    for(int ff = 0;ff<4;ff++) {
      ierr = VecSetValues(total_Force,3,forces_idx,&*faces_force[ff].data().begin(),ADD_VALUES); CHKERRQ(ierr);
    }

    ierr = OpStudentEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    ierr = VecGhostUpdateBegin(total_Force,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(total_Force,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(total_Force); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(total_Force); CHKERRQ(ierr);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      rval = moab.write_file("body_normals.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    }

    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);


    PetscFunctionReturn(0);
  }

};

struct PostProcPotentialFlowOnRefMesh: public PostProcDisplacemenysAndStarinOnRefMesh {

  Tag th_phi,th_p,th_u;
  PostProcPotentialFlowOnRefMesh(Interface& _moab): PostProcDisplacemenysAndStarinOnRefMesh(_moab) {
    double def_VAL2[3] = { 0.0, 0.0, 0.0 };
    rval = moab_post_proc.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("P",1,MB_TYPE_DOUBLE,th_p,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("U",3,MB_TYPE_DOUBLE,th_u,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
  }

    PetscErrorCode do_operator() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      Range ref_nodes;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
      if(4*ref_nodes.size()!=g_NTET.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(ref_nodes.size()!=coords_at_Gauss_nodes.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      Range::iterator nit = ref_nodes.begin();
      node_map.clear();
      for(int nn = 0;nit!=ref_nodes.end();nit++,nn++) {
	EntityHandle &node = node_map[*nit];
	rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nn]).data()[0],node); CHKERR_PETSC(rval);
      }
      Range ref_tets;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBTET,ref_tets); CHKERR_PETSC(rval);
      Range::iterator tit = ref_tets.begin();
      for(;tit!=ref_tets.end();tit++) {
	const EntityHandle *conn_ref;
        int num_nodes;
	rval = moab_ref.get_connectivity(*tit,conn_ref,num_nodes,true); CHKERR_PETSC(rval);
	EntityHandle conn_post_proc[num_nodes];
	for(int nn = 0;nn<num_nodes;nn++) {
	  conn_post_proc[nn] = node_map[conn_ref[nn]];
	}
	EntityHandle ref_tet;
	rval = moab_post_proc.create_element(MBTET,conn_post_proc,4,ref_tet); CHKERR_PETSC(rval);
      }

      PetscFunctionReturn(0);
    }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    //Loop over elements
    ierr = do_operator(); CHKERRQ(ierr);

    //Strains to Nodes in PostProc Mesh: create vector containing matrices
    vector< ublas::matrix< FieldData > > negativeVelocities;
    vector< ublas::vector< FieldData > > phi;
    vector< ublas::vector< FieldData > > p;

    ierr = GetGaussDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);
    ierr = GetGaussDataVector("PRESSURE_FIELD",p); CHKERRQ(ierr);
    ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",negativeVelocities); CHKERRQ(ierr);

    map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
    for(;mit!=node_map.end();mit++) {
      
      int gg = distance(node_map.begin(),mit);

      negativeVelocities[gg] *= -1;
      rval = moab_post_proc.tag_set_data(th_u,&mit->second,1,&(negativeVelocities[gg].data()[0])); CHKERR_PETSC(rval);
      rval = moab_post_proc.tag_set_data(th_phi,&mit->second,1,&(phi[gg].data()[0])); CHKERR_PETSC(rval);
      rval = moab_post_proc.tag_set_data(th_p,&mit->second,1,&(p[gg].data()[0])); CHKERR_PETSC(rval);

    }

    PetscFunctionReturn(0);

  }

};

};

#endif // __POTENTIALFLOWFEMETHOD_HPP__

