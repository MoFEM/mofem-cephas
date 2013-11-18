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


#ifndef __THERMALFLOWFEMETHOD_HPP__
#define __THERMALFLOWFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {

struct ThermalElem: public FEMethod_UpLevelStudent {

    FieldInterface& mField;
    Mat A;
    Vec F;
    ThermalElem(FieldInterface& _mField,Mat _A,Vec _F):
      FEMethod_UpLevelStudent(_mField.get_moab()),mField(_mField),A(_A),F(_F) {}; 

    const double *G_TET_W,*G_TRI_W;
    vector<double> g_NTET,g_NTRI;
    PetscErrorCode preProcess() {
      PetscFunctionBegin;    // PetSc start of a function
      g_NTET.resize(45*4);   // shape function for tetrahedron at each (45) gauss point, i.e. 45*4
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);  // G_TET_X45,G_TET_Y45,G_TET_Z45 these are coordinates of gauss points given in FEM.h
      G_TET_W = G_TET_W45; // weight of gauss points
      g_NTRI.resize(3*13); // shape function for triangular elements (faces of tetrahedron), i.e. 13*3 (13 gauss points are used here)
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); //G_TRI_X13, G_TRI_Y13 are coordinates of gauss points given in FEM.h
      G_TRI_W = G_TRI_W13; //weight of triangular guass points
      ierr = VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr); //The negative entries will be ignored
      PetscFunctionReturn(0);   // PetSc end of a function
    }

    vector<vector<DofIdx> > RowGlobDofs;    // DofIdx ?
    vector<vector< ublas::matrix<FieldData> > > diffRowNMatrix;
    
    
     PetscErrorCode operator()() {
        
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);  //function to run at the begining of the function when tet ele is evaluated

      RowGlobDofs.resize(1+6+4+1);
      diffRowNMatrix.resize(1+6+4+1);
      ierr = GetRowGlobalIndices("THERMAL_FIELD",RowGlobDofs[0]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("THERMAL_FIELD",diffRowNMatrix[0]); CHKERRQ(ierr);
         
         //cout<<"size [0]: "<<RowGlobDofs[0].size()<<"\n";
         //cout<<"GetGaussRowDiffNMatrix: "<<diffRowNMatrix[0][0]<<"\n";
         //cout<<"GetGaussRowDiffNMatrix: "<<diffRowNMatrix[0][1]<<"\n";
         
         //cout<<"GetGaussRowDiffNMatrix[0]: "<<diffRowNMatrix[0].size()<<"\n";
         //for(int zz=0; zz<diffRowNMatrix[0].size(); zz++){
         //    cout<<"GetGaussRowDiffNMatrix[0][0]: "<<diffRowNMatrix[0][zz]<<"\n";}
             
         //cout<<"GetGaussRowDiffNMatrix [0][0]: "<<GetGaussRowDiffNMatrix[0][0].size()<<"\n";
         
         
         //for(int zz=0; zz<RowGlobDofs[0].size(); zz++){
         //   cout<<RowGlobDofs[0][zz]<<"\n";}
         //cout<<"\n\n";
    
      //cout<<diffRowNMatrix[0][0]<<"\n";
         
         
      int ee = 0;
      for(;ee<6;ee++) {
	ierr = GetRowGlobalIndices("THERMAL_FIELD",MBEDGE,RowGlobDofs[1+ee],ee); CHKERRQ(ierr);
    
          //cout<<"size [1-6]: "<<RowGlobDofs[1+ee].size()<<"\n";
          
	if(RowGlobDofs[1+ee].size()>0) {
	  ierr = GetGaussRowDiffNMatrix("THERMAL_FIELD",MBEDGE,diffRowNMatrix[1+ee],ee); CHKERRQ(ierr);
	} 
      }
      int ff = 0;
      for(;ff<4;ff++) {
	ierr = GetRowGlobalIndices("THERMAL_FIELD",MBTRI,RowGlobDofs[1+6+ff],ff); CHKERRQ(ierr);
          
          //cout<<"size [7-10]: "<<RowGlobDofs[1+6+ff].size()<<"\n";
          
	if(RowGlobDofs[1+6+ff].size()>0) {
	  ierr = GetGaussRowDiffNMatrix("THERMAL_FIELD",MBTRI,diffRowNMatrix[1+6+ff],ff); CHKERRQ(ierr);
	}
      }
      ierr = GetRowGlobalIndices("THERMAL_FIELD",MBTET,RowGlobDofs[1+6+4]); CHKERRQ(ierr);
         
         //cout<<"size [11]: "<<RowGlobDofs[1+6+4].size()<<"\n";
         
      if(RowGlobDofs[1+6+4].size()>0) {
	ierr = GetGaussRowDiffNMatrix("THERMAL_FIELD",MBTET,diffRowNMatrix[1+6+4]); CHKERRQ(ierr);
      }


//To apply the homogenious Dirichlet boundary conditions       UnknownCubitName  TemperatureSet
    //cout<<"temp set in 1"<<"\n";
for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|TemperatureSet,it)) {
	//if(it->get_Cubit_name() != "ZeroTemp") continue;
    //cout<<it->get_Cubit_name()<<"\n";
    //cout<<"temp set in"<<"\n";

	SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
	SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX,0));
	SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBVERTEX,4));
    
    //cout<<"it->meshset: "<<it->meshset<<"\n";
    
    
	Range meshsets;
	rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
    //rval = moab.list_entities(meshsets); CHKERR_PETSC(rval);
    
  for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
      
      //cout<<"it->meshset: "<<it->meshset<<"\n";
      //cout<<"meshsets.begin(): "<<meshsets.begin()<<"\n";
      //cout<<"meshsets.end(): "<<meshsets.end()<<"\n";
      
	  
      //int zz=1;
      //cout<<"for loop"<<"\n";
      for(;siit!=hi_siit;siit++) {
        //cout<<"zz "<<zz<<"\n"; zz++;
        //cout<<"mit "<<*mit<<"\n";
        //cout<<"&siit->ent "<<siit->ent<<"\n";
          if(!moab.contains_entities(*mit,&siit->ent,1)) continue;
        
          //int aa=1;
          //cout<<"dofs loop"<<"\n";
	    for(_IT_GET_FEROW_BY_SIDE_DOFS_FOR_LOOP_(this,"THERMAL_FIELD",MBVERTEX,siit->side_number,dof)) {
            //cout<<"aa "<<aa<<"\t"; aa++;
            
	      vector<DofIdx>::iterator iit = find(RowGlobDofs[0].begin(),RowGlobDofs[0].end(),dof->get_petsc_gloabl_dof_idx());
	      if(iit!=RowGlobDofs[0].end()) {
              *iit = -1;
	      }
	    }
        //cout<<"\n";
	  }
	//cout<<"\n\n";
    }
  }
   //cout<<"\n\n";


   ublas::matrix<FieldData> K;
   for(int rr = 0;rr<(1+6+4+1); rr++) {
       //cout<<"RowGlobDofs[rr].size() "<<RowGlobDofs[rr].size()<<"\n";
	if(RowGlobDofs[rr].size()==0) continue;
	for(int cc = rr;cc<(1+6+4+1);cc++) {
        //cout<<"RowGlobDofs[cc].size() "<<RowGlobDofs[cc].size()<<"\n";
        //getch();
	  if(RowGlobDofs[cc].size()==0) continue;
	  K.resize(RowGlobDofs[rr].size(),RowGlobDofs[cc].size());
	  for(int gg = 0;gg<g_NTET.size()/4;gg++) {
	    ublas::noalias(K) = prod( trans(diffRowNMatrix[rr][gg]),diffRowNMatrix[cc][gg] );
          
        //cout<<"/n/V "<<V<<"\n";
	    K *= V*G_TET_W[gg];   // i think V is the jacobian
	    ierr = MatSetValues(A,     // PetSc function to add values to k matrix
	      RowGlobDofs[rr].size(),&(RowGlobDofs[rr])[0],
	      RowGlobDofs[cc].size(),&(RowGlobDofs[cc])[0],
	      &(K.data())[0],ADD_VALUES); CHKERRQ(ierr);
	    if(cc!=rr) continue;
	    ierr = MatSetValues(A,
	      RowGlobDofs[cc].size(),&(RowGlobDofs[cc])[0],
	      RowGlobDofs[rr].size(),&(RowGlobDofs[rr])[0],
	      &(K.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }      

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|HeatfluxSet,it)) {

	heatflux_cubit_bc_data mydata;
	ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
	double flux = mydata.data.value1;
        
    cout<<"flux: "<<flux;
        
	SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
	SideNumber_multiIndex::nth_index<1>::type::iterator sit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
	SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
	Range meshsets;
	rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
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
	  double normal[3];
	  ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,normal); CHKERRQ(ierr);
	  double area = cblas_dnrm2(3,normal,1);

	  //nodes
	  vector<DofIdx>& RowGlob_nodes = RowGlobDofs[0];
	  vector< ublas::matrix<FieldData> > FaceNMatrix_nodes;
	  ierr = GetGaussRowFaceNMatrix(siit->ent,"THERMAL_FIELD",FaceNMatrix_nodes,MBVERTEX); CHKERRQ(ierr);

	  //edges
	  vector<vector<ublas::matrix<FieldData> > > FaceNMatrix_edges(6);
	  SideNumber_multiIndex::nth_index<1>::type::iterator siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siiit = side_table.get<1>().upper_bound(boost::make_tuple(MBEDGE,6));
	  for(;siiit!=hi_siiit;siiit++) {
	    if(RowGlobDofs[1+siiit->side_number].size()>0) {
	      ierr = GetGaussRowFaceNMatrix(siit->ent,"THERMAL_FIELD",FaceNMatrix_edges[siiit->side_number],MBEDGE,siiit->ent); CHKERRQ(ierr);
	    }
	  }

	  //faces
	  vector<DofIdx>& RowGlob_face = RowGlobDofs[1+6+siit->side_number];
	  vector< ublas::matrix<FieldData> > FaceNMatrix_face;
	  if(RowGlob_face.size()>0) {
	    ierr = GetGaussRowFaceNMatrix(siit->ent,"THERMAL_FIELD",FaceNMatrix_face,MBTRI); CHKERRQ(ierr);
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

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|TemperatureSet,it)) {
	//if(it->get_Cubit_name() != "ZeroPressure") continue;
	Range nodes;
	rval = moab.get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
	for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
	  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
	    ierr = MatSetValue(A,dof->get_petsc_gloabl_dof_idx(),dof->get_petsc_gloabl_dof_idx(),1,INSERT_VALUES); CHKERRQ(ierr);
	}}
      }
      ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

}

#endif // __THERMALFLOWFEMETHOD_HPP_HPP__

