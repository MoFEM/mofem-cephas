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

#ifndef __MOABFEMETHOD_COMPLEXFORLAZY_HPP__
#define __MOABFEMETHOD_COMPLEXFORLAZY_HPP__

using namespace boost::numeric;

namespace ObosleteUsersModules {

struct FEMethod_ComplexForLazy_Data: public FEMethod_UpLevelStudent {

  FieldInterface& mField;
  FEMethod_ComplexForLazy_Data(FieldInterface& _mField,int _verbose = 0):
    FEMethod_UpLevelStudent(_mField.get_moab(),_verbose),mField(_mField) {}

};

/** \brief The user interface for NonLineae FE method (tangent is calculated
 * using complex direvatives)
 * 
 * This class give user some data structures and methods on those that
 * structures which could be useful. Method of calulating complex direvatives
 * is numerically inefficient, however it is simple to use.
 *
 * 
*/
struct FEMethod_ComplexForLazy: public virtual FEMethod_ComplexForLazy_Data {

  enum eRowGlob { i_nodes = 0, i_edge0=1+0, i_edge1=1+1, i_edge2=1+2, i_edge3=1+3, i_edge4=1+4, i_edge5=1+5, 
    i_face0=1+6+0, i_face1=1+6+1, i_face2=1+6+2, i_face3=1+6+3, i_volume=1+6+4, i_last=1+6+4+1 };
  enum analysis { 
    spatail_analysis = 1, material_analysis = 1<<1, mesh_quality_analysis = 1<<2, 
    scaled_themp_direvative_spatial = 1<<3, scaled_themp_direvative_material = 1<<4, analaysis_none = 1<<5 };
  analysis type_of_analysis;
  enum forces { conservative = 1, nonconservative = 2};
  forces type_of_forces;

  double lambda,mu,thermal_expansion;
  void *ptr_matctx;
  bool propeties_from_BLOCKSET_MAT_ELASTICSET;
  PetscErrorCode GetMatParameters(double *_lambda,double *_mu,double *_termal_epansion,void **ptr_matctx);
  ctx_EberleinHolzapfel1 EberleinHolzapfel1_mat_parameters;

  double eps;
  double thermal_load_factor;

  PetscBool flg_alpha2,flg_gamma;
  double alpha2,alpha22,gamma;

  string spatial_field_name;
  string material_field_name;	
  string termal_field_name;
  FEMethod_ComplexForLazy(FieldInterface& _mField,analysis _type,double _lambda,double _mu,double _thermal_expansion,int _verbose);
  ~FEMethod_ComplexForLazy();

  vector<double> g_NTET;
  const double *g_TET_W;
    
  ErrorCode rval;  
  PetscErrorCode ierr;

  vector<double*> diff_edgeNinvJac;
  vector<double*> diff_faceNinvJac;
  double *diff_volumeNinvJac;
  vector<double*> edgeN;
  vector<double*> faceN;
  double *volumeN;

  //Tangent_HH_hierachical
  ublas::matrix<double> KHH,KhH,KvolumeH;
  vector<ublas::matrix<double> > KedgeH_data,KfaceH_data;
  double* KedgeH[6];
  double* KfaceH[4];
  //Tangent_hh_hierachical
  ublas::matrix<double> Khh,KHh,Kvolumeh;
  vector<ublas::matrix<double> > Kedgeh_data,Kfaceh_data;
  double* Kedgeh[6];
  double* Kfaceh[4];
  //Tangent_hh_hierachical_edge
  vector<ublas::matrix<double> > Khedge_data,KHedge_data,Khh_volumeedge_data;
  double *Khedge[6],*KHedge[6];
  double *Khh_volumeedge[6];
  ublas::matrix<ublas::matrix<double> > Khh_edgeedge_data;
  ublas::matrix<ublas::matrix<double> > Khh_faceedge_data;
  double* Khh_edgeedge[6][6];
  double* Khh_faceedge[4][6];
  //Tangent_hh_hierachical_face
  vector<ublas::matrix<double> > Khface_data,KHface_data,Khh_volumeface_data;
  double *Khface[4],*KHface[4];
  double *Khh_volumeface[4];
  ublas::matrix<ublas::matrix<double> > Khh_faceface_data;
  ublas::matrix<ublas::matrix<double> > Khh_edgeface_data;
  double* Khh_faceface[4][4];
  double* Khh_edgeface[6][4];
  //Tangent_hh_hierachical_volume
  ublas::matrix<double> Khvolume,KHvolume;
  vector<ublas::matrix<double> > Khh_edgevolume_data,Khh_facevolume_data;
  double *Khh_edgevolume[6];
  double *Khh_facevolume[4];
  ublas::matrix<double> Khh_volumevolume;
  //Fint
  ublas::vector<double,ublas::bounded_array<double,12> > Fint_h,Fint_h_volume,Fint_H;
  vector<ublas::vector<double> > Fint_h_edge_data;
  vector<ublas::vector<double> > Fint_h_face_data;
  double* Fint_h_edge[6];
  double* Fint_h_face[4];
  ublas::vector<double,ublas::bounded_array<double,12> > iFint_h,iFint_h_volume,iFint_H;
  vector<ublas::vector<double> > iFint_h_edge_data;
  vector<ublas::vector<double> > iFint_h_face_data;
  double* iFint_h_edge[6];
  double* iFint_h_face[4];


  PetscErrorCode GetIndicesRow(
    vector<vector<DofIdx> >& RowGlob,
    string &field_name);
  PetscErrorCode GetIndicesCol(
    vector<vector<DofIdx> >& ColGlob,
    string &field_name);
  PetscErrorCode GetIndices(
    vector<vector<DofIdx> >& RowGlob,
    vector<vector<DofIdx> >& ColGlob,
    string &field_name);

  PetscErrorCode GetData(
    vector<int>& order_edges,vector<int>& order_faces,int& order_volume,
    vector<ublas::vector<double> >& dofs_edge_data,vector<double*>& dofs_edge,
    vector<ublas::vector<double> >& dofs_face_data,vector<double*>& dofs_face,
    ublas::vector<double>& dofs_volume,ublas::vector<double>& dofs_nodes,
    string &field_name);

  //space data
  vector<vector<DofIdx> > RowGlobSpatial;
  vector<vector<DofIdx> > ColGlobSpatial;
  ublas::vector<double> dofs_x,dofs_x_volume;
  vector<ublas::vector<double> > dofs_x_edge_data,dofs_x_face_data;
  vector<double*> dofs_x_edge,dofs_x_face;
  PetscErrorCode GetIndicesSpatial();

  //material data
  vector<vector<DofIdx> > RowGlobMaterial;
  vector<vector<DofIdx> > ColGlobMaterial;
  ublas::vector<double> dofs_X,dofs_X_volume;
  vector<ublas::vector<double> > dofs_X_edge_data,dofs_X_face_data;
  vector<double*> dofs_X_edge,dofs_X_face;
  PetscErrorCode GetIndicesMaterial();

  //temperature
  ublas::vector<double> dofs_temp;

  PetscErrorCode GetDofs_Termal_FromElementData();

  Tag th_quality0,th_quality,th_b;
  double *quality0,*quality,*b;
  PetscErrorCode GetTangent();
  PetscErrorCode GetFint();

  int face_order,face_order_material;
  EntityHandle GetFaceIndicesAndData_face;
  vector<int> FaceEdgeOrder,FaceEdgeOrder_Material;
  vector<DofIdx> FaceIndices;
  ublas::vector<double> FaceData,FaceData_Material;
  vector<double> N_face;
  vector<double> diffN_face;
  vector<DofIdx> FaceNodeIndices,FaceNodeIndices_Material;
  ublas::vector<double> FaceNodeData,FaceNodeData_Material;
  vector<vector<DofIdx> > FaceEdgeIndices_data;
  vector<ublas::vector<double> > FaceEdgeData_data,FaceEdgeData_data_Material;
  double *EdgeData[3],*EdgeData_Material[3];
  vector<vector<double> > N_edge_data;
  vector<vector<double> > diffN_edge_data;
  double* N_edge[3];
  double* diffN_edge[3];

  PetscErrorCode OpComplexForLazyStart();

  vector<int> order_x_edges;
  vector<int> order_x_faces;
  int order_x_volume;
  vector<int> order_X_edges;
  vector<int> order_X_faces;
  int order_X_volume;


};

}

#endif //__MOABFEMETHOD_COMPLEXFORLAZY_HPP__
