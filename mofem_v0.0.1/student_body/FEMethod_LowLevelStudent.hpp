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



#ifndef __MOABFEMETHOD_LOWLEVELSTUDENT_HPP__
#define __MOABFEMETHOD_LOWLEVELSTUDENT_HPP__

#include "FieldInterface.hpp"
#include "CoreDataStructures.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** 
 * \brief The core user interface for FE method
 * 
 * This class give user some data structures and methods on those that
 * structures which could be useful.
*/
struct FEMethod_LowLevelStudent: public FieldInterface::FEMethod {
  Interface& moab;

  //
  FEMethod_LowLevelStudent *ParentMethod;
  int verbose;

  double NTET[4],diffNTET[12],diffNTETinvJac[12];
  double NTRI[3],diffNTRI[6];

  const EntityHandle* conn;
  ublas::vector<double,ublas::bounded_array<double,18> > coords;

  FEMethod_LowLevelStudent(Interface& _moab,int verbose = 0);
  ~FEMethod_LowLevelStudent();

  PetscErrorCode GlobIndices();
  PetscErrorCode LocalIndices();
  PetscErrorCode DataOp();
  PetscErrorCode ParentData(const string &fe_name);

  typedef FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type dofs_by_Composite;
  typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type numered_dofs_by_Composite;

  template <typename T> static FieldData UnaryFunction_FieldData(const T *it) { return it->get_FieldData(); }
  template <typename T> static ApproximationRank UnaryFunction_ApproxRank(const T *it) { return it->get_dof_rank(); }
  template <typename T> static ApproximationOrder UnaryFunction_ApproxOrder(const T *it) { return it->get_dof_order(); }
  template <typename T> static DofIdx UnaryFunction_PetscGlobalIdx(const T *it) { return it->get_petsc_gloabl_dof_idx(); }
  template <typename T> static DofIdx UnaryFunction_PetscLocalIdx(const T *it) { return it->get_petsc_local_dof_idx(); }
  template <typename T> static DofIdx UnaryFunction_EntDofIdx(const T *it) { return it->get_EntDofIdx(); }

  typedef map<const MoFEMField*,vector<DofIdx> > Indices_Type;
  typedef map<const MoFEMEntity*,vector<DofIdx> > Indices_EntType;
  //Glocabl Indices
  Indices_Type row_nodesGlobIndices;
  Indices_EntType row_edgesGlobIndices;
  Indices_EntType row_facesGlobIndices;
  Indices_EntType row_elemGlobIndices;
  Indices_Type col_nodesGlobIndices;
  Indices_EntType col_edgesGlobIndices;
  Indices_EntType col_facesGlobIndices;
  Indices_EntType col_elemGlobIndices;
  //Local Indices
  Indices_Type row_nodesLocalIndices;
  Indices_EntType row_edgesLocalIndices;
  Indices_EntType row_facesLocalIndices;
  Indices_EntType row_elemLocalIndices;
  Indices_Type col_nodesLocalIndices;
  Indices_EntType col_edgesLocalIndices;
  Indices_EntType col_facesLocalIndices;
  Indices_EntType col_elemLocalIndices;

  typedef map<const MoFEMField*,ublas::vector<FieldData> > Data_Type;
  typedef map<const MoFEMEntity*,ublas::vector<FieldData> > Data_EntType;
  Data_Type data_nodes;
  Data_EntType data_edges;
  Data_EntType data_faces;
  Data_EntType data_elem;

  typedef map<string,vector< ublas::vector<FieldData> > > H1L2_Data_at_Gauss_pt;
  typedef map<string,vector< ublas::matrix<FieldData> > > H1_DiffData_at_Gauss_pt;
  typedef map<string,vector< ublas::matrix<FieldData> > > HcurlHdiv_Data_at_Gauss_pt;
  H1L2_Data_at_Gauss_pt h1l2_data_at_gauss_pt;
  H1_DiffData_at_Gauss_pt h1_diff_data_at_gauss_pt;
  HcurlHdiv_Data_at_Gauss_pt hcurl_hdiv_data_at_gauss_pt;
  PetscErrorCode Data_at_GaussPoints();
  PetscErrorCode DiffData_at_GaussPoints();

  typedef map<const MoFEMField*,vector< ublas::matrix<FieldData> > > N_Matrix_Type;
  typedef map<const MoFEMEntity*,vector< ublas::matrix<FieldData> > > N_Matrix_EntType;
  N_Matrix_Type row_N_Matrix_nodes;
  N_Matrix_EntType row_N_Matrix_edges;
  N_Matrix_EntType row_N_Matrix_faces;
  N_Matrix_EntType row_N_Matrix_elem;
  N_Matrix_Type col_N_Matrix_nodes;
  N_Matrix_EntType col_N_Matrix_edges;
  N_Matrix_EntType col_N_Matrix_faces;
  N_Matrix_EntType col_N_Matrix_elem;
  PetscErrorCode GetNMatrix_at_GaussPoint(
    Indices_Type& nodesGlobIndices, Indices_EntType& edgesGlobIndices,
    Indices_EntType& facesGlobIndices, Indices_EntType& volumeGlobIndices,
    N_Matrix_Type& N_Matrix_nodes,
    N_Matrix_EntType& N_Matrix_edges,
    N_Matrix_EntType& N_Matrix_faces,
    N_Matrix_EntType& N_Matrix_elem);
  PetscErrorCode GetRowNMatrix_at_GaussPoint();
  PetscErrorCode GetColNMatrix_at_GaussPoint();
  N_Matrix_Type row_diffN_Matrix_nodes;
  N_Matrix_EntType row_diffN_Matrix_edges;
  N_Matrix_EntType row_diffN_Matrix_faces;
  N_Matrix_EntType row_diffN_Matrix_elem;
  N_Matrix_Type col_diffN_Matrix_nodes;
  N_Matrix_EntType col_diffN_Matrix_edges;
  N_Matrix_EntType col_diffN_Matrix_faces;
  N_Matrix_EntType col_diffN_Matrix_elem;

  /**
   * \brief Calulate element GRADIENT matrices
   *
   * [ dN1/dx	....	.... 	]
   * [ dN1/dy	.... 	....	]
   * [ ....	....	....	]
   * [ 0	dN1/dx	....	]
   * [ 0	dN1/dy  ....	]
   * [ ....	....	....	]
   */
  PetscErrorCode GetDiffNMatrix_at_GaussPoint(
    Indices_Type& nodesGlobIndices, Indices_EntType& edgesGlobIndices,
    Indices_EntType& facesGlobIndices, Indices_EntType& volumeGlobIndices,
    N_Matrix_Type& diffNMatrix_nodes,
    N_Matrix_EntType& diffNMatrix_edges,
    N_Matrix_EntType& diffNMatrix_faces,
    N_Matrix_EntType& diffNMatrix_elem);
  PetscErrorCode GetRowDiffNMatrix_at_GaussPoint();
  PetscErrorCode GetColDiffNMatrix_at_GaussPoint();

  const EntMoFEMFiniteElement *fe_ent_ptr;

  /**
   * calulate element shape functions
   *
   * \param _gNTET_ vector of shape functions at Gauss pts.
   */
  PetscErrorCode ShapeFunctions_TET(vector<double>& _gNTET_);

  /**
   * calculate element faces shape functions
   *
   * \param _gNTRI_ vector of shape functions at Gauss pts.
   *
   */
  vector<double> gNTRIonPRISM;
  PetscErrorCode ShapeFunctions_PRISM(vector<double>& _gNTRI_);

  //H1,L2
  vector< vector<double> > H1edgeN,diffH1edgeN;
  vector< vector<double> > H1faceN,diffH1faceN;
  vector<double> H1elemN,diffH1elemN;
  vector<double> L2elemN,diffL2elemN;
  vector< vector<double> > diffH1edgeNinvJac;
  vector< vector<double> > diffH1faceNinvJac;
  vector<double> diffH1elemNinvJac;
  vector<double> diffL2elemNinvJac;
  //Hdiv Face
  vector< vector< ublas::vector<int> > > Hdiv_egde_faceOrder;
  vector< vector< ublas::vector<double> > > Hdiv_egde_faceN;
  vector< ublas::vector<int> > Hdiv_face_bubbleOrder;
  vector< ublas::vector<double> > Hdiv_face_bubbleN;
  //Hdiv Volume
  vector< ublas::vector<int> > Hdiv_edge_volumeOrder;
  vector< ublas::vector<double> > Hdiv_edge_volumeN;
  vector< ublas::vector<int> > Hdiv_face_volumeOrder;
  vector< ublas::vector<double> > Hdiv_face_volumeN;
  ublas::vector< int > Hdiv_volumeOrder;
  ublas::vector< double > Hdiv_volumeN;
  vector< vector<double> > Hdiv_faceN_byOrder;
  vector<double> Hdiv_volumeN_byOrder;

  /// get element shape functions
  PetscErrorCode get_ShapeFunction(
    vector<const double*> *base_functions_by_gauss_pt,
    vector<const double*> *diff_base_functions_by_gauss_pt,
    const MoFEMField* field_ptr,EntityType type,int side_number = -1);
  bool isH1,isHdiv,isHcurl,isL2;
  vector<int> maxOrderEdgeH1,maxOrderEdgeHcurl;
  vector<int> maxOrderFaceH1,maxOrderFaceHdiv,maxOrderFaceHcurl;
  int maxOrderElemH1,maxOrderElemHdiv,maxOrderElemHcurl,maxOrderElemL2;

  /**
   * calculate element faces shape functions
   *
   * \param _gNTRI_ vector of shape functions at Gauss pts.
   *
   */
  PetscErrorCode ShapeFunctions_TRI(EntityHandle ent,vector<double> &_gNTRI_);
  vector<vector<double> > gNTRIonElem;
  map<EntityHandle,map<EntityHandle,vector<double> > > H1edgeN_TRI,diffH1edgeN_TRI;
  map<EntityHandle,vector<double> > H1faceN_TRI,diffH1faceN_TRI;
  PetscErrorCode GetNMatrix_at_FaceGaussPoint(
    EntityHandle ent,
    const string& field_name,
    Indices_Type& nodesGlobIndices, 
    Indices_EntType& edgesGlobIndices,
    Indices_EntType& facesGlobIndices,
    N_Matrix_Type& N_Matrix_nodes,
    N_Matrix_EntType& N_Matrix_edges,
    N_Matrix_EntType& N_Matrix_faces,
    EntityType type = MBMAXTYPE,
    EntityHandle edge_handle = no_handle);

  PetscErrorCode InitDataStructures();

  PetscErrorCode ierr;
  ErrorCode rval;

  int get_dim_gNTET() const { return gNTET.size()/4; };
  const vector<double>& get_gNTET() const { return gNTET; };
  int get_dim_gNTRI() const { return gNTRI.size()/3; };
  const vector<double>& get_gNTRI() const { return gNTRI; };

  protected:
  vector<double> gNTET;
  vector<double> gNTRI;  

};

}

#endif // __MOABFEMETHOD_LOWLEVELSTUDENT_HPP__
