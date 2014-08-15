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

#ifndef __MOABFEMETHOD_UPLEVELSTUDENT_HPP__
#define __MOABFEMETHOD_UPLEVELSTUDENT_HPP__

namespace ObosleteUsersModules {

struct FEMethod_UpLevelStudent_ExceptionNegatvieTetVolume: public MofemException {
  FEMethod_UpLevelStudent_ExceptionNegatvieTetVolume(): 
    MofemException(MOFEM_DATA_INSONSISTENCY,"Negative volume") {}
};

/** 
 * \brief The student user interface for FE method
 * 
 * This class give user some data structures and methods on those that
 * structures which could be useful.
*/
struct FEMethod_UpLevelStudent: public FEMethod_LowLevelStudent {

  //TET
  double V;
  Tag th_volume;
  vector<ublas::vector<double,ublas::bounded_array<double, 3> > > coords_at_Gauss_nodes; ///< vector of coordinates at Gauss points

  //PRISM
  double area3,area4;
  double coords_face3[9];
  double coords_face4[9];
  double normal3[3];
  double normal4[3];
  EntityHandle conn_face3[3];
  EntityHandle conn_face4[3];

  FEMethod_UpLevelStudent(Interface& _moab,int verbose = 0);
  ~FEMethod_UpLevelStudent();

  /** \brief Initate data structures for running FE methods
   *
   * It has to be run at the begining of the function when tetrahedral element
   * is evaluated.
   *
   * \param _gNTET_ vector of shape of tetrahedral functions evaluated at Gauss
   * points
   */
  PetscErrorCode OpStudentStart_TET(vector<double>& _gNTET_);

  /** \brief Initate data structures for running FE methods
   *
   * It has to be run at the begining of the function when Interface PRISM
   * element is evaluated. The matrix element shape functions are calulated in
   * the form that degrees of shepe functon on one face (face4), have negative
   * value comparing to oposite face. For example 
   * gap = ShapeN_FunForPrism * nodal_displacements.
   *
   *
   * \param _gNTRI_ vector of shape of tetrahedral functions evaluated at Gauss
   * points
   */
  PetscErrorCode OpStudentStart_PRISM(vector<double>& _gNTRI_);

  /**
   * \brief Finalise data structures for running FE methods
   *
   * It has to be run at the end of the function
   */
  PetscErrorCode OpStudentEnd();

  /**
   * \brief Copy gloabl indices for dofs adjacent to nodes
   *
   * \param field_name name of the approx. field
   * \param vector on return stores global indices
   */
  PetscErrorCode GetRowGlobalIndices(const string &field_name,vector<DofIdx> &RowGlobDofs);
  PetscErrorCode GetRowLocalIndices(const string &field_name,vector<DofIdx> &RowLocalDofs);

  /**
   * \brief Copy gloabl indices for dofs adjacent to entities
   *
   * \param field_name name of the approx. field
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores global indices
   * \param side_number  do need to be specified for MBTET or MBPRISM
   */
  PetscErrorCode GetRowGlobalIndices(const string &field_name,EntityType type,vector<DofIdx> &RowGlobDofs,int side_number = -1);
  PetscErrorCode GetRowLocalIndices(const string &field_name,EntityType type,vector<DofIdx> &RowLocalDofs,int side_number = -1);

  /**
   * \brief Copy gloabl indices for dofs adjacent to nodes
   *
   * \param field_name name of the approx. field
   * \param vector on return stores global indices
   */
  PetscErrorCode GetColGlobalIndices(const string &field_name,vector<DofIdx> &ColGlobDofs);
  PetscErrorCode GetColLocalIndices(const string &field_name,vector<DofIdx> &ColLocalDofs);


  /**
   * \brief Copy gloabl indices for dofs adjacent to entities
   *
   * \param field_name name of the approx. field
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores global indices
   * \param side_number  do need to be specified for MBTET or MBPRISM   */
  PetscErrorCode GetColGlobalIndices(const string &field_name,EntityType type,vector<DofIdx> &ColGlobDofs,int side_number = -1);
  PetscErrorCode GetColLocalIndices(const string &field_name,EntityType type,vector<DofIdx> &ColGlobDofs,int side_number = -1);

  /**
   * \brief Copy dofs values for dofs adjacent to nodes
   *
   * \param field_name name of the approx. field
   * \param vector on return stores dofs values 
   */
  PetscErrorCode GetDataVector(const string &field_name,ublas::vector<FieldData> &Data);

  /**
   * \brief Copy dofs values for dofs adjacent to entities
   *
   * \param field_name name of the approx. field
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores dofs values 
   * \param side_number  do need to be specified for MBTET or MBPRISM
   */
  PetscErrorCode GetDataVector(const string &field_name,EntityType type,ublas::vector<FieldData> &Data,int side_number = -1);


  /** \brief Field data at Gauss points for L2 and H1 space
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param vector of vectors on return values of field at gauss points 
   &
   */
  PetscErrorCode GetGaussDataVector(const string &field_name,vector<ublas::vector<FieldData> > &Data);

  /** \brief Field data direvatives at Gauss points for L2 and H1 space
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param vector on return direvatives of the field at gauss points 
   */
  PetscErrorCode GetGaussDiffDataVector(const string &field_name,vector< ublas::matrix<FieldData> > &Data);

  /** \brief Copy shape functions for nodes
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param vector on return stores
   * values of field at gauss points 
   */
  PetscErrorCode GetGaussRowNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &NMatrix);

  /** \brief Copy shape functions for entities and given side number
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores
   * \param side_number  do need to be specified for MBTET or MBPRISM   
   */
  PetscErrorCode GetGaussRowNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &NMatrix,int side_number = -1);

  /** \brief Copy shape functions for nodes
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param vector on return stores
   * values of field at gauss points 
   */
  PetscErrorCode GetGaussColNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &NMatrix);

  /** \brief Copy shape functions for entities and given side number
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores
   * \param side_number  do need to be specified for MBTET or MBPRISM   
   */
  PetscErrorCode GetGaussColNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &NMatrix,int side_number = -1);

  /** \brief Copy derivative shape functions for nodes
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param vector on return stores
   * values of field at gauss points 
   */
  PetscErrorCode GetGaussRowDiffNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &diffNMatrix);

  /** \brief Copy derivative shape functions for entities and given side number
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores
   * \param side_number  do need to be specified for MBTET or MBPRISM   
   */
  PetscErrorCode GetGaussRowDiffNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &diffNMatrix,int side_number = -1);

  /** \brief Copy derivative shape functions for nodes
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param vector on return stores
   * values of field at gauss points 
   */
  PetscErrorCode GetGaussColDiffNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &diffNMatrix);

  /** \brief Copy derivative shape functions for entities and given side number
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores
   * \param side_number  do need to be specified for MBTET or MBPRISM   
   */
  PetscErrorCode GetGaussColDiffNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &diffNMatrix,int side_number = -1);

  /** 
   * \brief Make B matrix for 3D field
   *
   * For more detail look page 30 CHAPTER 6. DISPLACEMENT METHODS, FEAP Version 7.3 Theory Manual Robert L. Taylor
   *
   */ 
  PetscErrorCode MakeBMatrix3D(const string &field_name,
    vector<ublas::matrix<FieldData> > &diffNMatrix,vector<ublas::matrix<FieldData> > &BMatrix);

   /** \brief Get shape functions for integration on the face
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions on the face, face which is triangle g_dim =
   * g_NTRI.size()/3
   *
   * \param EntityHandle moab entity handle to face (MBTRI) on tetrahedral (MBTET)
   * \param field_name name of the approx. field 
   * \param vector on return stores
   * \param entity type on face (MBVERTEX, MBEDGE or MBTRI)
   * \param entity hande on face (it need to be given when entity type is MBEDEGE)
   * values of field at gauss points 
   */
  PetscErrorCode GetGaussRowFaceNMatrix(
    EntityHandle ent,const string &field_name,vector< ublas::matrix<FieldData> > &diffNMatrix,
    EntityType type = MBMAXTYPE,EntityHandle edge_handle = no_handle);

  /**
    * \brief hierarhical gemetry approximation, jacobian and determinant of jacobian
    */
  PetscErrorCode GetHierarchicalGeometryApproximation(vector< ublas::matrix<FieldData> > &invH,vector< FieldData > &detH);

  /**
    * \brief hierarhical gemetry approximation, diff shape functions
    */
  PetscErrorCode GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(int rank,vector< ublas::matrix<FieldData> > &invH,vector< ublas::matrix<FieldData> > &diffNMatrix);

};

}

#endif // __MOABFEMETHOD_UPLEVELSTUDENT_HPP__
