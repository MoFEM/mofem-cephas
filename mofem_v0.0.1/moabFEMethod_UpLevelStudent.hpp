/** \file moabField_Core.hpp
 * \brief Core moabField::FEMethod class for user interface
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __MOABFEMETHOD_STUDENT_HPP__
#define __MOABFEMETHOD_STUDENT_HPP__

#include "moabField.hpp"
#include "moabFEMethod_LowLevelStudent.hpp"
#include "Core_dataStructures.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** 
 * \brief The student user interface for FE method
 * 
 * This class give user some data structures and methods on those that
 * structures which could be useful.
*/
struct FEMethod_UpLevelStudent: public FEMethod_LowLevelStudent {

  double V;
  Tag th_volume;
  vector<ublas::vector<double,ublas::bounded_array<double, 3> > > coords_at_Gauss_nodes;

  FEMethod_UpLevelStudent(Interface& _moab,int verbose = 0);
  ~FEMethod_UpLevelStudent();

  /**
   * \brief Initate data structures for running FE methods
   *
   * It has to be run at the begining of the function
   *
   * \param _gNTET_ vector of shape of tetrahedral functions evaluated at Gauss points
   */
  PetscErrorCode OpStudentStart_TET(vector<double>& _gNTET_);

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
  PetscErrorCode GetRowIndices(const string &field_name,vector<DofIdx> &RowGlobDofs);

  /**
   * \brief Copy gloabl indices for dofs adjacent to entities
   *
   * \param field_name name of the approx. field
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores global indices
   * \param side_number  do need to be specified for MBTET or MBPRISM
   */
  PetscErrorCode GetRowIndices(const string &field_name,EntityType type,vector<DofIdx> &RowGlobDofs,int side_number = -1);

  /**
   * \brief Copy gloabl indices for dofs adjacent to nodes
   *
   * \param field_name name of the approx. field
   * \param vector on return stores global indices
   */
  PetscErrorCode GetColIndices(const string &field_name,vector<DofIdx> &ColGlobDofs);

  /**
   * \brief Copy gloabl indices for dofs adjacent to entities
   *
   * \param field_name name of the approx. field
   * \param type type of the entity (MBEDGE, MBTRI, MBTET or MBPRISM)
   * \param vector on return stores global indices
   * \param side_number  do need to be specified for MBTET or MBPRISM   */

  PetscErrorCode GetColIndices(const string &field_name,EntityType type,vector<DofIdx> &ColGlobDofs,int side_number = -1);

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


  /** \brief Copy Gauss pt. values for dofs
   *
   * Note that size of the global vector is determined by number of gauss
   * points in shape functions, f.e. in case of tetrahedral (MBTET) g_dim =
   * g_NTET.size()/4
   *
   * \param field_name name of the approx. field 
   * \param vector on return values of field at gauss points 
   */
  PetscErrorCode GetGaussDataVector(const string &field_name,vector<ublas::vector<FieldData> > &Data);

  /** \brief Copy Gauss pt. values direvatives
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

};

}

#endif //__MOABFEMETHOD_STUDENT_HPP__
