/* \file Dirichlet.hpp
 * \brief Implementation of Dirichlet boundary conditions
 *
 *
 * Structures and method in this file erase rows and column, set value on
 * matrix diagonal and on the right hand side vector to enforce boundary
 * condition.
 *
 * Current implementation is suboptimal, classes name too long. Need to
 * rethinking and improved, more elegant and more efficient implementation.
 *
 */

/* Notes:

 DirichletBCFromBlockSetFEMethodPreAndPostProc implemented by Zahur Ullah
 (Zahur.Ullah@glasgow.ac.uk)

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

#ifndef __DIRICHLETBC_HPP__
#define __DIRICHLETBC_HPP__

using namespace boost::numeric;

/** \brief Set Dirichlet boundary conditions on displacements
  * \ingroup Dirichlet_bc
  */
struct DisplacementBCFEMethodPreAndPostProc: public MoFEM::FEMethod {

  FieldInterface& mField;
  const std::string fieldName;			///< field name to set Dirichlet BC
  double dIag;					      ///< diagonal value set on zeroed column and rows

  DisplacementBCFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name,
    Mat Aij,Vec X,Vec F
  );
  DisplacementBCFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name
  );

  PetscErrorCode ierr;
  ErrorCode rval;

  std::map<DofIdx,FieldData> mapZeroRows;
  std::vector<int> dofsIndices;
  std::vector<double> dofsValues;
  std::vector<double> dofsXValues;
  virtual PetscErrorCode iNitalize();

  PetscErrorCode preProcess();
  PetscErrorCode postProcess();

  boost::ptr_vector<MethodForForceScaling> methodsOp;

};

/** \brief Set Dirichlet boundary conditions on spatial displacements
  * \ingroup Dirichlet_bc
  */
struct SpatialPositionsBCFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {

  SpatialPositionsBCFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name,Mat aij,Vec x,Vec f):
    DisplacementBCFEMethodPreAndPostProc(m_field,field_name,aij,x,f) {}

  SpatialPositionsBCFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name):
    DisplacementBCFEMethodPreAndPostProc(m_field,field_name) {}

  std::vector<std::string> fixFields;

  ublas::vector<double> cOords;
  PetscErrorCode iNitalize();

};

struct TemperatureBCFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {

  TemperatureBCFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name,Mat aij,Vec x,Vec f):
    DisplacementBCFEMethodPreAndPostProc(m_field,field_name,aij,x,f) {}

  TemperatureBCFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name):
    DisplacementBCFEMethodPreAndPostProc(m_field,field_name) {}

  PetscErrorCode iNitalize();

};

/** \brief Fix dofs on entities
  * \ingroup Dirichlet_bc
  */
struct FixBcAtEntities: public DisplacementBCFEMethodPreAndPostProc {

  Range &eNts;
  std::vector<std::string> fieldNames;
  FixBcAtEntities(
    FieldInterface& m_field,
    const std::string &field_name,
    Mat aij,Vec x,Vec f,
    Range &ents
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field_name,aij,x,f),eNts(ents) {
    fieldNames.push_back(fieldName);
  }

  FixBcAtEntities(
    FieldInterface& m_field,const std::string &field_name,Range &ents
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field_name),eNts(ents) {
    fieldNames.push_back(fieldName);
  }

  PetscErrorCode iNitalize();
  PetscErrorCode preProcess();
  PetscErrorCode postProcess();

};


/** \brief Blockset boundary conditions
  * \ingroup Dirichlet_bc
  *
  * Implementation of generalized Dirichlet Boundary Conditions from CUBIT Blockset
  * (or not using CUBIT building boundary conditions, e.g. Temperature or Displacements etc).
  * It can work for any Problem rank (1,2,3)
  *
  * FIXME Show example form journal file how to add block set
  *
**/
struct DirichletBCFromBlockSetFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {

  const std::string blocksetName;
  DirichletBCFromBlockSetFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name,const std::string &blockset_name,Mat aij,Vec x,Vec f
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field_name,aij,x,f),
  blocksetName(blockset_name) {
  }

  DirichletBCFromBlockSetFEMethodPreAndPostProc(
    FieldInterface& m_field,const std::string &field_name,const std::string &blockset_name
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field_name),
  blocksetName(blockset_name) {
  }

  PetscErrorCode iNitalize();

};

/**
 * \brief Add boundary conditions form block set having 6 attributes
 *
 * FIXME: Detailed description how to use it.
 * FIXME Show example form journal file how to add block set
 */
struct DirichletBCFromBlockSetFEMethodPreAndPostProcWithFlags: public DisplacementBCFEMethodPreAndPostProc {

  const std::string blocksetName;
  DirichletBCFromBlockSetFEMethodPreAndPostProcWithFlags(
    FieldInterface& m_field,const std::string &field_name,const std::string &blockset_name,Mat aij,Vec x,Vec f
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field_name,aij,x,f),
  blocksetName(blockset_name) {
  }

  DirichletBCFromBlockSetFEMethodPreAndPostProcWithFlags(
    FieldInterface& m_field,const std::string &field_name,const std::string &blockset_name
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field_name),
  blocksetName(blockset_name) {
  }

  PetscErrorCode iNitalize();

};

#endif //__DIRICHLETBC_HPP__

/***************************************************************************//**
 * \defgroup Dirichlet_bc Dirichlet boundary conditions
 * \ingroup user_modules
 ******************************************************************************/
