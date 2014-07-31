/** \file LoopMethods.hpp
 * \brief MoFEM interface 
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

#ifndef __LOOPMETHODS_HPP__
#define __LOOPMETHODS_HPP__

#include "FieldUnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMBasicMethod = MOFEMuuid( BitIntefaceId(BASIC_METHOD) );
static const MOFEMuuid IDD_MOFEMFEMethod = MOFEMuuid( BitIntefaceId(FE_METHOD) );
static const MOFEMuuid IDD_MOFEMEntMethod = MOFEMuuid( BitIntefaceId(ENT_METHOD) );

/**
 * \brief data structure for snes (nonlinear solver) context
 * \ingroup mofem_loop_methods
 *
 * Struture stores context data which are set in finctions run by PETSc SNES functions.
 *
 */
struct SnesMethod {

  enum SNESContext { CTX_SNESSETFUNCTION, CTX_SNESSETJACOBIAN, CTX_SNESNONE };
  //
  SNESContext snes_ctx;
  SnesMethod(): snes_ctx(CTX_SNESNONE) {};
  //
  PetscErrorCode set_snes_ctx(const SNESContext ctx_);
  //
  SNES snes;
  PetscErrorCode set_snes(SNES _snes);
  Vec snes_x,snes_f;
  Mat snes_A,snes_B;
  virtual ~SnesMethod() {};
};

/**
 * \brief data structure for ts (time stepping) context
 * \ingroup mofem_loop_methods
 *
 * Struture stores context data which are set in finctions run by PETSc Time Stepping functions.
 */
struct TSMethod {
  enum TSContext { CTX_TSSETRHSFUNCTION, CTX_TSSETRHSJACOBIAN, CTX_TSSETIFUNCTION, CTX_TSSETIJACOBIAN, CTX_TSTSMONITORSET, CTX_TSNONE };
  //
  TSContext ts_ctx;
  TSMethod(): ts_ctx(CTX_TSNONE),ts_a(0),ts_t(0) {};
  //
  PetscErrorCode set_ts_ctx(const TSContext ctx_);
  //
  TS ts;
  PetscErrorCode set_ts(TS _ts);
  Vec ts_u,ts_u_t,ts_F;
  Mat ts_A,ts_B;
  //
  PetscInt ts_step;
  PetscReal ts_a,ts_t;
  virtual ~TSMethod() {};
};

/**
 * \brief Data strutucture to exchange data between mofem and User Loop Methods.
 * \ingroup mofem_loop_methods
 *
 * It allows to exchange data between MoFEM and user functoions. It stores informaton about multi-indices.
 */
struct BasicMethod: public FieldUnknownInterface,SnesMethod,TSMethod {

  PetscErrorCode queryInterface (const MOFEMuuid& uuid, FieldUnknownInterface** iface) {
    if(uuid == IDD_MOFEMBasicMethod) {
      *iface = dynamic_cast<BasicMethod*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"unknown inteface");
  }

  BasicMethod();    
  //
  virtual PetscErrorCode preProcess() = 0;
  virtual PetscErrorCode operator()() = 0;
  virtual PetscErrorCode postProcess() = 0;
  //
  const RefMoFEMEntity_multiIndex *refinedEntitiesPtr;
  const RefMoFEMElement_multiIndex *refinedFiniteElementsPtr;
  const MoFEMProblem *problemPtr;
  const MoFEMField_multiIndex *fieldsPtr;
  const MoFEMEntity_multiIndex *entitiesPtr;
  const DofMoFEMEntity_multiIndex *dofsPtr;
  const MoFEMFiniteElement_multiIndex *finiteElementsPtr;
  const EntMoFEMFiniteElement_multiIndex *finiteElementsEntitiesPtr;
  const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex *adjacenciesPtr;
  virtual ~BasicMethod() {};
};

/**
  * \brief structure for User Loop Methods on finite elements
  * \ingroup mofem_loop_methods
  *
  * It can be used to calulate stiffnes matrices, residuals, load vectors etc.
  */  
struct FEMethod: public BasicMethod {

  PetscErrorCode queryInterface (const MOFEMuuid& uuid, FieldUnknownInterface** iface) {
    PetscFunctionBegin;
    if(uuid == IDD_MOFEMFEMethod) {
      *iface = dynamic_cast<FEMethod*>(this);
      PetscFunctionReturn(0);
    }
    PetscErrorCode ierr; 
    ierr = queryInterface(uuid,iface); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  FEMethod();

  /** \brief function is run at the beginig of looop
   *
   * It is used to zeroing matrices and vectors, calculation of shape
   * functions on reference element, preporocessing boundary conditions, etc.
   */
  PetscErrorCode preProcess();

  /** \brief function is run for every finite element 
   *
   * It is used to calulate element local matrices and assembly. It can be
   * used for post-processing.
   */
  PetscErrorCode operator()();

  /** \brief function is run at the end of looop
   *
   * It is used to assembly matrices and vectors, calulating global variables,
   * f.e. total internal energy, ect.
   * 
   * Iterating over dofs:
   * Example1 iterating over dofs in row by name of the field
   * for(_IT_GET_FEROW_BY_NAME_DOFS_FOR_LOOP_(this,"DISPLACEMENT",it)) { ... } 
   * 
   * 
   */
  PetscErrorCode postProcess();
  string feName;
  const NumeredMoFEMFiniteElement *fePtr;
  const FEDofMoFEMEntity_multiIndex *dataPtr;
  const FENumeredDofMoFEMEntity_multiIndex *rowPtr;
  const FENumeredDofMoFEMEntity_multiIndex *colPtr;

  /** \brief loop over all dofs which are on a particular FE row 
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEROW_DOFS_FOR_LOOP_(FE,IT) \
  FENumeredDofMoFEMEntity_multiIndex::iterator IT = FE->rowPtr->begin(); IT != FE->rowPtr->end();IT++ 

  /** \brief loop over all dofs which are on a particular FE column 
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FECOL_DOFS_FOR_LOOP_(FE,IT) \
  FENumeredDofMoFEMEntity_multiIndex::iterator IT = FE->colPtr->begin(); IT != FE->colPtr->end();IT++ 


  /** \brief loop over all dofs which are on a particular FE data 
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEDATA_DOFS_FOR_LOOP_(FE,IT) \
  FEDofMoFEMEntity_multiIndex::iterator IT = FE->dataPtr->begin(); IT != FE->dataPtr->end();IT++ 

  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,  
    const string &field_name,const EntityType type,const int side_number) const {
    return index.lower_bound(boost::make_tuple(field_name,type,side_number));
  } 
  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,  
    const string &field_name,const EntityType type,const int side_number) const {
    return index.upper_bound(boost::make_tuple(field_name,type,side_number));
  } 

    /** \brief loop over all dofs which are on a particular FE row, field, entity type and canonical side number
     * \ingroup mofem_loop_methods
     *
     * \param FE finite elements
     * \param Name field name
     * \param Type moab entity type (MBVERTEX, MBEDGE etc)
     * \param Side side canonical number
     * \param IT the interator in use
     */
  #define _IT_GET_FEROW_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->rowPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->rowPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

  /** \brief loop over all dofs which are on a particular FE column, field, entity type and canonical side number
    * \ingroup mofem_loop_methods
    */ 
  #define _IT_GET_FECOL_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->colPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->colPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

  /** \brief loop over all dofs which are on a particular FE data, field, entity type and canonical side number
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEDATA_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
  FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
    IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->dataPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
    IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->dataPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name,const EntityType type) const {
    return index.lower_bound(boost::make_tuple(field_name,type));
  } 
  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name,const EntityType type) const {
    return index.upper_bound(boost::make_tuple(field_name,type));
  } 

  /** \brief loop over all dofs which are on a particular FE row, field and entity type
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); IT++

  /** \brief loop over all dofs which are on a particular FE column, field and entity type
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FECOL_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); IT++

  /** \brief loop over all dofs which are on a particular FE data, field and entity type
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator \
    IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); \
    IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); IT++

  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name) const {
    return index.lower_bound(field_name);
  } 
  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name) const {
    return index.upper_bound(field_name);
  } 

  /** \brief loop over all dofs which are on a particular FE row and field
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEROW_BY_NAME_DOFS_FOR_LOOP_(FE,NAME,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->rowPtr->get<FieldName_mi_tag>(),NAME); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->rowPtr->get<FieldName_mi_tag>(),NAME); IT++

  /** \brief loop over all dofs which are on a particular FE column and field
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FECOL_BY_NAME_DOFS_FOR_LOOP_(FE,NAME,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->colPtr->get<FieldName_mi_tag>(),NAME); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->colPtr->get<FieldName_mi_tag>(),NAME); IT++

  /** \brief loop over all dofs which are on a particular FE data and field
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEDATA_BY_NAME_DOFS_FOR_LOOP_(FE,NAME,IT) \
  FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
    IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->dataPtr->get<FieldName_mi_tag>(),NAME); \
    IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->dataPtr->get<FieldName_mi_tag>(),NAME); IT++

  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const EntityHandle ent) const {
    return index.lower_bound(ent);
  } 
  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const EntityHandle ent) const {
    return index.upper_bound(ent);
  } 

  /** \brief loop over all dofs which are on a particular FE row and given element entity (handle from moab)
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEROW_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->rowPtr->get<MoABEnt_mi_tag>(),ENT); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->rowPtr->get<MoABEnt_mi_tag>(),ENT); IT++

  /** \brief loop over all dofs which are on a particular FE column and given element entity (handle from moab)
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FECOL_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->colPtr->get<MoABEnt_mi_tag>(),ENT); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->colPtr->get<MoABEnt_mi_tag>(),ENT); IT++

  /** \brief loop over all dofs which are on a particular FE data and given element entity (handle from moab)
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEDATA_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
  FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
    IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->dataPtr->get<MoABEnt_mi_tag>(),ENT); \
    IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->dataPtr->get<MoABEnt_mi_tag>(),ENT); IT++

  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name,const EntityHandle ent) const {
    return index.lower_bound(boost::make_tuple(field_name,ent));
  } 
  template<class MULTIINDEX>
  typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name,const EntityHandle ent) const {
    return index.upper_bound(boost::make_tuple(field_name,ent));
  } 

  /** \brief loop over all dofs which are on a particular FE row, field and given element entity (handle from moab)
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); IT++

  /** \brief loop over all dofs which are on a particular FE column, field and given element entity (handle from moab)
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FECOL_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator \
    IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); \
    IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); IT++

  /** \brief loop over all dofs which are on a particular FE data, field and given element entity (handle from moab)
    * \ingroup mofem_loop_methods
    */
  #define _IT_GET_FEDATA_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator \
    IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); \
    IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); IT++

};

/**
 * \brief Data strutucture to exchange data between mofem and User Loop Methods on Entirties.
 * \ingroup mofem_loop_methods
 *
 * It allows to exchange data between MoFEM and user functoions. It stores informaton about multi-indices.
 */
struct EntMethod: public BasicMethod {

  PetscErrorCode queryInterface (const MOFEMuuid& uuid, FieldUnknownInterface** iface) {
    PetscFunctionBegin;
    if(uuid == IDD_MOFEMEntMethod) {
      *iface = dynamic_cast<EntMethod*>(this);
      PetscFunctionReturn(0);
    }
    PetscErrorCode ierr;
    ierr = queryInterface(uuid,iface); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  EntMethod();
  
  PetscErrorCode preProcess();
  PetscErrorCode operator()();
  PetscErrorCode postProcess();
 
  const DofMoFEMEntity *dofPtr;
  const NumeredDofMoFEMEntity *dofNumeredPtr;
};

}

#endif // __LOOPMETHODS_HPP__

/***************************************************************************//**
 * \defgroup mofem_loop_methods Methods for Loops
 * \ingroup mofem
 ******************************************************************************/

