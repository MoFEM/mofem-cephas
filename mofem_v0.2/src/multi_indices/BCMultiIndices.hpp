/** \file BCMultiIndices.hpp
 * \brief Multi-index containers, data boundary data structures and other low-level functions
 *
 */ 

/*
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
 *
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

#ifndef __BCMULTIINDICES_HPP__
#define __BCMULTIINDICES_HPP__

namespace MoFEM {


/** 
 * \brief this struct keeps basic methods for moab meshset about material and boundary conditions
 */
struct CubitMeshSets {
  EntityHandle meshset;
  CubitBC_BitSet CubitBCType;
  vector<Tag> tag_handles;
  int *msId;
  char* tag_bc_data;
  int tag_bc_size;
  unsigned int *tag_block_header_data;
  double* tag_block_attributes;
  int tag_block_attributes_size;
  char* tag_name_data;
  const CubitBC_BitSet meshsets_mask;
  CubitMeshSets(Interface &moab,const EntityHandle _meshset);
  CubitMeshSets(Interface &moab,const CubitBC_BitSet _CubitBCType,const int _msId);
  inline int get_msId() const { return *msId; }
  inline CubitBC_BitSet get_CubitBCType() const { return CubitBCType; }

  inline EntityHandle get_meshset() const { return meshset; }
  inline unsigned long int get_CubitBCType_ulong() const { return CubitBCType.to_ulong(); }
  inline unsigned long int get_CubitBCType_mask_meshset_types_ulong() const { return (CubitBCType&meshsets_mask).to_ulong(); }
  inline unsigned long int get_CubitBCType_bc_data_types_ulong() const { return (CubitBCType&(~meshsets_mask)).to_ulong(); }

  PetscErrorCode get_Cubit_msId_entities_by_dimension(Interface &moab,const int dimension,Range &entities,const bool recursive = false) const;
  PetscErrorCode get_Cubit_msId_entities_by_dimension(Interface &moab,Range &entities,const bool recursive = false)  const;
  PetscErrorCode get_Cubit_msId_entities_by_type(Interface &moab,const EntityType type,Range &entities,const bool recursive = false) const;

  /** 
   *  \brief Function that returns the CubitBC_BitSet type of the contents of bc_data
   */
  PetscErrorCode get_type_from_bc_data(const vector<char> &bc_data,CubitBC_BitSet &type) const;

  /** 
   *  \brief Function that returns the CubitBC_BitSet type of the contents of bc_data
  */
  PetscErrorCode get_type_from_bc_data(CubitBC_BitSet &type) const;
    
  /**
   * \brief get bc_data vector from MoFEM database
   * 
   * \param bc_data is the in/out vector were bc_data will be stored
   */
  PetscErrorCode get_Cubit_bc_data(vector<char>& bc_data) const;
    
  /**
  * \brief get block_headers vector from MoFEM database
  *
  * \param material_data is the in/out vector were the material data will be stored
  */
  PetscErrorCode get_Cubit_block_header_data(vector<unsigned int>& material_data) const;

  /**
  * \brief print material_data int stream given by os
  *
  * f.e. it->print_Cubit_material_data(cout), i.e. printing to standard output
  * f.e. it->print_Cubit_material_data(cerr), i.e. printing to standard error output
  */
  PetscErrorCode print_Cubit_block_header_data(ostream& os) const;
    
  /**
   * \brief print bc_data int stream given by os
   *
   * f.e. it->print_Cubit_bc_data(cout), i.e. printing to standard output
   * f.e. it->print_Cubit_bc_data(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_bc_data(ostream& os) const;

  template<class _CUBIT_BC_DATA_TYPE_>
  PetscErrorCode get_cubit_bc_data_structure(_CUBIT_BC_DATA_TYPE_& data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((CubitBCType&data.type).none()) {
      SETERRQ(PETSC_COMM_SELF,1,"bc_data are not for _CUBIT_BC_DATA_TYPE_ structure");  
    }
    vector<char> bc_data;
    get_Cubit_bc_data(bc_data);
    ierr = data.fill_data(bc_data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   *  \brief Function that returns the CubitBC_BitSet type of the block name, sideset name etc.
   */
  PetscErrorCode get_type_from_Cubit_name(const string &name,CubitBC_BitSet &type) const;

  /**
   *  \brief Function that returns the CubitBC_BitSet type of the block name, sideset name etc.
   */
  PetscErrorCode get_type_from_Cubit_name(CubitBC_BitSet &type) const;
    
  /**
   * \brief get Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be stored
   */
  PetscErrorCode get_Cubit_attributes(vector<double> &attributes) const;

  /**
   * \brief print the attributes vector
   *
   * f.e. it->print_Cubit_attributes(cout), i.e. printing to standard output
   * f.e. it->print_Cubit_attributes(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_attributes(ostream& os) const;

  /**
   * \brief get name of block, sideset etc. (this is set in Cubit block properties)
   *
   * Block Name Conventions:
   * -----------------------
   * Materials are defined with block names starting with MAT_
   * e.g. MAT_ELASTIC_abcd, MAT_FRACTcdef etc.
   * Solution procedures are defined with block names starting with SOL_ e.g.
   * SOL_ELASTIC_xx, SOL_NLELASTICxx, SOL_FRACTabcd etc.
   *
   * List of materials/solution procedures
   * ---------------------------------------------------------------------------
   * Block name /  Number of attributes  / (1) Attribute 1, (2) Attribute 2 etc.
   * ---------------------------------------------------------------------------
   *
   * MAT_ELASTIC / 10 /  (1) Young's  modulus
   *                    (2) Poisson's ratio
   *                    (3) User attribute 8
   *                    ...
   *                    (10) User attribute 8
   *
   * MAT_ELASTIC_TRANSISO / 5 / (1) Young's modulus in xy plane (Ep)
   *                    (2) Young's modulus in z-direction (Ez)
   *                    (3) Poisson's ratio in xy plane (vp)
   *                    (4) Poisson's ratio in z-direction (vpz)
   *                    (5) Shear modulus in z-direction (Gzp)
   *
   * MAT_INTERF / 1 /   (1) Elastic modulus multiplier
   *
   * To be extended as appropriate
   */
   string get_Cubit_name() const;

  /**
   * \brief print name of block, sideset etc. (this is set in Cubit setting properties)
   *
   * e.g. it->print_Cubit_name(cout), i.e. printing to standard output
   * e.g it->print_Cubit_name(cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_Cubit_name(ostream& os) const;
    
  template<class _ATTRIBUTE_TYPE_>
  PetscErrorCode get_attribute_data_structure(_ATTRIBUTE_TYPE_ &data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((CubitBCType&data.type).none()) {
        SETERRQ(PETSC_COMM_SELF,1,"attributes are not for _ATTRIBUTE_TYPE_ structure");
    }
    vector<double> attributes;
    ierr = get_Cubit_attributes(attributes); CHKERRQ(ierr);
    ierr = data.fill_data(attributes); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  template<class _ATTRIBUTE_TYPE_>
  PetscErrorCode set_attribute_data_structure(_ATTRIBUTE_TYPE_ &data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((CubitBCType&data.type).none()) {
        SETERRQ(PETSC_COMM_SELF,1,"attributes are not for _ATTRIBUTE_TYPE_ structure");
    }
    double *ptr = const_cast<double*>(tag_block_attributes);
    ierr = data.set_data(ptr,8*tag_block_attributes_size); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
     
  friend ostream& operator<<(ostream& os,const CubitMeshSets& e);

  private:
  Tag nsTag,ssTag,nsTag_data,ssTag_data,bhTag,bhTag_header,block_attribs,entityNameTag;
  PetscErrorCode get_tags_hanlders(Interface &moab);
    
};
    
/**
 * @relates multi_index_container
 * \brief CubitMeshSet_multiIndex
 *
 * \param hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_bc_data_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_bc_data_types_ulong> >,
 * \param ordered_non_unique<
      tag<CubitMeshSets_name>, const_mem_fun<CubitMeshSets,string,&CubitMeshSets::get_Cubit_name> >,
 *
 * \param    hashed_unique<
      tag<Composite_mi_tag>,       
      composite_key<
	CubitMeshSets, <br>
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> > >
 *
 */
typedef multi_index_container<
  CubitMeshSets,
  indexed_by<
    hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
    ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_bc_data_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_bc_data_types_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_name>, const_mem_fun<CubitMeshSets,string,&CubitMeshSets::get_Cubit_name> >,
    hashed_unique<
      tag<Composite_Cubit_msId_and_MeshSetType_mi_tag>,       
      composite_key<
	CubitMeshSets,
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_CubitBCType_mask_meshset_types_ulong> > >
  > > CubitMeshSet_multiIndex;

}

#endif // __BCMULTIINDICES_HPP__
