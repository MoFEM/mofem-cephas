/** \file BCMultiIndices.hpp
 * \brief Multi-index containers, data boundary data structures and other low-level functions
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
 * \ingroup mofem_bc
 *
 */
struct CubitMeshSets {
  EntityHandle meshset;
  CubitBCType cubitBcType; 	///< type of meshset from cubit NodeSet, BlockSet, SideSet and more
  std::vector<Tag> tag_handles;	///< vector of tag handles to types of data passed from cubit
  int *msId;			///< cubit meshset ID
  char* tag_bc_data;
  int tag_bc_size;
  unsigned int *tag_block_header_data;
  double* tag_block_attributes;
  int tag_block_attributes_size;
  char* tag_name_data;
  const CubitBCType meshsets_mask;
  CubitMeshSets(Interface &moab,const EntityHandle _meshset);
  CubitMeshSets(Interface &moab,const CubitBCType _cubit_bc_type,const int _msId);

  inline int get_msId() const { return *msId; }
  inline CubitBCType get_cubit_bc_type() const { return cubitBcType; }
  inline EntityHandle getMeshSet() const { return meshset; }

  /** \deprecated Use getMeshSet() instead
  */
  DEPRECATED inline EntityHandle get_meshset() const { return getMeshSet(); }

  inline unsigned long int get_cubit_bc_type_ulong() const { return cubitBcType.to_ulong(); }
  inline unsigned long int get_cubit_bc_type_mask_meshset_types_ulong() const { return (cubitBcType&meshsets_mask).to_ulong(); }
  inline unsigned long int get_cubit_bc_type_bc_data_types_ulong() const { return (cubitBcType&(~meshsets_mask)).to_ulong(); }

  PetscErrorCode get_cubit_msId_entities_by_dimension(Interface &moab,const int dimension,Range &entities,const bool recursive = false) const;
  PetscErrorCode get_cubit_msId_entities_by_dimension(Interface &moab,Range &entities,const bool recursive = false)  const;
  PetscErrorCode get_cubit_msId_entities_by_type(Interface &moab,const EntityType type,Range &entities,const bool recursive = false) const;

  /**
   *  \brief Function that returns the CubitBCType type of the contents of bc_data
   */
  PetscErrorCode get_type_from_bc_data(const std::vector<char> &bc_data,CubitBCType &type) const;

  /**
   *  \brief Function that returns the CubitBCType type of the contents of bc_data
  */
  PetscErrorCode get_type_from_bc_data(CubitBCType &type) const;

  /**
   * \brief get bc_data vector from MoFEM database
   *
   * \param bc_data is the in/out vector were bc_data will be stored
   */
  PetscErrorCode get_bc_data(std::vector<char>& bc_data) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode get_Cubit_bc_data(std::vector<char>& bc_data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = get_bc_data(bc_data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  * \brief get block_headers vector from MoFEM database
  *
  * \param material_data is the in/out vector were the material data will be stored
  */
  PetscErrorCode get_block_header_data(std::vector<unsigned int>& material_data) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode get_Cubit_block_header_data(std::vector<unsigned int>& material_data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = get_block_header_data(material_data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
  * \brief print material_data int stream given by os
  *
  * f.e. it->print_Cubit_material_data(cout), i.e. printing to standard output
  * f.e. it->print_Cubit_material_data(std::cerr), i.e. printing to standard error output
  */
  PetscErrorCode print_block_header_data(std::ostream& os) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode print_Cubit_block_header_data(std::ostream& os) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = print_block_header_data(os); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * \brief print bc_data int stream given by os
   *
   * f.e. it->print_bc_data(cout), i.e. printing to standard output
   * f.e. it->print_bc_data(std::cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_bc_data(std::ostream& os) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode print_Cubit_bc_data(std::ostream& os) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = print_bc_data(os); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  template<class _CUBIT_BC_DATA_TYPE_>
  PetscErrorCode get_bc_data_structure(_CUBIT_BC_DATA_TYPE_& data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((cubitBcType&data.type).none()) {
      SETERRQ(PETSC_COMM_SELF,1,"bc_data are not for _CUBIT_BC_DATA_TYPE_ structure");
    }
    std::vector<char> bc_data;
    get_bc_data(bc_data);
    ierr = data.fill_data(bc_data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** deprecated \deprecated
  */
  template<class _CUBIT_BC_DATA_TYPE_>
  DEPRECATED
  PetscErrorCode get_cubit_bc_data_structure(_CUBIT_BC_DATA_TYPE_& data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = get_bc_data_structure(data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   *  \brief Function that returns the CubitBCType type of the block name, sideset name etc.
   */
  PetscErrorCode get_type_from_name(const std::string &name,CubitBCType &type) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode get_type_from_Cubit_name(const std::string &name,CubitBCType &type) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = get_type_from_name(name,type); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   *  \brief Function that returns the CubitBCType type of the block name, sideset name etc.
   */
  PetscErrorCode get_type_from_name(CubitBCType &type) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode get_type_from_Cubit_name(CubitBCType &type) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = get_type_from_name(type); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * \brief get Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be stored
   */
  PetscErrorCode get_attributes(std::vector<double> &attributes) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode get_Cubit_attributes(std::vector<double> &attributes) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = get_attributes(attributes); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * \brief print the attributes vector
   *
   * f.e. it->print_attributes(cout), i.e. printing to standard output
   * f.e. it->print_attributes(std::cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_attributes(std::ostream& os) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode print_Cubit_attributes(std::ostream& os) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = print_attributes(os); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /**
   * \brief get name of block, sideset etc. (this is set in Cubit block properties)
   *
   * Block Name Conventions:

   * Materials are defined with block names starting with MAT_
   * e.g. MAT_ELASTIC_abcd, MAT_FRACTcdef etc.
   * Solution procedures are defined with block names starting with SOL_ e.g.
   * SOL_ELASTIC_xx, SOL_NLELASTICxx, SOL_FRACTabcd etc.
   *
   * List of materials/solution procedures

   * Block name /  Number of attributes  / (1) Attribute 1, (2) Attribute 2 etc.
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
  std::string getName() const;

  /** deprecated \deprecated
  */
  DEPRECATED std::string get_name() const { return getName(); }

  /**
   * \brief print name of block, sideset etc. (this is set in Cubit setting properties)
   *
   * e.g. it->print_name(cout), i.e. printing to standard output
   * e.g it->print_name(std::cerr), i.e. printing to standard error output
   */
  PetscErrorCode print_name(std::ostream& os) const;

  /** deprecated \deprecated
  */
  DEPRECATED PetscErrorCode print_Cubit_name(std::ostream& os) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = print_name(os); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  template<class _ATTRIBUTE_TYPE_>
  PetscErrorCode get_attribute_data_structure(_ATTRIBUTE_TYPE_ &data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((cubitBcType&data.type).none()) {
        SETERRQ(PETSC_COMM_SELF,1,"attributes are not for _ATTRIBUTE_TYPE_ structure");
    }
    std::vector<double> attributes;
    ierr = get_attributes(attributes); CHKERRQ(ierr);
    ierr = data.fill_data(attributes); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  template<class _ATTRIBUTE_TYPE_>
  PetscErrorCode set_attribute_data_structure(_ATTRIBUTE_TYPE_ &data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((cubitBcType&data.type).none()) {
        SETERRQ(PETSC_COMM_SELF,1,"attributes are not for _ATTRIBUTE_TYPE_ structure");
    }
    double *ptr = const_cast<double*>(tag_block_attributes);
    ierr = data.set_data(ptr,8*tag_block_attributes_size); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  friend std::ostream& operator<<(std::ostream& os,const CubitMeshSets& e);

  Tag nsTag,ssTag,nsTag_data,ssTag_data,bhTag,bhTag_header,block_attribs,entityNameTag;
  PetscErrorCode get_tags_hanlders(Interface &moab);

};

/**
 * @relates multi_index_container
 * \brief CubitMeshSet_multiIndex
 */
typedef multi_index_container<
  CubitMeshSets,
  indexed_by<
    hashed_unique<
      tag<Meshset_mi_tag>, member<CubitMeshSets,EntityHandle,&CubitMeshSets::meshset> >,
    ordered_non_unique<
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_cubit_bc_type_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_cubit_bc_type_mask_meshset_types_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_bc_data_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_cubit_bc_type_bc_data_types_ulong> >,
    ordered_non_unique<
      tag<CubitMeshSets_name>, const_mem_fun<CubitMeshSets,std::string,&CubitMeshSets::getName> >,
    hashed_unique<
      tag<Composite_Cubit_msId_And_MeshSetType_mi_tag>,
      composite_key<
	CubitMeshSets,
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::get_msId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::get_cubit_bc_type_mask_meshset_types_ulong> > >
  > > CubitMeshSet_multiIndex;

  struct CubitMeshSets_change_add_bit_to_cubit_bc_type {
    CubitBCType bit;
    CubitMeshSets_change_add_bit_to_cubit_bc_type(const CubitBCType &_bit): bit(_bit) {};
    void operator()(CubitMeshSets &e) {
      e.cubitBcType |= bit;
    }
  };

  struct CubitMeshSets_change_name {
    Interface &mOab;
    std::string nAme;
    CubitMeshSets_change_name(Interface &moab,const std::string &name):
    mOab(moab),
    nAme(name) {
    };
    void operator()(CubitMeshSets &e);
  };

}

#endif // __BCMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup mofem_bc Boundary conditions
 * \ingroup mofem
 ******************************************************************************/
