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

  /**
   * \brief get meshset id as it set in preprocessing software
   * @return id of meshset
   */
  inline int getMeshSetId() const { return *msId; }

  /** \deprecated use getMeshSetId() instead
  */
  DEPRECATED inline int get_msId() const { return getMeshSetId(); }

  /**
   * \brief get type of meshset
   *
   * See CubitBC for set of types of meshsets.
   *
   * @return meshset type
   */
  inline CubitBCType getBcType() const { return cubitBcType; }

  /**
   * \brief get bc meshset
   * @return meshset entity handle
   */
  inline EntityHandle getMeshSet() const { return meshset; }

  /** \deprecated Use getMeshSet() instead
  */
  DEPRECATED inline EntityHandle get_meshset() const { return getMeshSet(); }

  /**
   * \brief get bc meshset type
   * @return return type as unsigned integer
   */
  inline unsigned long int getBcTypeULong() const {
    return cubitBcType.to_ulong();
  }

  /**
   * \brief get meshset type and mask
   * @return type is returned as unsigned integer
   */
  inline unsigned long int getMaksedBcTypeULong() const {
    return (cubitBcType&meshsets_mask).to_ulong();
  }

  /**
   * \brief get entities form meshset
   * @param  moab      moab instance
   * @param  dimension dimension of entities
   * @param  entities  range of returned entities
   * @param  recursive true if meshset should be searched recursively
   * @return           error code
   */
  PetscErrorCode getMeshSetIdEntitiesByDimension(
    Interface &moab,const int dimension,Range &entities,const bool recursive = false
  ) const;

  /** \deprecated Use getMeshSetIdEntitiesByDimension() instead
  */
  DEPRECATED inline PetscErrorCode get_cubit_msId_entities_by_dimension(
    Interface &moab,const int dimension,Range &entities,const bool recursive = false
  ) const {
    return getMeshSetIdEntitiesByDimension(moab,dimension,entities,recursive);
  }

  /**
   * \brief get entities form meshset
   *
   * Use if meshset have predefined dimension
   *
   * @param  moab      moab instance
   * @param  entities  range of returned entities
   * @param  recursive true if meshset should be searched recursively
   * @return           error code
   *
   */
  PetscErrorCode getMeshSetIdEntitiesByDimension(
    Interface &moab,Range &entities,const bool recursive = false
  )  const;

  /** \deprecated Use getMeshSetIdEntitiesByDimension() instead
  */
  DEPRECATED inline PetscErrorCode get_cubit_msId_entities_by_dimension(
    Interface &moab,Range &entities,const bool recursive = false
  )  const {
    return getMeshSetIdEntitiesByDimension(moab,entities,recursive);
  }

  /**
   * \brief get entities by type
   * @param  moab      moab instance
   * @param  type      type of entity
   * @param  entities  returned entities
   * @param  recursive true if meshset should be searched recursively
   * @return           error code
   */
  PetscErrorCode getMeshSetIdEntitiesByType(
    Interface &moab,const EntityType type,Range &entities,const bool recursive = false
  ) const;

  /** \deprecated Use getMeshSetIdEntitiesByType() instead
  */
  DEPRECATED inline PetscErrorCode get_cubit_msId_entities_by_type(
    Interface &moab,const EntityType type,Range &entities,const bool recursive = false
  ) const {
    return getMeshSetIdEntitiesByType(moab,type,entities,recursive);
  }

  /**
   *  \brief Function that returns the CubitBCType type of the contents of bc_data
   */
  PetscErrorCode getTypeFromBcData(const std::vector<char> &bc_data,CubitBCType &type) const;

  /**
   *  \brief Function that returns the CubitBCType type of the contents of bc_data
  */
  PetscErrorCode getTypeFromBcData(CubitBCType &type) const;

  /**
   * \brief get bc_data vector from MoFEM database
   *
   * \param bc_data is the in/out vector were bc_data will be stored
   */
  PetscErrorCode getBcData(std::vector<char>& bc_data) const;

  /**
  * \brief get block_headers vector from MoFEM database
  *
  * \param material_data is the in/out vector were the material data will be stored
  */
  PetscErrorCode getBlockHeaderData(std::vector<unsigned int>& material_data) const;

  /**
  * \brief print material_data int stream given by os
  *
  * f.e. it->print_Cubit_material_data(cout), i.e. printing to standard output
  * f.e. it->print_Cubit_material_data(std::cerr), i.e. printing to standard error output
  */
  PetscErrorCode printBlockHeaderData(std::ostream& os) const;

  /**
   * \brief print bc_data int stream given by os
   *
   * f.e. it->printBcData(cout), i.e. printing to standard output
   * f.e. it->printBcData(std::cerr), i.e. printing to standard error output
   */
  PetscErrorCode printBcData(std::ostream& os) const;

  /**
   *  \brief Function that returns the CubitBCType type of the block name, sideset name etc.
   */
  PetscErrorCode getTypeFromName(const std::string &name,CubitBCType &type) const;

  /**
   *  \brief Function that returns the CubitBCType type of the block name, sideset name etc.
   */
  PetscErrorCode getTypeFromName(CubitBCType &type) const;

  /**
   * \brief get Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be stored
   */
  PetscErrorCode getAttributes(std::vector<double> &attributes) const;

  /**
   * \brief cet Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be stored
   */
  PetscErrorCode setAttributes(moab::Interface &moab,const std::vector<double> &attributes);

  /** \deprecated Use getAttributes() instead
  */
  DEPRECATED inline PetscErrorCode get_attributes(std::vector<double> &attributes) const {
    return getAttributes(attributes);
  }

  /**
   * \brief print the attributes vector
   *
   * f.e. it->printAttributes(cout), i.e. printing to standard output
   * f.e. it->printAttributes(std::cerr), i.e. printing to standard error output
   */
  PetscErrorCode printAttributes(std::ostream& os) const;

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

  /** \deprecated Use getName() instead
  */
  DEPRECATED inline std::string get_name() const { return getName(); }

  /**
   * \brief print name of block, sideset etc. (this is set in Cubit setting properties)
   *
   * e.g. it->printName(cout), i.e. printing to standard output
   * e.g it->printName(std::cerr), i.e. printing to standard error output
   */
  PetscErrorCode printName(std::ostream& os) const;

  /**
   * \brief fill data structure with data saved on meshset
   */
  template<class ATTRIBUTE_TYPE>
  PetscErrorCode getAttributeDataStructure(ATTRIBUTE_TYPE &data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((cubitBcType&data.getType()).none()) {
      SETERRQ(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "attributes are not for ATTRIBUTE_TYPE structure"
      );
    }
    std::vector<double> attributes;
    ierr = getAttributes(attributes); CHKERRQ(ierr);
    ierr = data.fill_data(attributes); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** \deprecated Use getAttributeDataStructure() instead
  */
  template<class ATTRIBUTE_TYPE>
  DEPRECATED inline PetscErrorCode get_attribute_data_structure(ATTRIBUTE_TYPE &data) const {
    return getAttributeDataStructure(data);
  }

  /**
   * \brief fill meshset data with data on structure
   */
  template<class ATTRIBUTE_TYPE>
  PetscErrorCode setAttributeDataStructure(const ATTRIBUTE_TYPE &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((cubitBcType&data.getType()).none()) {
        SETERRQ(PETSC_COMM_SELF,1,"attributes are not for ATTRIBUTE_TYPE structure");
    }
    double *ptr = const_cast<double*>(tag_block_attributes);
    ierr = data.set_data(ptr,8*tag_block_attributes_size); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  template<class CUBIT_BC_DATA_TYPE>
  PetscErrorCode getBcDataStructure(CUBIT_BC_DATA_TYPE& data) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if((cubitBcType&data.tYpe).none()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"bc_data are not for CUBIT_BC_DATA_TYPE structure");
    }
    std::vector<char> bc_data;
    getBcData(bc_data);
    ierr = data.fill_data(bc_data); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** \deprecated Use getBcDataStructure() instead
  */
  template<class CUBIT_BC_DATA_TYPE>
  DEPRECATED inline PetscErrorCode get_bc_data_structure(CUBIT_BC_DATA_TYPE& data) const {
    return getBcDataStructure(data);
  }

  template<class CUBIT_BC_DATA_TYPE>
  PetscErrorCode setBcDataStructure(CUBIT_BC_DATA_TYPE& data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    char *ptr = const_cast<char*>(tag_bc_data);
    ierr = data.set_data(ptr,tag_bc_size); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  friend std::ostream& operator<<(std::ostream& os,const CubitMeshSets& e);

  Tag nsTag,ssTag,nsTag_data,ssTag_data,bhTag,bhTag_header,thBlockAttribs,entityNameTag;

  PetscErrorCode getTagsHanlders(Interface &moab);

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
      tag<CubitMeshSets_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::getBcTypeULong> >,
    ordered_non_unique<
      tag<CubitMeshSets_mask_meshset_mi_tag>, const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::getMaksedBcTypeULong> >,
    ordered_non_unique<
      tag<CubitMeshSets_name>, const_mem_fun<CubitMeshSets,std::string,&CubitMeshSets::getName> >,
    hashed_unique<
      tag<Composite_Cubit_msId_And_MeshSetType_mi_tag>,
      composite_key<
	CubitMeshSets,
	  const_mem_fun<CubitMeshSets,int,&CubitMeshSets::getMeshSetId>,
	  const_mem_fun<CubitMeshSets,unsigned long int,&CubitMeshSets::getMaksedBcTypeULong> > >
  > > CubitMeshSet_multiIndex;

/** \brief change meshset type
*/
struct CubitMeshSets_change_add_bit_to_cubit_bc_type {
  CubitBCType bIt;
  CubitMeshSets_change_add_bit_to_cubit_bc_type(const CubitBCType &bit):
  bIt(bit) {};
  void operator()(CubitMeshSets &e);
};

/**
 * \brief change meshset name
 */
struct CubitMeshSets_change_name {
  Interface &mOab;
  std::string nAme;
  CubitMeshSets_change_name(Interface &moab,const std::string &name):
  mOab(moab),
  nAme(name) {
  };
  void operator()(CubitMeshSets &e);
};

/**
 * change meshset attributes
 */
struct CubitMeshSets_change_attributes {
  Interface &mOab;
  const std::vector<double> &aTtr;
  CubitMeshSets_change_attributes(Interface &moab,const std::vector<double> &attr):
  mOab(moab),
  aTtr(attr) {}
  void operator()(CubitMeshSets &e);
};

/**
 * change meshset attributes for material data structure
 */
struct CubitMeshSets_change_attributes_data_structure {
  Interface &mOab;
  const GenericAttributeData &aTtr;
  CubitMeshSets_change_attributes_data_structure(
    Interface &moab,const GenericAttributeData &attr
  ):
  mOab(moab),
  aTtr(attr) {}
  void operator()(CubitMeshSets &e);
};

/**
 * change meshset attributes for material data structure
 */
struct CubitMeshSets_change_bc_data_structure {
  Interface &mOab;
  const GenericCubitBcData &bcData;
  CubitMeshSets_change_bc_data_structure(
    Interface &moab,const GenericCubitBcData &bc_data
  ):
  mOab(moab),
  bcData(bc_data) {}
  void operator()(CubitMeshSets &e);
};

}

#endif // __BCMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup mofem_bc Boundary conditions
 * \ingroup mofem
 ******************************************************************************/
