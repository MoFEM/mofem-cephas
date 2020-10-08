/** \file BCMultiIndices.hpp
 * \brief Multi-index containers, data boundary data structures and other
 * low-level functions
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
 * \brief this struct keeps basic methods for moab meshset about material and
 * boundary conditions \ingroup mofem_bc
 *
 */
struct CubitMeshSets {

  EntityHandle meshset;
  CubitBCType cubitBcType; ///< type of meshset from cubit NodeSet, BlockSet,
                           ///< SideSet and more
  std::vector<Tag>
      tag_handles; ///< vector of tag handles to types of data passed from cubit
  int *msId;       ///< cubit meshset ID
  char *tag_bc_data;
  int tag_bc_size;
  unsigned int *tag_block_header_data;
  double *tag_block_attributes;
  int tag_block_attributes_size;
  char *tagName;
  const CubitBCType meshsets_mask;

  CubitMeshSets(Interface &moab, const EntityHandle _meshset);
  CubitMeshSets(Interface &moab, const CubitBCType _cubit_bc_type,
                const int _msId);

  /**
   * \brief get meshset id as it set in preprocessing software
   * @return id of meshset
   */
  inline int getMeshsetId() const { return *msId; }

  // /** \deprecated use getMeshsetId() instead
  // */
  // DEPRECATED inline int get_msId() const { return getMeshsetId(); }

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
  inline EntityHandle getMeshset() const { return meshset; }

  /**
   * \brief get bc meshset type
   * @return return type as unsigned int
   */
  inline unsigned long int getBcTypeULong() const {
    return cubitBcType.to_ulong();
  }

  /**
   * \brief get meshset type and mask
   * @return type is returned as unsigned integer
   */
  inline unsigned long int getMaksedBcTypeULong() const {
    return (cubitBcType & meshsets_mask).to_ulong();
  }

  /**
   * @brief Get the meshset entities dimension
   *
   * \note If dimension is -1, then dimension for meshset ins undetermined.
   *
   * @return unsigned int
   */
  unsigned int getMeshsetEntitiesDimension() const {
    if (tag_block_header_data)
      return tag_block_header_data[2];
    else
      return -1;
  }

  /**
   * \brief get entities form meshset
   * @param  moab      moab instance
   * @param  dimension dimension of entities
   * @param  entities  range of returned entities
   * @param  recursive true if meshset should be searched recursively
   * @return           error code
   */
  MoFEMErrorCode
  getMeshsetIdEntitiesByDimension(Interface &moab, const int dimension,
                                  Range &entities,
                                  const bool recursive = false) const;

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
  MoFEMErrorCode
  getMeshsetIdEntitiesByDimension(Interface &moab, Range &entities,
                                  const bool recursive = false) const;

  /**
   * \brief get entities by type
   * @param  moab      moab instance
   * @param  type      type of entity
   * @param  entities  returned entities
   * @param  recursive true if meshset should be searched recursively
   * @return           error code
   */
  MoFEMErrorCode getMeshsetIdEntitiesByType(Interface &moab,
                                            const EntityType type,
                                            Range &entities,
                                            const bool recursive = false) const;

  /**
   *  \brief Function that returns the CubitBCType type of the contents of
   * bc_data
   */
  MoFEMErrorCode getTypeFromBcData(const std::vector<char> &bc_data,
                                   CubitBCType &type) const;

  /**
   *  \brief Function that returns the CubitBCType type of the contents of
   * bc_data
   */
  MoFEMErrorCode getTypeFromBcData(CubitBCType &type) const;

  /**
   * \brief get bc_data vector from MoFEM database
   *
   * \param bc_data is the in/out vector were bc_data will be stored
   */
  MoFEMErrorCode getBcData(std::vector<char> &bc_data) const;

  /**
   * \brief get block_headers vector from MoFEM database
   *
   * \param material_data is the in/out vector were the material data will be
   * stored
   */
  MoFEMErrorCode
  getBlockHeaderData(std::vector<unsigned int> &material_data) const;

  /**
   * \brief print material_data int stream given by os
   *
   * f.e. it->print_Cubit_material_data(cout), i.e. printing to standard output
   * f.e. it->print_Cubit_material_data(std::cerr), i.e. printing to standard
   * error output
   */
  MoFEMErrorCode printBlockHeaderData(std::ostream &os) const;

  /**
   * \brief print bc_data int stream given by os
   *
   * f.e. it->printBcData(cout), i.e. printing to standard output
   * f.e. it->printBcData(std::cerr), i.e. printing to standard error output
   */
  MoFEMErrorCode printBcData(std::ostream &os) const;

  /**
   *  \brief Function that returns the CubitBCType type of the block name,
   * sideset name etc.
   */
  MoFEMErrorCode getTypeFromName(const std::string &name,
                                 CubitBCType &type) const;

  /**
   *  \brief Function that returns the CubitBCType type of the block name,
   * sideset name etc.
   */
  MoFEMErrorCode getTypeFromName(CubitBCType &type) const;

  /**
   * \brief get Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be
   * stored
   */
  MoFEMErrorCode getAttributes(std::vector<double> &attributes) const;

  /**
   * \brief cet Cubit block attributes
   *
   * \param attributes is the vector where the block attribute data will be
   * stored
   */
  MoFEMErrorCode setAttributes(moab::Interface &moab,
                               const std::vector<double> &attributes);

  /**
   * \brief print the attributes vector
   *
   * f.e. it->printAttributes(cout), i.e. printing to standard output
   * f.e. it->printAttributes(std::cerr), i.e. printing to standard error output
   */
  MoFEMErrorCode printAttributes(std::ostream &os) const;

  /**
   * \brief get name of block, sideset etc. (this is set in Cubit block
   properties)
   *
   * Block Name Conventions:

   * Materials are defined with block names starting with MAT_
   * e.g. MAT_ELASTIC_abcd.
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

  /**
   * \brief print name of block, sideset etc. (this is set in Cubit setting
   * properties)
   *
   * e.g. it->printName(cout), i.e. printing to standard output
   * e.g it->printName(std::cerr), i.e. printing to standard error output
   */
  MoFEMErrorCode printName(std::ostream &os) const;

  /**
   * \brief fill data structure with data saved on meshset
   */
  template <class ATTRIBUTE_TYPE>
  MoFEMErrorCode getAttributeDataStructure(ATTRIBUTE_TYPE &data) const {
    MoFEMFunctionBegin;
    if ((cubitBcType & data.getType()).none()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "attributes are not for ATTRIBUTE_TYPE structure");
    }
    std::vector<double> attributes;
    CHKERR getAttributes(attributes);
    CHKERR data.fill_data(attributes);
    MoFEMFunctionReturn(0);
  }

  /**
   * \brief fill meshset data with data on structure
   */
  template <class ATTRIBUTE_TYPE>
  MoFEMErrorCode setAttributeDataStructure(const ATTRIBUTE_TYPE &data) {
    MoFEMFunctionBegin;
    if ((cubitBcType & data.getType()).none()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "attributes are not for ATTRIBUTE_TYPE structure");
    }
    double *ptr = const_cast<double *>(tag_block_attributes);
    CHKERR data.set_data(ptr, 8 * tag_block_attributes_size);
    MoFEMFunctionReturn(0);
  }

  template <class CUBIT_BC_DATA_TYPE>
  MoFEMErrorCode getBcDataStructure(CUBIT_BC_DATA_TYPE &data) const {
    MoFEMFunctionBeginHot;

    if ((cubitBcType & data.tYpe).none()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "bc_data are not for CUBIT_BC_DATA_TYPE structure");
    }
    std::vector<char> bc_data;
    getBcData(bc_data);
    ierr = data.fill_data(bc_data);
    CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
  }

  template <class CUBIT_BC_DATA_TYPE>
  MoFEMErrorCode setBcDataStructure(CUBIT_BC_DATA_TYPE &data) {
    MoFEMFunctionBeginHot;

    char *ptr = const_cast<char *>(tag_bc_data);
    ierr = data.set_data(ptr, tag_bc_size);
    CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
  }

  friend std::ostream &operator<<(std::ostream &os, const CubitMeshSets &e);

  Tag nsTag, ssTag, nsTag_data, ssTag_data, bhTag, bhTag_header, thBlockAttribs,
      entityNameTag;

  MoFEMErrorCode getTagsHandlers(Interface &moab);
};

/**
 * \typedef CubitMeshSet_multiIndex
 * \brief Stores data about meshsets (see CubitMeshSets) storing data about
 * boundary conditions, interfaces, sidesets, nodests, blocksets
 *
 * \param Meshset_mi_tag  index by meshset handle
 * \param CubitMeshSets_mi_tag index by bc type, see CubitBC
 * \param CubitMeshSets_mask_meshset_mi_tag index by NODESET, SIDESET, BLOCKSET
 * only
 *
 * \param CubitMeshSets_name index by meshset name
 *
 * \param Composite_Cubit_msId_And_MeshSetType_mi_tag index by meshset id and
 * type NODESET, SIDESET or BLOCKSET
 *
 *  Example:
 *  \code
 *   MeshsetsManager *m_mng;
 *   CHKERR m_field.getInterface(m_mng);
 *   auto &index = m_mng->getMeshsetsMultindex();
 *
 *
 *   auto mit =
 * index.get<CubitMeshSets_mask_meshset_mi_tag>().lower_bound(BLOCKSET); auto
 * hi_mit =
 * index.get<CubitMeshSets_mask_meshset_mi_tag>().upper_bound(BLOCKSET);
 *
 *   // Make a loop over all BLOCKSET
 *   for(;mit!=hi_mit;mit++) {
 *     int id = mit->getMeshsetId();            // get blockset id
 *     EntityHandle handle = mit->getMeshset(); // get block meshset
 *     std::vector< double > attributes;
 *     // get block attributes
 *     auto mit->getAttributes(attributes);
 *     // do something
 *   }
 *  \endcode
 *
 */
typedef multi_index_container<
    CubitMeshSets,
    indexed_by<
        hashed_unique<tag<Meshset_mi_tag>, member<CubitMeshSets, EntityHandle,
                                                  &CubitMeshSets::meshset>>,
        ordered_non_unique<tag<CubitMeshSets_mi_tag>,
                           const_mem_fun<CubitMeshSets, unsigned long int,
                                         &CubitMeshSets::getBcTypeULong>>,
        ordered_non_unique<tag<CubitMeshSets_mask_meshset_mi_tag>,
                           const_mem_fun<CubitMeshSets, unsigned long int,
                                         &CubitMeshSets::getMaksedBcTypeULong>>,
        ordered_non_unique<
            tag<CubitMeshSets_name>,
            const_mem_fun<CubitMeshSets, std::string, &CubitMeshSets::getName>>,
        hashed_unique<
            tag<Composite_Cubit_msId_And_MeshSetType_mi_tag>,
            composite_key<
                CubitMeshSets,
                const_mem_fun<CubitMeshSets, int, &CubitMeshSets::getMeshsetId>,
                const_mem_fun<CubitMeshSets, unsigned long int,
                              &CubitMeshSets::getMaksedBcTypeULong>>>>>
    CubitMeshSet_multiIndex;

/** \brief change meshset type
 */
struct CubitMeshSets_change_add_bit_to_cubit_bc_type {
  CubitBCType bIt;
  CubitMeshSets_change_add_bit_to_cubit_bc_type(const CubitBCType &bit)
      : bIt(bit){};
  void operator()(CubitMeshSets &e);
};

/**
 * \brief change meshset name
 */
struct CubitMeshSets_change_name {
  Interface &mOab;
  std::string nAme;
  CubitMeshSets_change_name(Interface &moab, const std::string &name)
      : mOab(moab), nAme(name){};
  void operator()(CubitMeshSets &e);
};

/**
 * change meshset attributes
 */
struct CubitMeshSets_change_attributes {
  Interface &mOab;
  const std::vector<double> &aTtr;
  CubitMeshSets_change_attributes(Interface &moab,
                                  const std::vector<double> &attr)
      : mOab(moab), aTtr(attr) {}
  void operator()(CubitMeshSets &e);
};

/**
 * change meshset attributes for material data structure
 */
struct CubitMeshSets_change_attributes_data_structure {
  Interface &mOab;
  const GenericAttributeData &aTtr;
  CubitMeshSets_change_attributes_data_structure(
      Interface &moab, const GenericAttributeData &attr)
      : mOab(moab), aTtr(attr) {}
  void operator()(CubitMeshSets &e);
};

/**
 * change meshset attributes for material data structure
 */
struct CubitMeshSets_change_bc_data_structure {
  Interface &mOab;
  const GenericCubitBcData &bcData;
  CubitMeshSets_change_bc_data_structure(Interface &moab,
                                         const GenericCubitBcData &bc_data)
      : mOab(moab), bcData(bc_data) {}
  void operator()(CubitMeshSets &e);
};

} // namespace MoFEM

#endif // __BCMULTIINDICES_HPP__

/**
 * \defgroup mofem_bc Boundary conditions
 * \ingroup mofem
 */
