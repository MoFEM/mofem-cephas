/** \file MeshsetsManager.hpp
 * \brief MeshsetsManager interface

 Interface to manage material and boundary sets

 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __MESHSETSMANAGER_HPP__
#define __MESHSETSMANAGER_HPP__

namespace MoFEM {

  /**
   * \brief Iterator that loops over all the Cubit MeshSets in a moFEM field
   * \ingroup mofem_bc

   *
   * \param MESHSET_MANAGER meshset manager (works as well with Interface)
   * \param iterator
   */
   #define _IT_CUBITMESHSETS_FOR_LOOP_(MESHSET_MANAGER,IT) \
   CubitMeshSet_multiIndex::iterator IT = MESHSET_MANAGER.get_meshsets_manager_ptr()->getBegin(); \
   IT!=MESHSET_MANAGER.get_meshsets_manager_ptr()->getEnd(); IT++

   /**
   * \brief Iterator that loops over a specific Cubit MeshSet in a moFEM field
   * \ingroup mofem_bc

   *
   * \param mField moFEM Field
   * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
   * \param iterator
   */
   #define _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MESHSET_MANAGER,CUBITBCTYPE,IT) \
   CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator IT = MESHSET_MANAGER.get_meshsets_manager_ptr()->getBegin(CUBITBCTYPE); \
   IT!=MESHSET_MANAGER.get_meshsets_manager_ptr()->getEnd(CUBITBCTYPE); IT++

   /**
   * \brief Iterator that loops over a specific Cubit MeshSet having a particular BC meshset in a moFEM field
   * \ingroup mofem_bc

   *
   * \param MESHSET_MANAGER meshset manager (works as well with Interface)
   * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
   * \param iterator
   *
   * Example: \code
   for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
   ...
   * } \endcode
   */
   #define _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MESHSET_MANAGER,CUBITBCTYPE,IT) \
   CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator IT = MESHSET_MANAGER.get_meshsets_manager_ptr()->getBySetTypeBegin(CUBITBCTYPE); \
   IT!=MESHSET_MANAGER.get_meshsets_manager_ptr()->getBySetTypeEnd(CUBITBCTYPE); IT++

   /**
   * \brief Iterator that loops over Cubit BlockSet having a particular name
   * \ingroup mofem_bc


   * \param MESHSET_MANAGER meshset manager (works as well with Interface)
   * \param NAME name
   * \param IT iterator
   *
   * Example: \code
   for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,"SOME_BLOCK_NANE",it) {
   ...
   * } \endcode
   */
   #define _IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(MESHSET_MANAGER,NAME,IT) \
   CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator IT = MESHSET_MANAGER.get_meshsets_manager_ptr()->getBegin(NAME); \
   IT!=MESHSET_MANAGER.get_meshsets_manager_ptr()->getEnd(NAME); IT++

  static const MOFEMuuid IDD_MOFEMMeshsetsManager = MOFEMuuid( BitIntefaceId(MESHSETSMANAGER_INTERFACE) );

  struct Core;

  /** \brief Interface for managing meshsets containing materials and boundary conditions
   * \ingroup mofem_bc
   */
  struct MeshsetsManager: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    MeshsetsManager(const MoFEM::Core& core);

    /**
     * \brief get tags handlers used on meshsets

     * On meshsets range of tages in set. Depending on tag type and data on that
     * tag type of meshset could be determined. This function get hanldes to
     * tags.
     *
     * Most of the tags are followinf convention used by MoAB or Cubit and other
     * meshing softwares, f.e. gmesh.

     */
    PetscErrorCode getTags(int verb = -1);

    /**
     * \brief get tag handle used to store "id" of NODESET
     */
    inline Tag get_nsTag() const { return nsTag; }

    /**
     * \brief get tag handle used to store "id" of SIDESET
     */
    inline Tag get_ssTag() const { return ssTag; }

    /**
     * \brief get tag handle used to store boundary data on NODESET
     */
    inline Tag get_nsTag_data() const { return nsTag_data; }

    /**
     * \brief get tag handle used to store boundary data on SIDESET
     */
    inline Tag get_ssTag_data() const { return ssTag_data; }

    /**
     * \brief get tag handle used to store "id" of BLOCKSET
     */
    inline Tag get_bhTag() const { return bhTag; }

    /**
     * \brief get tag handle used to store of block set header (Used by Cubit)
     */
    inline Tag get_bhTag_header() const { return bhTag_header; }

    /**
     * \brief return pointer to meshset manager
     */
    MeshsetsManager* get_meshsets_manager_ptr() { return this; }

    /**
     * \brief return pointer to meshset manager
     */
    const MeshsetsManager* get_meshsets_manager_ptr() const { return this; }

    /**
     * \brief clear multi-index container
     * @return error code
     */
    PetscErrorCode clearMap();

    /**
     * \brier initialize container form data on mesh
     * @return [description]
     */
    PetscErrorCode initialiseDatabseInformationFromMesh(int verb = 0);

    template<class CUBIT_BC_DATA_TYPE>
    PetscErrorCode printBcSet(CUBIT_BC_DATA_TYPE& data,unsigned long int type) const {
      PetscErrorCode ierr;
      MoABErrorCode rval;
      PetscFunctionBegin;
      try {
        const MoFEM::Interface& m_field = cOre;
        const moab::Interface& moab = m_field.get_moab();
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_((*this),type,it)) {
          ierr = it->getBcDataStructure(data); CHKERRQ(ierr);
          std::ostringstream ss;
          ss << *it << std::endl;
          ss << data << std::endl;
          Range tets,tris,edges,nodes;
          rval = moab.get_entities_by_type(it->meshset,MBTET,tets,true); CHKERRQ_MOAB(rval);
          rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
          rval = moab.get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERRQ_MOAB(rval);
          rval = moab.get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERRQ_MOAB(rval);
          ss << "name "<< it->getName() << std::endl;
          ss << "msId "<< it->getMeshsetId() << " nb. tets " << tets.size() << std::endl;
          ss << "msId "<< it->getMeshsetId() << " nb. tris " << tris.size() << std::endl;
          ss << "msId "<< it->getMeshsetId() << " nb. edges " << edges.size() << std::endl;
          ss << "msId "<< it->getMeshsetId() << " nb. nodes " << nodes.size() << std::endl;
          ss << std::endl;
          PetscPrintf(m_field.get_comm(),ss.str().c_str());
        }
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      }
      PetscFunctionReturn(0);
    }

    /**
     * \brief print meshsets with displacement boundary conditions data structure
     */
    PetscErrorCode printDisplacementSet() const;

    /**
     * \brief print meshsets with pressure boundary conditions data structure
     */
    PetscErrorCode printPressureSet() const;

    /**
     * \brief print meshsets with force boundary conditions data structure
     */
    PetscErrorCode printForceSet() const;

    /**
     * \brief print meshsets with temperature boundary conditions data structure
     */
    PetscErrorCode printTemperatureSet() const;

    /**
     * \brief print meshsets with heat flux boundary conditions data structure
     */
    PetscErrorCode printHeatFluxSet() const;

    /**
     * \brief print meshsets with material data structure set on it
     */
    PetscErrorCode printMaterialsSet() const;

    inline CubitMeshSet_multiIndex& getMeshsetsMultindex() {
      return cubitMeshsets;
    }

    /**
     * \ingroup mofem_bc
     * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
     *
     * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
     * 	...
     * }
     *
     */
    inline CubitMeshSet_multiIndex::iterator getBegin() const { return cubitMeshsets.begin(); }

    /**
     * \ingroup mofem_bc
     * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
     *
     * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
     * 	...
     * }
     *
     */
    CubitMeshSet_multiIndex::iterator getEnd() const { return cubitMeshsets.end(); }

    /**
      * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
      * \ingroup mofem_bc

      *
      * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
      * 	...
      * }
      *
      * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
      */
    inline CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator
    getBegin(const unsigned int cubit_bc_type) const {
      return cubitMeshsets.get<CubitMeshSets_mi_tag>().lower_bound(cubit_bc_type);
    }

    /**
      * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
      * \ingroup mofem_bc

      *
      * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
      * 	...
      * }
      *
      * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
      */
    inline CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator
    getEnd(const unsigned int cubit_bc_type) const {
      return cubitMeshsets.get<CubitMeshSets_mi_tag>().upper_bound(cubit_bc_type);
    }

    /**
      * \brief get end iterator of cubit meshset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
      * \ingroup mofem_bc

      *
      * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET,it) {
      * 	...
      * }
      *
      * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
      */
    inline CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator
    getBySetTypeBegin(const unsigned int cubit_bc_type) const {
      return cubitMeshsets.get<CubitMeshSets_mask_meshset_mi_tag>().lower_bound(cubit_bc_type);
    }

    /**
      * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
      * \ingroup mofem_bc

      *
      * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET,it) {
      * 	...
      * }
      *
      * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
      */
    inline CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator
    getBySetTypeEnd(const unsigned int cubit_bc_type) const {
      return cubitMeshsets.get<CubitMeshSets_mask_meshset_mi_tag>().upper_bound(cubit_bc_type);
    }

    /**
      * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
      * \ingroup mofem_bc

      *
      * for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
      * 	...
      * }
      *
      * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
      */
    inline CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator
    getBegin(const std::string& name) const {
      return cubitMeshsets.get<CubitMeshSets_name>().lower_bound(name);
    }

    /**
      * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
      * \ingroup mofem_bc

      *
      * for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
      * 	...
      * }
      *
      * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
      */
    inline CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator
    getEnd(const std::string& name) const {
      return cubitMeshsets.get<CubitMeshSets_name>().upper_bound(name);
    }

    /**
      * \brief check for CUBIT Id and CUBIT type
      * \ingroup mofem_bc

      \todo All cubit interface functions should be outsourced to dedicated interface

      * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
      * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
      */
    bool checkMeshset(const int ms_id,const CubitBCType cubit_bc_type);

    /**
      * \brief add cubit meshset
      * \ingroup mofem_bc

      *
      * \param see CubitBC (NODESET, SIDESET or BLOCKSET and more)
      * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET
      * \param name of set

      */
    PetscErrorCode addMeshset(const CubitBCType cubit_bc_tyep,const int ms_id,const std::string name = "");

    /**
     * \brief add entities to cubit meshset
     * @param  cubit_bc_tyep type of meshset, f.e. NODESET, SIDESET or BLOCKSET
     * @param  ms_id         id of meshset
     * @param  ents          entities to add
     * @return               error code
     */
    PetscErrorCode addEntitiesToMeshset(const CubitBCType cubit_bc_tyep,const int ms_id,Range &ents);

    /**
     * \brief set attributes to cubit meshset
     * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET, SIDESET
     * @param  ms_id         id of meshset
     * @param  attributes    attributes
     * @return               error code
     */
    PetscErrorCode setAttribites(
      const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name = ""
    );

    /**
     * \brief set (material) data structure to cubit meshset
     * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET, SIDESET
     * @param  ms_id         id of meshset
     * @param  attributes    attributes
     * @return               error code
     */
    PetscErrorCode setAttribitesByDataStructure(
      const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name = ""
    );

    /**
     * \brief set boundary data structure to meshset
     * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET, NODESET, SIDESET
     * @param  ms_id         id of meshset
     * @param  data          data structure
     * @return               error code
     */
    PetscErrorCode setBcData(
      const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
    );

    /**
      * \brief delete cubit meshset
      * \ingroup mopfem_bc

      *
      * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
      * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
      *
      */
    PetscErrorCode deleteMeshset(const CubitBCType cubit_bc_type,const int ms_id);

    /**
      * \brief get cubit meshset
      * \ingroup mofem_bc

      */
    PetscErrorCode getCubitMeshsetPtr(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr);

    /**
      * \brief get entities from CUBIT/meshset of a particular entity dimension
      * \ingroup mofem_bc

      * Nodeset can contain nodes, edges, triangles and tets. This applies to other  meshsets too.
      * The nodeset's meshset contain the nodes in the MIDDLE of the surface or volume which is done by default in Cubit,
      * Hence if all nodes on a particular nodeset are required,
      * one should get all triangles or tetrahedrons for which the nodeset was create in Cubit,
      * and get all the connectivities of tris/tets.

      * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
      * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
      * \param dimensions (0 - Nodes, 1 - Edges, 2 - Faces, 3 - Volume(tetrahedral))
      * \param Range containing the retreived entities
      * \param recursive If true, meshsets containing meshsets are queried recursively. Returns the contents of meshsets, but not the meshsets themselves if true.
      */
    PetscErrorCode getEntitiesByDimension(
      const int ms_id,const unsigned int cubit_bc_type, const int dimension,Range &entities,const bool recursive = false
    );

    /**
      * \brief get entities related to CUBIT/meshset,
      * \ingroup mofem_bc

      * NODESET will get Vertices only, even if the NODESET contains edges, tris and tets
      * SIDESET will get Tris, BLOCKSET will get Tets, DISPLACEMENTSET and FORCESET are stored in NODESET, PRESSURESET is stored in Sideset.

      * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
      * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
      * \param Range containing the retreived entities related to the
      * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
      */
    PetscErrorCode getEntitiesByDimension(
      const int ms_id,const unsigned int cubit_bc_type, Range &entities,const bool recursive = false
    );

    /**
      * \ingroup mofem_bc
      * \brief get meshset from CUBIT Id and CUBIT type
      *
      * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
      * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
      * \param meshset where to store the retrieved entities
      */
    PetscErrorCode getMeshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset);

    /**
      * \ingroup mofem_bc
      * \brief get all CUBIT meshsets by CUBIT type
      *
      * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more).
      * \param meshsets is range of meshsets
      */
    PetscErrorCode getMeshsetsByType(const unsigned int cubit_bc_type,Range &meshsets);

  protected:

    Tag nsTag;
    Tag ssTag;
    Tag nsTag_data;
    Tag ssTag_data;
    Tag bhTag;
    Tag bhTag_header;

    //cubit
    CubitMeshSet_multiIndex cubitMeshsets;	   ///< cubit meshsets

  };

}

#endif //__MESHSETSMANAGER_HPP__

/***************************************************************************//**
 * \defgroup mofem_bc Managing meshsets which boundary conditions and materials
 * \ingroup mofem
 ******************************************************************************/
