/** \file MeshsetsManager.hpp
 * \brief MeshsetsManager interface

 Interface to manage material and boundary sets

 * \ingroup mofem_meshset_mng
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

using Sev = MoFEM::LogManager::SeverityLevel;

typedef CubitMeshSet_multiIndex::index<CubitMeshsetType_mi_tag>::type
    CubitMeshsetByType;

typedef CubitMeshSet_multiIndex::index<CubitMeshsetMaskedType_mi_tag>::type
    CubitMeshsetByMask;

typedef CubitMeshSet_multiIndex::index<CubitMeshsets_name>::type
    CubitMeshsetByName;

typedef CubitMeshSet_multiIndex::index<CubitMeshsetType_mi_tag>::type
    CubitMeshsetById;

/**
 * \brief Iterator that loops over all the Cubit MeshSets in a moFEM field
 * \ingroup mofem_meshset_mng

 *
 * \param MESHSET_MANAGER meshset manager (works as well with Interface)
 * \param iterator
 */
#define _IT_CUBITMESHSETS_FOR_LOOP_(MESHSET_MANAGER, IT)                       \
  CubitMeshSet_multiIndex::iterator IT =                                       \
      MESHSET_MANAGER.get_meshsets_manager_ptr()->getBegin();                  \
  IT != MESHSET_MANAGER.get_meshsets_manager_ptr()->getEnd();                  \
  IT++

/**
* \brief Iterator that loops over a specific Cubit MeshSet in a moFEM field
* \ingroup mofem_meshset_mng

*
* \param mField moFEM Field
* \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
* \param iterator
*/
#define _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MESHSET_MANAGER,            \
                                                   CUBITBCTYPE, IT)            \
  CubitMeshsetByType::iterator IT =                                            \
      MESHSET_MANAGER.get_meshsets_manager_ptr()->getBegin(CUBITBCTYPE);       \
  IT != MESHSET_MANAGER.get_meshsets_manager_ptr()->getEnd(CUBITBCTYPE);       \
  IT++

/**
* \brief Iterator that loops over a specific Cubit MeshSet having a particular
BC meshset in a moFEM field
* \ingroup mofem_meshset_mng

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
#define _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MESHSET_MANAGER, CUBITBCTYPE,  \
                                                IT)                            \
  CubitMeshsetByMask::iterator IT =                                            \
      MESHSET_MANAGER.get_meshsets_manager_ptr()->getBySetTypeBegin(           \
          CUBITBCTYPE);                                                        \
  IT != MESHSET_MANAGER.get_meshsets_manager_ptr()->getBySetTypeEnd(           \
            CUBITBCTYPE);                                                      \
  IT++

/**
* \brief Iterator that loops over Cubit BlockSet having a particular name
* \ingroup mofem_meshset_mng


* \param MESHSET_MANAGER meshset manager (works as well with Interface)
* \param NAME name
* \param IT iterator
*
* Example: \code
for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,"SOME_BLOCK_NAME",it) {
...
* } \endcode
*/
#define _IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(MESHSET_MANAGER, NAME, IT)         \
  CubitMeshsetByName::iterator IT =                                            \
      MESHSET_MANAGER.get_meshsets_manager_ptr()->getBegin(NAME);              \
  IT != MESHSET_MANAGER.get_meshsets_manager_ptr()->getEnd(NAME);              \
  IT++

/** \brief Interface for managing meshsets containing materials and boundary
 * conditions
 * \ingroup mofem_meshset_mng
 */
struct MeshsetsManager : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  MeshsetsManager(const MoFEM::Core &core);
  virtual ~MeshsetsManager() = default;

  /**
   * \brief get tags handlers used on meshsets

   * On meshsets range of tags in set. Depending on tag type and data on that
   * tag type of meshset could be determined. This function get hanldes to
   * tags.
   *
   * Most of the tags are followinf convention used by MoAB or Cubit and other
   * meshing softwares, f.e. gmesh.

   */
  MoFEMErrorCode getTags(int verb = -1);

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
  MeshsetsManager *get_meshsets_manager_ptr() { return this; }

  /**
   * \brief return pointer to meshset manager
   */
  const MeshsetsManager *get_meshsets_manager_ptr() const { return this; }

  /**
   * \brief clear multi-index container
   * @return error code
   */
  MoFEMErrorCode clearMap();

  /**
   * \brier initialize container form data on mesh
   * @return [description]
   */
  MoFEMErrorCode initialiseDatabaseFromMesh(int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Boradcats meshsets
   *
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode readMeshsets(int verb = DEFAULT_VERBOSITY);

  /**
   * @brief Boradcats meshsets
   *
   * @param verb
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode broadcastMeshsets(int verb = DEFAULT_VERBOSITY);

  template <class CUBIT_BC_DATA_TYPE>
  MoFEMErrorCode printBcSet(CUBIT_BC_DATA_TYPE &data,
                            unsigned long int type) const;

  /**
   * \brief print meshsets with displacement boundary conditions data
   * structure
   */
  MoFEMErrorCode printDisplacementSet() const;

  /**
   * \brief print meshsets with pressure boundary conditions data structure
   */
  MoFEMErrorCode printPressureSet() const;

  /**
   * \brief print meshsets with force boundary conditions data structure
   */
  MoFEMErrorCode printForceSet() const;

  /**
   * \brief print meshsets with temperature boundary conditions data structure
   */
  MoFEMErrorCode printTemperatureSet() const;

  /**
   * \brief print meshsets with heat flux boundary conditions data structure
   */
  MoFEMErrorCode printHeatFluxSet() const;

  /**
   * \brief print meshsets with material data structure set on it
   */
  MoFEMErrorCode printMaterialsSet() const;

  inline CubitMeshSet_multiIndex &getMeshsetsMultindex() {
    return cubitMeshsets;
  }

  /**
   * \ingroup mofem_meshset_mng
   * \brief get begin iterator of cubit mehset of given type (instead you can
   * use _IT_CUBITMESHSETS_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
   *
   * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
   * 	...
   * }
   *
   */
  inline CubitMeshSet_multiIndex::iterator getBegin() const {
    return cubitMeshsets.begin();
  }

  /**
   * \ingroup mofem_meshset_mng
   * \brief get begin iterator of cubit mehset of given type (instead you can
   * use _IT_CUBITMESHSETS_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
   *
   * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
   * 	...
   * }
   *
   */
  CubitMeshSet_multiIndex::iterator getEnd() const {
    return cubitMeshsets.end();
  }

  /**
    * \brief get begin iterator of cubit mehset of given type (instead you can
    use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    * \ingroup mofem_meshset_mng

    *
    *
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)
    {
    * 	...
    * }
    *
    * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  inline CubitMeshsetByType::iterator
  getBegin(const unsigned int cubit_bc_type) const {
    return cubitMeshsets.get<CubitMeshsetType_mi_tag>().lower_bound(
        cubit_bc_type);
  }

  /**
    * \brief get begin iterator of cubit mehset of given type (instead you can
    use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    * \ingroup mofem_meshset_mng

    *
    *
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)
    {
    * 	...
    * }
    *
    * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  inline CubitMeshsetByType::iterator
  getEnd(const unsigned int cubit_bc_type) const {
    return cubitMeshsets.get<CubitMeshsetType_mi_tag>().upper_bound(
        cubit_bc_type);
  }

  /**
    * \brief get end iterator of cubit meshset of given type (instead you can
    use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    * \ingroup mofem_meshset_mng

    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET,it) {
    * 	...
    * }
    *
    * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  inline CubitMeshsetByMask::iterator
  getBySetTypeBegin(const unsigned int cubit_bc_type) const {
    return cubitMeshsets.get<CubitMeshsetMaskedType_mi_tag>().lower_bound(
        cubit_bc_type);
  }

  /**
    * \brief get end iterator of cubit mehset of given type (instead you can
    use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    * \ingroup mofem_meshset_mng

    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET,it) {
    * 	...
    * }
    *
    * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  inline CubitMeshsetByMask::iterator
  getBySetTypeEnd(const unsigned int cubit_bc_type) const {
    return cubitMeshsets.get<CubitMeshsetMaskedType_mi_tag>().upper_bound(
        cubit_bc_type);
  }

  /**
    * \brief get begin iterator of cubit mehset of given type (instead you can
    use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    * \ingroup mofem_meshset_mng

    *
    *
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)
    {
    * 	...
    * }
    *
    * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  inline CubitMeshsetByName::iterator getBegin(const std::string &name) const {
    return cubitMeshsets.get<CubitMeshsets_name>().lower_bound(name);
  }

  /**
    * \brief get begin iterator of cubit mehset of given type (instead you can
    use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    * \ingroup mofem_meshset_mng

    *
    *
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)
    {
    * 	...
    * }
    *
    * \param  type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  inline CubitMeshsetByName::iterator getEnd(const std::string &name) const {
    return cubitMeshsets.get<CubitMeshsets_name>().upper_bound(name);
  }

  /**
    * \brief check for CUBIT Id and CUBIT type
    * \ingroup mofem_meshset_mng

    \todo All cubit interface functions should be outsourced to dedicated
    interface

    * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    */
  bool checkMeshset(const int ms_id, const CubitBCType cubit_bc_type) const;

  /**
   * \brief check if meshset of given name exist
   * @param  name name of meshset
   * @return      error code
   */
  bool checkMeshset(const string name,
                    int *const number_of_meshsets_ptr = NULL) const;

  /**
    * \brief add cubit meshset
    * \ingroup mofem_meshset_mng

    *
    * \param see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET
    * \param name of set

    */
  MoFEMErrorCode addMeshset(const CubitBCType cubit_bc_type, const int ms_id,
                            const std::string name = "");

  /**
   * \brief add entities to cubit meshset
   * @param  cubit_bc_type type of meshset, f.e. NODESET, SIDESET or BLOCKSET
   * @param  ms_id         id of meshset
   * @param  ents          entities to add
   * @return               error code
   */
  MoFEMErrorCode addEntitiesToMeshset(const CubitBCType cubit_bc_type,
                                      const int ms_id, const Range &ents);

  /**
   * \brief add entities to cubit meshset
   * @param  cubit_bc_type type of meshset, f.e. NODESET, SIDESET or BLOCKSET
   * @param  ms_id         id of meshset
   * @param  ents          pointer to entities array
   * @param  nb_ents       number of entities in array
   * @return               error code
   */
  MoFEMErrorCode addEntitiesToMeshset(const CubitBCType cubit_bc_type,
                                      const int ms_id, const EntityHandle *ents,
                                      const int nb_ents);

  /**
   * \brief set attributes to cubit meshset
   * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET,
   * NODESET, SIDESET
   * @param  ms_id         id of meshset
   * @param  attributes    attributes
   * @param  name          set name to blockset
   * @return               error code
   */
  MoFEMErrorCode setAtributes(const CubitBCType cubit_bc_type, const int ms_id,
                              const std::vector<double> &attributes,
                              const std::string name = "");

  /**
   * \brief set (material) data structure to cubit meshset
   * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET,
   * NODESET, SIDESET
   * @param  ms_id         id of meshset
   * @param  attributes    attributes
   * @return               error code
   */
  MoFEMErrorCode setAtributesByDataStructure(const CubitBCType cubit_bc_type,
                                             const int ms_id,
                                             const GenericAttributeData &data,
                                             const std::string name = "");

  /**
   * \brief set boundary data structure to meshset
   * @param  cubit_bc_type type of meshset, see CubitBC, i.e. BLOCKSET,
   * NODESET, SIDESET
   * @param  ms_id         id of meshset
   * @param  data          data structure
   * @return               error code
   */
  MoFEMErrorCode setBcData(const CubitBCType cubit_bc_type, const int ms_id,
                           const GenericCubitBcData &data);

  /**
    * \brief delete cubit meshset
    * \ingroup mopfem_bc

    *
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    *
    */
  MoFEMErrorCode deleteMeshset(const CubitBCType cubit_bc_type, const int ms_id,
                               const MoFEMTypes bh = MF_EXIST);

  /**
   * \brief get cubit meshset
   * \ingroup mofem_meshset_mng
   *
   *
   */
  MoFEMErrorCode
  getCubitMeshsetPtr(const int ms_id, const CubitBCType cubit_bc_type,
                     const CubitMeshSets **cubit_meshset_ptr) const;

  /**
   * \brief get cubit meshset
   *
   * \ingroup mofem_meshset_mng
   */
  MoFEMErrorCode
  getCubitMeshsetPtr(const string name,
                     const CubitMeshSets **cubit_meshset_ptr) const;

  /**
   * @brief Get vector of poointer to blocksets with name satisfying regular
   * expression.
   *
   * \ingroup mofem_meshset_mng
   *
   * @param reg_exp_name
   * @param std::vector<const CubitMeshSets *>
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  getCubitMeshsetPtr(const std::regex reg_exp_name,
                     std::vector<const CubitMeshSets *> &vec_ptr) const;

  /**
    * \brief get entities from CUBIT/meshset of a particular entity dimension
    * \ingroup mofem_meshset_mng

    * Nodeset can contain nodes, edges, triangles and tets. This applies to
    other  meshsets too.
    * The nodeset's meshset contain the nodes in the MIDDLE of the surface or
    volume which is done by default in Cubit,
    * Hence if all nodes on a particular nodeset are required,
    * one should get all triangles or tetrahedrons for which the nodeset was
    create in Cubit,
    * and get all the connectivities of tris/tets.

    * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param dimensions (0 - Nodes, 1 - Edges, 2 - Faces, 3 -
    Volume(tetrahedral))
    * \param Range containing the retrieved entities
    * \param recursive If true, meshsets containing meshsets are queried
    recursively. Returns the contents of meshsets, but not the meshsets
    themselves if true.
    */
  MoFEMErrorCode getEntitiesByDimension(const int ms_id,
                                        const unsigned int cubit_bc_type,
                                        const int dimension, Range &entities,
                                        const bool recursive = true) const;

  /**
    * \brief get entities related to CUBIT/meshset,
    * \ingroup mofem_meshset_mng

    * NODESET will get Vertices only, even if the NODESET contains edges, tris
    and tets
    * SIDESET will get Tris, BLOCKSET will get Tets, DISPLACEMENTSET and
    FORCESET are stored in NODESET, PRESSURESET is stored in Sideset.

    * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
    * \param Range containing the retrieved entities related to the
    * \param recursive If true, meshsets containing meshsets are queried
    recursively.  Returns the contents of meshsets, but not the meshsets
    themselves if true.
    */
  MoFEMErrorCode getEntitiesByDimension(const int ms_id,
                                        const unsigned int cubit_bc_type,
                                        Range &entities,
                                        const bool recursive = true) const;

  /**
   * \ingroup mofem_meshset_mng
   * \brief get meshset from CUBIT Id and CUBIT type
   *
   * \param ms_id id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
   * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more)
   * \param meshset where to store the retrieved entities
   */
  MoFEMErrorCode getMeshset(const int ms_id, const unsigned int cubit_bc_type,
                            EntityHandle &meshset) const;

  /**
   * @brief Check if meshset constains entities
   *
   * @param ms_id
   * @param cubit_bc_type
   * @param entities
   * @param num_entities
   * @param operation_type
   * @return true
   * @return false
   */
  bool checkIfMeshsetContainsEntities(
      const int ms_id, const unsigned int cubit_bc_type,
      const EntityHandle *entities, int num_entities,
      const int operation_type = moab::Interface::INTERSECT);

  /**
   * \ingroup mofem_meshset_mng
   * \brief get all CUBIT meshsets by CUBIT type
   *
   * \param  see CubitBC (NODESET, SIDESET or BLOCKSET and more).
   * \param meshsets is range of meshsets
   */
  MoFEMErrorCode getMeshsetsByType(const unsigned int cubit_bc_type,
                                   Range &meshsets) const;

  /**
   * \brief add blocksets reading config file

   Example of config file
   \code

   [block_1001]

   # Example applying attributes to blockset

   id=2001
   add=BLOCKSET
   user1=1.0  # attribute value 1
   user2=2.0  # you can set up to 10 attributes (if needed could be easily
   extended to more, let us know)
   user3=3.0

   [block_1002]

   # Example applying material block (isotropic elastic material)

   id=2002
   add=BLOCKSET
   name=MAT_ELASTIC
   young=10
   poisson=0.25
   thermalexpansion=0

   [block_1003]

   # Example applying displacement constrains

   id=2003
   add=NODESET

   # Each flag means that boundary consition on displacements is set.
   disp_flag1=1 # Setting constrains in x- direction
   disp_flag2=1 # Setting constrains in y- direction
   disp_flag3=1 # Setting constrains in z- direction
   disp_flag4=1 # Setting constrains on rotation over x- axis
   disp_flag5=1 # Setting constrains on rotation over y- axis
   disp_flag6=1 # Setting constrains on rotation over z-axis
   disp_ux=1 # value of disp in x- direction
   disp_uy=2
   disp_uz=3
   disp_rx=4 # value of rotation in y-direction
   disp_ry=5
   disp_rz=6

   # Note above values could be interpreted differently if needed.

   [block_1004]

   # Example applying force boundary conditions

   id=2004
   add=NODESET
   force_magnitude=1
   moment_magnitude=1
   force_fx=1
   force_fy=1
   force_fz=1
   moment_mx=1
   moment_my=1
   moment_mz=1

   [block_1005]

   # Example applying pressure boundary conditions

   id=2005
   add=SIDESET
   pressure_flag2=1       # 0: Pressure is interpreted as pure pressure 1:
   pressure is interpreted as total force
   pressure_magnitude=1

   # Example applying temperature boundary conditions

   [block_1006]

   id=2006
   add=NODESET
   temperature_flag1=1        # 0: N/A, 1: temperature value applied
   temperature_t=1

   [block_1007]

   id=2007
   add=SIDESET
   heatflux_flag1=1        # 0: N/A, 1: heat flux applied
   heatflux_magnitude=1

   [block_1008]

   # Example applying material block (isotropic thermal material)

   id=2008
   add=BLOCKSET
   name=MAT_THERMAL # Hast to be set for Thermal Mat
   conductivity=1
   capacity=1

   [block_1009]

   # Example applying interface
   id=2009
   add=SIDESET
   interface_type=1

   [block_1010]

   # Example applying material block for interface element

   id=2010
   add=BLOCKSET
   name=MAT_INTERF
   interface_alpha = 1
   interface_beta = 0
   interface_ft = 1
   interface_Gf = 1


   [block_1009]

   # Example applying material block (isotropic trans iso material)

   id=2011
   add=BLOCKSET
   name=MAT_ELASTIC_TRANS_ISO
   Youngp=1
   Youngz=2
   Poissonp=3
   Poissonpz=4
   Shearzp=5

   [SET_ATTR_foo]

   # Example set atttributes to block name "foo"
   number_of_attributes=3
   user1=1
   user2=2
   user3=3

   \endcode

   * @param  file_name config file
   * @return           error code

   */
  MoFEMErrorCode setMeshsetFromFile(const string file_name,
                                    const bool clean_file_options = true);

  /**
   * \brief get name of config file from line command '-meshsets_config'
   * @return error code

   Option is "-meshsets_config file_name.cfg"

   */
  MoFEMErrorCode setMeshsetFromFile();

  /**
   * @brief save cubit meshset entities on the moab mesh
   *
   * @param ms_id id of the cubit meshset (NODESET  SIDESET  BLOCKSET)
   * @param cubit_bc_type type of a cubit mesheset
   * @param file_name optional name for the file
   * @param file_type optional file type for moab (VTK MOAB)
   * @param options optional parameters for moab writer (PARALLEL=WRITE_PART)
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode
  saveMeshsetToFile(const int ms_id, const unsigned int cubit_bc_type,
                    const std::string file_name = "out_meshset.vtk",
                    const std::string file_type = "VTK",
                    const std::string options = "") const;

  /**
   * @brief save cubit meshset entities on the moab mesh
   *
   * @param ms_id id of the cubit meshset
   * @param cubit_bc_type type of a cubit mesheset (NODESET  SIDESET  BLOCKSET)
   * @param dim dimension of the entities
   * @param file_name optional name for the file
   * @param file_type optional file type for moab (VTK MOAB)
   * @param options optional parameters for moab writer (PARALLEL=WRITE_PART)
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode saveMeshsetToFile(
      const int ms_id, const unsigned int cubit_bc_type, const int dim,
      const std::string file_name = "out_meshset.vtk",
      const bool recursive = false, const std::string file_type = "VTK",
      const std::string options = "") const;

  /**
   * \brief Get config file options, use with care
   * @return error code
   */
  inline boost::shared_ptr<boost::program_options::options_description> &
  getConfigFileOptionsPtr() {
    return configFileOptionsPtr;
  }

  MoFEMErrorCode updateAllMeshsetsByEntitiesChildren(const BitRefLevel &bit);

  static bool brodcastMeshsets; ///< if true meshsets are synchrinised between
                                ///< processors

protected:
  Tag nsTag;
  Tag ssTag;
  Tag nsTag_data;
  Tag ssTag_data;
  Tag bhTag;
  Tag bhTag_header;

  // cubit
  CubitMeshSet_multiIndex cubitMeshsets; ///< cubit meshsets
  boost::shared_ptr<boost::program_options::options_description>
      configFileOptionsPtr; ///< config file options
};

template <class CUBIT_BC_DATA_TYPE>
MoFEMErrorCode MeshsetsManager::printBcSet(CUBIT_BC_DATA_TYPE &data,
                                           unsigned long int type) const {
  MoFEMFunctionBegin;
  const MoFEM::Interface &m_field = cOre;
  const moab::Interface &moab = m_field.get_moab();
  for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_((*this), type, it)) {
    CHKERR it->getBcDataStructure(data);
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << *it;
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << data;
    MOFEM_LOG("MeshsetMngWorld", Sev::inform) << "name " << it->getName();
    for (EntityType t = MBVERTEX; t != MBENTITYSET; ++t) {
      int nb;
      CHKERR moab.get_number_entities_by_type(it->meshset, t, nb, true);
      MOFEM_LOG("MeshsetMngWorld", Sev::inform)
          << "msId " << it->getMeshsetId() << " number of "
          << moab::CN::EntityTypeName(t) << " " << nb;
    }
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM

#endif //__MESHSETSMANAGER_HPP__

/**
 * \defgroup mofem_meshset_mng MeshsetsManager
 * \brief Interface for meshsets with entities with data and boundary conditions
 *
 * \ingroup mofem
 **/
