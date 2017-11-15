/** \file UnknownInterface.hpp
 * \brief MoFEM interface
 *
 * Low level data structures not used directly by user
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

#ifndef __MOFEMUNKNOWNINTERFACE_HPP__
#define __MOFEMUNKNOWNINTERFACE_HPP__

namespace MoFEM {

/**
* \brief MoFEM interface unique ID
* \ingroup mofem
*/
struct MOFEMuuid {

  MOFEMuuid() { memset(this, 0, sizeof(MOFEMuuid)); }
  MOFEMuuid(const BitIntefaceId &uuid) { uUId = uuid; }

  /** \brief returns whether two uuid's are equal
  **/
  inline bool operator==(const MOFEMuuid &orig) const {
    return (uUId & orig.uUId) == orig.uUId;
  }

  // uuid
  BitIntefaceId uUId;
};

/** uuid for an unknown interface
* this can be used to either return a default interface
* or a NULL interface
**/
static const MOFEMuuid IDD_MOFEMUnknown =
    MOFEMuuid(BitIntefaceId(UNKNOWNINTERFACE));

struct Version {
  int majorVersion;
  int minorVersion;
  int buildVersion;
  Version()
      : majorVersion(MoFEM_VERSION_MAJOR), minorVersion(MoFEM_VERSION_MINOR),
        buildVersion(MoFEM_VERSION_BUILD) {}
  Version(const int v[3])
      : majorVersion(v[0]), minorVersion(v[1]), buildVersion(v[2]) {}
  Version(const int minor, const int major, const int build)
      : majorVersion(minor), minorVersion(major), buildVersion(build) {}
};

/** \brief base class for all interface classes
* \ingroup mofem
**/
struct UnknownInterface {

  virtual MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                         UnknownInterface **iface) const = 0;

  /**
   * @brief Register interface
   *
   * Example:
   * \code
   * ierr = regSubInterface<Simple>(IDD_MOFEMSimple);
   * CHKERRABORT(PETSC_COMM_SELF, ierr);
   * \endcode
   *
   * @param uuid
   * @param true
   * @return MoFEMErrorCode
   */
  template <class IFACE>
  MoFEMErrorCode registerInterface(const MOFEMuuid &uuid,
                                   bool error_if_registration_failed = true) {
    MoFEMFunctionBeginHot;
    std::pair<iFaceTypeMap_multiIndex::iterator, bool> p;
    p = iFaceTypeMap.insert(
        UIdTypeMap(uuid, boost::typeindex::type_id<IFACE>()));
    if (error_if_registration_failed && (!p.second)) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
               "Registration of interface typeid(IFACE).name() = %s failed",
               typeid(IFACE).name());
    }
    MoFEMFunctionReturnHot(0);
  }

  /**
   * @brief Get interface by uuid and return reference to pointer of interface
   *
   * \note uuid of interface and interface are verified, if second template
   * parameter is true. Function throw error if both do not match.
   *
   * \note Do not use this function directly, it is called by other overload
   * getInterface methods.
   *
   * @param uuid
   * @param iface reference to a interface pointer
   * @return MoFEMErrorCode
   */
  template <class IFACE, bool VERIFY /* =false C++11 needed to have this */ >
  inline MoFEMErrorCode getInterface(const MOFEMuuid &uuid,
                                     IFACE *&iface) const {
    MoFEMFunctionBeginHot;
    if (VERIFY) {
      if (boost::typeindex::type_id<IFACE>()!=getClassIdx(uuid)) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
                 "Inconsistency between interface Id and type");
      }
    }
    UnknownInterface *ptr;
    ierr = getInterface<UnknownInterface, false>(uuid, ptr);
    CHKERRG(ierr);
    iface = static_cast<IFACE *>(ptr);
    MoFEMFunctionReturnHot(0);
  }

  /**
   * @brief Get interface refernce to pointer of interface
   *
   * \code
   * // Create moab database
   * moab::Core mb_instance;
   * // Access moab database by interface
   * moab::Interface &moab = mb_instance;
   *
   * // Create MoFEM database
   * MoFEM::Core core(moab);
   * // Acces MoFEM database by Interface
   * MoFEM::Interface &m_field = core;
   *
   * // Get interface
   * // Get simple interface is simplified version enabling quick and
   * // easy construction of problem.
   * Simple *simple_interface;
   * // Query interface and get pointer to Simple interface
   *  CHKERR m_field.getInterface(simple_interface);
   *
   * \endcode
   *
   * @param iface reference to a interface pointer
   * @return MoFEMErrorCode
   */
  template <class IFACE>
  inline MoFEMErrorCode getInterface(IFACE *&iface) const {
    return getInterface<IFACE, false>(
        getUId(boost::typeindex::type_id<IFACE>()), iface);
  }

  /**
   * @brief Get interface pointer to pointer of interface
   *
   * \code
   * // Create moab database
   * moab::Core mb_instance;
   * // Access moab database by interface
   * moab::Interface &moab = mb_instance;
   *
   * // Create MoFEM database
   * MoFEM::Core core(moab);
   * // Acces MoFEM database by Interface
   * MoFEM::Interface &m_field = core;
   *
   * // Get interface
   * // Get simple interface is simplified version enabling quick and
   * // easy construction of problem.
   * Simple *simple_interface;
   * // Query interface and get pointer to Simple interface
   *  CHKERR m_field.getInterface(&simple_interface);
   *
   * \endcode
   *
   *
   * @param iface const pointer to a interface pointer
   * @return MoFEMErrorCode
   */
  template <class IFACE>
  inline MoFEMErrorCode getInterface(IFACE **const iface) const {
    return getInterface<IFACE, false>(boost::typeindex::type_id<IFACE>(),
                                      *iface);
  }

  /**
   * @brief Get interface pointer to pointer of interface
   *
   * \code
   * // Create moab database
   * moab::Core mb_instance;
   * // Access moab database by interface
   * moab::Interface &moab = mb_instance;
   *
   * // Create MoFEM database
   * MoFEM::Core core(moab);
   * // Acces MoFEM database by Interface
   * MoFEM::Interface &m_field = core;
   *
   * // Get interface
   * // Get simple interface is simplified version enabling quick and
   * // easy construction of problem.
   * Simple *simple_interface = m_field.getInterface<Simple*,0>();
   *
   * \endcode
   *
   * @return IFACE*
   */
  template <
    class IFACE, 
    typename boost::enable_if<boost::is_pointer<IFACE>, int>::type /* =0 C++11 needed to have this */  
  >
  inline IFACE getInterface() const {
    typedef typename boost::remove_pointer<IFACE>::type IFaceType;
    IFaceType* iface = NULL;
    ierr = getInterface<IFaceType, false>(
        getUId(boost::typeindex::type_id<IFaceType>()), iface);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return iface;
  }

 /**
   * @brief Get reference to interface 
   *
   * \code
   * // Create moab database
   * moab::Core mb_instance;
   * // Access moab database by interface
   * moab::Interface &moab = mb_instance;
   *
   * // Create MoFEM database
   * MoFEM::Core core(moab);
   * // Acces MoFEM database by Interface
   * MoFEM::Interface &m_field = core;
   *
   * // Get interface
   * // Get simple interface is simplified version enabling quick and
   * // easy construction of problem.
   * Simple &simple_interface = m_field.getInterface<Simple&,0>();
   *
   * \endcode
   *
   * @return IFACE&
   */
  template <
    class IFACE, 
    typename boost::enable_if<boost::is_reference<IFACE>, int>::type /* =0 C++11 needed to have this */  
  >
  inline IFACE getInterface() const {
    typedef typename boost::remove_reference<IFACE>::type IFaceType;
    IFaceType* iface = NULL;
    ierr = getInterface<IFaceType, false>(
        getUId(boost::typeindex::type_id<IFaceType>()), iface);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return *iface;
  }

  /**
   * @brief Function returning pointer to interface
   *
   * \code
   * // Create moab database
   * moab::Core mb_instance;
   * // Access moab database by interface
   * moab::Interface &moab = mb_instance;
   *
   * // Create MoFEM database
   * MoFEM::Core core(moab);
   * // Acces MoFEM database by Interface
   * MoFEM::Interface &m_field = core;
   *
   * // Get interface
   * // Get simple interface is simplified version enabling quick and
   * // easy construction of problem.
   * Simple *simple_interface = m_field.getInterface<Simple,0>();
   *
   * \endcode
   *
   * @return IFACE*
   */
  template <class IFACE>
  inline IFACE* getInterface() const {
    return getInterface<IFACE*,0>();
  }

  virtual ~UnknownInterface() {}

  /**
  * \brief Get library version
  *
  * This is library version.
  *
  * @return error code
  */
  virtual MoFEMErrorCode getLibVersion(Version &version) const {
    MoFEMFunctionBeginHot;
    version =
        Version(MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR, MoFEM_VERSION_BUILD);
    MoFEMFunctionReturnHot(0);
  }

  /**
  * \brief Get database major version
  *
  * This is database version. MoFEM can read DataBase from file created by older
  * version. Then library version and database version could be different.
  *
  * @return error code
  */
  virtual const MoFEMErrorCode getFileVersion(moab::Interface &moab,
                                              Version &version) const {
    MoFEMFunctionBegin;
    const EntityHandle root_meshset = 0;
    const int def_version[] = {-1, -1, -1};
    Tag th;
    rval = moab.tag_get_handle("MOFEM_VERSION", 3, MB_TYPE_INTEGER, th,
                               MB_TAG_CREAT | MB_TAG_MESH, &def_version);
    int *version_ptr;
    if (rval == MB_ALREADY_ALLOCATED) {
      const void *tag_data[1];
      CHKERR moab.tag_get_by_ptr(th, &root_meshset, 1, tag_data);
      version_ptr = (int *)tag_data[0];
    } else {
      const void *tag_data[1];
      CHKERR moab.tag_get_by_ptr(th, &root_meshset, 1, tag_data);
      version_ptr = (int *)tag_data[0];
      version_ptr[0] = MoFEM_VERSION_MAJOR;
      version_ptr[1] = MoFEM_VERSION_MINOR;
      version_ptr[2] = MoFEM_VERSION_BUILD;
    }
    version = Version(version_ptr);
    MoFEMFunctionReturn(0);
  }

  /**
  * \brief Get database major version
  *
  * Implementation of particular interface could be different than main lib. For
  * example user could use older interface, to keep back compatibility.
  *
  * @return error code
  */
  virtual MoFEMErrorCode getInterfaceVersion(Version &version) const {
    MoFEMFunctionBeginHot;
    version =
        Version(MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR, MoFEM_VERSION_BUILD);
    MoFEMFunctionReturnHot(0);
  }

protected:

  struct NotKnownClass {};

  /**
   * \brief Get type name for interface Id
   * @param  uid interface Id
   * @return     class name
   */
  inline boost::typeindex::type_index getClassIdx(const MOFEMuuid &uid) const {
    iFaceTypeMap_multiIndex::nth_index<0>::type::iterator it;
    it = iFaceTypeMap.get<0>().find(uid);
    if (it != iFaceTypeMap.get<0>().end()) {
      return it->classIdx;
    }
    return boost::typeindex::type_id<NotKnownClass>();
  }

  /**
   * \brief Get interface Id for class name
   * @param  class_name
   * @return            Id
   */
  inline MOFEMuuid getUId(const boost::typeindex::type_index &class_idx) const {
    iFaceTypeMap_multiIndex::nth_index<1>::type::iterator it;
    it = iFaceTypeMap.get<1>().find(class_idx);
    if (it != iFaceTypeMap.get<1>().end()) {
      return it->uID;
    }
    return IDD_MOFEMUnknown;
  }

private:

  struct UIdTypeMap {
    MOFEMuuid uID;
    boost::typeindex::type_index classIdx;
    UIdTypeMap(const MOFEMuuid &uid, const boost::typeindex::type_index &idx)
        : uID(uid), classIdx(idx) {}
  };

  struct HashMOFEMuuid {
    inline unsigned int operator()(const MOFEMuuid &value) const {
      return value.uUId.to_ulong();
    }
  };

  /// Data structure for interfaces Id and class names
  typedef multi_index_container<
    UIdTypeMap,
    indexed_by<
     hashed_unique<
      member<UIdTypeMap, MOFEMuuid, &UIdTypeMap::uID>, 
      HashMOFEMuuid
     >,
     hashed_unique<
      member<UIdTypeMap, boost::typeindex::type_index, &UIdTypeMap::classIdx>
     >
    >
  > iFaceTypeMap_multiIndex;

  mutable iFaceTypeMap_multiIndex
      iFaceTypeMap; ///< Maps MOFEMuuid to interface type name
};

template <>
inline MoFEMErrorCode UnknownInterface::getInterface<UnknownInterface, false>(
    const MOFEMuuid &uuid, UnknownInterface *&iface) const {
  return query_interface(uuid, &iface);
}
}

#endif // __MOFEMUNKNOWNINTERFACE_HPP__
