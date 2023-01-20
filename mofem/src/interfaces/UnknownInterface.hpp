/** \file UnknownInterface.hpp
 * \brief MoFEM interface
 *
 * Low level data structures not used directly by user
 */

#ifndef __MOFEMUNKNOWNINTERFACE_HPP__
#define __MOFEMUNKNOWNINTERFACE_HPP__

namespace MoFEM {

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

  std::string strVersion() {
    auto str = [](auto v) { return boost::lexical_cast<std::string>(v); };
    return str(majorVersion) + "." + str(minorVersion) + "." +
           str(buildVersion);
  }
};

/** \brief base class for all interface classes
 * \ingroup mofem
 **/
struct UnknownInterface {

  virtual MoFEMErrorCode
  query_interface(boost::typeindex::type_index type_index,
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
  MoFEMErrorCode registerInterface(bool error_if_registration_failed = true) {
    MoFEMFunctionBeginHot;
    auto p =
        iFaceTypeMap.insert(UIdTypeMap(boost::typeindex::type_id<IFACE>()));
    if (error_if_registration_failed && (!p.second)) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
               "Registration of interface typeid(IFACE).name() = %s failed",
               typeid(IFACE).name());
    }
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
   * CHKERR m_field.getInterface(simple_interface);
   *
   * \endcode
   *
   * @param iface reference to a interface pointer
   * @return MoFEMErrorCode
   */
  template <class IFACE>
  inline MoFEMErrorCode getInterface(IFACE *&iface) const {
    MoFEMFunctionBegin;
    UnknownInterface *ptr;
    CHKERR query_interface(boost::typeindex::type_id<IFACE>(), &ptr);
    if (!(iface = static_cast<IFACE *>(ptr)))
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "Cast Impossible");
    MoFEMFunctionReturn(0);
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
    return getInterface<IFACE>(*iface);
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
   * Simple *simple_interface = m_field.getInterface<Simple*>();
   *
   * \endcode
   *
   * @return IFACE*
   */
  template <class IFACE,
            typename boost::enable_if<boost::is_pointer<IFACE>, int>::type = 0>
  inline IFACE getInterface() const {
    typedef typename boost::remove_pointer<IFACE>::type IFaceType;
    IFaceType *iface = NULL;
    ierr = getInterface<IFACE>(iface);
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
   * Simple &simple_interface = m_field.getInterface<Simple&>();
   *
   * \endcode
   *
   * @return IFACE&
   */
  template <class IFACE, typename boost::enable_if<boost::is_reference<IFACE>,
                                                   int>::type = 0>
  inline IFACE getInterface() const {
    typedef typename boost::remove_reference<IFACE>::type IFaceType;
    IFaceType *iface = NULL;
    ierr = getInterface<IFaceType>(iface);
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
   * Simple *simple_interface = m_field.getInterface<Simple>();
   *
   * \endcode
   *
   * @return IFACE*
   */
  template <class IFACE> inline IFACE *getInterface() const {
    IFACE *iface = NULL;
    ierr = getInterface<IFACE>(iface);
    CHKERRABORT(PETSC_COMM_SELF, ierr);
    return iface;
  }

  virtual ~UnknownInterface() = default;

  /**
   * \brief Get library version
   *
   * This is library version.
   *
   * @return error code
   */
  static MoFEMErrorCode getLibVersion(Version &version);

  /**
   * \brief Get database major version
   *
   * This is database version. MoFEM can read DataBase from file created by
   * older version. Then library version and database version could be
   * different.
   *
   * @return error code
   */
  static MoFEMErrorCode getFileVersion(moab::Interface &moab, Version &version);

  /**
   * \brief Get database major version
   *
   * This is database version. MoFEM can read DataBase from file created by
   * older version. Then library version and database version could be
   * different.
   *
   * @return error code
   */
  static MoFEMErrorCode setFileVersion(
      moab::Interface &moab,
      Version version = Version(MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR,
                                MoFEM_VERSION_BUILD));

  /**
   * \brief Get database major version
   *
   * Implementation of particular interface could be different than main lib.
   * For example user could use older interface, to keep back compatibility.
   *
   * @return error code
   */
  static MoFEMErrorCode getInterfaceVersion(Version &version);

protected:
  struct NotKnownClass {};

private:
  struct UIdTypeMap {
    boost::typeindex::type_index classIdx;
    UIdTypeMap(const boost::typeindex::type_index &idx) : classIdx(idx) {}
  };

  /// Data structure for interfaces Id and class names
  typedef multi_index_container<
      UIdTypeMap,
      indexed_by<

          hashed_unique<member<UIdTypeMap, boost::typeindex::type_index,
                               &UIdTypeMap::classIdx>>

          >>
      iFaceTypeMap_multiIndex;

  mutable iFaceTypeMap_multiIndex
      iFaceTypeMap; ///< Maps implementation to interface type name
};

} // namespace MoFEM

#endif // __MOFEMUNKNOWNINTERFACE_HPP__
