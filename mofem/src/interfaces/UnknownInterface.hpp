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

    MOFEMuuid() { memset( this, 0, sizeof(MOFEMuuid)); }
    MOFEMuuid(const BitIntefaceId& uuid) { uUId = uuid; }

    /** \brief returns whether two uuid's are equal
    **/
    inline bool operator==(const MOFEMuuid& orig) const {
      return (uUId&orig.uUId) == orig.uUId;
    }

    //uuid
    BitIntefaceId uUId;

  };

  /** uuid for an unknown interface
  * this can be used to either return a default interface
  * or a NULL interface
  **/
  static const MOFEMuuid IDD_MOFEMUnknown = MOFEMuuid( BitIntefaceId(UNKNOWNINTERFACE) );

  struct Version {
    int majorVersion;
    int minorVersion;
    int buildVersion;
    Version():
    majorVersion(MoFEM_VERSION_MAJOR),
    minorVersion(MoFEM_VERSION_MINOR),
    buildVersion(MoFEM_VERSION_BUILD) {
    }
    Version(const int v[3]):
    majorVersion(v[0]),
    minorVersion(v[1]),
    buildVersion(v[2]) {
    }
    Version(const int minor,const int major,const int build):
    majorVersion(minor),
    minorVersion(major),
    buildVersion(build) {
    }
  };

  /** \brief base class for all interface classes
  * \ingroup mofem
  **/
  struct UnknownInterface {

    virtual PetscErrorCode queryInterface(
      const MOFEMuuid& uuid,UnknownInterface** iface
    ) = 0;

    template <class IFACE>
    inline PetscErrorCode getInterface(
      const MOFEMuuid& uuid,IFACE*& iface
    ) {
      MoFEMFunctionBegin;
      if(typeid(IFACE).name()==getClassName(uuid)) {
        UnknownInterface *ptr = static_cast<UnknownInterface*>(iface);
        ierr = getInterface(uuid,ptr); CHKERRQ(ierr);
        iface = static_cast<IFACE*>(ptr);
        MoFEMFunctionReturnHot(0);
      }
      SETERRQ2(
        PETSC_COMM_SELF,MOFEM_INVALID_DATA,
        "Inconsistency between interface Id and type "
        "typeid(IFACE).name() != iFaceTypeMap.at(uuid) "
        "%s != %s",
        typeid(IFACE).name(),getClassName(uuid).c_str()
      );
      MoFEMFunctionReturn(0);
    }

    template <class IFACE>
    inline PetscErrorCode getInterface(IFACE*& iface) {
      return getInterface(getUId(typeid(IFACE).name()),iface);
    }

    template <class IFACE>
    inline PetscErrorCode getInterface(IFACE**const iface) {
      return getInterface(getUId(typeid(IFACE).name()),*iface);
    }

    template <class IFACE>
    inline IFACE* getInterface() {
      IFACE* iface;
      ierr = getInterface(getUId(typeid(IFACE).name()),iface); CHKERRQ(ierr);
      return iface;
    }

    virtual ~UnknownInterface() {}

    /**
    * \brief Get library version
    *
    * This is library version.
    *
    * @return error code
    */
    virtual PetscErrorCode getLibVersion(Version &version) const {
      MoFEMFunctionBeginHot;
      version = Version(MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD);
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
    virtual const PetscErrorCode getFileVersion(moab::Interface &moab,Version &version) const {
      MoFEMFunctionBeginHot;
      const EntityHandle root_meshset = 0;
      const int def_version[] = {-1,-1,-1};
      Tag th;
      rval = moab.tag_get_handle(
        "MOFEM_VERSION",3,MB_TYPE_INTEGER,th,MB_TAG_CREAT|MB_TAG_MESH,&def_version
      );
      int *version_ptr;
      if(rval==MB_ALREADY_ALLOCATED) {
        const void* tag_data[1];
        rval = moab.tag_get_by_ptr(th,&root_meshset,1,tag_data); CHKERRQ_MOAB(rval);
        version_ptr = (int*)tag_data[0];
      } else {
        CHKERRQ_MOAB(rval);
        const void* tag_data[1];
        rval = moab.tag_get_by_ptr(th,&root_meshset,1,tag_data); CHKERRQ_MOAB(rval);
        version_ptr = (int*)tag_data[0];
        version_ptr[0] = MoFEM_VERSION_MAJOR;
        version_ptr[1] = MoFEM_VERSION_MINOR;
        version_ptr[2] = MoFEM_VERSION_BUILD;
      }
      version = Version(version_ptr);
      MoFEMFunctionReturnHot(0);
    }

    /**
    * \brief Get database major version
    *
    * Implementation of particular interface could be different than main lib. For
    * example user could use older interface, to keep back compatibility.
    *
    * @return error code
    */
    virtual PetscErrorCode getInterfaceVersion(Version &version) const {
      MoFEMFunctionBeginHot;
      version = Version(MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD);
      MoFEMFunctionReturnHot(0);
    }

  protected:

    struct HashMap {
      MOFEMuuid uID;
      std::string className;
      HashMap(const MOFEMuuid &uid,const std::string& name):
      uID(uid),
      className(name) {
      }
    };

    struct HashMOFEMuuid {
      inline unsigned int operator()(const MOFEMuuid& value) const {
        return value.uUId.to_ulong();
      }
    };

    typedef multi_index_container<
      HashMap,
      indexed_by<
        hashed_unique< member<HashMap,MOFEMuuid,&HashMap::uID>,HashMOFEMuuid >,
        hashed_unique< member<HashMap,std::string,&HashMap::className> >
      >
    > iFaceTypeMap_multiIndex;

    static iFaceTypeMap_multiIndex iFaceTypeMap;

    static std::string getClassName(const MOFEMuuid& uid) {
      if(iFaceTypeMap.get<0>().find(uid)!=iFaceTypeMap.get<0>().end()) {
        return iFaceTypeMap.get<0>().find(uid)->className;
      }
      return string("No class for UId found");
    }

    static MOFEMuuid getUId(const std::string& class_name) {
      if(iFaceTypeMap.get<1>().find(class_name)!=iFaceTypeMap.get<1>().end()) {
        return iFaceTypeMap.get<1>().find(class_name)->uID;
      }
      return IDD_MOFEMUnknown;
    }


  };

  template<>
  inline PetscErrorCode UnknownInterface::getInterface(
    const MOFEMuuid& uuid,UnknownInterface*& iface
  ) {
    return queryInterface(uuid,&iface);
  }

}

#endif // __MOFEMUNKNOWNINTERFACE_HPP__
