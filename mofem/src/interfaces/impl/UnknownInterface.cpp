/** \file UnknownInterface.cpp
 * \brief Unknown interfce implementation
 */



namespace MoFEM {

MoFEMErrorCode UnknownInterface::getLibVersion(Version &version) {
  MoFEMFunctionBeginHot;
  version =
      Version(MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR, MoFEM_VERSION_BUILD);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode UnknownInterface::getFileVersion(moab::Interface &moab,
                                                Version &version) {
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

MoFEMErrorCode UnknownInterface::setFileVersion(moab::Interface &moab,
                                                Version version) {
  MoFEMFunctionBegin;
  const EntityHandle root_meshset = 0;
  const int def_version[] = {-1, -1, -1};
  Tag th;
  rval = moab.tag_get_handle("MOFEM_VERSION", 3, MB_TYPE_INTEGER, th,
                             MB_TAG_CREAT | MB_TAG_MESH, &def_version);
  int *version_ptr;
  const void *tag_data[1];
  CHKERR moab.tag_get_by_ptr(th, &root_meshset, 1, tag_data);
  version_ptr = (int *)tag_data[0];
  version_ptr[0] = version.majorVersion;
  version_ptr[1] = version.minorVersion;
  version_ptr[2] = version.buildVersion;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode UnknownInterface::getInterfaceVersion(Version &version) {
  MoFEMFunctionBeginHot;
  version =
      Version(MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR, MoFEM_VERSION_BUILD);
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM