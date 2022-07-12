/** \file UnknownInterface.cpp
 * \brief Unknown interfce implementation
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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