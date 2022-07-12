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

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Read parameters from line command
    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    ierr = PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    CHKERRG(ierr);
#else
    ierr = PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
    CHKERRG(ierr);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    // Read mesh to MOAB
    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // Open mesh_file_name.txt for writing
    std::ofstream myfile;
    myfile.open((std::string(mesh_file_name) + ".txt").c_str());

    std::cout << "<<<< All BLOCKSETs, SIDESETs and NODESETs >>>>>" << std::endl;
    for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, it)) {
      std::cout << it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout << "<<<< BLOCKSETs >>>>>" << std::endl;
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {
      std::cout << it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout << "<<<< NODESETs >>>>>" << std::endl;
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, NODESET, it)) {
      std::cout << it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout << "<<<< SIDESETs >>>>>" << std::endl;
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, SIDESET, it)) {
      std::cout << it->getName() << std::endl;
      myfile << it->getName() << std::endl;
    }
    std::cout << "<<<< MeshSet of Name Moon >>>>" << std::endl;
    for (_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field, "Moon", it)) {
      std::cout << it->getName() << std::endl;
      myfile << it->getName() << std::endl;
      if (it->getBcTypeULong() & BLOCKSET) {
        std::cout << "BLOCKSET" << std::endl;
        myfile << "BLOCKSET" << std::endl;
      }
      if (it->getBcTypeULong() & SIDESET) {
        std::cout << "SIDESET" << std::endl;
        myfile << "SIDESET" << std::endl;
      }
      if (it->getBcTypeULong() & NODESET) {
        std::cout << "NODESET" << std::endl;
        myfile << "NODESET" << std::endl;
      }
    }

    // Close mesh_file_name.txt
    myfile.close();
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}
