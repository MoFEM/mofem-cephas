/** \file reading_med_file.cpp

  \brief Testing interface for reading and writing med files

*/

/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */


 #include <MoFEM.hpp>
 using namespace MoFEM;

 static char help[] = "...\n\n";

 int main(int argc, char *argv[]) {

   ErrorCode rval;
   PetscErrorCode ierr;

   PetscInitialize(&argc,&argv,(char *)0,help);

   moab::Core mb_instance;
   moab::Interface& moab = mb_instance;
   ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
   if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

   //Create MoFEM (Joseph) database
   MoFEM::Core core(moab);
   MoFEM::Interface& m_field = core;

   MedInterface *med_interface_ptr;
   ierr = m_field.query_interface(med_interface_ptr); CHKERRQ(ierr);

   ierr = med_interface_ptr->readMed(); CHKERRQ(ierr);
   ierr = med_interface_ptr->medGetFieldNames(); CHKERRQ(ierr);

   int ii = 0;
   const int check_list[] = { 2163, 624, 65, 104};
   for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,cit)) {
     EntityHandle meshset = cit->getMeshset();
     int nb_ents;
     rval = moab.get_number_entities_by_handle(meshset,nb_ents,true); CHKERRQ_MOAB(rval);
     ierr = PetscPrintf(PETSC_COMM_WORLD,"Nb of ents in %s %d\n",cit->getName().c_str(),nb_ents); CHKERRQ(ierr);
     if(nb_ents!=check_list[ii]) {
       SETERRQ2(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Wrong numbers of entities in meshset %d != %d",nb_ents,check_list[ii]);
     }
     ii++;
   }


   rval = moab.write_file("out.vtk","VTK",""); CHKERRQ_MOAB(rval);

   ierr = PetscFinalize(); CHKERRQ(ierr);

   return 0;
 }
