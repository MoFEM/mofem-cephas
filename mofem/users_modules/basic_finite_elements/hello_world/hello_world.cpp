/**
 * \file hello_world.cpp
 * \ingroup mofem_simple_interface
 * \example hello_world.cpp
 *
 * Prints basic information about users data operator.
 * See more details in \ref hello_world_tut1
 *
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

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

static map<EntityType,std::string> type_name;

struct OpRow: public ForcesAndSourcesCore::UserDataOperator {
 OpRow(const std::string& field_name):
 ForcesAndSourcesCore::UserDataOperator(field_name,field_name,OPROW) {
 }
 MoFEMErrorCode doWork(int side,EntityType type,DataForcesAndSourcesCore::EntData &data) {
   MoFEMFunctionBeginHot;
   if(type == MBVERTEX) {
     // get number of evaluated element in the loop
     std::cout << std::endl << "**** " << 	getNinTheLoop() << " **** " << std::endl;
     std::cout <<"**** Operators **** " << std::endl;
   }
   std::cout << "Hello Operator OpRow:"
   << " field name " << rowFieldName << " side " << side << " type " << type_name[type]
   << " nb dofs on entity " << data.getIndices().size()
   << std::endl;
   MoFEMFunctionReturnHot(0);
 }
};

struct OpRowCol: public ForcesAndSourcesCore::UserDataOperator {
 OpRowCol(
   const std::string row_field,
   const std::string col_field,
   const bool symm
 ):
 ForcesAndSourcesCore::UserDataOperator(row_field,col_field,OPROWCOL,symm) {
 }
 virtual MoFEMErrorCode doWork(
   int row_side,int col_side,
   EntityType row_type,EntityType col_type,
   DataForcesAndSourcesCore::EntData &row_data,
   DataForcesAndSourcesCore::EntData &col_data
 ) {
   MoFEMFunctionBeginHot;
   std::cout << "Hello Operator OpRowCol:"
   << " row field name " << rowFieldName
   << " row side " << row_side << " row type " << type_name[row_type]
   << " nb dofs on row entity" << row_data.getIndices().size() << " : "
   << " col field name " << colFieldName
   << " col side " << col_side << " col type " << type_name[col_type]
   << " nb dofs on col entity" << col_data.getIndices().size() << std::endl;
   MoFEMFunctionReturnHot(0);
 }
};

struct OpVolume: public VolumeElementForcesAndSourcesCore::UserDataOperator {
 OpVolume(const std::string& field_name):
 VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,field_name,OPROW) {
 }
 MoFEMErrorCode doWork(int side,EntityType type,DataForcesAndSourcesCore::EntData &data) {
   MoFEMFunctionBeginHot;
   if(type == MBVERTEX) {
    std::cout << "Hello Operator OpVolume:" << " volume " << getVolume() << endl;
   }
   MoFEMFunctionReturnHot(0);
 }
};

struct OpFace: public FaceElementForcesAndSourcesCore::UserDataOperator {
 OpFace(const std::string& field_name):
 FaceElementForcesAndSourcesCore::UserDataOperator(field_name,OPROW) {
 }
 MoFEMErrorCode doWork(int side,EntityType type,DataForcesAndSourcesCore::EntData &data) {
   MoFEMFunctionBeginHot;
   if(type == MBVERTEX) {
     std::cout << "Hello Operator OpFace:" << " normal " << getNormal() << endl;
   }
   MoFEMFunctionReturnHot(0);
 }
};

struct OpFaceSide: public FaceElementForcesAndSourcesCore::UserDataOperator {
 boost::shared_ptr<VolumeElementForcesAndSourcesCoreOnSide>& feSidePtr;
 OpFaceSide(
   const std::string& field_name,
   boost::shared_ptr<VolumeElementForcesAndSourcesCoreOnSide>& fe_side_ptr
 ):
 FaceElementForcesAndSourcesCore::UserDataOperator(field_name,OPROW),
 feSidePtr(fe_side_ptr) {
 }
 MoFEMErrorCode doWork(int side,EntityType type,DataForcesAndSourcesCore::EntData &data) {
   
   MoFEMFunctionBeginHot;
   if(type == MBVERTEX) {
     std::cout << "Hello Operator OpSideFace" << endl;
     ierr = loopSideVolumes("dFE",*feSidePtr); CHKERRQ(ierr);
   }
   MoFEMFunctionReturnHot(0);
 }
};

struct OpVolumeSide: public VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator {
 OpVolumeSide(const std::string& field_name):
 VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator(field_name,field_name,OPROW) {
 }
 MoFEMErrorCode doWork(int side,EntityType type,DataForcesAndSourcesCore::EntData &data) {
   MoFEMFunctionBeginHot;
   if(type == MBVERTEX) {
     std::cout << "Hello Operator OpVolumeSide:"
     << " volume " << getVolume()
     << " normal " << getNormal()
     << endl;
   }
   MoFEMFunctionReturnHot(0);
 }
};

int main(int argc, char *argv[]) {

 

 type_name[MBVERTEX] = "Vertex";
 type_name[MBEDGE] = "Edge";
 type_name[MBTRI] = "Triangle";
 type_name[MBTET] = "Tetrahedra";

 // initialize petsc
 PetscInitialize(&argc,&argv,(char *)0,help);
 // Register DM Manager
 DMType dm_name = "DMMOFEM";
 ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);

 try {

   // Create MoAB database
   moab::Core moab_core;
   moab::Interface& moab = moab_core;

   // Create MoFEM database and link it to MoAB
   MoFEM::Core mofem_core(moab);
   MoFEM::Interface& m_field = mofem_core;


   // Simple interface
   Simple *simple_interface;
   ierr = m_field.getInterface(simple_interface); CHKERRQ(ierr);

   // get options from command line
   ierr = simple_interface->getOptions(); CHKERRQ(ierr);
   // load mesh file
   ierr = simple_interface->loadFile(); CHKERRQ(ierr);
   // add fields
   ierr = simple_interface->addDomainField("U",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
   ierr = simple_interface->addBoundaryField("L",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
   ierr = simple_interface->addSkeletonField("S",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
   // set fields order
   ierr = simple_interface->setFieldOrder("U",4); CHKERRQ(ierr);
   ierr = simple_interface->setFieldOrder("L",3); CHKERRQ(ierr);
   ierr = simple_interface->setFieldOrder("S",3); CHKERRQ(ierr);
   // setup problem
   ierr = simple_interface->setUp(); CHKERRQ(ierr);
   // Create elements
   boost::shared_ptr<ForcesAndSourcesCore> domain_fe(new VolumeElementForcesAndSourcesCore(m_field));
   boost::shared_ptr<ForcesAndSourcesCore> boundary_fe(new FaceElementForcesAndSourcesCore(m_field));
   boost::shared_ptr<ForcesAndSourcesCore> skeleton_fe(new FaceElementForcesAndSourcesCore(m_field));
   boost::shared_ptr<VolumeElementForcesAndSourcesCoreOnSide> side_fe(new VolumeElementForcesAndSourcesCoreOnSide(m_field));
   // create distributed vector to accumulate values from processors.
   // set operator to the volume element
   domain_fe->getOpPtrVector().push_back(new OpRow("U"));
   domain_fe->getOpPtrVector().push_back(new OpRowCol("U","U",true));
   domain_fe->getOpPtrVector().push_back(new OpVolume("U"));
   // set operator to the face element
   boundary_fe->getOpPtrVector().push_back(new OpRow("L"));
   boundary_fe->getOpPtrVector().push_back(new OpRowCol("U","L",false));
   boundary_fe->getOpPtrVector().push_back(new OpFace("L"));
   // set operator to the face element on skeleton
   skeleton_fe->getOpPtrVector().push_back(new OpRow("S"));
   skeleton_fe->getOpPtrVector().push_back(new OpFaceSide("S",side_fe));
   // set operator to the volume on side of the Skeleton face
   side_fe->getOpPtrVector().push_back(new OpVolumeSide("U"));
   DM dm;
   // get dm
   ierr = simple_interface->getDM(&dm); CHKERRQ(ierr);
   // iterate domain elements
   ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),domain_fe); CHKERRQ(ierr);
   // iterate boundary elements
   ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getBoundaryFEName(),boundary_fe); CHKERRQ(ierr);
   // iterate skeleton element
   ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getSkeletonFEName(),skeleton_fe); CHKERRQ(ierr);
   // destroy dm
   ierr = DMDestroy(&dm); CHKERRQ(ierr);

 } catch (MoFEMException const &e) {
   SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
 }

 // finish work cleaning memory, getting statistics, etc.
 ierr = PetscFinalize(); CHKERRQ(ierr);

 return 0;
}
