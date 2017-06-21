/** \file simple_interface.cpp
  * \example simple_interface.hpp
  * \brief Simple interface
  *
  * This is simple test, calculate volume and apply divergence theorem to surface
  * integral to calculate volume. Integration on surface and in volume should be
  * equal to each other.
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

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

struct OpVolume: public VolumeElementForcesAndSourcesCore::UserDataOperator {
  Vec vOl;
  OpVolume(const std::string &field_name,Vec vol):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,OPROW),
  vOl(vol) {
  }
  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    if(type!=MBVERTEX) PetscFunctionReturn(0);
    const int nb_int_pts = getGaussPts().size2();
    FTensor::Tensor0<double*> t_w = getIntegrationWeight();
    double v = getMeasure();
    double vol = 0;
    for(int gg = 0;gg!=nb_int_pts;gg++) {
      vol += t_w*v;
      ++t_w;
    }
    ierr = VecSetValue(vOl,0,vol,ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode doWork(
    int row_side,int col_side,
    EntityType row_type, EntityType col_type,
    DataForcesAndSurcesCore::EntData &row_data,
    DataForcesAndSurcesCore::EntData &col_data
  )	{
    PetscFunctionBegin;
    // PetscPrintf(PETSC_COMM_WORLD,"domain: calculate matrix\n");
    PetscFunctionReturn(0);
  }
};

struct OpFace: public FaceElementForcesAndSourcesCore::UserDataOperator {
  Vec vOl;
  OpFace(const std::string &field_name,Vec vol):
  FaceElementForcesAndSourcesCore::UserDataOperator(field_name,OPROW),
  vOl(vol) {
  }
  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    if(type!=MBVERTEX) PetscFunctionReturn(0);
    const int nb_int_pts = getGaussPts().size2();
    FTensor::Tensor1<double*,3> t_normal = getTensor1Normal();
    FTensor::Tensor0<double*> t_w = getIntegrationWeight();
    FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
    FTensor::Index<'i',3> i;
    double vol = 0;
    for(int gg = 0;gg!=nb_int_pts;gg++) {
      vol -= (t_coords(i)*t_normal(i))*t_w;
      ++t_w;
      ++t_coords;
    }
    vol /= 6;
    ierr = VecSetValue(vOl,0,vol,ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode doWork(
    int row_side,int col_side,
    EntityType row_type, EntityType col_type,
    DataForcesAndSurcesCore::EntData &row_data,
    DataForcesAndSurcesCore::EntData &col_data
  )	{
    PetscFunctionBegin;
    // PetscPrintf(PETSC_COMM_WORLD,"boundary: calculate matrix\n");
    PetscFunctionReturn(0);
  }
};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  // initialize petsc
  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    // Create MoAB database
    moab::Core moab_core;
    moab::Interface& moab = moab_core;

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);
    MoFEM::Interface& m_field = mofem_core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);

    // Simple interface
    Simple *simple_interface;
    ierr = m_field.query_interface(simple_interface); CHKERRQ(ierr);
    {
      // get options from command line
      ierr = simple_interface->getOptions(); CHKERRQ(ierr);
      // load mesh file
      ierr = simple_interface->loadFile(); CHKERRQ(ierr);
      // add fields
      ierr = simple_interface->addDomainField("U",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      ierr = simple_interface->addBoundaryField("L",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
      // set fields order
      ierr = simple_interface->setFieldOrder("U",1); CHKERRQ(ierr);
      ierr = simple_interface->setFieldOrder("L",1); CHKERRQ(ierr);
      // setup problem
      ierr = simple_interface->setUp(); CHKERRQ(ierr);
      DM dm;
      // get dm
      ierr = simple_interface->getDM(&dm); CHKERRQ(ierr);
      // Declare finite element class and set integration rule
      MAKE_MY_FE_WITH_RULE(MyVol,VolumeElementForcesAndSourcesCore,4);
      MAKE_MY_FE_WITH_RULE(MyFace,FaceElementForcesAndSourcesCore,4);
      // create elements
      boost::shared_ptr<FEMethod> domainFE =
      boost::shared_ptr<ForcesAndSurcesCore>(new MyVol(m_field));
      boost::shared_ptr<FEMethod> boundaryFE =
      boost::shared_ptr<ForcesAndSurcesCore>(new MyFace(m_field));
      // set operators to the elements
      int ghosts[] = { 0 };
      Vec vol;
      ierr = VecCreateGhost(
        PETSC_COMM_WORLD,m_field.get_comm_rank()==0?1:0,1,m_field.get_comm_rank()==0?0:1,ghosts,&vol
      ); CHKERRQ(ierr);
      boost::static_pointer_cast<ForcesAndSurcesCore>(domainFE)->getOpPtrVector().push_back(new OpVolume("U",vol));
      boost::static_pointer_cast<ForcesAndSurcesCore>(boundaryFE)->getOpPtrVector().push_back(new OpFace("U",vol));
      boost::shared_ptr<FEMethod> null_fe;
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),domainFE); CHKERRQ(ierr);
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getBoundaryFEName(),boundaryFE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(vol); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(vol); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(vol,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(vol,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      if(m_field.get_comm_rank()==0) {
        double *array;
        ierr = VecGetArray(vol,&array); CHKERRQ(ierr);
        cout << "Volume = " << array[0] << endl;
        if(fabs(array[0])>1e-12) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"Should be zero");
        }
        ierr = VecRestoreArray(vol,&array); CHKERRQ(ierr);
      }
      ierr = VecDestroy(&vol); CHKERRQ(ierr);
      ierr = DMDestroy(&dm); CHKERRQ(ierr);
   }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  // finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
