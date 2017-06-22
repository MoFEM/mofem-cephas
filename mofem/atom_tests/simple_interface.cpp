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
    // cerr << nb_int_pts << endl;
    FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
    FTensor::Tensor0<double*> t_ho_det = getFTenosr0HoMeasure();
    double v = getMeasure();
    double vol = 0;
    for(int gg = 0;gg!=nb_int_pts;gg++) {
      vol += t_w*t_ho_det*v;
      // cerr << t_ho_det << endl;
      ++t_w;
      ++t_ho_det;
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
    FTensor::Tensor1<double*,3> t_normal = getTensor1NormalsAtGaussPt();
    FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
    FTensor::Tensor1<double*,3> t_coords = getTensor1HoCoordsAtGaussPts();
    FTensor::Index<'i',3> i;
    double vol = 0;
    for(int gg = 0;gg!=nb_int_pts;gg++) {
      vol += (t_coords(i)*t_normal(i))*t_w;
      ++t_normal;
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

struct VolRule { int operator()(int,int,int) const { return 2; } };
struct FaceRule { int operator()(int,int,int) const { return 4; } };

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
      ierr = simple_interface->addDomainField("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
      ierr = simple_interface->addBoundaryField("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
      // set fields order
      ierr = simple_interface->setFieldOrder("MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
      // setup problem
      ierr = simple_interface->setUp(); CHKERRQ(ierr);
      // Project mesh coordinate on mesh
      Projection10NodeCoordsOnField ent_method(m_field,"MESH_NODE_POSITIONS");
      ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRQ(ierr);
      DM dm;
      // get dm
      ierr = simple_interface->getDM(&dm); CHKERRQ(ierr);
      // create elements
      boost::shared_ptr<FEMethod> domainFE =
      boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(m_field));
      boost::shared_ptr<FEMethod> boundaryFE =
      boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(m_field));
      // set integration rule
      boost::static_pointer_cast<ForcesAndSurcesCore>(domainFE)->getRuleHook = VolRule();
      boost::static_pointer_cast<ForcesAndSurcesCore>(boundaryFE)->getRuleHook = FaceRule();
      // create distributed vector to accumulate values from processors.
      int ghosts[] = { 0 };
      Vec vol,surf_vol;
      ierr = VecCreateGhost(
        PETSC_COMM_WORLD,m_field.get_comm_rank()==0?1:0,1,m_field.get_comm_rank()==0?0:1,ghosts,&vol
      ); CHKERRQ(ierr);
      ierr = VecDuplicate(vol,&surf_vol); CHKERRQ(ierr);
      // set operator to the volume element
      boost::static_pointer_cast<ForcesAndSurcesCore>(domainFE)->getOpPtrVector().push_back(
        new OpVolume("MESH_NODE_POSITIONS",vol)
      );
      // set operator to the face element
      boost::static_pointer_cast<ForcesAndSurcesCore>(boundaryFE)->getOpPtrVector().push_back(
        new OpFace("MESH_NODE_POSITIONS",surf_vol)
      );
      // make integration in volume (here real calculations starts)
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getDomainFEName(),domainFE); CHKERRQ(ierr);
      // make integration on boundary
      ierr = DMoFEMLoopFiniteElements(dm,simple_interface->getBoundaryFEName(),boundaryFE); CHKERRQ(ierr);
      // assemble volumes from processors and accumulate on processor of rank 0
      ierr = VecAssemblyBegin(vol); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(vol); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(vol,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(vol,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(surf_vol); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(surf_vol); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(surf_vol,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(surf_vol,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      if(m_field.get_comm_rank()==0) {
        double *a_vol;
        ierr = VecGetArray(vol,&a_vol); CHKERRQ(ierr);
        double *a_surf_vol;
        ierr = VecGetArray(surf_vol,&a_surf_vol); CHKERRQ(ierr);
        cout << "Volume = " << a_vol[0] << endl;
        cout << "Surf Volume = " << a_surf_vol[0] << endl;
        if(fabs(a_vol[0]-a_surf_vol[0])>1e-12) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"Should be zero");
        }
        ierr = VecRestoreArray(vol,&a_vol); CHKERRQ(ierr);
        ierr = VecRestoreArray(vol,&a_surf_vol); CHKERRQ(ierr);
      }
      // destroy vector
      ierr = VecDestroy(&vol); CHKERRQ(ierr);
      ierr = VecDestroy(&surf_vol); CHKERRQ(ierr);
      // destroy dm
      ierr = DMDestroy(&dm); CHKERRQ(ierr);
   }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  // finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
