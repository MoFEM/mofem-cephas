/* \fiele SurfacePressure.cpp
  \brief Implementation of pressure and forces on triangles surface

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
#include <MethodForForceScaling.hpp>
#include <SurfacePressure.hpp>
#include <NodalForce.hpp>

using namespace boost::numeric;

NeummanForcesSurface::MyTriangleFE::MyTriangleFE(FieldInterface &m_field):
FaceElementForcesAndSourcesCore(m_field) {
}

NeummanForcesSurface::OpNeumannForce::OpNeumannForce(
  const std::string field_name,Vec _F,bCForce &data,
  boost::ptr_vector<MethodForForceScaling> &methods_op,
  bool ho_geometry
):
FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
F(_F),
dAta(data),
methodsOp(methods_op),
hoGeometry(ho_geometry) {

}

PetscErrorCode NeummanForcesSurface::OpNeumannForce::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(data.getIndices().size()==0) PetscFunctionReturn(0);
  EntityHandle ent = getNumeredEntFiniteElementPtr()->get_ent();
  if(dAta.tRis.find(ent)==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofEntity *dof_ptr;
  ierr = getNumeredEntFiniteElementPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(
    data.getIndices()[0],&dof_ptr
  ); CHKERRQ(ierr);
  int rank = dof_ptr->get_nb_of_coeffs();
  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size(),false);
  Nf.clear();

  for (unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getGaussPts()(2,gg);
    if(hoGeometry) {
      val *= 0.5*cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
    } else {
      val *= getArea();
    }

    for (int rr = 0;rr<rank;rr++) {

      double force;
      if(rr == 0) {
        force = dAta.data.data.value3;
      } else if(rr == 1) {
        force = dAta.data.data.value4;
      } else if(rr == 2) {
        force = dAta.data.data.value5;
      } else {
        SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      force *= dAta.data.data.value1;
      cblas_daxpy(nb_row_dofs,val*force,&data.getN()(gg,0),1,&Nf[rr],rank);

    }

  }

  ierr = MethodForForceScaling::applyScale(getFEMethod(), methodsOp, Nf); CHKERRQ(ierr);
  {
    Vec my_f;
    if(F == PETSC_NULL) {
      my_f = getFEMethod()->snes_f;
    } else {
      my_f = F;
    }
    ierr = VecSetValues(
      my_f,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES
    ); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

NeummanForcesSurface::OpNeumannPreassure::OpNeumannPreassure(
  const std::string field_name, Vec _F,bCPreassure &data,boost::ptr_vector<MethodForForceScaling> &methods_op,bool ho_geometry
):
FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
F(_F),
dAta(data),
methodsOp(methods_op),
hoGeometry(ho_geometry) {}

PetscErrorCode NeummanForcesSurface::OpNeumannPreassure::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(data.getIndices().size()==0) PetscFunctionReturn(0);
  if(dAta.tRis.find(getNumeredEntFiniteElementPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofEntity *dof_ptr;
  ierr = getNumeredEntFiniteElementPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->get_nb_of_coeffs();

  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size(),false);
  Nf.clear();

  //std::cerr << getNormal() << std::endl;
  //std::cerr << getNormals_at_GaussPt() << std::endl;

  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getGaussPts()(2,gg);
    for(int rr = 0;rr<rank;rr++) {

      double force;
      if(hoGeometry) {
        force = dAta.data.data.value1*getNormals_at_GaussPt()(gg,rr);
      } else {
        force = dAta.data.data.value1*getNormal()[rr];
      }
      cblas_daxpy(nb_row_dofs,0.5*val*force,&data.getN()(gg,0),1,&Nf[rr],rank);

    }

  }

  // if(type == MBTRI) {
  //   std::cerr << "Tri " << getNumeredEntFiniteElementPtr()->get_ent() << " getN " << data.getN() << std::endl;
  //   std::cerr << "Tri " << getNumeredEntFiniteElementPtr()->get_ent() << " getDiffN " << data.getDiffN() << std::endl;
  //   std::cerr << "Tri " << getNumeredEntFiniteElementPtr()->get_ent() << " Indices " << data.getIndices() << std::endl;
  // }

  /*std::cerr << "VecSetValues\n";
  std::cerr << Nf << std::endl;
  std::cerr << data.getIndices() << std::endl;*/
  ierr = MethodForForceScaling::applyScale(getFEMethod(),methodsOp,Nf); CHKERRQ(ierr);
  {
    Vec my_f;
    if(F == PETSC_NULL) {
      my_f = getFEMethod()->snes_f;
    } else {
      my_f = F;
    }
    ierr = VecSetValues(
      my_f,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES
    ); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

NeummanForcesSurface::OpNeumannFlux::OpNeumannFlux(
  const std::string field_name,Vec _F,
  bCPreassure &data,boost::ptr_vector<MethodForForceScaling> &methods_op,
  bool ho_geometry
):
FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
F(_F),
dAta(data),
methodsOp(methods_op),
hoGeometry(ho_geometry) {}

PetscErrorCode NeummanForcesSurface::OpNeumannFlux::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(data.getIndices().size()==0) PetscFunctionReturn(0);
  if(dAta.tRis.find(getNumeredEntFiniteElementPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofEntity *dof_ptr;
  ierr = getNumeredEntFiniteElementPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->get_nb_of_coeffs();

  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size(),false);
  Nf.clear();
  //std::cerr << getNormal() << std::endl;
  //std::cerr << getNormals_at_GaussPt() << std::endl;

  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getGaussPts()(2,gg);
    double flux;
    if(hoGeometry) {
      double area = 0.5*cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
      flux = dAta.data.data.value1*area;
    } else {
      flux = dAta.data.data.value1*getArea();
    }
    cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);

  }

  //std::cerr << "VecSetValues\n";
  //std::cerr << Nf << std::endl;
  //std::cerr << data.getIndices() << std::endl;
  ierr = MethodForForceScaling::applyScale(getFEMethod(), methodsOp, Nf); CHKERRQ(ierr);
  {
    Vec my_f;
    if(F == PETSC_NULL) {
      my_f = getFEMethod()->snes_f;
    } else {
      my_f = F;
    }
    ierr = VecSetValues(my_f,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


PetscErrorCode NeummanForcesSurface::addForce(const std::string field_name,Vec F,int ms_id,bool ho_geometry) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_cubit_msId(ms_id,NODESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_bc_data_structure(mapForce[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapForce[ms_id].tRis,true); CHKERRQ_MOAB(rval);
  fe.getOpPtrVector().push_back(new OpNeumannForce(field_name,F,mapForce[ms_id],methodsOp,ho_geometry));
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurface::addPreassure(const std::string field_name,Vec F,int ms_id,bool ho_geometry) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_bc_data_structure(mapPreassure[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERRQ_MOAB(rval);
  fe.getOpPtrVector().push_back(new OpNeumannPreassure(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurface::addFlux(const std::string field_name,Vec F,int ms_id,bool ho_geometry) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_bc_data_structure(mapPreassure[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERRQ_MOAB(rval);
  fe.getOpPtrVector().push_back(new OpNeumannFlux(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
  PetscFunctionReturn(0);
}
