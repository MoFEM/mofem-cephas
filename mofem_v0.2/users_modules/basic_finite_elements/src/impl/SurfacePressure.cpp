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
#include <SurfacePressure.hpp>
#include <NodalForce.hpp>

using namespace boost::numeric;

NeummanForcesSurface::MyTriangleFE::MyTriangleFE(FieldInterface &m_field):
FaceElementForcesAndSourcesCore(m_field) {
}

NeummanForcesSurface::OpNeumannForce::OpNeumannForce(
  const string field_name,Vec &_F,bCForce &data,
  boost::ptr_vector<MethodsForOp> &methods_op
):
FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
F(_F),
dAta(data),
methodsOp(methods_op) {

}

PetscErrorCode NeummanForcesSurface::OpNeumannForce::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(data.getIndices().size()==0) PetscFunctionReturn(0);
  EntityHandle ent = getMoFEMFEPtr()->get_ent();
  if(dAta.tRis.find(ent)==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofMoFEMEntity *dof_ptr;
  ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->get_max_rank();
  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size());
  Nf.clear();

  for (unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getArea()*getGaussPts()(2,gg);
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

  ierr = MethodsForOp::applyScale(getFEMethod(), methodsOp, Nf); CHKERRQ(ierr);
  ierr = VecSetValues(F,data.getIndices().size(), &data.getIndices()[0], &Nf[0], ADD_VALUES); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

NeummanForcesSurface::OpNeumannPreassure::OpNeumannPreassure(
  const string field_name, Vec &_F,bCPreassure &data,boost::ptr_vector<MethodsForOp> &methods_op,bool ho_geometry
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
  if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofMoFEMEntity *dof_ptr;
  ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->get_max_rank();

  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size());
  bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

  //cerr << getNormal() << endl;
  //cerr << getNormals_at_GaussPt() << endl;

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

  /*cerr << "VecSetValues\n";
  cerr << Nf << endl;
  cerr << data.getIndices() << endl;*/
  ierr = MethodsForOp::applyScale(getFEMethod(),methodsOp,Nf); CHKERRQ(ierr);
  ierr = VecSetValues(F,data.getIndices().size(),
  &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

NeummanForcesSurface::OpNeumannFlux::OpNeumannFlux(
  const string field_name,Vec &_F,
  bCPreassure &data,boost::ptr_vector<MethodsForOp> &methods_op,
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
  if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofMoFEMEntity *dof_ptr;
  ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->get_max_rank();

  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size());
  Nf.clear();
  //cerr << getNormal() << endl;
  //cerr << getNormals_at_GaussPt() << endl;

  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getGaussPts()(2,gg);
    double flux;
    if(hoGeometry) {
      double area = cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
      flux = dAta.data.data.value1*area;
    } else {
      flux = dAta.data.data.value1*getArea();
    }
    cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);

  }

  //cerr << "VecSetValues\n";
  //cerr << Nf << endl;
  //cerr << data.getIndices() << endl;
  ierr = MethodsForOp::applyScale(getFEMethod(), methodsOp, Nf); CHKERRQ(ierr);
  ierr = VecSetValues(F, data.getIndices().size(), &data.getIndices()[0], &Nf[0], ADD_VALUES); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode NeummanForcesSurface::addForce(const string field_name,Vec &F,int ms_id) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_cubit_msId(ms_id,NODESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_bc_data_structure(mapForce[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapForce[ms_id].tRis,true); CHKERR_PETSC(rval);
  fe.getOpPtrVector().push_back(new OpNeumannForce(field_name,F,mapForce[ms_id],methodsOp));
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurface::addPreassure(const string field_name,Vec &F,int ms_id,bool ho_geometry) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_bc_data_structure(mapPreassure[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERR_PETSC(rval);
  fe.getOpPtrVector().push_back(new OpNeumannPreassure(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurface::addFlux(const string field_name,Vec &F,int ms_id,bool ho_geometry) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_bc_data_structure(mapPreassure[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERR_PETSC(rval);
  fe.getOpPtrVector().push_back(new OpNeumannFlux(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
  PetscFunctionReturn(0);
}
