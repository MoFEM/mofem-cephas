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

NeummanForcesSurface::MyTriangleFE::MyTriangleFE(MoFEM::Interface &m_field):
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
  EntityHandle ent = getNumeredEntFiniteElementPtr()->getEnt();
  if(dAta.tRis.find(ent)==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofEntity *dof_ptr;
  ierr = getNumeredEntFiniteElementPtr()->getRowDofsByPetscGlobalDofIdx(
    data.getIndices()[0],&dof_ptr
  ); CHKERRQ(ierr);
  int rank = dof_ptr->getNbOfCoeffs();
  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size(),false);
  Nf.clear();

  for (unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getGaussPts()(2,gg);
    if(hoGeometry) {
      val *= 0.5*cblas_dnrm2(3,&getNormalsAtGaussPt()(gg,0),1);
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
      switch (getFEMethod()->ts_ctx) {
        case FEMethod::CTX_TSSETIFUNCTION: {
          const_cast<FEMethod*>(getFEMethod())->snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
          const_cast<FEMethod*>(getFEMethod())->snes_x = getFEMethod()->ts_u;
          const_cast<FEMethod*>(getFEMethod())->snes_f = getFEMethod()->ts_F;
          break;
        }
        default:
        break;
      }
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
  if(dAta.tRis.find(getNumeredEntFiniteElementPtr()->getEnt())==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofEntity *dof_ptr;
  ierr = getNumeredEntFiniteElementPtr()->getRowDofsByPetscGlobalDofIdx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->getNbOfCoeffs();

  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size(),false);
  Nf.clear();

  //std::cerr << getNormal() << std::endl;
  //std::cerr << getNormalsAtGaussPt() << std::endl;

  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getGaussPts()(2,gg);
    for(int rr = 0;rr<rank;rr++) {

      double force;
      if(hoGeometry) {
        force = dAta.data.data.value1*getNormalsAtGaussPt()(gg,rr);
      } else {
        force = dAta.data.data.value1*getNormal()[rr];
      }
      cblas_daxpy(nb_row_dofs,0.5*val*force,&data.getN()(gg,0),1,&Nf[rr],rank);

    }

  }

  // if(type == MBTRI) {
  //   std::cerr << "Tri " << getNumeredEntFiniteElementPtr()->getEnt() << " getN " << data.getN() << std::endl;
  //   std::cerr << "Tri " << getNumeredEntFiniteElementPtr()->getEnt() << " getDiffN " << data.getDiffN() << std::endl;
  //   std::cerr << "Tri " << getNumeredEntFiniteElementPtr()->getEnt() << " Indices " << data.getIndices() << std::endl;
  // }

  /*std::cerr << "VecSetValues\n";
  std::cerr << Nf << std::endl;
  std::cerr << data.getIndices() << std::endl;*/
  ierr = MethodForForceScaling::applyScale(getFEMethod(),methodsOp,Nf); CHKERRQ(ierr);
  {
    Vec my_f;
    if(F == PETSC_NULL) {
      switch (getFEMethod()->ts_ctx) {
        case FEMethod::CTX_TSSETIFUNCTION: {
          const_cast<FEMethod*>(getFEMethod())->snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
          const_cast<FEMethod*>(getFEMethod())->snes_x = getFEMethod()->ts_u;
          const_cast<FEMethod*>(getFEMethod())->snes_f = getFEMethod()->ts_F;
          break;
        }
        default:
        break;
      }
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
  if(dAta.tRis.find(getNumeredEntFiniteElementPtr()->getEnt())==dAta.tRis.end()) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofEntity *dof_ptr;
  ierr = getNumeredEntFiniteElementPtr()->getColDofsByPetscGlobalDofIdx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->getNbOfCoeffs();

  int nb_row_dofs = data.getIndices().size()/rank;

  Nf.resize(data.getIndices().size(),false);
  Nf.clear();
  //std::cerr << getNormal() << std::endl;
  //std::cerr << getNormalsAtGaussPt() << std::endl;

  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

    double val = getGaussPts()(2,gg);
    double flux;
    if(hoGeometry) {
      double area = 0.5*cblas_dnrm2(3,&getNormalsAtGaussPt()(gg,0),1);
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
      switch (getFEMethod()->ts_ctx) {
        case FEMethod::CTX_TSSETIFUNCTION: {
          const_cast<FEMethod*>(getFEMethod())->snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
          const_cast<FEMethod*>(getFEMethod())->snes_x = getFEMethod()->ts_u;
          const_cast<FEMethod*>(getFEMethod())->snes_f = getFEMethod()->ts_F;
          break;
        }
        default:
        break;
      }
      my_f = getFEMethod()->snes_f;
    } else {
      my_f = F;
    }
    ierr = VecSetValues(my_f,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


PetscErrorCode NeummanForcesSurface::addForce(const std::string field_name,Vec F,int ms_id,bool ho_geometry,bool block_set) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  if(block_set) {
    //Add data from block set.
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_cubit_msId(ms_id,BLOCKSET,&cubit_meshset_ptr); CHKERRQ(ierr);
    std::vector<double> mydata;
    ierr = cubit_meshset_ptr->getAttributes(mydata); CHKERRQ(ierr);
    ublas::vector<double> force(mydata.size());
    for(unsigned int ii = 0;ii<mydata.size();ii++) {
      force[ii] = mydata[ii];
    }
    //Read forces from BLOCKSET Force (if exists)
    if(force.empty()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Force not given");
    }
    //Assign values from BLOCKSET FORCE to RHS vector. Info about native Cubit BC data structure can be found in BCData.hpp
    const string name = "Force";
    strncpy(mapForce[ms_id].data.data.name,name.c_str(),name.size()>5?5:name.size());
    double magnitude = sqrt(force[0]*force[0] + force[1]*force[1] + force[2]*force[2]);
    mapForce[ms_id].data.data.value1 = -magnitude; //< Force magnitude
    mapForce[ms_id].data.data.value2 = 0;
    mapForce[ms_id].data.data.value3 = force[0] / magnitude; //< X-component of force vector
    mapForce[ms_id].data.data.value4 = force[1] / magnitude; //< Y-component of force vector
    mapForce[ms_id].data.data.value5 = force[2] / magnitude; //< Z-component of force vector
    mapForce[ms_id].data.data.value6 = 0;
    mapForce[ms_id].data.data.value7 = 0;
    mapForce[ms_id].data.data.value8 = 0;
    mapForce[ms_id].data.data.zero[0] = 0;
    mapForce[ms_id].data.data.zero[1] = 0;
    mapForce[ms_id].data.data.zero[2] = 0;
    mapForce[ms_id].data.data.zero2 = 0;
    // std::cout << "TETSING ONLY:" << std::endl;
    // std::cout << mapForce[ms_id].data << std::endl;

    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapForce[ms_id].tRis,true); CHKERRQ_MOAB(rval);
    fe.getOpPtrVector().push_back(new OpNeumannForce(field_name,F,mapForce[ms_id],methodsOp,ho_geometry));

    // SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not implemented");
  } else {
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_cubit_msId(ms_id,NODESET,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->getBcDataStructure(mapForce[ms_id].data); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapForce[ms_id].tRis,true); CHKERRQ_MOAB(rval);
    fe.getOpPtrVector().push_back(new OpNeumannForce(field_name,F,mapForce[ms_id],methodsOp,ho_geometry));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurface::addPreassure(const std::string field_name,Vec F,int ms_id,bool ho_geometry,bool block_set) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  if(block_set) {
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_cubit_msId(ms_id,BLOCKSET,&cubit_meshset_ptr); CHKERRQ(ierr);
    std::vector<double> mydata;
    ierr = cubit_meshset_ptr->getAttributes(mydata); CHKERRQ(ierr);
    ublas::vector<double> pressure(mydata.size());
    for(unsigned int ii = 0;ii<mydata.size();ii++) {
      pressure[ii] = mydata[ii];
    }
    if(pressure.empty()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Pressure not given");
    }
    const string name = "Pressure";
    strncpy(mapPreassure[ms_id].data.data.name,name.c_str(),name.size()>8?8:name.size());
    mapPreassure[ms_id].data.data.flag1 = 0;
    mapPreassure[ms_id].data.data.flag2 = 1;
    mapPreassure[ms_id].data.data.value1 = pressure[0];
    mapPreassure[ms_id].data.data.zero = 0;
    // std::cerr << "TETSING ONLY:" << std::endl;
    // std::cerr << mapPreassure[ms_id].data << std::endl;
    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERRQ_MOAB(rval);
    fe.getOpPtrVector().push_back(new OpNeumannPreassure(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
  } else {
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->getBcDataStructure(mapPreassure[ms_id].data); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERRQ_MOAB(rval);
    fe.getOpPtrVector().push_back(new OpNeumannPreassure(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurface::addFlux(const std::string field_name,Vec F,int ms_id,bool ho_geometry) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->getBcDataStructure(mapPreassure[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERRQ_MOAB(rval);
  fe.getOpPtrVector().push_back(new OpNeumannFlux(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
  PetscFunctionReturn(0);
}
