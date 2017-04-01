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
#include <DirichletBC.hpp>

using namespace boost::numeric;

DisplacementBCFEMethodPreAndPostProc::DisplacementBCFEMethodPreAndPostProc(
  MoFEM::Interface& m_field,const std::string &field_name,Mat Aij,Vec X,Vec F
):
mField(m_field),
fieldName(field_name),
dIag(1) {
  snes_B = Aij;
  snes_x = X;
  snes_f = F;
  ts_B = Aij;
  ts_u = X;
  ts_F = F;
};

DisplacementBCFEMethodPreAndPostProc::DisplacementBCFEMethodPreAndPostProc(
  MoFEM::Interface& m_field,const std::string &field_name
):
mField(m_field),
fieldName(field_name),
dIag(1) {
  snes_B = PETSC_NULL;
  snes_x = PETSC_NULL;
  snes_f = PETSC_NULL;
  ts_B = PETSC_NULL;
  ts_u = PETSC_NULL;
  ts_F = PETSC_NULL;
};

PetscErrorCode DisplacementBCFEMethodPreAndPostProc::iNitalize() {
  PetscFunctionBegin;
  if(mapZeroRows.empty() || !methodsOp.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)) {
      DisplacementCubitBcData mydata;
      ierr = it->getBcDataStructure(mydata); CHKERRQ(ierr);
      ublas::vector<double> scaled_values(3);
      scaled_values[0] = mydata.data.value1;
      scaled_values[1] = mydata.data.value2;
      scaled_values[2] = mydata.data.value3;
      ierr = MethodForForceScaling::applyScale(this,methodsOp,scaled_values); CHKERRQ(ierr);
      for(int dim = 0;dim<3;dim++) {
        Range ents;
        ierr = it->getMeshsetIdEntitiesByDimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
        if(dim>1) {
          Range _edges;
          ierr = mField.get_moab().get_adjacencies(ents,1,false,_edges,moab::Interface::UNION); CHKERRQ(ierr);
          ents.insert(_edges.begin(),_edges.end());
        }
        if(dim>0) {
          Range _nodes;
          rval = mField.get_moab().get_connectivity(ents,_nodes,true); CHKERRQ_MOAB(rval);
          ents.insert(_nodes.begin(),_nodes.end());
        }
        for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
          for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof_ptr)) {
            NumeredDofEntity *dof = dof_ptr->get();
            std::bitset<8> pstatus(dof->getPStatus());
            if(pstatus.test(0)) continue; //only local
            if(dof->getEntType() == MBVERTEX) {
              if(dof->getDofCoeffIdx() == 0 && mydata.data.flag1) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[0];
              }
              if(dof->getDofCoeffIdx() == 1 && mydata.data.flag2) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[1];
              }
              if(dof->getDofCoeffIdx() == 2 && mydata.data.flag3) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[2];
              }
            } else {
              if(dof->getDofCoeffIdx() == 0 && mydata.data.flag1) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
              }
              if(dof->getDofCoeffIdx() == 1 && mydata.data.flag2) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
              }
              if(dof->getDofCoeffIdx() == 2 && mydata.data.flag3) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
              }
            }
          }
        }
      }
    }
    dofsIndices.resize(mapZeroRows.size());
    dofsValues.resize(mapZeroRows.size());
    int ii = 0;
    std::map<DofIdx,FieldData>::iterator mit = mapZeroRows.begin();
    for(;mit!=mapZeroRows.end();mit++,ii++) {
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
      //std::cerr << dofsIndices[ii] << " " << dofsValues[ii] << std::endl;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DisplacementBCFEMethodPreAndPostProc::preProcess() {
  PetscFunctionBegin;

  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      snes_ctx = CTX_SNESSETFUNCTION;
      snes_x = ts_u;
      snes_f = ts_F;
      break;
    }
    case CTX_TSSETIJACOBIAN: {
      snes_ctx = CTX_SNESSETJACOBIAN;
      snes_B = ts_B;
      break;
    }
    default:
    break;
  }

  ierr = iNitalize(); CHKERRQ(ierr);

  if(snes_ctx == CTX_SNESNONE && ts_ctx == CTX_TSNONE) {
    if(dofsIndices.size()>0) {
      ierr = VecSetValues(
        snes_x,dofsIndices.size(),&dofsIndices[0],&dofsValues[0],INSERT_VALUES
      ); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(snes_x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(snes_x); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DisplacementBCFEMethodPreAndPostProc::postProcess() {
  PetscFunctionBegin;

  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      snes_ctx = CTX_SNESSETFUNCTION;
      snes_x = ts_u;
      snes_f = ts_F;
      break;
    }
    case CTX_TSSETIJACOBIAN: {
      snes_ctx = CTX_SNESSETJACOBIAN;
      snes_B = ts_B;
      break;
    }
    default:
    break;
  }

  if(snes_ctx == CTX_SNESNONE && ts_ctx == CTX_TSNONE) {
    if(snes_B) {
      ierr = MatAssemblyBegin(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatZeroRowsColumns(
        snes_B,
        dofsIndices.size(),
        dofsIndices.empty()?PETSC_NULL:&dofsIndices[0],
        dIag,PETSC_NULL,PETSC_NULL
      ); CHKERRQ(ierr);
    }
    if(snes_f) {
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      for(std::vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++) {
        ierr = VecSetValue(snes_f,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
    }
  }

  switch(snes_ctx) {
    case CTX_SNESNONE:
    break;
    case CTX_SNESSETFUNCTION: {
      if(!dofsIndices.empty()) {
        dofsXValues.resize(dofsIndices.size());
        ierr = VecGetValues(
          snes_x,dofsIndices.size(),
          dofsIndices.empty()?PETSC_NULL:&*dofsIndices.begin(),
          dofsXValues.empty()?PETSC_NULL:&*dofsXValues.begin()
        ); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      if(!dofsIndices.empty()) {
        int ii = 0;
        for(std::vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++,ii++) {
          double val = 0;
          if(!dofsXValues.empty()) {
            val += dofsXValues[ii];
            val += -mapZeroRows[*vit]; // in snes it is on the left hand side, that way -1
            dofsXValues[ii] = val;
          }
        }
        ierr = VecSetValues(
          snes_f,dofsIndices.size(),
          dofsIndices.empty()?PETSC_NULL:&*dofsIndices.begin(),
          dofsXValues.empty()?PETSC_NULL:&*dofsXValues.begin(),
          INSERT_VALUES
        ); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      ierr = MatAssemblyBegin(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatZeroRowsColumns(
        snes_B,
        dofsIndices.size(),
        dofsIndices.empty()?PETSC_NULL:&*dofsIndices.begin(),
        dIag,
        PETSC_NULL,
        PETSC_NULL
      ); CHKERRQ(ierr);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,1,"unknown snes stage");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SpatialPositionsBCFEMethodPreAndPostProc::iNitalize() {
  PetscFunctionBegin;
  if(mapZeroRows.empty() || !methodsOp.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)) {
      DisplacementCubitBcData mydata;
      ierr = it->getBcDataStructure(mydata); CHKERRQ(ierr);
      ublas::vector<double> scaled_values(3);
      scaled_values[0] = mydata.data.value1;
      scaled_values[1] = mydata.data.value2;
      scaled_values[2] = mydata.data.value3;
      ierr = MethodForForceScaling::applyScale(this,methodsOp,scaled_values); CHKERRQ(ierr);
      for(int dim = 0;dim<3;dim++) {
        Range ents;
        ierr = it->getMeshsetIdEntitiesByDimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
        if(dim>1) {
          Range _edges;
          ierr = mField.get_moab().get_adjacencies(ents,1,false,_edges,moab::Interface::UNION); CHKERRQ(ierr);
          ents.insert(_edges.begin(),_edges.end());
        }
        if(dim>0) {
          Range _nodes;
          rval = mField.get_moab().get_connectivity(ents,_nodes,true); CHKERRQ_MOAB(rval);
          ents.insert(_nodes.begin(),_nodes.end());
        }
        for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
          for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof_ptr)) {
            NumeredDofEntity *dof = dof_ptr->get();
            if(dof->getEntType() == MBVERTEX) {
              EntityHandle node = dof->getEnt();
              cOords.resize(3);
              rval = mField.get_moab().get_coords(&node,1,&*cOords.data().begin()); CHKERRQ_MOAB(rval);
              if(dof->getDofCoeffIdx() == 0 && mydata.data.flag1) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = cOords[0]+scaled_values[0];
              }
              if(dof->getDofCoeffIdx() == 1 && mydata.data.flag2) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = cOords[1]+scaled_values[1];
              }
              if(dof->getDofCoeffIdx() == 2 && mydata.data.flag3) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = cOords[2]+scaled_values[2];
              }
            } else {
              if(dof->getDofCoeffIdx() == 0 && mydata.data.flag1) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = dof->getFieldData();
              }
              if(dof->getDofCoeffIdx() == 1 && mydata.data.flag2) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = dof->getFieldData();
              }
              if(dof->getDofCoeffIdx() == 2 && mydata.data.flag3) {
                mapZeroRows[dof->getPetscGlobalDofIdx()] = dof->getFieldData();
              }
            }
          }
          for(std::vector<std::string>::iterator fit = fixFields.begin();fit!=fixFields.end();fit++) {
            for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,*fit,*eit,pcomm->rank(),dof_ptr)) {
              NumeredDofEntity *dof = dof_ptr->get();
              mapZeroRows[dof->getPetscGlobalDofIdx()] = dof->getFieldData();
            }
          }
        }
      }
    }
    dofsIndices.resize(mapZeroRows.size());
    dofsValues.resize(mapZeroRows.size());
    int ii = 0;
    std::map<DofIdx,FieldData>::iterator mit = mapZeroRows.begin();
    for(;mit!=mapZeroRows.end();mit++,ii++) {
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TemperatureBCFEMethodPreAndPostProc::iNitalize() {
  PetscFunctionBegin;
  if(mapZeroRows.empty() || !methodsOp.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|TEMPERATURESET,it)) {
      TemperatureCubitBcData mydata;
      ierr = it->getBcDataStructure(mydata); CHKERRQ(ierr);
      ublas::vector<double> scaled_values(1);
      scaled_values[0] = mydata.data.value1;
      ierr = MethodForForceScaling::applyScale(this,methodsOp,scaled_values); CHKERRQ(ierr);
      for(int dim = 0;dim<3;dim++) {
        Range ents;
        ierr = it->getMeshsetIdEntitiesByDimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
        if(dim>1) {
          Range _edges;
          ierr = mField.get_moab().get_adjacencies(ents,1,false,_edges,moab::Interface::UNION); CHKERRQ(ierr);
          ents.insert(_edges.begin(),_edges.end());
        }
        if(dim>0) {
          Range _nodes;
          rval = mField.get_moab().get_connectivity(ents,_nodes,true); CHKERRQ_MOAB(rval);
          ents.insert(_nodes.begin(),_nodes.end());
        }
        for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
          for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof_ptr)) {
            NumeredDofEntity *dof = dof_ptr->get();
            if(dof->getEntType() == MBVERTEX) {
              mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[0];
            } else {
              mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
            }
          }
        }
      }
    }
    dofsIndices.resize(mapZeroRows.size());
    dofsValues.resize(mapZeroRows.size());
    int ii = 0;
    std::map<DofIdx,FieldData>::iterator mit = mapZeroRows.begin();
    for(;mit!=mapZeroRows.end();mit++,ii++) {
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FixBcAtEntities::iNitalize() {
  PetscFunctionBegin;
  if(mapZeroRows.empty()) {
    for(std::vector<std::string>::iterator fit = fieldNames.begin();fit!=fieldNames.end();fit++) {
      for(Range::iterator eit = eNts.begin();eit!=eNts.end();eit++) {
        for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,*fit,*eit,mField.get_comm_rank(),dof_ptr)) {
          NumeredDofEntity *dof = dof_ptr->get();
          mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
        }
      }
    }
    dofsIndices.resize(mapZeroRows.size());
    dofsValues.resize(mapZeroRows.size());
    int ii = 0;
    std::map<DofIdx,FieldData>::iterator mit = mapZeroRows.begin();
    for(;mit!=mapZeroRows.end();mit++,ii++) {
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FixBcAtEntities::preProcess() {
  PetscFunctionBegin;

  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      snes_ctx = CTX_SNESSETFUNCTION;
      snes_x = ts_u;
      snes_f = ts_F;
      break;
    }
    case CTX_TSSETIJACOBIAN: {
      snes_ctx = CTX_SNESSETJACOBIAN;
      snes_B = ts_B;
      break;
    }
    default:
    break;
  }

  ierr = iNitalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FixBcAtEntities::postProcess() {
  PetscFunctionBegin;
  if(snes_ctx == CTX_SNESNONE && ts_ctx == CTX_TSNONE) {
    if(snes_B) {
      ierr = MatAssemblyBegin(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatZeroRowsColumns(
        snes_B,dofsIndices.size(),dofsIndices.empty()?PETSC_NULL:&dofsIndices[0],dIag,PETSC_NULL,PETSC_NULL
      ); CHKERRQ(ierr);
    }
    if(snes_f) {
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      int ii = 0;
      for(std::vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++,ii++) {
        ierr = VecSetValue(snes_f,*vit,dofsValues[ii],INSERT_VALUES); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
    }
  }

  switch(snes_ctx) {
    case CTX_SNESNONE: {}
    break;
    case CTX_SNESSETFUNCTION: {
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      int ii = 0;
      for(std::vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++,ii++) {
        ierr = VecSetValue(snes_f,*vit,dofsValues[ii],INSERT_VALUES); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      ierr = MatAssemblyBegin(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatZeroRowsColumns(
        snes_B,
        dofsIndices.size(),
        dofsIndices.empty()?PETSC_NULL:&*dofsIndices.begin(),
        dIag,
        PETSC_NULL,
        PETSC_NULL
      ); CHKERRQ(ierr);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,1,"unknown snes stage");
  }

  PetscFunctionReturn(0);
}


PetscErrorCode DirichletBCFromBlockSetFEMethodPreAndPostProc::iNitalize() {
  PetscFunctionBegin;
  if(mapZeroRows.empty() || !methodsOp.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->getName().compare(0,blocksetName.length(),blocksetName) == 0) {
        std::vector<double> mydata;
        ierr = it->getAttributes(mydata); CHKERRQ(ierr);
        ublas::vector<double> scaled_values(mydata.size());
        for(unsigned int ii = 0;ii<mydata.size();ii++) {
          scaled_values[ii] = mydata[ii];
        }
        ierr = MethodForForceScaling::applyScale(this,methodsOp,scaled_values); CHKERRQ(ierr);
        for(int dim = 0;dim<3;dim++) {
          Range ents;
          ierr = it->getMeshsetIdEntitiesByDimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
          if(dim>1) {
            Range edges;
            ierr = mField.get_moab().get_adjacencies(ents,1,false,edges,moab::Interface::UNION); CHKERRQ(ierr);
            ents.insert(edges.begin(),edges.end());
          }
          if(dim>0) {
            Range nodes;
            rval = mField.get_moab().get_connectivity(ents,nodes,true); CHKERRQ_MOAB(rval);
            ents.insert(nodes.begin(),nodes.end());
          }
          for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
            for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof_ptr)) {
              NumeredDofEntity *dof = dof_ptr->get();
              if(dof->getEntType() == MBVERTEX) {
                if(dof->getDofCoeffIdx() == 0) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[0];
                }
                if(dof->getDofCoeffIdx() == 1) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[1];
                }
                if(dof->getDofCoeffIdx() == 2) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[2];
                }
              } else {
                if(dof->getDofCoeffIdx() == 0) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
                }
                if(dof->getDofCoeffIdx() == 1) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
                }
                if(dof->getDofCoeffIdx() == 2) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
                }
              }
            }
          }
        }
      }
    }
    dofsIndices.resize(mapZeroRows.size());
    dofsValues.resize(mapZeroRows.size());
    int ii = 0;
    std::map<DofIdx,FieldData>::iterator mit = mapZeroRows.begin();
    for(;mit!=mapZeroRows.end();mit++,ii++) {
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}


PetscErrorCode DirichletBCFromBlockSetFEMethodPreAndPostProcWithFlags::iNitalize() {
  PetscFunctionBegin;
  if(mapZeroRows.empty() || !methodsOp.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->getName().compare(0,blocksetName.length(),blocksetName) == 0) {
        std::vector<double> mydata;
        ierr = it->getAttributes(mydata); CHKERRQ(ierr);
        ublas::vector<double> scaled_values(mydata.size());
        for(unsigned int ii = 0;ii<mydata.size();ii++) {
          scaled_values[ii] = mydata[ii];
        }

        ierr = MethodForForceScaling::applyScale(this,methodsOp,scaled_values); CHKERRQ(ierr);
        for(int dim = 0;dim<3;dim++) {
          Range ents;
          ierr = it->getMeshsetIdEntitiesByDimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
          if(dim>1) {
            Range edges;
            ierr = mField.get_moab().get_adjacencies(ents,1,false,edges,moab::Interface::UNION); CHKERRQ(ierr);
            ents.insert(edges.begin(),edges.end());
          }
          if(dim>0) {
            Range nodes;
            rval = mField.get_moab().get_connectivity(ents,nodes,true); CHKERRQ_MOAB(rval);
            ents.insert(nodes.begin(),nodes.end());
          }
          for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
            for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof_ptr)) {
              NumeredDofEntity *dof = dof_ptr->get();
              if(dof->getEntType() == MBVERTEX) {
                if(dof->getDofCoeffIdx() == 0) {
                  if(mydata[3] == 1)
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[0];
                }
                if(dof->getDofCoeffIdx() == 1) {
                  if(mydata[4] == 1)
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[1];
                }
                if(dof->getDofCoeffIdx() == 2) {
                  if(mydata[5] == 1)
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = scaled_values[2];
                }
              } else {
                if(dof->getDofCoeffIdx() == 0) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
                }
                if(dof->getDofCoeffIdx() == 1) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
                }
                if(dof->getDofCoeffIdx() == 2) {
                  mapZeroRows[dof->getPetscGlobalDofIdx()] = 0;
                }
              }
            }
          }
        }
      }
    }
    dofsIndices.resize(mapZeroRows.size());
    dofsValues.resize(mapZeroRows.size());
    int ii = 0;
    std::map<DofIdx,FieldData>::iterator mit = mapZeroRows.begin();
    for(;mit!=mapZeroRows.end();mit++,ii++) {
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}
