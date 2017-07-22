/** \file Simple.cpp
 * \brief Implementation of simple interface
 * \ingroup mofem_simple_interface
 */


/* MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <Core.hpp>

#include <KspCtx.hpp>
#include <SnesCtx.hpp>
#include <TsCtx.hpp>
#include <DMMoFEM.hpp>

#include <Simple.hpp>

namespace MoFEM {

  PetscErrorCode Simple::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMSimple) {
      *iface = dynamic_cast<Simple*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  Simple::Simple(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)),
  bitLevel(BitRefLevel().set(0)),
  meshSet(0),
  boundaryMeshset(0),
  domainFE("dFE"),
  boundaryFE("bFE"),
  skeletonFE("sFE"),
  dIm(-1),
  dM(PETSC_NULL) {
    PetscLogEventRegister("LoadMesh",0,&USER_EVENT_SimpleLoadMesh);
    PetscLogEventRegister("buildFields",0,&USER_EVENT_SimpleBuildFields);
    PetscLogEventRegister("buildFiniteElements",0,&USER_EVENT_SimpleBuildFiniteElements);
    PetscLogEventRegister("SimpleSetUp",0,&USER_EVENT_SimpleBuildProblem);
    PetscLogEventRegister("SimpleKSPSolve",0,&USER_EVENT_SimpleKSPSolve);
    strcpy(meshFileName,"mesh.h5m");
  }
  Simple::~Simple() {
    
    if(dM!=PETSC_NULL) {
      ierr = DMDestroy(&dM); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
  }

  PetscErrorCode Simple::getOptions() {
    
    PetscBool flg = PETSC_TRUE;
    PetscFunctionBegin;
    ierr = PetscOptionsBegin(
      PETSC_COMM_WORLD,"",
      "Simple interface options","none"
    ); CHKERRQ(ierr);
    ierr = PetscOptionsString(
      "-file_name",
      "file name","", "mesh.h5m",meshFileName,255,&flg
    ); CHKERRQ(ierr);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::loadFile() {
    
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    PetscLogEventBegin(USER_EVENT_SimpleLoadMesh,0,0,0,0);
    // This is a case of distributed mesh and algebra. In that case each processor
    // keep only part of the problem.
    const char *option;
    option = "PARALLEL=READ_PART;"
    "PARALLEL_RESOLVE_SHARED_ENTS;"
    "PARTITION=PARALLEL_PARTITION;";
    rval = m_field.get_moab().load_file(meshFileName,0,option); CHKERRQ_MOAB(rval);
    ierr = m_field.rebuild_database(); CHKERRQ(ierr);
    // determine problem dimension
    if(dIm==-1) {
      int nb_ents_3d;
      rval = m_field.get_moab().get_number_entities_by_dimension(meshSet,3,nb_ents_3d,true); CHKERRQ_MOAB(rval);
      if(nb_ents_3d>0) {
        dIm = 3;
      } else {
        int nb_ents_2d;
        rval = m_field.get_moab().get_number_entities_by_dimension(meshSet,2,nb_ents_2d,true); CHKERRQ_MOAB(rval);
        if(nb_ents_2d>0) {
          dIm = 2;
        } else {
          dIm = 1;
        }
      }
    }
    Range ents;
    ierr = m_field.get_moab().get_entities_by_dimension(meshSet,dIm,ents,true);
    ierr = m_field.seed_ref_level(ents,bitLevel,false); CHKERRQ(ierr);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&m_field.get_moab(),m_field.get_comm());
    PetscLogEventEnd(USER_EVENT_SimpleLoadMesh,0,0,0,0);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::addDomainField(
    const std::string& name,
    const FieldSpace space,
    const FieldApproximationBase base,
    const FieldCoefficientsNumber nb_of_cooficients,
    const TagType tag_type,
    const enum MoFEMTypes bh,
    int verb
  ) {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    ierr = m_field.add_field(
      name, space, base, nb_of_cooficients, tag_type, bh, verb
    ); CHKERRQ(ierr);
    domainFields.push_back(name);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::addBoundaryField(
    const std::string& name,
    const FieldSpace space,
    const FieldApproximationBase base,
    const FieldCoefficientsNumber nb_of_cooficients,
    const TagType tag_type,
    const enum MoFEMTypes bh,
    int verb
  ) {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    ierr = m_field.add_field(
      name, space, base, nb_of_cooficients, tag_type, bh, verb
    ); CHKERRQ(ierr);
    boundaryFields.push_back(name);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::addSkeletionField(
    const std::string& name,
    const FieldSpace space,
    const FieldApproximationBase base,
    const FieldCoefficientsNumber nb_of_cooficients,
    const TagType tag_type,
    const enum MoFEMTypes bh,
    int verb
  ) {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    ierr = m_field.add_field(
      name, space, base, nb_of_cooficients, tag_type, bh, verb
    ); CHKERRQ(ierr);
    skeletonFields.push_back(name);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::addDataField(
    const std::string& name,
    const FieldSpace space,
    const FieldApproximationBase base,
    const FieldCoefficientsNumber nb_of_cooficients,
    const TagType tag_type,
    const enum MoFEMTypes bh,
    int verb
  ) {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    ierr = m_field.add_field(
      name, space, base, nb_of_cooficients, tag_type, bh, verb
    ); CHKERRQ(ierr);
    dataFields.push_back(name);
    PetscFunctionReturn(0);
  }


  PetscErrorCode Simple::defineFiniteElements() {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    // Define finite elements
    ierr = m_field.add_finite_element(domainFE); CHKERRQ(ierr);
    for(unsigned int ff = 0;ff!=domainFields.size();ff++) {
      ierr = m_field.modify_finite_element_add_field_row(domainFE,domainFields[ff]); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col(domainFE,domainFields[ff]); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data(domainFE,domainFields[ff]); CHKERRQ(ierr);
    }
    for(unsigned int ff = 0;ff!=dataFields.size();ff++) {
      ierr = m_field.modify_finite_element_add_field_data(domainFE,dataFields[ff]); CHKERRQ(ierr);
    }
    if(!boundaryFields.empty()) {
      ierr = m_field.add_finite_element(boundaryFE); CHKERRQ(ierr);
      for(unsigned int ff = 0;ff!=domainFields.size();ff++) {
        ierr = m_field.modify_finite_element_add_field_row(boundaryFE,domainFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_col(boundaryFE,domainFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data(boundaryFE,domainFields[ff]); CHKERRQ(ierr);
      }
      for(unsigned int ff = 0;ff!=boundaryFields.size();ff++) {
        ierr = m_field.modify_finite_element_add_field_row(boundaryFE,boundaryFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_col(boundaryFE,boundaryFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data(boundaryFE,boundaryFields[ff]); CHKERRQ(ierr);
      }
      for(unsigned int ff = 0;ff!=skeletonFields.size();ff++) {
        ierr = m_field.modify_finite_element_add_field_row(boundaryFE,skeletonFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_col(boundaryFE,skeletonFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data(boundaryFE,skeletonFields[ff]); CHKERRQ(ierr);
      }
    }
    if(!skeletonFields.empty()) {
      ierr = m_field.add_finite_element(skeletonFE); CHKERRQ(ierr);
      for(unsigned int ff = 0;ff!=domainFields.size();ff++) {
        ierr = m_field.modify_finite_element_add_field_row(skeletonFE,domainFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_col(skeletonFE,domainFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data(skeletonFE,domainFields[ff]); CHKERRQ(ierr);
      }
      for(unsigned int ff = 0;ff!=boundaryFields.size();ff++) {
        ierr = m_field.modify_finite_element_add_field_row(skeletonFE,boundaryFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_col(skeletonFE,boundaryFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data(skeletonFE,boundaryFields[ff]); CHKERRQ(ierr);
      }
      for(unsigned int ff = 0;ff!=skeletonFields.size();ff++) {
        ierr = m_field.modify_finite_element_add_field_row(skeletonFE,skeletonFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_col(skeletonFE,skeletonFields[ff]); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data(skeletonFE,skeletonFields[ff]); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::defineProblem() {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    if(dM!=PETSC_NULL) {
      ierr = DMDestroy(&dM); CHKERRQ(ierr);
    }
    // Create dm instance
    ierr = DMCreate(m_field.get_comm(),&dM);CHKERRQ(ierr);
    ierr = DMSetType(dM,"DMMOFEM");CHKERRQ(ierr);
    //set dm data structure which created mofem data structures
    ierr = DMMoFEMCreateMoFEM(dM,&m_field,"SimpleProblem",bitLevel); CHKERRQ(ierr);
    ierr = DMSetFromOptions(dM); CHKERRQ(ierr);
    ierr = DMMoFEMAddElement(dM,domainFE.c_str()); CHKERRQ(ierr);
    if(!boundaryFields.empty()) {
      ierr = DMMoFEMAddElement(dM,boundaryFE.c_str()); CHKERRQ(ierr);
    }
    if(!skeletonFields.empty()) {
      ierr = DMMoFEMAddElement(dM,skeletonFE.c_str()); CHKERRQ(ierr);
    }
    ierr = DMMoFEMSetIsPartitioned(dM,PETSC_TRUE); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::setFieldOrder(const std::string field_name,const int order,const Range* ents) {
    PetscFunctionBegin;
    fieldsOrder[field_name] = std::pair<int,Range>(order,ents==NULL?Range():Range(*ents));
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::buildFields() {
    
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    PetscLogEventBegin(USER_EVENT_SimpleBuildFields,0,0,0,0);
    // take skin
    {
      Range domain_ents;
      rval = m_field.get_moab().get_entities_by_dimension(meshSet,dIm,domain_ents,true); CHKERRQ_MOAB(rval);
      Skinner skin(&m_field.get_moab());
      Range domain_skin;
      rval = skin.find_skin(0,domain_ents,false,domain_skin); CHKERRQ_MOAB(rval);
      // filter not owned entities, those are not on boundary
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);
      Range proc_domain_skin;
      rval = pcomm->filter_pstatus(
        domain_skin,PSTATUS_SHARED|PSTATUS_MULTISHARED,PSTATUS_NOT,-1,&proc_domain_skin
      ); CHKERRQ_MOAB(rval);
      // cerr << proc_domain_skin << endl;
      // create boundary meshset
      if(boundaryMeshset!=0) {
        rval = m_field.get_moab().delete_entities(&boundaryMeshset,1); CHKERRQ_MOAB(rval);
      }
      rval = m_field.get_moab().create_meshset(MESHSET_SET,boundaryMeshset); CHKERRQ_MOAB(rval);
      rval = m_field.get_moab().add_entities(boundaryMeshset,proc_domain_skin); CHKERRQ_MOAB(rval);
      for(int dd = 0;dd!=dIm-1;dd++) {
        Range adj;
        rval = m_field.get_moab().get_adjacencies(
          proc_domain_skin,dd,false,adj,moab::Interface::UNION
        ); CHKERRQ_MOAB(rval);
        rval = m_field.get_moab().add_entities(boundaryMeshset,adj); CHKERRQ_MOAB(rval);
      }
    }
    // Add entities to the fields
    for(unsigned int ff = 0;ff!=domainFields.size();ff++) {
      ierr = m_field.add_ents_to_field_by_dim(meshSet,dIm,domainFields[ff]); CHKERRQ(ierr);
      ierr = m_field.synchronise_field_entities(domainFields[ff],0); CHKERRQ(ierr);
    }
    for(unsigned int ff = 0;ff!=dataFields.size();ff++) {
      ierr = m_field.add_ents_to_field_by_dim(meshSet,dIm,dataFields[ff]); CHKERRQ(ierr);
      ierr = m_field.synchronise_field_entities(dataFields[ff],0); CHKERRQ(ierr);
    }
    for(unsigned int ff = 0;ff!=boundaryFields.size();ff++) {
      ierr = m_field.add_ents_to_field_by_dim(boundaryMeshset,dIm-1,boundaryFields[ff]); CHKERRQ(ierr);
      ierr = m_field.synchronise_field_entities(boundaryFields[ff],0); CHKERRQ(ierr);
    }
    for(unsigned int ff = 0;ff!=skeletonFields.size();ff++) {
      ierr = m_field.add_ents_to_field_by_dim(meshSet,dIm-1,skeletonFields[ff]); CHKERRQ(ierr);
      ierr = m_field.synchronise_field_entities(skeletonFields[ff],0); CHKERRQ(ierr);
    }
    // Set order
    for(unsigned int ff = 0;ff!=domainFields.size();ff++) {
      if(fieldsOrder.find(domainFields[ff])==fieldsOrder.end()) {
        SETERRQ1(
          PETSC_COMM_WORLD,MOFEM_INVALID_DATA,
          "Order for field not set %s",domainFields[ff].c_str()
        );
      }
      int dds = 0;
      const Field *field = m_field.get_field_structure(domainFields[ff]);
      switch (field->getSpace()) {
        case L2: dds = dIm; break;
        case HDIV: dds = 2; break;
        case HCURL: dds = 1; break;
        case H1: dds = 1; break;
        default:
        SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Huston we have a problem");
      }
      if(field->getSpace()==H1) {
        ierr = m_field.set_field_order(meshSet,MBVERTEX,domainFields[ff],1); CHKERRQ(ierr);
      }
      for(int dd = dds;dd<=dIm;dd++) {
        Range ents;
        ierr = m_field.get_field_entities_by_dimension(domainFields[ff],dd,ents); CHKERRQ(ierr);
        if(!fieldsOrder.at(domainFields[ff]).second.empty()) { ents = intersect(ents,fieldsOrder.at(domainFields[ff]).second); }
        ierr = m_field.set_field_order(ents,domainFields[ff],fieldsOrder.at(domainFields[ff]).first); CHKERRQ(ierr);
      }
    }
    // Set order to data fiels
    for(unsigned int ff = 0;ff!=dataFields.size();ff++) {
      if(fieldsOrder.find(dataFields[ff])==fieldsOrder.end()) {
        SETERRQ1(
          PETSC_COMM_WORLD,MOFEM_INVALID_DATA,
          "Order for field not set %s",dataFields[ff].c_str()
        );
      }
      int dds = 0;
      const Field *field = m_field.get_field_structure(dataFields[ff]);
      switch (field->getSpace()) {
        case L2: dds = dIm; break;
        case HDIV: dds = 2; break;
        case HCURL: dds = 1; break;
        case H1: dds = 1; break;
        default:
        SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Huston we have a problem");
      }
      if(field->getSpace()==H1) {
        ierr = m_field.set_field_order(meshSet,MBVERTEX,dataFields[ff],1); CHKERRQ(ierr);
      }
      for(int dd = dds;dd<=dIm;dd++) {
        Range ents;
        ierr = m_field.get_field_entities_by_dimension(dataFields[ff],dd,ents); CHKERRQ(ierr);
        ierr = m_field.set_field_order(ents,dataFields[ff],fieldsOrder.at(dataFields[ff]).first); CHKERRQ(ierr);
      }
    }
    // Set order to boundary
    for(unsigned int ff = 0;ff!=boundaryFields.size();ff++) {
      if(fieldsOrder.find(boundaryFields[ff])==fieldsOrder.end()) {
        SETERRQ1(
          PETSC_COMM_WORLD,MOFEM_INVALID_DATA,
          "Order for field not set %s",boundaryFields[ff].c_str()
        );
      }
      int dds = 0;
      const Field *field = m_field.get_field_structure(boundaryFields[ff]);
      switch (field->getSpace()) {
        case L2: dds = dIm-1; break;
        case HDIV: dds = 2; break;
        case HCURL: dds = 1; break;
        case H1: dds = 1; break;
        default:
        SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Huston we have a problem");
      }
      if(field->getSpace()==H1) {
        ierr = m_field.set_field_order(meshSet,MBVERTEX,boundaryFields[ff],1); CHKERRQ(ierr);
      }
      for(int dd = dds;dd<=dIm-1;dd++) {
        Range ents;
        ierr = m_field.get_field_entities_by_dimension(boundaryFields[ff],dd,ents); CHKERRQ(ierr);
        ierr = m_field.set_field_order(ents,boundaryFields[ff],fieldsOrder.at(boundaryFields[ff]).first); CHKERRQ(ierr);
      }
    }
    // Set order to skeleton
    for(unsigned int ff = 0;ff!=skeletonFields.size();ff++) {
      if(fieldsOrder.find(skeletonFields[ff])==fieldsOrder.end()) {
        SETERRQ1(
          PETSC_COMM_WORLD,MOFEM_INVALID_DATA,
          "Order for field not set %s",skeletonFields[ff].c_str()
        );
      }
      int dds = 0;
      const Field *field = m_field.get_field_structure(skeletonFields[ff]);
      switch (field->getSpace()) {
        case L2: dds = dIm-1; break;
        case HDIV: dds = 2; break;
        case HCURL: dds = 1; break;
        case H1: dds = 1; break;
        default:
        SETERRQ(PETSC_COMM_WORLD,MOFEM_DATA_INCONSISTENCY,"Huston we have a problem");
      }
      if(field->getSpace()==H1) {
        ierr = m_field.set_field_order(meshSet,MBVERTEX,skeletonFields[ff],1); CHKERRQ(ierr);
      }
      for(int dd = dds;dd<=dIm-1;dd++) {
        Range ents;
        ierr = m_field.get_field_entities_by_dimension(skeletonFields[ff],dd,ents); CHKERRQ(ierr);
        ierr = m_field.set_field_order(ents,skeletonFields[ff],fieldsOrder.at(skeletonFields[ff]).first); CHKERRQ(ierr);
      }
    }
    // Build fields
    ierr = m_field.build_fields(); CHKERRQ(ierr);
    PetscLogEventEnd(USER_EVENT_SimpleBuildFields,0,0,0,0);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::buildFiniteElements() {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    PetscLogEventBegin(USER_EVENT_SimpleBuildFiniteElements,0,0,0,0);
    // Add finite elements
    ierr = m_field.add_ents_to_finite_element_by_dim(meshSet,dIm,domainFE,true); CHKERRQ(ierr);
    ierr = m_field.build_finite_elements(domainFE); CHKERRQ(ierr);
    if(!boundaryFields.empty()) {
        ierr = m_field.add_ents_to_finite_element_by_dim(boundaryMeshset,dIm-1,boundaryFE,true); CHKERRQ(ierr);
        ierr = m_field.build_finite_elements(boundaryFE); CHKERRQ(ierr);
    }
    if(!skeletonFields.empty()) {
      ierr = m_field.add_ents_to_finite_element_by_dim(meshSet,dIm-1,skeletonFE,true); CHKERRQ(ierr);
      ierr = m_field.build_finite_elements(skeletonFE); CHKERRQ(ierr);
    }
    PetscLogEventEnd(USER_EVENT_SimpleBuildFiniteElements,0,0,0,0);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::buildProblem() {
    
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    PetscLogEventBegin(USER_EVENT_SimpleBuildProblem,0,0,0,0);
    ierr = m_field.build_adjacencies(bitLevel); CHKERRQ(ierr);
    ierr = DMSetUp(dM); CHKERRQ(ierr);
    PetscLogEventEnd(USER_EVENT_SimpleBuildProblem,0,0,0,0);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::setUp() {
    
    PetscFunctionBegin;
    ierr = defineFiniteElements(); CHKERRQ(ierr);
    ierr = defineProblem(); CHKERRQ(ierr);
    ierr = buildFields(); CHKERRQ(ierr);
    ierr = buildFiniteElements(); CHKERRQ(ierr);
    ierr = buildProblem(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode Simple::getDM(DM *dm) {
    PetscFunctionBegin;
    PetscObjectReference((PetscObject)dM);
    *dm = dM;
    PetscFunctionReturn(0);
  }


}
