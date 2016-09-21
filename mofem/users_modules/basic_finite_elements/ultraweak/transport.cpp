/** \fi;e transport.cpp
\brief Example implementation of transport problem using ultra-week formulation

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

struct BcFluxData {
  Range eNts;
  double fLux;
};
typedef map<int,BcFluxData> BcFluxMap;

/** thefine sources and other stuff
  *
  * UltraWeakTransportElement is a class collecting functons, opertors and
  * data for ultra week implementation of transport element. See there to
  * learn how elements are created or how operators look like.
  *
  * Some methods in UltraWeakTransportElement are abstract, f.e. user need to
  * implement own surce therm.
  *
  */
struct MyUltraWeakFE: public UltraWeakTransportElement {

  BcFluxMap &bcFluxMap;
  EntityHandle lastEnt;
  double lastFlux;

  MyUltraWeakFE(MoFEM::Interface &m_field,BcFluxMap &bc_flux_map):
  UltraWeakTransportElement(m_field),
  bcFluxMap(bc_flux_map),
  lastEnt(0),
  lastFlux(0) {
  }

  /**
   * \brief set source term
   * @param  ent  handle to entity on which function is evaluated
   * @param  x    coord
   * @param  y    coord
   * @param  z    coord
   * @param  flux reference to source term set by function
   * @return      error code
   */
  PetscErrorCode getFlux(EntityHandle ent,const double x,const double y,const double z,double &flux) {
    PetscFunctionBegin;
    flux = 0;
    PetscFunctionReturn(0);
  }

  /**
   * \brief natural (Dirihlet) boundary conditions (set values)
   * @param  ent   handle to finite element entity
   * @param  x     coord
   * @param  y     coord
   * @param  z     coord
   * @param  value reference to value set by function
   * @return       error code
   */
  PetscErrorCode getBcOnValues(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &value) {
    PetscFunctionBegin;
    value = 0;
    PetscFunctionReturn(0);
  }

  /**
   * \brief essential (Neumann) boundary condition (set fluxes)
   * @param  ent  handle to finite element entity
   * @param  x    coord
   * @param  y    coord
   * @param  z    coord
   * @param  flux reference to flux which is set by function
   * @return      [description]
   */
  PetscErrorCode getBcOnFluxes(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &flux) {
    PetscFunctionBegin;
    if(lastEnt==ent) {
      flux = lastFlux;
    } else {
      flux = 0;
      for(BcFluxMap::iterator mit = bcFluxMap.begin();mit!=bcFluxMap.end();mit++) {
        Range &tris = mit->second.eNts;
        if(tris.find(ent)!=tris.end()) {
          flux = mit->second.fLux;
        }
      }
      lastEnt = ent;
      lastFlux = flux;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode addBoundaryElements(BitRefLevel &ref_level) {
    MoABErrorCode rval;
    PetscErrorCode ierr;
    PetscFunctionBegin;
    Range tets;
    EntityHandle fe_meshset = mField.get_finite_element_meshset("ULTRAWEAK"); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(fe_meshset,MBTET,tets); CHKERR_MOAB(rval);
    Skinner skin(&mField.get_moab());
    Range skin_faces; // skin faces from 3d ents
    rval = skin.find_skin(0,tets,false,skin_faces); CHKERR_MOAB(rval);
    // note: what is essential (dirichlet) is natural (neumann) for ultra weak compared to classical FE
    Range natural_bc;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|TEMPERATURESET,it)) {
      Range tris;
      ierr = it->getMeshsetIdEntitiesByDimension(mField.get_moab(),2,tris,true); CHKERRQ(ierr);
      natural_bc.insert(tris.begin(),tris.end());
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|HEATFLUXSET,it)) {
      HeatFluxCubitBcData mydata;
      ierr = it->getBcDataStructure(mydata); CHKERRQ(ierr);
      if(mydata.data.flag1==1) {
        Range tris;
        ierr = it->getMeshsetIdEntitiesByDimension(mField.get_moab(),2,tris,true); CHKERRQ(ierr);
        bcFluxMap[it->getMeshsetId()].eNts = tris;
        bcFluxMap[it->getMeshsetId()].fLux = mydata.data.value1;
        // cerr << bcFluxMap[it->getMeshsetId()].eNts << endl;
        // cerr << bcFluxMap[it->getMeshsetId()].fLux << endl;
      }
    }
    Range essential_bc = subtract(skin_faces,natural_bc);
    Range bit_tris;
    ierr = mField.get_entities_by_type_and_ref_level(ref_level,BitRefLevel().set(),MBTRI,bit_tris);
    essential_bc = intersect(bit_tris,essential_bc);
    natural_bc = intersect(bit_tris,natural_bc);
    ierr = mField.add_ents_to_finite_element_by_TRIs(essential_bc,"ULTRAWEAK_BCFLUX"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_finite_element_by_TRIs(natural_bc,"ULTRAWEAK_BCVALUE"); CHKERRQ(ierr);
    // ierr = mField.add_ents_to_finite_element_by_TRIs(skin_faces,"ULTRAWEAK_BCVALUE"); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode squashBits() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    BitRefLevel all_but_0;
    all_but_0.set(0);
    all_but_0.flip();
    BitRefLevel new_bit;
    new_bit.set(BITREFLEVEL_SIZE-1); // Garbage level
    const RefEntity_multiIndex *refined_ents_ptr;
    ierr = mField.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
    RefEntity_multiIndex::iterator mit = refined_ents_ptr->begin();
    for(;mit!=refined_ents_ptr->end();mit++) {
      if(mit->get()->getEntType() == MBENTITYSET) continue;
      BitRefLevel bit = mit->get()->getBitRefLevel();
      if((all_but_0&bit)==bit) {
        const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(mit,RefEntity_change_set_bit(new_bit));
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode refienMesh(
    UltraWeakTransportElement &ufe,const int nb_levels
  ) {
    PetscErrorCode ierr;
    MoABErrorCode rval;
    MeshRefinment *refine_ptr;
    PetscFunctionBegin;
    Range refined_edges;
    BitRefLevel all_but_0;
    all_but_0.set(0);
    all_but_0.flip();
    ierr = mField.get_entities_by_type_and_ref_level(
      all_but_0,all_but_0,MBEDGE,refined_edges
    ); CHKERRQ(ierr);
    Range tets_to_refine;
    int size = ((double)2/3)*ufe.errorMap.size();
    for(
      map<double,EntityHandle>::iterator mit = ufe.errorMap.begin();
      mit!=ufe.errorMap.end();
      mit++
    ) {
      // cerr << mit->first << " " << mit->second << endl;
      if((size--)>0) continue;
      tets_to_refine.insert(mit->second);
    }
    Range tets_to_refine_edges;
    rval = mField.get_moab().get_adjacencies(
      tets_to_refine,1,false,tets_to_refine_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    // cerr << tets_to_refine_edges << endl;
    ierr = mField.query_interface(refine_ptr); CHKERRQ(ierr);
    refined_edges.merge(tets_to_refine_edges);
    // cerr << refined_edges << endl;
    for(int ll = 0;ll!=nb_levels;ll++) {
      ierr = refine_ptr->add_verices_in_the_middel_of_edges(refined_edges,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      Range tets;
      ierr = mField.get_entities_by_type_and_ref_level(
        BitRefLevel().set(ll),BitRefLevel().set(),MBTET,tets
      ); CHKERRQ(ierr);
      ierr = refine_ptr->refine_TET(tets,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      // cerr << BitRefLevel().set(ll+1) << endl;
    }

    // {
    //   EntityHandle out_meshset_tet;
    //   rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset_tet); CHKERRQ_MOAB(rval);
    //   // cerr << BitRefLevel().set(nb_levels) << endl;
    //   ierr = mField.get_entities_by_type_and_ref_level(
    //     BitRefLevel().set(nb_levels),BitRefLevel().set(),MBTET,out_meshset_tet
    //   ); CHKERRQ(ierr);
    //   if(mField.getCommRank()==0) {
    //     rval = mField.get_moab().write_file("ref_mesh.vtk","VTK","",&out_meshset_tet,1); CHKERRQ_MOAB(rval);
    //   }
    // }

    PetscFunctionReturn(0);
  }

  PetscErrorCode updateMeshsetsFieldsAndElements(const int nb_levels,const int order) {
    MoABErrorCode rval;
    PetscErrorCode ierr;
    BitRefLevel ref_level;
    PetscFunctionBegin;
    ref_level.set(nb_levels);
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it)) {
      EntityHandle meshset = it->meshset;
      ierr = mField.update_meshset_by_entities_children(meshset,ref_level,meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(meshset,ref_level,meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(meshset,ref_level,meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(meshset,ref_level,meshset,MBTET,true); CHKERRQ(ierr);
    }
    // update fields and elements
    EntityHandle out_meshset_tet;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset_tet); CHKERRQ_MOAB(rval);
    // cerr << BitRefLevel().set(nb_levels) << endl;
    ierr = mField.get_entities_by_type_and_ref_level(
      BitRefLevel().set(nb_levels),BitRefLevel().set(),MBTET,out_meshset_tet
    ); CHKERRQ(ierr);
    //add entities to field
    ierr = mField.add_ents_to_field_by_TETs(out_meshset_tet,"FLUXES"); CHKERRQ(ierr);
    ierr = mField.add_ents_to_field_by_TETs(out_meshset_tet,"VALUES"); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTET,"FLUXES",order+1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"FLUXES",order+1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTET,"VALUES",order); CHKERRQ(ierr);
    Range ref_tets;
    rval = mField.get_moab().get_entities_by_type(out_meshset_tet,MBTET,ref_tets); CHKERRQ_MOAB(rval);
    //add entities to finite elements
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_THERMALSET,it)) {
      Mat_Thermal temp_data;
      ierr = it->getAttributeDataStructure(temp_data); CHKERRQ(ierr);
      setOfBlocks[it->getMeshsetId()].cOnductivity = temp_data.data.Conductivity;
      setOfBlocks[it->getMeshsetId()].cApacity = temp_data.data.HeatCapacity;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->getMeshsetId()].tEts,true); CHKERRQ_MOAB(rval);
      setOfBlocks[it->getMeshsetId()].tEts = intersect(ref_tets,setOfBlocks[it->getMeshsetId()].tEts);
      ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->getMeshsetId()].tEts,"ULTRAWEAK"); CHKERRQ(ierr);
      // {
      //   if(mField.getCommRank()==0) {
      //     EntityHandle out_meshset_tet;
      //     rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset_tet); CHKERRQ_MOAB(rval);
      //     rval = mField.get_moab().add_entities(out_meshset_tet,setOfBlocks[it->getMeshsetId()].tEts); CHKERRQ_MOAB(rval);
      //     rval = mField.get_moab().write_file("block_mesh.vtk","VTK","",&out_meshset_tet,1); CHKERRQ_MOAB(rval);
      //   }
      // }
    }
    rval = mField.get_moab().delete_entities(&out_meshset_tet,1); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  //Create mofem interface
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  // Add meshsets with material and boundary conditions
  MeshsetsManager *meshsets_manager_ptr;
  ierr = m_field.query_interface(meshsets_manager_ptr); CHKERRQ(ierr);
  ierr = meshsets_manager_ptr->setMeshsetFromFile(); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,it)) {
    cerr << *it << endl;
  }

  //set entities bit level
  BitRefLevel ref_level;
  ref_level.set(0);
  ierr = m_field.seed_ref_level_3D(0,ref_level); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 0;
  }

  //finite elements

  BcFluxMap bc_flux_map;
  MyUltraWeakFE ufe(m_field,bc_flux_map);

  ierr = ufe.addFields("VALUES","FLUXES",order); CHKERRQ(ierr);
  ierr = ufe.addFiniteElements("FLUXES","VALUES"); CHKERRQ(ierr);
  // Set boundary conditions
  ierr = ufe.addBoundaryElements(ref_level);
  ierr = ufe.buildProblem(ref_level); CHKERRQ(ierr);
  ierr = ufe.createMatrices(); CHKERRQ(ierr);
  ierr = ufe.solveProblem(); CHKERRQ(ierr);
  ierr = ufe.calculateResidual(); CHKERRQ(ierr);
  ierr = ufe.evaluateError(); CHKERRQ(ierr);
  ierr = ufe.destroyMatrices(); CHKERRQ(ierr);
  ierr = ufe.postProc("out_0.h5m"); CHKERRQ(ierr);

  {
    const int nb_levels = 1;
    ierr = ufe.squashBits(); CHKERRQ(ierr);
    ierr = ufe.refienMesh(ufe,nb_levels); CHKERRQ(ierr);
    ierr = ufe.updateMeshsetsFieldsAndElements(nb_levels,order); CHKERRQ(ierr);
    ref_level = BitRefLevel().set(nb_levels);
    bc_flux_map.clear();
    ierr = ufe.addBoundaryElements(ref_level);
    ierr = ufe.buildProblem(ref_level); CHKERRQ(ierr);
    // {
    //   EntityHandle out_meshset;
    //   rval = m_field.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERRQ_MOAB(rval);
    //   ierr = m_field.get_problem_finite_elements_entities("ULTRAWEAK","ULTRAWEAK",out_meshset); CHKERRQ(ierr);
    //   if(m_field.getCommRank()==0) {
    //     rval = m_field.get_moab().write_file("vol_mesh.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    //   }
    // }
    // {
    //   EntityHandle out_meshset;
    //   rval = m_field.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERRQ_MOAB(rval);
    //   ierr = m_field.get_problem_finite_elements_entities("ULTRAWEAK","ULTRAWEAK_BCVALUE",out_meshset); CHKERRQ(ierr);
    //   if(m_field.getCommRank()==0) {
    //     rval = m_field.get_moab().write_file("bc_nat_mesh.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    //   }
    // }
    // {
    //   EntityHandle out_meshset;
    //   rval = m_field.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERRQ_MOAB(rval);
    //   ierr = m_field.get_problem_finite_elements_entities("ULTRAWEAK","ULTRAWEAK_BCFLUX",out_meshset); CHKERRQ(ierr);
    //   if(m_field.getCommRank()==0) {
    //     rval = m_field.get_moab().write_file("bc_ess_mesh.vtk","VTK","",&out_meshset,1); CHKERRQ_MOAB(rval);
    //   }
    // }
    ierr = ufe.createMatrices(); CHKERRQ(ierr);
    ierr = ufe.solveProblem(); CHKERRQ(ierr);
    ierr = ufe.calculateResidual(); CHKERRQ(ierr);
    ierr = ufe.evaluateError(); CHKERRQ(ierr);
    ierr = ufe.destroyMatrices(); CHKERRQ(ierr);
    ierr = ufe.postProc(
      static_cast<std::ostringstream&>
      (std::ostringstream().seekp(0) << "out_" << nb_levels << ".h5m").str()
    ); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
