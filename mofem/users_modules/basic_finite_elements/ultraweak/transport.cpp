/** \file transport.cpp
\brief Example implementation of transport problem using ultra-week formulation

\todo Should be implemented and tested problem from this article
Demkowicz, Leszek, and Jayadeep Gopalakrishnan. "Analysis of the DPG method for
the Poisson equation." SIAM Journal on Numerical Analysis 49.5 (2011):
1788-1809.

\ingroup mofem_ultra_weak_transport_elem
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

/**
 * Data structure to pass information between function evaluating boundary
 * values and fluxes and generic data structures for boundary conditions on
 * meshsets.
 */
struct BcFluxData {
  Range eNts;
  double fLux;
};
typedef map<int,BcFluxData> BcFluxMap;

/** \brief Application of ultraweak data structure
  *
  * UltraWeakTransportElement is a class collecting functions, operators and
  * data for ultra week implementation of transport element. See there to
  * learn how elements are created or how operators look like.
  *
  * Some methods in UltraWeakTransportElement are abstract, f.e. user need to
  * implement own source therm.

  * \ingroup mofem_ultra_weak_transport_elem
  */
struct ExampleUltraWeak: public UltraWeakTransportElement {

  BcFluxMap &bcFluxMap;
  EntityHandle lastEnt;
  double lastFlux;

  ExampleUltraWeak(MoFEM::Interface &m_field,BcFluxMap &bc_flux_map):
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
  PetscErrorCode getSource(EntityHandle ent,const double x,const double y,const double z,double &flux) {
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

  /**
   * \bried set-up boundary conditions
   * @param  ref_level mesh refinement level

   \note It is assumed that user would like to something non-standard with boundary
   conditions, have a own type of data structures to pass to functions calculating
   values and fluxes on boundary. For example BcFluxMap. That way this function
   is implemented here not in generic class UltraWeakTransportElement.

   * @return           error code
   */
  PetscErrorCode addBoundaryElements(BitRefLevel &ref_level) {
    MoABErrorCode rval;
    PetscErrorCode ierr;
    PetscFunctionBegin;
    Range tets;
    ierr = mField.get_entities_by_type_and_ref_level(ref_level,BitRefLevel().set(),MBTET,tets);
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

  /**
   * \brief Refine mesh
   * @param  ufe       general data structure
   * @param  nb_levels number of refinement levels
   * @param  order     set order of approximation
   * @return           errpr code

   Refinement of could result in distorted mesh, for example, imagine when you
   have two levels of non-uniform refinement. Some tetrahedra on the mesh at
   first refinement instance are only refined by splitting subset of edges on
   it. Then refined child tetrahedra usually will have worse quality than
   quality of parent element. Refining such element in subsequent mesh
   refinement, potentially will deteriorate elements quality even worse. To
   prevent that adding new refinement level, recreate whole hierarchy of meshes.

   Note on subsequent improvement could include refinement of
   tetrahedra from different levels, including initial mesh. So refinement two
   could split elements created during refinement one and also split elements
   from an initial mesh.

   That adding the new refinement level creates refinement hierarchy of meshes
   from a scratch,  not adding to existing one.

   Entities from previous hierarchy are used in that process, but bit levels on
   those entities are squashed.

   */
  PetscErrorCode refineMesh(
    UltraWeakTransportElement &ufe,const int nb_levels,const int order
  ) {
    PetscErrorCode ierr;
    MoABErrorCode rval;
    MeshRefinement *refine_ptr;
    PetscFunctionBegin;
    // get refined edges having child vertex
    const RefEntity_multiIndex *ref_ents_ptr;
    ierr = mField.get_ref_ents(&ref_ents_ptr); CHKERRQ(ierr);
    typedef RefEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type RefEntsByComposite;
    const RefEntsByComposite &ref_ents = ref_ents_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>();
    RefEntsByComposite::iterator rit,hi_rit;
    rit = ref_ents.lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
    hi_rit = ref_ents.upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
    Range refined_edges;
    // thist loop is over vertices which parent is edge
    for(;rit!=hi_rit;rit++) {
      refined_edges.insert((*rit)->getParentEnt()); // get parent edge
    }
    // get tets which has large error
    Range tets_to_refine;
    const double max_error = ufe.errorMap.rbegin()->first;
    // int size = ((double)5/6)*ufe.errorMap.size();
    for(
      map<double,EntityHandle>::iterator mit = ufe.errorMap.begin();
      mit!=ufe.errorMap.end();
      mit++
    ) {
      // cerr << mit->first << " " << mit->second << endl;
      // if((size--)>0) continue;
      if(mit->first<0.25*max_error) continue;
      tets_to_refine.insert(mit->second);
    }
    Range tets_to_refine_edges;
    rval = mField.get_moab().get_adjacencies(
      tets_to_refine,1,false,tets_to_refine_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    refined_edges.merge(tets_to_refine_edges);
    ierr = mField.query_interface(refine_ptr); CHKERRQ(ierr);
    for(int ll = 0;ll!=nb_levels;ll++) {
      Range edges;
      ierr = mField.get_entities_by_type_and_ref_level(
        BitRefLevel().set(ll),BitRefLevel().set(),MBEDGE,edges
      ); CHKERRQ(ierr);
      edges = intersect(edges,refined_edges);
      // add edges to refine at current level edges (some of the where refined before)
      ierr = refine_ptr->add_verices_in_the_middel_of_edges(edges,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      //  get tets at current level
      Range tets;
      ierr = mField.get_entities_by_type_and_ref_level(
        BitRefLevel().set(ll),BitRefLevel().set(),MBTET,tets
      ); CHKERRQ(ierr);
      ierr = refine_ptr->refine_TET(tets,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      ierr = updateMeshsetsFieldsAndElements(ll+1); CHKERRQ(ierr);
    }

    // update fields and elements
    EntityHandle ref_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,ref_meshset); CHKERRQ_MOAB(rval);
    {
      // cerr << BitRefLevel().set(nb_levels) << endl;
      ierr = mField.get_entities_by_type_and_ref_level(
        BitRefLevel().set(nb_levels),BitRefLevel().set(),MBTET,ref_meshset
      ); CHKERRQ(ierr);

      Range ref_tets;
      rval = mField.get_moab().get_entities_by_type(ref_meshset,MBTET,ref_tets); CHKERRQ_MOAB(rval);

      //add entities to field
      ierr = mField.add_ents_to_field_by_TETs(ref_meshset,"FLUXES"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_field_by_TETs(ref_meshset,"VALUES"); CHKERRQ(ierr);
      ierr = mField.set_field_order(0,MBTET,"FLUXES",order+1); CHKERRQ(ierr);
      ierr = mField.set_field_order(0,MBTRI,"FLUXES",order+1); CHKERRQ(ierr);
      ierr = mField.set_field_order(0,MBTET,"VALUES",order); CHKERRQ(ierr);

      // add entities to skeleton
      Range ref_tris;
      ierr = mField.get_entities_by_type_and_ref_level(
        BitRefLevel().set(nb_levels),BitRefLevel().set(),MBTRI,ref_tris
      ); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TRIs(
        ref_tris,"ULTRAWEAK_SKELETON"
      ); CHKERRQ(ierr);

      //add entities to finite elements
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_THERMALSET,it)) {
        Mat_Thermal temp_data;
        ierr = it->getAttributeDataStructure(temp_data); CHKERRQ(ierr);
        setOfBlocks[it->getMeshsetId()].cOnductivity = temp_data.data.Conductivity;
        setOfBlocks[it->getMeshsetId()].cApacity = temp_data.data.HeatCapacity;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->getMeshsetId()].tEts,true); CHKERRQ_MOAB(rval);
        setOfBlocks[it->getMeshsetId()].tEts = intersect(ref_tets,setOfBlocks[it->getMeshsetId()].tEts);
        ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->getMeshsetId()].tEts,"ULTRAWEAK"); CHKERRQ(ierr);
      }
    }
    rval = mField.get_moab().delete_entities(&ref_meshset,1); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  /**
   * \brief Squash bits of entities

   Information about hierarchy of meshses is lost, but entities are not deleted
   from the mesh. After squshing entities bits, new hierarchy can be created.

   * @return error code
   */
  PetscErrorCode squashBits() {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    BitRefLevel all_but_0;
    all_but_0.set(0);
    all_but_0.flip();
    BitRefLevel garbage_bit;
    garbage_bit.set(BITREFLEVEL_SIZE-1); // Garbage level
    const RefEntity_multiIndex *refined_ents_ptr;
    ierr = mField.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
    RefEntity_multiIndex::iterator mit = refined_ents_ptr->begin();
    for(;mit!=refined_ents_ptr->end();mit++) {
      if(mit->get()->getEntType() == MBENTITYSET) continue;
      BitRefLevel bit = mit->get()->getBitRefLevel();
      if((all_but_0&bit)==bit) {
        const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
          mit,RefEntity_change_set_bit(garbage_bit)
        );
      } else {
        const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->modify(
          mit,RefEntity_change_set_bit(BitRefLevel().set(0))
        );
      }
    }
    PetscFunctionReturn(0);
  }

  /**
   * \brief update meshsets with new entities after mesh refinement
   * @param  nb_levels nb_levels
   * @param  order     appropriate order
   * @return           error code
   */
  PetscErrorCode updateMeshsetsFieldsAndElements(const int nb_levels) {
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

  // get file name form command line
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  // create MOAB communicator
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

  PetscPrintf(PETSC_COMM_WORLD,"Read meshsets add added meshsets for bc.cfg\n");
  for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,it)) {
    PetscPrintf(
      PETSC_COMM_WORLD,
      "%s",static_cast<std::ostringstream&>(std::ostringstream().seekp(0) << *it << endl).str().c_str()
    );
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
  ExampleUltraWeak ufe(m_field,bc_flux_map);

  // Initially calculate problem on coarse mesh

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

  int nb_levels = 5; // default number of refinement levels
  // get number of refinement levels form command line
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-nb_levels",&nb_levels,PETSC_NULL); CHKERRQ(ierr);

  // refine mesh, solve problem and do it again until number of refinement levels are exceeded.
  for(int ll = 1;ll!=nb_levels;ll++) {
    const int nb_levels = ll;
    ierr = ufe.squashBits(); CHKERRQ(ierr);
    ierr = ufe.refineMesh(ufe,nb_levels,order); CHKERRQ(ierr);
    ref_level = BitRefLevel().set(nb_levels);
    bc_flux_map.clear();
    ierr = ufe.addBoundaryElements(ref_level);
    ierr = ufe.buildProblem(ref_level); CHKERRQ(ierr);
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
