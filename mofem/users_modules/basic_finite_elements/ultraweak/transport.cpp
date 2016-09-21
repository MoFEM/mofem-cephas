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

  // for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,it)) {
  //   cerr << *it << endl;
  // }

  //set entities bit level
  BitRefLevel bit_ref_level;
  bit_ref_level.set(0);
  ierr = m_field.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

  //finite elements

  BcFluxMap bc_flux_map;

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


  };

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }

  MyUltraWeakFE ufe(m_field,bc_flux_map);
  ierr = ufe.addFields("VALUES","FLUXES",order); CHKERRQ(ierr);
  ierr = ufe.addFiniteElements("FLUXES","VALUES"); CHKERRQ(ierr);

  // Set boundary conditions
  {
    Range tets;
    EntityHandle fe_meshset = m_field.get_finite_element_meshset("ULTRAWEAK"); CHKERRQ(ierr);
    rval = m_field.get_moab().get_entities_by_type(fe_meshset,MBTET,tets); CHKERR_MOAB(rval);
    Skinner skin(&moab);
    Range skin_faces; // skin faces from 3d ents
    rval = skin.find_skin(0,tets,false,skin_faces); CHKERR_MOAB(rval);
    // note: what is essential (dirichlet) is natural (neumann) for ultra weak compared to classical FE
    Range natural_bc;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|TEMPERATURESET,it)) {

      Range tris;
      ierr = it->getMeshsetIdEntitiesByDimension(m_field.get_moab(),2,tris,true); CHKERRQ(ierr);
      natural_bc.insert(tris.begin(),tris.end());

    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|HEATFLUXSET,it)) {
      HeatFluxCubitBcData mydata;
      ierr = it->getBcDataStructure(mydata); CHKERRQ(ierr);
      if(mydata.data.flag1==1) {
        Range tris;
        ierr = it->getMeshsetIdEntitiesByDimension(m_field.get_moab(),2,tris,true); CHKERRQ(ierr);
        bc_flux_map[it->getMeshsetId()].eNts = tris;
        bc_flux_map[it->getMeshsetId()].fLux = mydata.data.value1;
        // cerr << bc_flux_map[it->getMeshsetId()].eNts << endl;
        // cerr << bc_flux_map[it->getMeshsetId()].fLux << endl;
      }
    }
    Range essential_bc = subtract(skin_faces,natural_bc);
    Range bit_tris;
    ierr = m_field.get_entities_by_type_and_ref_level(bit_ref_level,BitRefLevel().set(),MBTRI,bit_tris);
    essential_bc = intersect(bit_tris,essential_bc);
    natural_bc = intersect(bit_tris,natural_bc);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(essential_bc,"ULTRAWEAK_BCFLUX"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(natural_bc,"ULTRAWEAK_BCVALUE"); CHKERRQ(ierr);
    // ierr = m_field.add_ents_to_finite_element_by_TRIs(skin_faces,"ULTRAWEAK_BCVALUE"); CHKERRQ(ierr);
  }

  ierr = ufe.buildProblem(bit_ref_level); CHKERRQ(ierr);
  ierr = ufe.createMatrices(); CHKERRQ(ierr);
  ierr = ufe.solveProblem(); CHKERRQ(ierr);
  ierr = ufe.calculateResidual(); CHKERRQ(ierr);
  ierr = ufe.destroyMatrices(); CHKERRQ(ierr);
  ierr = ufe.postProc(); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
