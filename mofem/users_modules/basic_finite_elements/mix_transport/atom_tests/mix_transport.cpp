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
#include <MixTransportElement.hpp>
using namespace MoFEM;
using namespace MixTransport;

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

static char help[] = "...\n\n";

/** define sources and other stuff
  *
  * MixTransportElement is a class collecting functions, operators and
  * data for ultra week implementation of transport element. See there to
  * learn how elements are created or how operators look like.
  *
  * Some methods in MixTransportElement are abstract, f.e. user need to
  * implement own source therm.
  *
  */
struct MyTransport: public MixTransportElement {

  MyTransport(MoFEM::Interface &m_field): MixTransportElement(m_field) {};

  PetscErrorCode getSource(EntityHandle ent,const double x,const double y,const double z,double &flux) {
    PetscFunctionBegin;
    //double d = sqrt(x*x+y*y+z*z);
    flux = 1;//-pow(d,5./4.);
    PetscFunctionReturn(0);
  }

  PetscErrorCode getBcOnValues(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &value) {
    PetscFunctionBegin;
    value = 1;
    PetscFunctionReturn(0);
  }

  PetscErrorCode getBcOnFluxes(
    const EntityHandle ent,
    const double x,const double y,const double z,
    double &flux) {
    PetscFunctionBegin;
    flux = 0;
    PetscFunctionReturn(0);
  }

};

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

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
    option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
    BARRIER_RANK_START(pcomm)
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
    BARRIER_RANK_END(pcomm)

    //Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    //set entities bit level
    BitRefLevel ref_level;
    ref_level.set(0);
    ierr = m_field.seed_ref_level_3D(0,ref_level); CHKERRQ(ierr);

    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    PetscInt order;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      order = 2;
    }

    //finite elements
    MyTransport ufe(m_field);

    ierr = ufe.addFields("VALUES","FLUXES",order); CHKERRQ(ierr);
    ierr = ufe.addFiniteElements("FLUXES","VALUES"); CHKERRQ(ierr);

    // Set boundary conditions
    Range tets;
    ierr = m_field.get_entities_by_type_and_ref_level(BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
    Skinner skin(&moab);
    Range skin_faces; // skin faces from 3d ents
    rval = skin.find_skin(0,tets,false,skin_faces); CHKERRQ_MOAB(rval);
    ierr = m_field.add_ents_to_finite_element_by_type(skin_faces,MBTRI,"MIX_BCVALUE"); CHKERRQ(ierr);

    ierr = ufe.buildProblem(ref_level); CHKERRQ(ierr);
    ierr = ufe.createMatrices(); CHKERRQ(ierr);
    ierr = ufe.solveLinearProblem(); CHKERRQ(ierr);
    ierr = ufe.calculateResidual(); CHKERRQ(ierr);
    ierr = ufe.evaluateError(); CHKERRQ(ierr);

    double nrm2_F;
    ierr = VecNorm(ufe.F,NORM_2,&nrm2_F); CHKERRQ(ierr);
    // PetscPrintf(PETSC_COMM_WORLD,"nrm2_F = %6.4e\n",nrm2_F);
    const double eps = 1e-8;
    if(nrm2_F > eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"problem with residual");
    }

    ierr = ufe.destroyMatrices(); CHKERRQ(ierr);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
