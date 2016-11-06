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

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

static const double eps = 1e-6;
static const double eps_diff = 1e-4;

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  enum bases {
    H1TRI,
    HCURLTRI,
    LASTOP
  };

  const char *list[] = {
    "h1tri",
  };

  PetscBool flg;
  PetscInt choise_value = H1TRI;
  ierr = PetscOptionsGetEList(
    PETSC_NULL,NULL,"-base",list,LASTOP,&choise_value,&flg
  ); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"base not set");
  }

  FieldSpace space = LASTSPACE;
  if(choise_value == H1TRI) {
    space = H1;
  } else if(choise_value == HCURLTRI) {
    space = HCURL;
  }

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  // create one tet
  double tri_coords[] = {
    0,0,0,
    .5,0,0,
    0,1.,0
  };
  EntityHandle nodes[3];
  for(int nn = 0;nn<3;nn++) {
    rval = moab.create_vertex(&tri_coords[3*nn],nodes[nn]); CHKERRQ_MOAB(rval);
  }
  EntityHandle tri;
  rval = moab.create_element(MBTRI,nodes,3,tri); CHKERRQ_MOAB(rval);
  // Create adjacencies entitities
  Range adj;
  rval = moab.get_adjacencies(&tri,1,1,true,adj); CHKERRQ(ierr);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_2D(0,bit_level0,10); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FIELD",space,1); CHKERRQ(ierr);

  //FE TET
  ierr = m_field.add_finite_element("TRI_FE"); CHKERRQ(ierr);
  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TRI_FE","FIELD"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TRI_FE","FIELD"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TRI_FE","FIELD"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TRI_FE"); CHKERRQ(ierr);
  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TRIs(root_set,"FIELD"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TRIs(root_set,"TRI_FE"); CHKERRQ(ierr);

  //set app. order
  int order = 4;
  if(space == H1) {
    ierr = m_field.set_field_order(root_set,MBTRI,"FIELD",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD",1); CHKERRQ(ierr);
  }
  if(space == HCURL) {
    ierr = m_field.set_field_order(root_set,MBTRI,"FIELD",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD",order); CHKERRQ(ierr);
  }

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  typedef tee_device<std::ostream, std::ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  std::ofstream ofs("forces_and_sources_checking_direvatives.txt");
  TeeDevice my_tee(std::cout, ofs);
  TeeStream my_split(my_tee);

  struct OpCheckingDirevatives: public FaceElementForcesAndSourcesCore::UserDataOperator {

    TeeStream &mySplit;
    OpCheckingDirevatives(TeeStream &my_split):
    FaceElementForcesAndSourcesCore::UserDataOperator("FIELD",UserDataOperator::OPROW),
    mySplit(my_split) {
    }

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      if(data.getFieldDofs()[0]->getSpace()==H1) {

        // mySplit << std::fixed << data.getN() << std::endl;
        // mySplit << std::fixed << data.getDiffN() << std::endl;

        const int nb_dofs = data.getN().size2();
        for(int dd = 0;dd!=nb_dofs;dd++) {
          const double dksi = (data.getN()(1,dd)-data.getN()(0,dd))/eps;
          const double deta = (data.getN()(3,dd)-data.getN()(2,dd))/(2*eps);
          mySplit << "DKsi " << dksi << std::endl;
          mySplit << "DEta " << deta << std::endl;
          mySplit << "diffN " << data.getDiffN()(4,2*dd+0) << std::endl;
          mySplit << "diffN " << data.getDiffN()(4,2*dd+1) << std::endl;
          if(fabs(dksi-data.getDiffN()(4,2*dd+0))>eps_diff) {
            SETERRQ(
              PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"H1 inconsistent dKsi derivative"
            );
          }
          if(fabs(deta-data.getDiffN()(4,2*dd+1))>eps_diff) {
            SETERRQ(
              PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"H1 inconsistent dEta derivative"
            );
          }
        }

      }

      if(data.getFieldDofs()[0]->getSpace()==HCURL) {

        SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"Test not implemented yet");

      }

      PetscFunctionReturn(0);
    }

  };

  struct MyFE: public FaceElementForcesAndSourcesCore {

    MyFE(MoFEM::Interface &m_field): FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return -1; };

    PetscErrorCode setGaussPts(int order) {
      PetscFunctionBegin;

      const double ksi = G_TRI_X1[0];
      const double eta = G_TRI_Y1[0];

      gaussPts.resize(3,5);
      gaussPts.clear();

      gaussPts(0,0) = ksi-eps;
      gaussPts(1,0) = eta;
      gaussPts(0,1) = ksi+eps;
      gaussPts(1,1) = eta;

      gaussPts(0,2) = ksi;
      gaussPts(1,2) = eta-eps;
      gaussPts(0,3) = ksi;
      gaussPts(1,3) = eta+eps;

      gaussPts(0,4) = ksi;
      gaussPts(1,4) = eta;

      for(int ii = 0;ii!=gaussPts.size2();ii++) {
        gaussPts(2,ii) = 1;
      }

      // cerr << gaussPts << endl;

      PetscFunctionReturn(0);
    }

  };

  MyFE tri_fe(m_field);

  MatrixDouble inv_jac;
  tri_fe.getOpPtrVector().push_back(new OpCalculateInvJacForFace("FIELD",inv_jac));
  tri_fe.getOpPtrVector().push_back(new OpSetInvJacH1ForFace("FIELD",inv_jac));
  tri_fe.getOpPtrVector().push_back(new OpCheckingDirevatives(my_split));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TRI_FE",tri_fe);  CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);
}
