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

  try {

  enum bases {
    H1TRI,
    HCURLTRI,
    LASTOP
  };

  const char *list[] = {
    "h1tri",
    "hcurltri"
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
    // 0,0,0,
    // 1,0,0,
    // 0,1,0
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
  ierr = m_field.seed_ref_level_2D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FIELD",space,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

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
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTRI,"FIELD"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTRI,"TRI_FE"); CHKERRQ(ierr);

  //set app. order
  int order = 5;
  if(space == H1) {
    ierr = m_field.set_field_order(root_set,MBTRI,"FIELD",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD",1); CHKERRQ(ierr);
  }
  if(space == HCURL) {
    ierr = m_field.set_field_order(root_set,MBTRI,"FIELD",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD",order); CHKERRQ(ierr);
  }

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRQ(ierr);

  //partition
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRQ(ierr);

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
      mySplit << "type " << type << " side " << side << endl;

      if(data.getFieldDofs()[0]->getSpace()==H1) {

        // mySplit << std::fixed << data.getN() << std::endl;
        // mySplit << std::fixed << data.getDiffN() << std::endl;

        const int nb_dofs = data.getN().size2();
        for(int dd = 0;dd!=nb_dofs;dd++) {
          const double dksi = (data.getN()(1,dd)-data.getN()(0,dd))/(eps);
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

        FTensor::Tensor1<double*,3> base_ksi_m(
          &data.getN()(0,HCURL0),&data.getN()(0,HCURL1),&data.getN()(0,HCURL2),3
        );
        FTensor::Tensor1<double*,3> base_ksi_p(
          &data.getN()(1,HCURL0),&data.getN()(1,HCURL1),&data.getN()(1,HCURL2),3
        );
        FTensor::Tensor1<double*,3> base_eta_m(
          &data.getN()(2,HCURL0),&data.getN()(2,HCURL1),&data.getN()(2,HCURL2),3
        );
        FTensor::Tensor1<double*,3> base_eta_p(
          &data.getN()(3,HCURL0),&data.getN()(3,HCURL1),&data.getN()(3,HCURL2),3
        );

        // cerr << data.getN() << endl;
        // cerr << data.getDiffN() << endl;

        FTensor::Tensor2<double*,3,2> diff_base(
          &data.getDiffN()(4,HCURL0_0),&data.getDiffN()(4,HCURL0_1),
          &data.getDiffN()(4,HCURL1_0),&data.getDiffN()(4,HCURL1_1),
          &data.getDiffN()(4,HCURL2_0),&data.getDiffN()(4,HCURL2_1),6
        );

        FTensor::Index<'i',3> i;
        FTensor::Number<0> N0;
        FTensor::Number<1> N1;

        const int nb_dofs = data.getN().size2()/3;
        for(int dd = 0;dd!=nb_dofs;dd++) {

          mySplit << "MoFEM " << diff_base(0,0) << " " << diff_base(1,0) << " " << diff_base(2,0) << endl;
          mySplit << "MoFEM " << diff_base(0,1) << " " << diff_base(1,1) << " " << diff_base(2,1) << endl;


          FTensor::Tensor1<double,3> dksi;
          dksi(i) = (base_ksi_p(i)-base_ksi_m(i))/(eps);
          FTensor::Tensor1<double,3> deta;
          deta(i) = (base_eta_p(i)-base_eta_m(i))/(2*eps);
          mySplit << "Finite difference dKsi " << dksi(0) << "  " << dksi(1) << " " << dksi(2) << endl;
          mySplit << "Finite diffrence dEta " << deta(0) << "  " << deta(1) << " " << deta(2) << endl;

          dksi(i) -= diff_base(i,N0);
          deta(i) -= diff_base(i,N1);
          if(sqrt(dksi(i)*dksi(i))>eps_diff) {
            // mySplit << "KSI ERROR\n";
            SETERRQ2(
              PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,
              "%s inconsistent dKsi derivative  for type %d",
              FieldSpaceNames[data.getFieldDofs()[0]->getSpace()],
              type
            );
          } else {
            mySplit << "OK" << std::endl;
          }
          if(sqrt(deta(i)*deta(i))>eps_diff) {
            // mySplit << "ETA ERROR\n";
            SETERRQ2(
              PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,
              "%s inconsistent dEta derivative for type %d",
              FieldSpaceNames[data.getFieldDofs()[0]->getSpace()],
              type
            );
          } else {
            mySplit << "OK" << std::endl;
          }

          ++base_ksi_m;
          ++base_ksi_p;
          ++base_eta_m;
          ++base_eta_p;
          ++diff_base;
        }
      }

      mySplit << endl;

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
  tri_fe.getOpPtrVector().push_back(new OpCalculateInvJacForFace(inv_jac));
  if(space==H1) {
    tri_fe.getOpPtrVector().push_back(new OpSetInvJacH1ForFace(inv_jac));
  }
  if(space==HCURL) {
    tri_fe.getOpPtrVector().push_back(new OpSetInvJacHcurlFace(inv_jac));
  }
  tri_fe.getOpPtrVector().push_back(new OpCheckingDirevatives(my_split));
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TRI_FE",tri_fe);  CHKERRQ(ierr);
  // cerr << inv_jac << endl;


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
}
