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

// namespace bio = boost::iostreams;
// using bio::tee_device;
// using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //create one tet
  double tet_coords[] = {
    0,0,0,
    1.,0,0,
    0,1.,0,
    0,0,1.
  };

  EntityHandle nodes[4];
  for(int nn = 0;nn<4;nn++) {
    rval = moab.create_vertex(&tet_coords[3*nn],nodes[nn]); CHKERRQ_MOAB(rval);
  }

  EntityHandle tet;
  rval = moab.create_element(MBTET,nodes,4,tet); CHKERRQ_MOAB(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;
  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //fields
  ierr = m_field.add_field("HCURL",HCURL,1); CHKERRQ(ierr);
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"HCURL"); CHKERRQ(ierr);
  //set app. order
  int order = 4;
  ierr = m_field.set_field_order(root_set,MBTET,"HCURL",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"HCURL",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"HCURL",order); CHKERRQ(ierr);
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //finite elements
  ierr = m_field.add_finite_element("TET_FE"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("SKIN_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TET_FE","HCURL"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TET_FE","HCURL"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TET_FE","HCURL"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("SKIN_FE","HCURL"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("SKIN_FE","HCURL"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("SKIN_FE","HCURL"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"TET_FE"); CHKERRQ(ierr);
  Range tets;
  ierr = m_field.get_entities_by_type_and_ref_level(
    BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets
  ); CHKERRQ(ierr);
  Skinner skin(&moab);
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0,tets,false,skin_faces); CHKERR_MOAB(rval);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(skin_faces,"SKIN_FE"); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TET_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","SKIN_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  struct OpTetCurl: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    FTensor::Tensor1<double*,3> &cUrl;
    OpTetCurl(FTensor::Tensor1<double*,3> &curl):
    VolumeElementForcesAndSourcesCore::UserDataOperator("HCURL",UserDataOperator::OPROW),
    cUrl(curl) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      int nb_gauss_pts = data.getDiffHcurlN().size1();
      int nb_dofs = data.getFieldData().size();

      MatrixDouble curl_mat;
      FTensor::Index<'i',3> i;

      int gg = 0;
      for(;gg<nb_gauss_pts;gg++) {
        double w = getGaussPts()(3,gg)*getVolume();
        ierr = getCurlOfHCurlBaseFunctions(side,type,data,gg,curl_mat); CHKERRQ(ierr);
        FTensor::Tensor1<double*,3> t_curl(&curl_mat(0,0),&curl_mat(0,1),&curl_mat(0,2),3);
        for(int dd = 0;dd!=nb_dofs;dd++) {
          cUrl(i) += w*t_curl(i);
          ++t_curl;
        }
      }

      PetscFunctionReturn(0);
    }

  };

  struct MyFE: public VolumeElementForcesAndSourcesCore {

    MyFE(MoFEM::Interface &m_field):
    VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return order; }; //order/2; };

  };

  struct MyTriFE: public FaceElementForcesAndSourcesCore {

    MyTriFE(MoFEM::Interface &m_field):
    FaceElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return order; };//2*order; }; //order/2; };

  };

  struct OpFacesRot: public FaceElementForcesAndSourcesCore::UserDataOperator {

    FTensor::Tensor1<double*,3> &cUrl;
    OpFacesRot(FTensor::Tensor1<double*,3> &curl):
    FaceElementForcesAndSourcesCore::UserDataOperator("HCURL",UserDataOperator::OPROW),
    cUrl(curl) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;

      int nb_dofs = data.getFieldData().size();
      if(nb_dofs==0) PetscFunctionReturn(0);
      int nb_gauss_pts = data.getHcurlN().size1();

      FTensor::Tensor1<double*,3> t_curl_base = data.getFTensor1HcurlN<3>();
      // double area = getArea();
      double n0 = getNormal()[0]*0.5;
      double n1 = getNormal()[1]*0.5;
      double n2 = getNormal()[2]*0.5;

      // cout << area << " " << norm_2(getNormal()) << endl;

      // for(int ii = 0;ii!=3;ii++) {
      //   for(int jj = 0;jj!=3;jj++) {
      //     cout << t_spin_normal(ii,jj) << " ";
      //   }
      //   cout << endl;
      // }
      // cout << endl;

      FTensor::Index<'i',3> i;
      FTensor::Index<'j',3> j;

      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        for(int dd = 0;dd<nb_dofs;dd++) {
          double w = getGaussPts()(2,gg);
          cUrl(0) += (n1*t_curl_base(2)-n2*t_curl_base(1))*w;
          cUrl(1) += (n2*t_curl_base(0)-n0*t_curl_base(2))*w;
          cUrl(2) += (n0*t_curl_base(1)-n1*t_curl_base(0))*w;
          ++t_curl_base;
        }
      }

      PetscFunctionReturn(0);
    }

  };

  VectorDouble curl_vol(3);
  VectorDouble curl_skin(3);
  FTensor::Tensor1<double*,3> t_curl_vol(&curl_vol[0],&curl_vol[1],&curl_vol[2]);
  FTensor::Tensor1<double*,3> t_curl_skin(&curl_skin[0],&curl_skin[1],&curl_skin[2]);
  curl_vol.clear();
  curl_skin.clear();

  MyFE tet_fe(m_field);
  tet_fe.getOpPtrVector().push_back(new OpTetCurl(t_curl_vol));
  MyTriFE skin_fe(m_field);
  skin_fe.getOpPtrVector().push_back(new OpFacesRot(t_curl_skin));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TET_FE",tet_fe);  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","SKIN_FE",skin_fe);  CHKERRQ(ierr);

  std::cout.precision(12);

  std::cout << "curl_vol " << curl_vol << std::endl;
  std::cout << "curl_skin " << curl_skin << std::endl;

  FTensor::Index<'i',3> i;
  t_curl_vol(i)-=t_curl_skin(i);
  double nrm2 = sqrt(t_curl_vol(i)*t_curl_vol(i));

  const double eps = 1e-8;
  if(fabs(nrm2)>eps) {
     SETERRQ(
       PETSC_COMM_SELF,
       MOFEM_ATOM_TEST_INVALID,
       "Curl operator not passed test\n"
     );
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
}
