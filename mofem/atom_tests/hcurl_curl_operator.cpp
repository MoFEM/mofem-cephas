/**
 * \file hcurl_curl_operator.cpp
 * \brief Testich curl-curl operator by applying Stokes theorem
 * \example hcurl_curl_operator.cpp
 *
 * Using Basic interface calculate the curl of base functions, and integral of
 * the vector tangent vector with normal on the boundary. Since the h-curl space
 * is used, volume integral and boundary integral should give the same result,
 * as a result, as we are applying Stokes theorem on h-curl space.
 *
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

static char help[] = "...\n\n";

struct OpTetCurl : public VolumeElementForcesAndSourcesCore::UserDataOperator {

  FTensor::Tensor1<double, 3> &tCurl;
  OpTetCurl(FTensor::Tensor1<double, 3> &t_curl)
      : VolumeElementForcesAndSourcesCore::UserDataOperator(
            "HCURL", UserDataOperator::OPROW),
        tCurl(t_curl) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

struct OpFacesRot : public FaceElementForcesAndSourcesCore::UserDataOperator {

  FTensor::Tensor1<double, 3> &tCurl;
  OpFacesRot(FTensor::Tensor1<double, 3> &t_curl)
      : FaceElementForcesAndSourcesCore::UserDataOperator(
            "HCURL", UserDataOperator::OPROW),
        tCurl(t_curl) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    enum bases { AINSWORTH, DEMKOWICZ, LASTOP };

    const char *list[] = {"ainsworth", "demkowicz"};

    PetscBool flg;
    PetscInt choise_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list, LASTOP,
                                &choise_value, &flg);
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    }

    PetscBool ho_geometry = PETSC_FALSE;
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-ho_geometry", &ho_geometry,
                               PETSC_NULL);

    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Create MoAB
    moab::Core mb_instance;              ///< database
    moab::Interface &moab = mb_instance; ///< interface

    // Create MoFEM
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    auto basic_interface = m_field.getInterface<Basic>();
    CHKERR basic_interface->getOptions();
    CHKERR basic_interface->loadFile("");

    // fields
    switch (choise_value) {
    case AINSWORTH:
      CHKERR basic_interface->addDomainField("HCURL", HCURL,
                                             AINSWORTH_LEGENDRE_BASE, 1);
      CHKERR basic_interface->addBoundaryField("HCURL", HCURL,
                                               AINSWORTH_LEGENDRE_BASE, 1);
      break;
    case DEMKOWICZ:
      CHKERR basic_interface->addDomainField("HCURL", HCURL,
                                             DEMKOWICZ_JACOBI_BASE, 1);
      CHKERR basic_interface->addBoundaryField("HCURL", HCURL,
                                               DEMKOWICZ_JACOBI_BASE, 1);
      break;
    }

    if (ho_geometry == PETSC_TRUE)
      CHKERR basic_interface->addDataField("MESH_NODE_POSITIONS", H1,
                                           AINSWORTH_LEGENDRE_BASE, 3);

    constexpr int order = 5;
    CHKERR basic_interface->setFieldOrder("HCURL", order);
    if (ho_geometry == PETSC_TRUE)
      CHKERR basic_interface->setFieldOrder("MESH_NODE_POSITIONS", 2);
    CHKERR basic_interface->setUp();

    auto integration_rule = [](int, int, int p_data) { return 2 * p_data; };
    basic_interface->getOpDomainRhsRuleHook() = integration_rule;
    basic_interface->getOpBoundaryRhsRuleHook() = integration_rule;

    FTensor::Tensor1<double, 3> t_curl_vol;
    FTensor::Tensor1<double, 3> t_curl_skin;

    basic_interface->getOpDomainRhsPipeline().push_back(
        new OpTetCurl(t_curl_vol));
    basic_interface->getOpBoundaryRhsPipeline().push_back(
        new OpFacesRot(t_curl_skin));

    FTensor::Index<'i', 3> i;
    t_curl_vol(i) = 0;
    t_curl_skin(i) = 0;

    // project geometry form 10 node tets on higher order approx. functions
    if (ho_geometry == PETSC_TRUE) {
      Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
      CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);
    }

    // Run pipelines on mesh
    CHKERR basic_interface->loopFiniteElements();

    std::cout.precision(12);

    std::cout << "curl_vol " << t_curl_vol << std::endl;
    std::cout << "curl_skin " << t_curl_skin << std::endl;

    t_curl_vol(i) -= t_curl_skin(i);
    double nrm2 = sqrt(t_curl_vol(i) * t_curl_vol(i));

    constexpr double eps = 1e-8;
    if (fabs(nrm2) > eps) 
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "Curl operator not passed test\n");
    
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}

MoFEMErrorCode OpTetCurl::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (data.getFieldData().size() == 0)
    MoFEMFunctionReturnHot(0);

  const unsigned int nb_gauss_pts = data.getDiffN().size1();
  const unsigned int nb_dofs = data.getFieldData().size();

  MatrixDouble curl_mat;
  FTensor::Index<'i', 3> i;

  unsigned int gg = 0;
  for (; gg < nb_gauss_pts; gg++) {
    double w = getGaussPts()(3, gg) * getVolume();
    if (getHoGaussPtsDetJac().size() == nb_gauss_pts) {
      // if ho geometry is given
      w *= getHoGaussPtsDetJac()(gg);
    }
    CHKERR getCurlOfHCurlBaseFunctions(side, type, data, gg, curl_mat);
    FTensor::Tensor1<double *, 3> t_curl(&curl_mat(0, 0), &curl_mat(0, 1),
                                         &curl_mat(0, 2), 3);
    for (unsigned int dd = 0; dd != nb_dofs; dd++) {
      tCurl(i) += w * t_curl(i);
      ++t_curl;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpFacesRot::doWork(int side, EntityType type,
                                  DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);
  int nb_gauss_pts = data.getN().size1();

  auto t_curl_base = data.getFTensor1N<3>();
  // double area = getArea();
  double n0 = getNormal()[0] * 0.5;
  double n1 = getNormal()[1] * 0.5;
  double n2 = getNormal()[2] * 0.5;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  for (int gg = 0; gg < nb_gauss_pts; gg++) {
    for (int dd = 0; dd < nb_dofs; dd++) {
      double w = getGaussPts()(2, gg);
      if (getNormalsAtGaussPts().size1() == (unsigned int)nb_gauss_pts) {
        n0 = getNormalsAtGaussPts(gg)[0] * 0.5;
        n1 = getNormalsAtGaussPts(gg)[1] * 0.5;
        n2 = getNormalsAtGaussPts(gg)[2] * 0.5;
      }

      tCurl(0) += (n1 * t_curl_base(2) - n2 * t_curl_base(1)) * w;
      tCurl(1) += (n2 * t_curl_base(0) - n0 * t_curl_base(2)) * w;
      tCurl(2) += (n0 * t_curl_base(1) - n1 * t_curl_base(0)) * w;

      ++t_curl_base;
    }
  }

  MoFEMFunctionReturn(0);
}