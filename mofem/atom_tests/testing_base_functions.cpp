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
#include <quad.h>

using namespace MoFEM;

static char help[] = "testing interface inserting algorithm\n\n";

static double sum_matrix(MatrixDouble &m) {
  double s = 0;
  for (unsigned int ii = 0; ii < m.size1(); ii++) {
    for (unsigned int jj = 0; jj < m.size2(); jj++) {
      s += m(ii, jj);
    }
  }
  return s;
}

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    MoFEM::Core core(moab);
    // MoFEM::Interface& m_field = core;

    enum bases {
      LEGENDREPOLYNOMIAL,
      LOBATTOPOLYNOMIAL,
      JACOBIPOLYNOMIAL,
      INTEGRATEDJACOBIPOLYNOMIAL,
      H1TET_AINSWORTH,
      H1TET_BERNSTEIN_BEZIER,
      HDIVTET_AINSWORTH,
      HDIVTET_DEMKOWICZ,
      HCURLTET_AINSWORTH,
      HCURLTET_DEMKOWICZ,
      L2TET,
      H1TRI_AINSWORTH,
      H1TRI_BERNSTEIN_BEZIER,
      H1QUAD,
      HDIVTRI_AINSWORTH,
      HDIVTRI_DEMKOWICZ,
      HCURLTRI_AINSWORTH,
      HCURLTRI_DEMKOWICZ,
      L2TRI,
      H1EDGE_AINSWORTH,
      H1EDGE_BERNSTEIN_BEZIER,
      HCURLEDGE_AINSWORTH,
      HCURLEDGE_DEMKOWICZ,
      H1FLATPRIS,
      H1FATPRISM,
      LASTOP
    };

    const char *list[] = {"legendre",
                          "lobatto",
                          "jacobi",
                          "integrated_jacobi",
                          "h1tet_ainsworth",
                          "h1tet_bernstein_bezier",
                          "hdivtet_ainsworth",
                          "hdivtet_demkowicz",
                          "hcurltet_ainsworth",
                          "hcurltet_demkowicz",
                          "l2tet",
                          "h1tri_ainsworth",
                          "h1tri_bernstein_bezier",
                          "h1quad",
                          "hdivtri_ainsworth",
                          "hdivtri_demkowicz",
                          "hcurltri_ainsworth",
                          "hcurltri_demkowicz",
                          "l2tri",
                          "h1edge_ainsworth",
                          "h1edge_bernstein_bezier",
                          "hcurledge_ainsworth",
                          "hcurledge_demkowicz",
                          "h1flatprism",
                          "h1fatprism"};

    PetscBool flg;
    PetscInt choice_value = LEGENDREPOLYNOMIAL;
    ierr = PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list, LASTOP,
                                &choice_value, &flg);
    CHKERRG(ierr);
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    }

    MatrixDouble pts_1d(1, 3);
    pts_1d(0, 0) = -0.5;
    pts_1d(0, 1) = 0.;
    pts_1d(0, 2) = +0.5;

    boost::shared_ptr<MatrixDouble> base_ptr(new MatrixDouble);
    boost::shared_ptr<MatrixDouble> diff_base_ptr(new MatrixDouble);

    const double eps = 1e-3;

    if (choice_value == LEGENDREPOLYNOMIAL) {
      double diff_s = 1;
      ierr = LegendrePolynomial().getValue(
          pts_1d, boost::shared_ptr<BaseFunctionCtx>(new LegendrePolynomialCtx(
                      4, &diff_s, 1, base_ptr, diff_base_ptr)));
      CHKERRG(ierr);

      std::cout << "LegendrePolynomial\n";
      std::cout << pts_1d << std::endl;
      std::cout << *base_ptr << std::endl;
      std::cout << *diff_base_ptr << std::endl;
      double sum = sum_matrix(*base_ptr);
      double diff_sum = sum_matrix(*diff_base_ptr);
      std::cout << sum << std::endl;
      std::cout << diff_sum << std::endl;
      if (fabs(2.04688 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs(2.25 - diff_sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    pts_1d.resize(1, 11, false);
    for (int ii = 0; ii != 11; ii++) {
      pts_1d(0, ii) = 2 * ((double)ii / 10) - 1;
    }

    boost::shared_ptr<MatrixDouble> kernel_base_ptr(new MatrixDouble);
    boost::shared_ptr<MatrixDouble> diff_kernel_base_ptr(new MatrixDouble);

    if (choice_value == LOBATTOPOLYNOMIAL) {
      double diff_s = 1;
      ierr = LobattoPolynomial().getValue(
          pts_1d, boost::shared_ptr<BaseFunctionCtx>(new LobattoPolynomialCtx(
                      7 + 2, &diff_s, 1, base_ptr, diff_base_ptr)));
      CHKERRG(ierr);
      ierr = KernelLobattoPolynomial().getValue(
          pts_1d,
          boost::shared_ptr<BaseFunctionCtx>(new KernelLobattoPolynomialCtx(
              7, &diff_s, 1, kernel_base_ptr, diff_kernel_base_ptr)));
      CHKERRG(ierr);
      for (int ii = 0; ii != 11; ii++) {
        double s = pts_1d(0, ii);
        std::cerr << "lobatto_plot " << s << " " << (*base_ptr)(ii, 1) << " "
                  << (*diff_base_ptr)(ii, 1) << " " << (*kernel_base_ptr)(ii, 1)
                  << " " << (*diff_kernel_base_ptr)(ii, 1) << " "
                  << (*kernel_base_ptr)(ii, 1) * (1 - s * s) << " "
                  << (*kernel_base_ptr)(ii, 1) * (-2 * s) +
                         (*diff_kernel_base_ptr)(ii, 1) * (1 - s * s)
                  << " " << std::endl;
      }
      std::cout << "LobattoPolynomial\n";
      std::cout << pts_1d << std::endl;
      std::cout << *base_ptr << std::endl;
      std::cout << *diff_base_ptr << std::endl;
      {
        double sum = sum_matrix(*base_ptr);
        double diff_sum = sum_matrix(*diff_base_ptr);
        std::cout << sum << std::endl;
        std::cout << diff_sum << std::endl;
        if (fabs(-1.6053 - sum) > eps) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
        }
        if (fabs(3.07745 - diff_sum) > eps) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
        }
      }
      {
        double sum = sum_matrix(*kernel_base_ptr);
        double diff_sum = sum_matrix(*diff_kernel_base_ptr);
        std::cout << sum << std::endl;
        std::cout << diff_sum << std::endl;
        if (fabs(-13.9906 * 4 - sum) > eps) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
        }
        if (fabs(-101.678 * 4 - diff_sum) > eps) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
        }
      }
    }

    if (choice_value == JACOBIPOLYNOMIAL) {

      int n = 21;
      MatrixDouble pts_1d(1, n);
      pts_1d.resize(1, n, false);
      MatrixDouble pts_1d_t(1, n);
      for (int ii = 0; ii != n; ii++) {
        pts_1d(0, ii) = (double)ii / 20.;
        pts_1d_t(0, ii) = 1;
      }

      base_ptr->clear();
      diff_base_ptr->clear();

      double diff_x = 1;
      double diff_t = 0;
      ierr = JacobiPolynomial().getValue(
          pts_1d, pts_1d_t,
          boost::shared_ptr<BaseFunctionCtx>(new JacobiPolynomialCtx(
              5, &diff_x, &diff_t, 1, 0, base_ptr, diff_base_ptr)));
      CHKERRG(ierr);
      for (int ii = 0; ii != n; ii++) {
        double s = pts_1d(0, ii);
        std::cerr << "jacobi_plot " << s << " " << (*base_ptr)(ii, 4) << " "
                  << (*diff_base_ptr)(ii, 4) << std::endl;
      }
      std::cout << "JacobiPolynomial\n";
      std::cout << pts_1d << std::endl;
      std::cout << *base_ptr << std::endl;
      std::cout << *diff_base_ptr << std::endl;
      {
        double sum = sum_matrix(*base_ptr);
        double diff_sum = sum_matrix(*diff_base_ptr);
        std::cout << sum << std::endl;
        std::cout << diff_sum << std::endl;
        if (fabs(23.2164 - sum) > eps) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
        }
        if (fabs(167.995 - diff_sum) > eps) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
        }
      }
    }

    if (choice_value == INTEGRATEDJACOBIPOLYNOMIAL) {

      int n = 21;
      MatrixDouble pts_1d(1, n);
      pts_1d.resize(1, n, false);
      MatrixDouble pts_1d_t(1, n);
      for (int ii = 0; ii != n; ii++) {
        pts_1d(0, ii) = (double)ii / 20.;
        pts_1d_t(0, ii) = 1 - pts_1d(0, ii);
      }

      base_ptr->clear();
      diff_base_ptr->clear();

      double diff_x = 1;
      double diff_t = 0;
      ierr = IntegratedJacobiPolynomial().getValue(
          pts_1d, pts_1d_t,
          boost::shared_ptr<BaseFunctionCtx>(new IntegratedJacobiPolynomialCtx(
              6, &diff_x, &diff_t, 1, 1, base_ptr, diff_base_ptr)));
      CHKERRG(ierr);
      for (int ii = 0; ii != n; ii++) {
        double s = pts_1d(0, ii);
        std::cerr << "integrated_jacobi_plot " << s << " " << (double)ii / 20.
                  << " " << (*base_ptr)(ii, 1) << " " << (*base_ptr)(ii, 2)
                  << " " << (*base_ptr)(ii, 3) << " " << (*base_ptr)(ii, 4)
                  << " " << (*base_ptr)(ii, 5) << endl;
        ;
        // << " " << (*diff_base_ptr)(ii, 4) << std::endl;
      }
      std::cout << "IntegratedJacobiPolynomial\n";
      std::cout << pts_1d << std::endl;
      std::cout << *base_ptr << std::endl;
      std::cout << *diff_base_ptr << std::endl;
      {
        double sum = sum_matrix(*base_ptr);
        double diff_sum = sum_matrix(*diff_base_ptr);
        std::cout << sum << std::endl;
        std::cout << diff_sum << std::endl;
        // if(fabs(7.1915-sum)>eps) {
        //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
        // }
        // if(fabs(23.2164-diff_sum)>eps) {
        //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
        // }
      }
    }

    DataForcesAndSourcesCore tet_data(MBTET);
    for (int type = MBVERTEX; type != MBMAXTYPE; type++) {
      tet_data.spacesOnEntities[type].set(L2);
      tet_data.spacesOnEntities[type].set(H1);
      tet_data.spacesOnEntities[type].set(HDIV);
      tet_data.spacesOnEntities[type].set(HCURL);
    }
    tet_data.dataOnEntities[MBVERTEX].resize(1);
    tet_data.dataOnEntities[MBVERTEX][0].getDataOrder() = 1;
    tet_data.dataOnEntities[MBEDGE].resize(6);
    for (int ee = 0; ee < 6; ee++) {
      tet_data.dataOnEntities[MBEDGE][ee].getDataOrder() = 3;
      tet_data.dataOnEntities[MBEDGE][ee].getSense() = 1;
    }
    tet_data.dataOnEntities[MBTRI].resize(4);
    for (int ff = 0; ff < 4; ff++) {
      tet_data.dataOnEntities[MBTRI][ff].getDataOrder() = 4;
      tet_data.dataOnEntities[MBTRI][ff].getSense() = 1;
    }
    tet_data.dataOnEntities[MBTET].resize(1);
    tet_data.dataOnEntities[MBTET][0].getDataOrder() = 5;
    tet_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().resize(4, false);
    std::fill(tet_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().begin(),
              tet_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().end(), 5);

    MatrixDouble pts_tet;
    int tet_rule = 2;
    int nb_gauss_pts = QUAD_3D_TABLE[tet_rule]->npoints;
    pts_tet.resize(3, nb_gauss_pts, false);
    cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[tet_rule]->points[1], 4,
                &pts_tet(0, 0), 1);
    cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[tet_rule]->points[2], 4,
                &pts_tet(1, 0), 1);
    cblas_dcopy(nb_gauss_pts, &QUAD_3D_TABLE[tet_rule]->points[3], 4,
                &pts_tet(2, 0), 1);
    tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 4,
                                                             false);
    {
      double *shape_ptr =
          &*tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(4 * nb_gauss_pts, QUAD_3D_TABLE[tet_rule]->points, 1,
                  shape_ptr, 1);
    }
    tet_data.facesNodes.resize(4, 3);
    const int cannonical_tet_face[4][3] = {
        {0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1}};
    for (int ff = 0; ff < 4; ff++) {
      for (int nn = 0; nn < 3; nn++) {
        tet_data.facesNodes(ff, nn) = cannonical_tet_face[ff][nn];
      }
    }

    if (choice_value == H1TET_AINSWORTH) {

      tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 4,
                                                               false);
      ierr = ShapeMBTET(
          &*tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &pts_tet(0, 0), &pts_tet(1, 0), &pts_tet(2, 0), nb_gauss_pts);
      CHKERRG(ierr);

      ierr = TetPolynomialBase().getValue(
          pts_tet, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tet_data, H1, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      CHKERRG(ierr);
      if (tet_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          tet_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }

      double sum = 0, diff_sum = 0;
      std::cout << "Edges\n";
      for (int ee = 0; ee < 6; ee++) {
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getDiffN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_LEGENDRE_BASE));
        diff_sum += sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getDiffN(
            AINSWORTH_LEGENDRE_BASE));
      }
      std::cout << "Faces\n";
      for (int ff = 0; ff < 4; ff++) {
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_LEGENDRE_BASE));
        diff_sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getDiffN(
            AINSWORTH_LEGENDRE_BASE));
      }
      std::cout << "Tets\n";
      std::cout << tet_data.dataOnEntities[MBTET][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBTET][0].getDiffN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_LEGENDRE_BASE));
      diff_sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      if (fabs(1.3509 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs(0.233313 - diff_sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == H1TET_BERNSTEIN_BEZIER) {

      tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 4,
                                                               false);
      CHKERR Tools::shapeFunMBTET(
          &*tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &pts_tet(0, 0), &pts_tet(1, 0), &pts_tet(2, 0), nb_gauss_pts);
      CHKERR TetPolynomialBase().getValue(
          pts_tet, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tet_data, "TEST_FIELD", H1,
                       AINSWORTH_BERNSTEIN_BEZIER_BASE, NOBASE)));
      if (tet_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() ==
          tet_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_BERNSTEIN_BEZIER_BASE)
              .get())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "The same pointers");

      double sum = 0, diff_sum = 0;
      std::cout << "Vertices\n";
      std::cout << tet_data.dataOnEntities[MBVERTEX][0].getN("TEST_FIELD")
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBVERTEX][0].getDiffN("TEST_FIELD")
                << std::endl;
      sum +=
          sum_matrix(tet_data.dataOnEntities[MBVERTEX][0].getN("TEST_FIELD"));
      diff_sum += sum_matrix(
          tet_data.dataOnEntities[MBVERTEX][0].getDiffN("TEST_FIELD"));
      std::cout << "Edges\n";
      for (int ee = 0; ee < 6; ee++) {
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getN("TEST_FIELD")
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getDiffN("TEST_FIELD")
                  << std::endl;
        sum +=
            sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getN("TEST_FIELD"));
        diff_sum += sum_matrix(
            tet_data.dataOnEntities[MBEDGE][ee].getDiffN("TEST_FIELD"));
      }
      std::cout << "Faces\n";
      for (int ff = 0; ff < 4; ff++) {
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getN("TEST_FIELD")
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN("TEST_FIELD")
                  << std::endl;
        sum +=
            sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getN("TEST_FIELD"));
        diff_sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getDiffN("TEST_FIELD"));
      }
      std::cout << "Tets\n";
      std::cout << tet_data.dataOnEntities[MBTET][0].getN("TEST_FIELD")
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBTET][0].getDiffN("TEST_FIELD")
                << std::endl;
      sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getN("TEST_FIELD"));
      diff_sum +=
          sum_matrix(tet_data.dataOnEntities[MBTET][0].getDiffN("TEST_FIELD"));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      if (fabs(4.38395 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs(diff_sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == HDIVTET_AINSWORTH) {
      ierr = TetPolynomialBase().getValue(
          pts_tet, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tet_data, HDIV, AINSWORTH_LEGENDRE_BASE)));
      CHKERRG(ierr);
      double sum = 0, diff_sum = 0;
      std::cout << "Faces\n";
      for (int ff = 0; ff < 4; ff++) {
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_LEGENDRE_BASE));
        diff_sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getDiffN(
            AINSWORTH_LEGENDRE_BASE));
      }
      std::cout << "Tets\n";
      std::cout << tet_data.dataOnEntities[MBTET][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBTET][0].getDiffN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_LEGENDRE_BASE));
      diff_sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << 1e8 * sum << std::endl;
      std::cout << "diff_sum " << 1e8 * diff_sum << std::endl;
      if (fabs(0.188636 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs(32.9562 - diff_sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == HDIVTET_DEMKOWICZ) {
      ierr = TetPolynomialBase().getValue(
          pts_tet, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tet_data, HDIV, DEMKOWICZ_JACOBI_BASE)));
      CHKERRG(ierr);
      double sum = 0, diff_sum = 0;
      std::cout << "Faces\n";
      for (int ff = 0; ff < 4; ff++) {
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getN(
                         DEMKOWICZ_JACOBI_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(
                         DEMKOWICZ_JACOBI_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getN(DEMKOWICZ_JACOBI_BASE));
        diff_sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getDiffN(DEMKOWICZ_JACOBI_BASE));
      }
      std::cout << "Tets\n";
      std::cout << tet_data.dataOnEntities[MBTET][0].getN(DEMKOWICZ_JACOBI_BASE)
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBTET][0].getDiffN(
                       DEMKOWICZ_JACOBI_BASE)
                << std::endl;
      sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getN(DEMKOWICZ_JACOBI_BASE));
      diff_sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getDiffN(DEMKOWICZ_JACOBI_BASE));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      const double expected_result = -2.70651;
      const double expected_diff_result = 289.421;
      if (fabs((expected_result - sum) / expected_result) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs((expected_diff_result - diff_sum) / expected_diff_result) >
          eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == HCURLTET_AINSWORTH) {
      ierr = TetPolynomialBase().getValue(
          pts_tet, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tet_data, HCURL, AINSWORTH_LEGENDRE_BASE)));
      CHKERRG(ierr);
      double sum = 0, diff_sum = 0;
      std::cout << "Edges\n";
      for (int ee = 0; ee < 6; ee++) {
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getDiffN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_LEGENDRE_BASE));
        diff_sum += sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getDiffN(
            AINSWORTH_LEGENDRE_BASE));
      }
      std::cout << "Faces\n";
      for (int ff = 0; ff < 4; ff++) {
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_LEGENDRE_BASE));
        diff_sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getDiffN(
            AINSWORTH_LEGENDRE_BASE));
      }
      std::cout << "Tets\n";
      std::cout << tet_data.dataOnEntities[MBTET][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBTET][0].getDiffN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_LEGENDRE_BASE));
      diff_sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      if (fabs(-1.7798 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs(-67.1793 - diff_sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == HCURLTET_DEMKOWICZ) {
      CHKERR TetPolynomialBase().getValue(
          pts_tet, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tet_data, HCURL, DEMKOWICZ_JACOBI_BASE)));
      double sum = 0, diff_sum = 0;
      std::cout << "Edges\n";
      for (int ee = 0; ee < 6; ee++) {
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getN(
                         DEMKOWICZ_JACOBI_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBEDGE][ee].getDiffN(
                         DEMKOWICZ_JACOBI_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBEDGE][ee].getN(DEMKOWICZ_JACOBI_BASE));
        diff_sum += sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getDiffN(
            DEMKOWICZ_JACOBI_BASE));
      }
      std::cout << "Faces\n";
      for (int ff = 0; ff < 4; ff++) {
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getN(
                         DEMKOWICZ_JACOBI_BASE)
                  << std::endl;
        std::cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(
                         DEMKOWICZ_JACOBI_BASE)
                  << std::endl;
        sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getN(DEMKOWICZ_JACOBI_BASE));
        diff_sum += sum_matrix(
            tet_data.dataOnEntities[MBTRI][ff].getDiffN(DEMKOWICZ_JACOBI_BASE));
      }
      std::cout << "Tets\n";
      std::cout << tet_data.dataOnEntities[MBTET][0].getN(DEMKOWICZ_JACOBI_BASE)
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBTET][0].getDiffN(
                       DEMKOWICZ_JACOBI_BASE)
                << std::endl;
      sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getN(DEMKOWICZ_JACOBI_BASE));
      diff_sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getDiffN(DEMKOWICZ_JACOBI_BASE));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      if (fabs(7.35513 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs(62.4549 - diff_sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == L2TET) {
      ierr = TetPolynomialBase().getValue(
          pts_tet, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tet_data, L2, AINSWORTH_LEGENDRE_BASE)));
      CHKERRG(ierr);
      double sum = 0, diff_sum = 0;
      std::cout << "Tets\n";
      std::cout << tet_data.dataOnEntities[MBTET][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      std::cout << tet_data.dataOnEntities[MBTET][0].getDiffN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_LEGENDRE_BASE));
      diff_sum += sum_matrix(
          tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      if (fabs(3.60352 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
      if (fabs(-36.9994 - diff_sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    DataForcesAndSourcesCore tri_data(MBTRI);
    for (int type = MBVERTEX; type != MBMAXTYPE; type++) {
      tri_data.spacesOnEntities[type].set(L2);
      tri_data.spacesOnEntities[type].set(H1);
      tri_data.spacesOnEntities[type].set(HDIV);
      tri_data.spacesOnEntities[type].set(HCURL);
    }
    tri_data.dataOnEntities[MBVERTEX].resize(1);
    tri_data.dataOnEntities[MBVERTEX][0].getDataOrder() = 1;
    tri_data.dataOnEntities[MBEDGE].resize(3);
    for (int ee = 0; ee < 3; ee++) {
      tri_data.dataOnEntities[MBEDGE][ee].getDataOrder() = 3;
      tri_data.dataOnEntities[MBEDGE][ee].getSense() = 1;
    }
    tri_data.dataOnEntities[MBTRI].resize(1);
    tri_data.dataOnEntities[MBTRI][0].getDataOrder() = 4;
    tri_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().resize(3, false);
    std::fill(tri_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().begin(),
              tri_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().end(), 4);

    MatrixDouble pts_tri;
    int tri_rule = 2;
    nb_gauss_pts = QUAD_2D_TABLE[tri_rule]->npoints;
    pts_tri.resize(2, nb_gauss_pts, false);
    cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[tri_rule]->points[1], 3,
                &pts_tri(0, 0), 1);
    cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[tri_rule]->points[2], 3,
                &pts_tri(1, 0), 1);
    tri_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 3,
                                                             false);
    {
      double *shape_ptr =
          &*tri_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(3 * nb_gauss_pts, QUAD_2D_TABLE[tri_rule]->points, 1,
                  shape_ptr, 1);
    }
    tri_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(3, 2, false);
    ierr = ShapeDiffMBTRI(
        &*tri_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).data().begin());
    CHKERRG(ierr);

    if (choice_value == H1TRI_AINSWORTH) {
      CHKERR TriPolynomialBase().getValue(
          pts_tri, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tri_data, H1, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      if (tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          tri_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");

      double sum = 0, diff_sum = 0;
      std::cout << "Edges\n";
      for (int ee = 0; ee < 3; ee++) {
        std::cout << tri_data.dataOnEntities[MBEDGE][ee].getN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        std::cout << tri_data.dataOnEntities[MBEDGE][ee].getDiffN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        sum += sum_matrix(
            tri_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_LEGENDRE_BASE));
        diff_sum += sum_matrix(tri_data.dataOnEntities[MBEDGE][ee].getDiffN(
            AINSWORTH_LEGENDRE_BASE));
      }
      std::cout << "Face\n";
      std::cout << tri_data.dataOnEntities[MBTRI][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      std::cout << tri_data.dataOnEntities[MBTRI][0].getDiffN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_LEGENDRE_BASE));
      diff_sum += sum_matrix(
          tri_data.dataOnEntities[MBTRI][0].getDiffN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      if (fabs(0.805556 - sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");

      if (fabs(0.0833333 - diff_sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
    }

    if (choice_value == H1TRI_BERNSTEIN_BEZIER) {
      CHKERR TriPolynomialBase().getValue(
          pts_tri, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tri_data, "TET_FIELD", H1,
                       AINSWORTH_BERNSTEIN_BEZIER_BASE, NOBASE)));
      if (tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() ==
          tri_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_BERNSTEIN_BEZIER_BASE)
              .get())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Pointers should be diffrent");

      double sum = 0, diff_sum = 0;
      std::cout << "Vertex\n";
      std::cout << tri_data.dataOnEntities[MBVERTEX][0].getN("TET_FIELD")
                << std::endl;
      std::cout << tri_data.dataOnEntities[MBVERTEX][0].getDiffN("TET_FIELD")
                << std::endl;
      sum += sum_matrix(tri_data.dataOnEntities[MBVERTEX][0].getN("TET_FIELD"));
      diff_sum += sum_matrix(
          tri_data.dataOnEntities[MBVERTEX][0].getDiffN("TET_FIELD"));
      std::cout << "Edges\n";
      for (int ee = 0; ee < 3; ee++) {
        std::cout << tri_data.dataOnEntities[MBEDGE][ee].getN("TET_FIELD")
                  << std::endl;
        std::cout << tri_data.dataOnEntities[MBEDGE][ee].getDiffN("TET_FIELD")
                  << std::endl;
        sum +=
            sum_matrix(tri_data.dataOnEntities[MBEDGE][ee].getN("TET_FIELD"));
        diff_sum += sum_matrix(
            tri_data.dataOnEntities[MBEDGE][ee].getDiffN("TET_FIELD"));
      }
      std::cout << "Face\n";
      std::cout << tri_data.dataOnEntities[MBTRI][0].getN("TET_FIELD")
                << std::endl;
      std::cout << tri_data.dataOnEntities[MBTRI][0].getDiffN("TET_FIELD")
                << std::endl;
      sum += sum_matrix(tri_data.dataOnEntities[MBTRI][0].getN("TET_FIELD"));
      diff_sum +=
          sum_matrix(tri_data.dataOnEntities[MBTRI][0].getDiffN("TET_FIELD"));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff_sum " << diff_sum << std::endl;
      if (std::abs(3.01389 - sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");

      if (std::abs(diff_sum) > eps)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong result %3.4e != $3.4e", 0, diff_sum);
    }

    if (choice_value == HDIVTRI_AINSWORTH) {
      CHKERR TriPolynomialBase().getValue(
          pts_tri, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tri_data, HDIV, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      if (tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          tri_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");

      double sum = 0;
      std::cout << "Face\n";
      std::cout << tri_data.dataOnEntities[MBTRI][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      if (fabs(1.93056 - sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
    }

    if (choice_value == HDIVTRI_DEMKOWICZ) {
      ierr = TriPolynomialBase().getValue(
          pts_tri, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tri_data, HDIV, DEMKOWICZ_JACOBI_BASE, NOBASE)));
      CHKERRG(ierr);
      if (tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          tri_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(DEMKOWICZ_JACOBI_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      double sum = 0;
      std::cout << "Face\n";
      std::cout << tri_data.dataOnEntities[MBTRI][0].getN(DEMKOWICZ_JACOBI_BASE)
                << std::endl;
      sum += sum_matrix(
          tri_data.dataOnEntities[MBTRI][0].getN(DEMKOWICZ_JACOBI_BASE));
      std::cout << "sum  " << sum << std::endl;
      const double expected_result = 28.25;
      if (fabs((expected_result - sum) / expected_result) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == HCURLTRI_AINSWORTH) {
      ierr = TriPolynomialBase().getValue(
          pts_tri, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tri_data, HCURL, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      CHKERRG(ierr);
      if (tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          tri_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      double sum = 0; //,diff_sum = 0;
      std::cout << "Edges\n";
      for (int ee = 0; ee < 3; ee++) {
        std::cout << tri_data.dataOnEntities[MBEDGE][ee].getN(
                         AINSWORTH_LEGENDRE_BASE)
                  << std::endl;
        sum += sum_matrix(
            tri_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_LEGENDRE_BASE));
      }
      std::cout << "Face\n";
      std::cout << tri_data.dataOnEntities[MBTRI][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      if (fabs(0.333333 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == HCURLTRI_DEMKOWICZ) {
      CHKERR TriPolynomialBase().getValue(
          pts_tri, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tri_data, HCURL, DEMKOWICZ_JACOBI_BASE, NOBASE)));
      if (tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          tri_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(DEMKOWICZ_JACOBI_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      double sum = 0; //,diff_sum = 0;
      std::cout << "Edges\n";
      for (int ee = 0; ee < 3; ee++) {
        std::cout << tri_data.dataOnEntities[MBEDGE][ee].getN(
                         DEMKOWICZ_JACOBI_BASE)
                  << std::endl;
        sum += sum_matrix(
            tri_data.dataOnEntities[MBEDGE][ee].getN(DEMKOWICZ_JACOBI_BASE));
      }
      std::cout << "Face\n";
      std::cout << tri_data.dataOnEntities[MBTRI][0].getN(DEMKOWICZ_JACOBI_BASE)
                << std::endl;
      sum += sum_matrix(
          tri_data.dataOnEntities[MBTRI][0].getN(DEMKOWICZ_JACOBI_BASE));
      std::cout << "sum  " << sum << std::endl;
      if (fabs(1.04591 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == L2TRI) {
      ierr = TriPolynomialBase().getValue(
          pts_tri, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                       tri_data, L2, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      CHKERRG(ierr);
      if (tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          tri_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      double sum = 0;
      std::cout << "Face\n";
      std::cout << tri_data.dataOnEntities[MBTRI][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      if (fabs(0.671875 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    DataForcesAndSourcesCore edge_data(MBTRI);
    for (int type = MBVERTEX; type != MBMAXTYPE; type++) {
      edge_data.spacesOnEntities[type].set(L2);
      edge_data.spacesOnEntities[type].set(H1);
      edge_data.spacesOnEntities[type].set(HDIV);
      edge_data.spacesOnEntities[type].set(HCURL);
    }
    edge_data.dataOnEntities[MBVERTEX].resize(1);
    edge_data.dataOnEntities[MBVERTEX][0].getDataOrder() = 1;
    edge_data.dataOnEntities[MBEDGE].resize(1);
    edge_data.dataOnEntities[MBEDGE][0].getDataOrder() = 4;
    edge_data.dataOnEntities[MBEDGE][0].getSense() = 1;
    edge_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().resize(2, false);
    std::fill(edge_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().begin(),
              edge_data.dataOnEntities[MBVERTEX][0].getBBNodeOrder().end(), 4);

    MatrixDouble pts_edge;
    int edge_rule = 6;
    nb_gauss_pts = QUAD_1D_TABLE[edge_rule]->npoints;
    pts_edge.resize(1, nb_gauss_pts, false);
    cblas_dcopy(nb_gauss_pts, &QUAD_1D_TABLE[edge_rule]->points[1], 2,
                &pts_edge(0, 0), 1);
    edge_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 2,
                                                              false);
    {
      double *shape_ptr =
          &*edge_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(2 * nb_gauss_pts, QUAD_1D_TABLE[edge_rule]->points, 1,
                  shape_ptr, 1);
    }

    if (choice_value == H1EDGE_AINSWORTH) {
      CHKERR EdgePolynomialBase().getValue(
          pts_edge, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                        edge_data, H1, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      if (edge_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          edge_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      double sum = 0, diff_sum = 0;
      std::cout << "Edge\n";
      std::cout << edge_data.dataOnEntities[MBEDGE][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      std::cout << edge_data.dataOnEntities[MBEDGE][0].getDiffN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      sum += sum_matrix(
          edge_data.dataOnEntities[MBEDGE][0].getN(AINSWORTH_LEGENDRE_BASE));
      diff_sum += sum_matrix(edge_data.dataOnEntities[MBEDGE][0].getDiffN(
          AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff sum " << diff_sum << std::endl;
      if (std::abs(0.506122 - sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");

      if (std::abs(2.85714 - diff_sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
    }

    if (choice_value == H1EDGE_BERNSTEIN_BEZIER) {
      CHKERR EdgePolynomialBase().getValue(
          pts_edge, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                        edge_data, "TET_FIELD", H1,
                        AINSWORTH_BERNSTEIN_BEZIER_BASE, NOBASE)));
      if (edge_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() ==
          edge_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_BERNSTEIN_BEZIER_BASE)
              .get())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Should be diffrent pointers");

      double sum = 0, diff_sum = 0;
      std::cout << "Vertex\n";
      std::cout << edge_data.dataOnEntities[MBVERTEX][0].getN("TET_FIELD")
                << std::endl;
      std::cout << edge_data.dataOnEntities[MBVERTEX][0].getDiffN("TET_FIELD")
                << std::endl;
      std::cout << "Edge\n";
      std::cout << edge_data.dataOnEntities[MBEDGE][0].getN("TET_FIELD")
                << std::endl;
      std::cout << edge_data.dataOnEntities[MBEDGE][0].getDiffN("TET_FIELD")
                << std::endl;
      sum +=
          sum_matrix(edge_data.dataOnEntities[MBVERTEX][0].getN("TET_FIELD"));
      diff_sum += sum_matrix(
          edge_data.dataOnEntities[MBVERTEX][0].getDiffN("TET_FIELD"));
      sum += sum_matrix(edge_data.dataOnEntities[MBEDGE][0].getN("TET_FIELD"));
      diff_sum +=
          sum_matrix(edge_data.dataOnEntities[MBEDGE][0].getDiffN("TET_FIELD"));
      std::cout << "sum  " << sum << std::endl;
      std::cout << "diff sum " << diff_sum << std::endl;
      if (std::abs(4 - sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");

      if (std::abs(diff_sum) > eps)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
    }

    if (choice_value == HCURLEDGE_AINSWORTH) {
      ierr = EdgePolynomialBase().getValue(
          pts_edge, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                        edge_data, HCURL, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      CHKERRG(ierr);
      if (edge_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          edge_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      std::cout << edge_data.dataOnEntities[MBEDGE][0].getN(
                       AINSWORTH_LEGENDRE_BASE)
                << std::endl;
      int sum = 0;
      sum += sum_matrix(
          edge_data.dataOnEntities[MBEDGE][0].getN(AINSWORTH_LEGENDRE_BASE));
      std::cout << "sum  " << sum << std::endl;
      if (fabs(-4 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    if (choice_value == HCURLEDGE_DEMKOWICZ) {
      CHKERR EdgePolynomialBase().getValue(
          pts_edge, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                        edge_data, HCURL, DEMKOWICZ_JACOBI_BASE, NOBASE)));
      if (edge_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          edge_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(DEMKOWICZ_JACOBI_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      std::cout << edge_data.dataOnEntities[MBEDGE][0].getN(
                       DEMKOWICZ_JACOBI_BASE)
                << std::endl;
      int sum = 0;
      sum += sum_matrix(
          edge_data.dataOnEntities[MBEDGE][0].getN(DEMKOWICZ_JACOBI_BASE));
      std::cout << "sum  " << sum << std::endl;
      if (fabs(4 - sum) > eps) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      }
    }

    DataForcesAndSourcesCore quad_data(MBQUAD);
    for (int type = MBVERTEX; type != MBMAXTYPE; type++) {
      quad_data.spacesOnEntities[type].set(H1);
    }
    quad_data.dataOnEntities[MBVERTEX].resize(1);
    quad_data.dataOnEntities[MBVERTEX][0].getDataOrder() = 1;
    quad_data.dataOnEntities[MBEDGE].resize(4);
    for (int ee = 0; ee < 4; ee++) {
      quad_data.dataOnEntities[MBEDGE][ee].getDataOrder() = 4;
      quad_data.dataOnEntities[MBEDGE][ee].getSense() = 1;
    }
    quad_data.dataOnEntities[MBQUAD].resize(1);
    quad_data.dataOnEntities[MBQUAD][0].getDataOrder() = 6;

    MatrixDouble pts_quad;
    int rule_ksi = 6;
    int rule_eta = 4;
    CHKERR Tools::outerProductOfEdgeIntegrationPtsForQuad(pts_quad, rule_ksi,
                                                          rule_eta);
    nb_gauss_pts = pts_quad.size2();
    quad_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 4,
                                                              false);
    quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(nb_gauss_pts,
                                                                  8, false);
    for (int i = 0; i < nb_gauss_pts; ++i) {
      quad_data.dataOnEntities[MBVERTEX][0].getN(NOBASE)(i, 0) =
          N_MBQUAD0(pts_quad(0, i), pts_quad(1, i));
      quad_data.dataOnEntities[MBVERTEX][0].getN(NOBASE)(i, 1) =
          N_MBQUAD1(pts_quad(0, i), pts_quad(1, i));
      quad_data.dataOnEntities[MBVERTEX][0].getN(NOBASE)(i, 2) =
          N_MBQUAD2(pts_quad(0, i), pts_quad(1, i));
      quad_data.dataOnEntities[MBVERTEX][0].getN(NOBASE)(i, 3) =
          N_MBQUAD3(pts_quad(0, i), pts_quad(1, i));

      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 0) =
          diffN_MBQUAD0x(pts_quad(1, i));
      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 1) =
          diffN_MBQUAD0y(pts_quad(0, i));
      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 2) =
          diffN_MBQUAD1x(pts_quad(1, i));
      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 3) =
          diffN_MBQUAD1y(pts_quad(0, i));
      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 4) =
          diffN_MBQUAD2x(pts_quad(1, i));
      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 5) =
          diffN_MBQUAD2y(pts_quad(0, i));
      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 6) =
          diffN_MBQUAD3x(pts_quad(1, i));
      quad_data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(i, 7) =
          diffN_MBQUAD3y(pts_quad(0, i));
    }

    if (choice_value == H1QUAD) {
      CHKERR QuadPolynomialBase().getValue(
          pts_quad, boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                        quad_data, H1, AINSWORTH_LEGENDRE_BASE, NOBASE)));
      if (quad_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get() !=
          quad_data.dataOnEntities[MBVERTEX][0]
              .getNSharedPtr(AINSWORTH_LEGENDRE_BASE)
              .get()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Different pointers");
      }
      double sum = 0, diff_sum = 0;
      for (int ee = 0; ee < 4; ee++) {
        sum += sum_matrix(
            quad_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_LEGENDRE_BASE));
        diff_sum += sum_matrix(quad_data.dataOnEntities[MBEDGE][ee].getDiffN(
            AINSWORTH_LEGENDRE_BASE));
      }
      sum += sum_matrix(
          quad_data.dataOnEntities[MBQUAD][0].getN(AINSWORTH_LEGENDRE_BASE));
      diff_sum += sum_matrix(quad_data.dataOnEntities[MBQUAD][0].getDiffN(
          AINSWORTH_LEGENDRE_BASE));

      std::cout << sum << " " << diff_sum << endl;

      if (std::abs(3.58249 - sum) > eps) 
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      

      if (std::abs(-0.134694 - diff_sum) > eps) 
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
      
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}
