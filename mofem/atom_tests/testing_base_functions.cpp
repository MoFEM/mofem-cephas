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

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "testing interface inserting algorithm\n\n";

static double sum_matrix(ublas::matrix<double> &m) {
  double s = 0;
  for(int ii = 0;ii<m.size1();ii++) {
    for(int jj = 0;jj<m.size2();jj++) {
      s +=m (ii,jj);
    }
  }
  return s;
}

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  enum bases {
    LEGENDREPOLYNOMIAL,
    LOBATTOPOLYNOMIAL,
    H1TET,
    HDIVTET,
    HCURLTET,
    L2TET,
    H1TRI,
    HDIVTRI,
    HCURLTRI,
    L2TRI,
    H1EDGE,
    H1FLATPRIS,
    H1FATPRISM,
    LASTOP
  };

  const char *list[] = {
    "legendre",
    "lobatto",
    "h1tet",
    "hdivtet",
    "hcurltet",
    "l2tet",
    "h1tri",
    "hdivtri",
    "hcurltri",
    "l2tri",
    "h1edge",
    "h1flatprism",
    "h1fatprism"
  };

  PetscErrorCode ierr;

  PetscBool flg;
  PetscInt choise_value = LEGENDREPOLYNOMIAL;
  ierr = PetscOptionsGetEList(
    NULL,"-base",list,LASTOP,&choise_value,&flg
  ); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"base not set");
  }

  ublas::matrix<double> pts_1d(1,3);
  pts_1d(0,0)=-0.5;
  pts_1d(0,1)=0.;
  pts_1d(0,2)=+0.5;

  boost::shared_ptr<ublas::matrix<double> > base_ptr(new ublas::matrix<double>);
  boost::shared_ptr<ublas::matrix<double> > diff_base_ptr(new ublas::matrix<double>);

  const double eps = 1e-4;

  if(choise_value==LEGENDREPOLYNOMIAL) {
    double diff_s = 1;
    ierr = LegendrePolynomial().getValue(
      pts_1d,
      boost::shared_ptr<BaseFunctionCtx>(
        new LegendrePolynomialCtx(4,&diff_s,1,base_ptr,diff_base_ptr)
      )
    ); CHKERRQ(ierr);

    cout << "LegendrePolynomial\n";
    cout << *base_ptr << endl;
    cout << *diff_base_ptr << endl;
    double sum = sum_matrix(*base_ptr);
    double diff_sum = sum_matrix(*diff_base_ptr);
    cout << sum << endl;
    cout << diff_sum << endl;
    if(fabs(0.648438-sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
    if(fabs(1.3125-diff_sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==LOBATTOPOLYNOMIAL) {
    double diff_s = 1;
    ierr = LobattoPolynomial().getValue(
      pts_1d,
      boost::shared_ptr<BaseFunctionCtx>(
        new LobattoPolynomialCtx(4,&diff_s,1,base_ptr,diff_base_ptr)
      )
    ); CHKERRQ(ierr);
    cout << "LobattoPolynomial\n";
    cout << *base_ptr << endl;
    cout << *diff_base_ptr << endl;
    if(1) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  DataForcesAndSurcesCore tet_data(MBTET);
  for(int type = MBVERTEX;type!=MBMAXTYPE;type++) {
    tet_data.spacesOnEntities[type].set(L2);
    tet_data.spacesOnEntities[type].set(H1);
    tet_data.spacesOnEntities[type].set(HDIV);
  }
  tet_data.dataOnEntities[MBVERTEX].resize(1);
  tet_data.dataOnEntities[MBVERTEX][0].getDataOrder() = 1;
  tet_data.dataOnEntities[MBEDGE].resize(6);
  for(int ee = 0;ee<6;ee++) {
    tet_data.dataOnEntities[MBEDGE][ee].getDataOrder() = 3;
    tet_data.dataOnEntities[MBEDGE][ee].getSense() = 1;
  }
  tet_data.dataOnEntities[MBTRI].resize(4);
  for(int ff = 0;ff<4;ff++) {
    tet_data.dataOnEntities[MBTRI][ff].getDataOrder() = 4;
    tet_data.dataOnEntities[MBTRI][ff].getSense() = 1;
  }
  tet_data.dataOnEntities[MBTET].resize(1);
  tet_data.dataOnEntities[MBTET][0].getDataOrder() = 5;

  ublas::matrix<double> pts_tet;
  int tet_rule = 2;
  int nb_gauss_pts = QUAD_3D_TABLE[tet_rule]->npoints;
  pts_tet.resize(3,nb_gauss_pts,false);
  cblas_dcopy(
    nb_gauss_pts,&QUAD_3D_TABLE[tet_rule]->points[1],4,&pts_tet(0,0),1
  );
  cblas_dcopy(
    nb_gauss_pts,&QUAD_3D_TABLE[tet_rule]->points[2],4,&pts_tet(1,0),1
  );
  cblas_dcopy(
    nb_gauss_pts,&QUAD_3D_TABLE[tet_rule]->points[3],4,&pts_tet(2,0),1
  );
  tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,4,false);
  {
    double *shape_ptr = tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
    cblas_dcopy(
      4*nb_gauss_pts,QUAD_3D_TABLE[tet_rule]->points,1,shape_ptr,1
    );
  }
  tet_data.facesNodes.resize(4,3);
  const int cannonical_tet_face[4][3] = { {0,1,3}, {1,2,3}, {0,3,2}, {0,2,1} };
  for(int ff = 0;ff<4;ff++) {
    for(int nn = 0;nn<3;nn++) {
      tet_data.facesNodes(ff,nn) = cannonical_tet_face[ff][nn];
    }
  }

  if(choise_value==H1TET) {

    tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,4,false);
    ierr = ShapeMBTET(
      &*tet_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
      &pts_tet(0,0),
      &pts_tet(1,0),
      &pts_tet(2,0),
      nb_gauss_pts
    ); CHKERRQ(ierr);

    ierr = TetPolynomialBase().getValue(
      pts_tet,
      boost::shared_ptr<BaseFunctionCtx>(
        new EntPolynomialBaseCtx(tet_data,H1,AINSWORTH_COLE_BASE,NOBASE)
      )
    ); CHKERRQ(ierr);
    if(
      tet_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get()!=
      tet_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(AINSWORTH_COLE_BASE).get()
    ) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Different pointers");
    }

    double sum = 0,diff_sum = 0;
    cout << "Edges\n";
    for(int ee = 0;ee<6;ee++) {
      cout << tet_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_COLE_BASE) << endl;
      cout << tet_data.dataOnEntities[MBEDGE][ee].getDiffN(AINSWORTH_COLE_BASE) << endl;
      sum += sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_COLE_BASE));
      diff_sum += sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getDiffN(AINSWORTH_COLE_BASE));
    }
    cout << "Faces\n";
    for(int ff = 0;ff<4;ff++) {
      cout << tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_COLE_BASE) << endl;
      cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(AINSWORTH_COLE_BASE) << endl;
      sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_COLE_BASE));
      diff_sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getDiffN(AINSWORTH_COLE_BASE));
    }
    cout << "Tets\n";
    cout << tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE) << endl;
    cout << tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE) << endl;
    sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE));
    diff_sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE));
    cout << "sum  " << sum << endl;
    cout << "diff_sum " << diff_sum << endl;
    if(fabs(1.3509-sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
    if(fabs(0.233313-diff_sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==HDIVTET) {
    ierr = TetPolynomialBase().getValue(
      pts_tet,
      boost::shared_ptr<BaseFunctionCtx>(
        new EntPolynomialBaseCtx(tet_data,HDIV,AINSWORTH_COLE_BASE)
      )
    ); CHKERRQ(ierr);
    double sum = 0,diff_sum = 0;
    cout << "Faces\n";
    for(int ff = 0;ff<4;ff++) {
      cout << tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_COLE_BASE) << endl;
      cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(AINSWORTH_COLE_BASE) << endl;
      sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_COLE_BASE));
      diff_sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getDiffN(AINSWORTH_COLE_BASE));
    }
    cout << "Tets\n";
    cout << tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE) << endl;
    cout << tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE) << endl;
    sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE));
    diff_sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE));
    cout << "sum  " << sum << endl;
    cout << "diff_sum " << diff_sum << endl;
    if(fabs(0.188636-sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
    if(fabs(32.9562-diff_sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==HCURLTET) {
    ierr = TetPolynomialBase().getValue(
      pts_tet,
      boost::shared_ptr<BaseFunctionCtx>(
        new EntPolynomialBaseCtx(tet_data,HCURL,AINSWORTH_COLE_BASE)
      )
    ); CHKERRQ(ierr);
    double sum = 0,diff_sum = 0;
    cout << "Edges\n";
    for(int ee = 0;ee<6;ee++) {
      cout << tet_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_COLE_BASE) << endl;
      cout << tet_data.dataOnEntities[MBEDGE][ee].getDiffN(AINSWORTH_COLE_BASE) << endl;
      sum += sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_COLE_BASE));
      diff_sum += sum_matrix(tet_data.dataOnEntities[MBEDGE][ee].getDiffN(AINSWORTH_COLE_BASE));
    }
    cout << "Faces\n";
    for(int ff = 0;ff<4;ff++) {
      cout << tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_COLE_BASE) << endl;
      cout << tet_data.dataOnEntities[MBTRI][ff].getDiffN(AINSWORTH_COLE_BASE) << endl;
      sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getN(AINSWORTH_COLE_BASE));
      diff_sum += sum_matrix(tet_data.dataOnEntities[MBTRI][ff].getDiffN(AINSWORTH_COLE_BASE));
    }
    cout << "Tets\n";
    cout << tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE) << endl;
    cout << tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE) << endl;
    sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE));
    diff_sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE));
    cout << "sum  " << sum << endl;
    cout << "diff_sum " << diff_sum << endl;
    // if(fabs(1.3509-sum)>eps) {
    //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    // }
    // if(fabs(0.233313-diff_sum)>eps) {
    //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    // }
  }

  if(choise_value==L2TET) {
    ierr = TetPolynomialBase().getValue(
      pts_tet,
      boost::shared_ptr<BaseFunctionCtx>(
        new EntPolynomialBaseCtx(tet_data,L2,AINSWORTH_COLE_BASE)
      )
    ); CHKERRQ(ierr);
    double sum = 0,diff_sum = 0;
    cout << "Tets\n";
    cout << tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE) << endl;
    cout << tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE) << endl;
    sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getN(AINSWORTH_COLE_BASE));
    diff_sum += sum_matrix(tet_data.dataOnEntities[MBTET][0].getDiffN(AINSWORTH_COLE_BASE));
    cout << "sum  " << sum << endl;
    cout << "diff_sum " << diff_sum << endl;
    if(fabs(3.60352-sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
    if(fabs(-36.9994-diff_sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  DataForcesAndSurcesCore tri_data(MBTRI);
  for(int type = MBVERTEX;type!=MBMAXTYPE;type++) {
    tri_data.spacesOnEntities[type].set(L2);
    tri_data.spacesOnEntities[type].set(H1);
    tri_data.spacesOnEntities[type].set(HDIV);
  }
  tri_data.dataOnEntities[MBVERTEX].resize(1);
  tri_data.dataOnEntities[MBVERTEX][0].getDataOrder() = 1;
  tri_data.dataOnEntities[MBEDGE].resize(3);
  for(int ee = 0;ee<3;ee++) {
    tri_data.dataOnEntities[MBEDGE][ee].getDataOrder() = 3;
    tri_data.dataOnEntities[MBEDGE][ee].getSense() = 1;
  }
  tri_data.dataOnEntities[MBTRI].resize(1);
  tri_data.dataOnEntities[MBTRI][0].getDataOrder() = 4;

  ublas::matrix<double> pts_tri;
  int tri_rule = 2;
  nb_gauss_pts = QUAD_2D_TABLE[tri_rule]->npoints;
  pts_tri.resize(2,nb_gauss_pts,false);
  cblas_dcopy(
    nb_gauss_pts,&QUAD_2D_TABLE[tri_rule]->points[1],3,&pts_tri(0,0),1
  );
  cblas_dcopy(
    nb_gauss_pts,&QUAD_2D_TABLE[tri_rule]->points[2],3,&pts_tri(1,0),1
  );
  tri_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,3,false);
  {
    double *shape_ptr = tri_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
    cblas_dcopy(
      3*nb_gauss_pts,QUAD_2D_TABLE[tri_rule]->points,1,shape_ptr,1
    );
  }

  if(choise_value==H1TRI) {
    ierr = TriPolynomialBase().getValue(
      pts_tri,
      boost::shared_ptr<BaseFunctionCtx>(
        new EntPolynomialBaseCtx(tri_data,H1,AINSWORTH_COLE_BASE,NOBASE)
      )
    ); CHKERRQ(ierr);
    if(
      tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get()!=
      tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(AINSWORTH_COLE_BASE).get()
    ) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Different pointers");
    }
    double sum = 0,diff_sum = 0;
    cout << "Edges\n";
    for(int ee = 0;ee<3;ee++) {
      cout << tri_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_COLE_BASE) << endl;
      cout << tri_data.dataOnEntities[MBEDGE][ee].getDiffN(AINSWORTH_COLE_BASE) << endl;
      sum += sum_matrix(tri_data.dataOnEntities[MBEDGE][ee].getN(AINSWORTH_COLE_BASE));
      diff_sum += sum_matrix(tri_data.dataOnEntities[MBEDGE][ee].getDiffN(AINSWORTH_COLE_BASE));
    }
    cout << "Face\n";
    cout << tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_COLE_BASE) << endl;
    cout << tri_data.dataOnEntities[MBTRI][0].getDiffN(AINSWORTH_COLE_BASE) << endl;
    sum += sum_matrix(tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_COLE_BASE));
    diff_sum += sum_matrix(tri_data.dataOnEntities[MBTRI][0].getDiffN(AINSWORTH_COLE_BASE));
    cout << "sum  " << sum << endl;
    cout << "diff_sum " << diff_sum << endl;
    if(fabs(0.805556-sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
    if(fabs(0.0833333-diff_sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==HDIVTRI) {
    ierr = TriPolynomialBase().getValue(
      pts_tri,
      boost::shared_ptr<BaseFunctionCtx>(
        new EntPolynomialBaseCtx(tri_data,HDIV,AINSWORTH_COLE_BASE,NOBASE)
      )
    ); CHKERRQ(ierr);
    if(
      tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get()!=
      tri_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(AINSWORTH_COLE_BASE).get()
    ) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Different pointers");
    }
    double sum = 0;
    cout << "Face\n";
    cout << tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_COLE_BASE) << endl;
    sum += sum_matrix(tri_data.dataOnEntities[MBTRI][0].getN(AINSWORTH_COLE_BASE));
    cout << "sum  " << sum << endl;
    if(fabs(1.93056-sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==HCURLTRI) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==L2TRI) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  DataForcesAndSurcesCore edge_data(MBTRI);
  for(int type = MBVERTEX;type!=MBMAXTYPE;type++) {
    edge_data.spacesOnEntities[type].set(L2);
    edge_data.spacesOnEntities[type].set(H1);
    edge_data.spacesOnEntities[type].set(HDIV);
  }
  edge_data.dataOnEntities[MBVERTEX].resize(1);
  edge_data.dataOnEntities[MBVERTEX][0].getDataOrder() = 1;
  edge_data.dataOnEntities[MBEDGE].resize(1);
  edge_data.dataOnEntities[MBEDGE][0].getDataOrder() = 4;
  edge_data.dataOnEntities[MBEDGE][0].getSense() = 1;

  ublas::matrix<double> pts_edge;
  int edge_rule = 6;
  nb_gauss_pts = QUAD_1D_TABLE[edge_rule]->npoints;
  pts_edge.resize(1,nb_gauss_pts,false);
  cblas_dcopy(
    nb_gauss_pts,&QUAD_1D_TABLE[edge_rule]->points[1],2,&pts_edge(0,0),1
  );
  edge_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts,2,false);
  {
    double *shape_ptr = edge_data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
    cblas_dcopy(
      2*nb_gauss_pts,QUAD_1D_TABLE[edge_rule]->points,1,shape_ptr,1
    );
  }

  if(choise_value==H1EDGE) {
    ierr = EdgePolynomialBase().getValue(
      pts_edge,
      boost::shared_ptr<BaseFunctionCtx>(
        new EntPolynomialBaseCtx(edge_data,H1,AINSWORTH_COLE_BASE,NOBASE)
      )
    ); CHKERRQ(ierr);
    if(
      edge_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(NOBASE).get()!=
      edge_data.dataOnEntities[MBVERTEX][0].getNSharedPtr(AINSWORTH_COLE_BASE).get()
    ) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Different pointers");
    }
    double sum = 0,diff_sum = 0;
    cout << "Edge\n";
    cout << edge_data.dataOnEntities[MBEDGE][0].getN(AINSWORTH_COLE_BASE) << endl;
    cout << edge_data.dataOnEntities[MBEDGE][0].getDiffN(AINSWORTH_COLE_BASE) << endl;
    sum += sum_matrix(edge_data.dataOnEntities[MBEDGE][0].getN(AINSWORTH_COLE_BASE));
    diff_sum += sum_matrix(edge_data.dataOnEntities[MBEDGE][0].getDiffN(AINSWORTH_COLE_BASE));
    cout << "sum  " << sum << endl;
    cout << "diff sum " << diff_sum << endl;
    if(fabs(0.506122-sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
    if(fabs(2.85714-diff_sum)>eps) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==H1FLATPRIS) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==H1FATPRISM) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  PetscFinalize();

}
