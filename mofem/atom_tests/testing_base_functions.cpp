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

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "teting interface inserting algorithm\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  enum bases {
    LEGENDREPOLYNOMIAL,
    LOBATTOPOLYNOMIAL,
    H1TET,
    HDIVTET,
    L2TET,
    H1TRI,
    HDIVTRI,
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
    "l2tet",
    "h1tri",
    "hdiftri",
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

  const double eps = 1e-3;

  if(choise_value==LEGENDREPOLYNOMIAL) {
    double diff_s = 1;
    ierr = LegendrePolynomial().getValue(
      pts_1d,
      base_ptr,
      diff_base_ptr,
      boost::shared_ptr<BaseFunctionCtx>(new LegendrePolynomialCtx(4,&diff_s,1))
    ); CHKERRQ(ierr);

    cout << "LegendrePolynomial\n";
    cout << *base_ptr << endl;
    cout << *diff_base_ptr << endl;
    if(fabs(-0.289062-(*base_ptr)(0,4))>1e-6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
    if(fabs(-1.5625-(*diff_base_ptr)(0,4))>1e-6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==LOBATTOPOLYNOMIAL) {
    double diff_s = 1;
    ierr = LobattoPolynomial().getValue(
      pts_1d,
      base_ptr,
      diff_base_ptr,
      boost::shared_ptr<BaseFunctionCtx>(new LobattoPolynomialCtx(4,&diff_s,1))
    ); CHKERRQ(ierr);
    cout << "LobattoPolynomial\n";
    cout << *base_ptr << endl;
    cout << *diff_base_ptr << endl;
    if(1) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
    }
  }

  if(choise_value==H1TET) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==HDIVTET) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==L2TET) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==H1TRI) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==HDIVTRI) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==H1EDGE) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==H1FLATPRIS) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  if(choise_value==H1FATPRISM) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"wrong result");
  }

  PetscFinalize();

}
