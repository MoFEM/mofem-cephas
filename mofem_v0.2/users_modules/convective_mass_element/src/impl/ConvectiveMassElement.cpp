/** \file ConvectiveMassElement.cpp
 * \brief Operators and data structures for mass and convective mass element
 * \ingroup convective_mass_elem
 *
 */

/* Implementation of convective mass element
 *
 * This file is part of MoFEM.
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

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <adolc/adolc.h>
#include <DirichletBC.hpp>
#include <ConvectiveMassElement.hpp>

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

ConvectiveMassElement::MyVolumeFE::MyVolumeFE(FieldInterface &_mField):
VolumeElementForcesAndSourcesCore(_mField),
A(PETSC_NULL),
F(PETSC_NULL),
initV(false) {
  meshPositionsFieldName = "NoNE";
}


int ConvectiveMassElement::MyVolumeFE::getRule(int order) { return order; };

PetscErrorCode ConvectiveMassElement::MyVolumeFE::preProcess() {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  ierr = VolumeElementForcesAndSourcesCore::preProcess(); CHKERRQ(ierr);

  if(A != PETSC_NULL) {
    ts_B = A;
  }

  if(F != PETSC_NULL) {
    ts_F = F;
  }

  int ghosts[] = { 0 };
  int rank;
  MPI_Comm_rank(mField.get_comm(),&rank);

  switch (ts_ctx) {
    case CTX_TSNONE:
    if(!initV) {
      if(rank == 0) {
        ierr = VecCreateGhost(mField.get_comm(),1,1,1,ghosts,&V); CHKERRQ(ierr);
      } else {
        ierr = VecCreateGhost(mField.get_comm(),0,1,1,ghosts,&V); CHKERRQ(ierr);
      }
      initV = true;
    }
    ierr = VecZeroEntries(V); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    break;
    default:
    break;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConvectiveMassElement::MyVolumeFE::postProcess() {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  ierr = VolumeElementForcesAndSourcesCore::postProcess(); CHKERRQ(ierr);

  double *array;
  switch (ts_ctx) {
    case CTX_TSNONE:
    ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(V,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(V,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGetArray(V,&array); CHKERRQ(ierr);
    eNergy = array[0];
    ierr = VecRestoreArray(V,&array); CHKERRQ(ierr);
    if(initV) {
      ierr = VecDestroy(&V); CHKERRQ(ierr);
      initV = false;
    }
    break;
    default:
    break;
  }

  PetscFunctionReturn(0);
}

ConvectiveMassElement::ConvectiveMassElement(
  FieldInterface &m_field,short int tag
):
feMassRhs(m_field),
feMassLhs(m_field),
feMassAuxLhs(m_field),
feVelRhs(m_field),
feVelLhs(m_field),
feTRhs(m_field),
feTLhs(m_field),
feEnergy(m_field),
mField(m_field),tAg(tag) {

}

ConvectiveMassElement::OpGetDataAtGaussPts::OpGetDataAtGaussPts(const string field_name,
  vector<ublas::vector<double> > &values_at_gauss_pts,
  vector<ublas::matrix<double> > &gardient_at_gauss_pts
):
VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
valuesAtGaussPts(values_at_gauss_pts),gradientAtGaussPts(gardient_at_gauss_pts),
zeroAtType(MBVERTEX) {

}

PetscErrorCode ConvectiveMassElement::OpGetDataAtGaussPts::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  try {

    int nb_dofs = data.getFieldData().size();
    if(nb_dofs == 0) {
      PetscFunctionReturn(0);
    }
    int nb_gauss_pts = data.getN().size1();

    //initialize
    ublas::vector<double>& values = data.getFieldData();
    valuesAtGaussPts.resize(nb_gauss_pts);
    gradientAtGaussPts.resize(nb_gauss_pts);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      valuesAtGaussPts[gg].resize(3);
      gradientAtGaussPts[gg].resize(3,3);
    }

    if(type == zeroAtType) {
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        valuesAtGaussPts[gg].clear();
        gradientAtGaussPts[gg].clear();
      }
    }

    //cerr << valuesAtGaussPts[0] << " : ";

    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      ublas::vector<double> N = data.getN(gg,nb_dofs/3);
      ublas::matrix<double> diffN = data.getDiffN(gg,nb_dofs/3);
      for(int dd = 0;dd<nb_dofs/3;dd++) {
        for(int rr1 = 0;rr1<3;rr1++) {
          valuesAtGaussPts[gg][rr1] += N[dd]*values[3*dd+rr1];
          for(int rr2 = 0;rr2<3;rr2++) {
            gradientAtGaussPts[gg](rr1,rr2) += diffN(dd,rr2)*values[3*dd+rr1];
          }
        }
      }
    }
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

ConvectiveMassElement::OpGetCommonDataAtGaussPts::OpGetCommonDataAtGaussPts(const string field_name,CommonData &common_data):
OpGetDataAtGaussPts(field_name,
  common_data.dataAtGaussPts[field_name],
  common_data.gradAtGaussPts[field_name]
) {

}

ConvectiveMassElement::OpMassJacobian::OpMassJacobian(
  const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian,bool linear
):
VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
dAta(data),
commonData(common_data),
tAg(tag),
jAcobian(jacobian),
lInear(linear),
fieldDisp(false) {
}

PetscErrorCode ConvectiveMassElement::OpMassJacobian::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
    PetscFunctionReturn(0);
  }

  //do it only once, no need to repeat this for edges,faces or tets
  if(row_type != MBVERTEX) PetscFunctionReturn(0);

  int nb_dofs = row_data.getIndices().size();
  if(nb_dofs==0) PetscFunctionReturn(0);

  try {

    a.resize(3);
    dot_W.resize(3);
    dp_dt.resize(3);
    a_res.resize(3);

    g.resize(3,3);
    G.resize(3,3);
    h.resize(3,3);
    H.resize(3,3);
    invH.resize(3,3);
    F.resize(3,3);

    dot_W.clear();
    H.clear();
    invH.clear();
    for(int dd = 0;dd<3;dd++) {
      H(dd,dd) = 1;
      invH(dd,dd) = 1;
    }

    a_res.resize(3);
    int nb_gauss_pts = row_data.getN().size1();
    commonData.valMass.resize(nb_gauss_pts);
    commonData.jacMassRowPtr.resize(nb_gauss_pts);
    commonData.jacMass.resize(nb_gauss_pts);

    int nb_active_vars = 0;
    for(int gg = 0;gg<nb_gauss_pts;gg++) {

      if(gg == 0) {

        trace_on(tAg);

        for(int nn1 = 0;nn1<3;nn1++) { //0
          a[nn1] <<= (commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg])[nn1];
          nb_active_vars++;
        }
        for(int nn1 = 0;nn1<3;nn1++) { //3
          for(int nn2 = 0;nn2<3;nn2++) {
            h(nn1,nn2) <<= (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(nn1,nn2);
            if(fieldDisp) {
              if(nn1==nn2) {
                h(nn1,nn2) += 1;
              }
            }
            nb_active_vars++;
          }
        }
        if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
          for(int nn1 = 0;nn1<3;nn1++) { //3+9=12
            for(int nn2 = 0;nn2<3;nn2++) {
              g(nn1,nn2) <<= (commonData.gradAtGaussPts[commonData.spatialVelocities][gg])(nn1,nn2);
              nb_active_vars++;
            }
          }
          for(int nn1 = 0;nn1<3;nn1++) { //3+9+9=21
            dot_W(nn1) <<= (commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg])[nn1];
            nb_active_vars++;
          }
          for(int nn1 = 0;nn1<3;nn1++) { //3+9+9+3=24
            for(int nn2 = 0;nn2<3;nn2++) {
              H(nn1,nn2) <<= (commonData.gradAtGaussPts[commonData.meshPositions][gg])(nn1,nn2);
              nb_active_vars++;
            }
          }
        }
        adouble detH;
        ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
        ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
        noalias(G) = prod(g,invH);
        double rho0 = dAta.rho0;
        ublas::vector<double>& a0 = dAta.a0;
        if(!lInear) {
          noalias(F) = prod(h,invH);
          adouble detF;
          ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
          //calulate current density
          adouble rho = rho0*detF;
          //momentum rate
          noalias(dp_dt) = rho*(a0 + a + prod(G,dot_W));
        } else {
          noalias(dp_dt) = rho0*(a0 + a + prod(G,dot_W));
        }
        noalias(a_res) = dp_dt*detH;
        //dependant
        ublas::vector<double>& res = commonData.valMass[gg];
        res.resize(3);
        for(int rr = 0;rr<3;rr++) {
          a_res[rr] >>= res[rr];
        }

        trace_off();

      }

      active.resize(nb_active_vars);
      int aa = 0;
      for(int nn1 = 0;nn1<3;nn1++) { //0
        active[aa++] = (commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg])[nn1];
      }
      for(int nn1 = 0;nn1<3;nn1++) { //3
        for(int nn2 = 0;nn2<3;nn2++) {
          if(fieldDisp&&nn1 == nn2) {
            active[aa++] = (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(nn1,nn2)+1;
          } else {
            active[aa++] = (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(nn1,nn2);
          }
        }
      }
      if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
        for(int nn1 = 0;nn1<3;nn1++) { //3+9=12
          for(int nn2 = 0;nn2<3;nn2++) {
            active[aa++] = (commonData.gradAtGaussPts[commonData.spatialVelocities][gg])(nn1,nn2);
          }
        }
        for(int nn1 = 0;nn1<3;nn1++) { //3+9+9=21
          active[aa++] = (commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg])[nn1];
        }
        for(int nn1 = 0;nn1<3;nn1++) { //3+9+9+3=24
          for(int nn2 = 0;nn2<3;nn2++) {
            active[aa++] = (commonData.gradAtGaussPts[commonData.meshPositions][gg])(nn1,nn2);
          }
        }
      }

      if(!jAcobian) {
        ublas::vector<double>& res = commonData.valMass[gg];
        if(gg>0) {
          res.resize(3);
          int r;
          r = ::function(tAg,3,nb_active_vars,&active[0],&res[0]);
          if(r!=3) { // function is locally analytic
            SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
          }
        }
        double val = getVolume()*getGaussPts()(3,gg);
        if(getHoGaussPtsDetJac().size()>0) {
          val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
        }
        res *= val;
      } else {
        commonData.jacMassRowPtr[gg].resize(3);
        commonData.jacMass[gg].resize(3,nb_active_vars);
        for(int nn1 = 0;nn1<3;nn1++) {
          (commonData.jacMassRowPtr[gg])[nn1] = &(commonData.jacMass[gg](nn1,0));
        }
        int r;
        r = jacobian(
          tAg,3,nb_active_vars,
          &active[0],&(commonData.jacMassRowPtr[gg])[0]);
          if(r!=3) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
          }
          double val = getVolume()*getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
          }
          commonData.jacMass[gg] *= val;
        }

      }

    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

  ConvectiveMassElement::OpMassRhs::OpMassRhs(
    const string field_name,BlockData &data,CommonData &common_data
  ):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
  dAta(data),
  commonData(common_data) {
  }

  PetscErrorCode ConvectiveMassElement::OpMassRhs::doWork(
    int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
  ) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
      PetscFunctionReturn(0);
    }
    if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
    int nb_dofs = row_data.getIndices().size();

    try {

      nf.resize(nb_dofs);
      nf.clear();

      for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
        ublas::vector<double>& res = commonData.valMass[gg];
        //cerr << res << endl;
        for(int dd = 0;dd<nb_dofs/3;dd++) {
          for(int rr = 0;rr<3;rr++) {
            nf[3*dd+rr] += row_data.getN()(gg,dd)*res[rr];
          }
        }
      }

      if((unsigned int)nb_dofs > 3*row_data.getN().size2()) {
        SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = VecSetValues(getFEMethod()->ts_F,nb_dofs,
      &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);

    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

  ConvectiveMassElement::OpMassLhs_dM_dv::OpMassLhs_dM_dv(
    const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
  ):
  VolumeElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name,ForcesAndSurcesCore::UserDataOperator::OPROWCOL),
  dAta(data),
  commonData(common_data) {
    sYmm = false;
    if(forcesonlyonentities_ptr!=NULL) {
      forcesOnlyOnEntities = *forcesonlyonentities_ptr;
    }
  }

  PetscErrorCode ConvectiveMassElement::OpMassLhs_dM_dv::getJac(
    DataForcesAndSurcesCore::EntData &col_data,int gg
  ) {
    PetscFunctionBegin;
    try {
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacMass[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        for(int nn = 0;nn<3;nn++) {
          jac(0,3*dd+nn) = commonData.jacMass[gg](0,nn)*N(dd)*getFEMethod()->ts_a;
          jac(1,3*dd+nn) = commonData.jacMass[gg](1,nn)*N(dd)*getFEMethod()->ts_a;
          jac(2,3*dd+nn) = commonData.jacMass[gg](2,nn)*N(dd)*getFEMethod()->ts_a;
        }
      }
      if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
        ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
        for(int dd = 0;dd<nb_col/3;dd++) {
          //h00 //h01 //h02
          jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+3*0+0)*diffN(dd,0);
          jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+3*0+1)*diffN(dd,1);
          jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+3*0+2)*diffN(dd,2);
          //h10 //h11 //h12
          jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+3*1+0)*diffN(dd,0);
          jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+3*1+1)*diffN(dd,1);
          jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+3*1+2)*diffN(dd,2);
          //h20 //h21 //h22
          jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+3*2+0)*diffN(dd,0);
          jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+3*2+1)*diffN(dd,1);
          jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+3*2+2)*diffN(dd,2);
        }
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::OpMassLhs_dM_dv::doWork(
    int row_side,int col_side,
    EntityType row_type,EntityType col_type,
    DataForcesAndSurcesCore::EntData &row_data,
    DataForcesAndSurcesCore::EntData &col_data
  ) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
      PetscFunctionReturn(0);
    }

    int nb_row = row_data.getIndices().size();
    int nb_col = col_data.getIndices().size();
    if(nb_row==0) PetscFunctionReturn(0);
    if(nb_col==0) PetscFunctionReturn(0);

    try {

      k.resize(nb_row,nb_col);
      k.clear();
      jac.resize(3,nb_col);

      for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

        try {
          ierr = getJac(col_data,gg); CHKERRQ(ierr);
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }

        { //integrate element stiffnes matrix
          for(int dd1 = 0;dd1<nb_row/3;dd1++) {
            for(int rr1 = 0;rr1<3;rr1++) {
              for(int dd2 = 0;dd2<nb_col;dd2++) {
                k(3*dd1+rr1,dd2) += row_data.getN()(gg,dd1)*jac(rr1,dd2);
              }
            }
          }
        }

      }

      if(!forcesOnlyOnEntities.empty()) {
        ublas::vector<DofIdx> indices = row_data.getIndices();
        ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
        ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
        for(int ii = 0;dit!=dofs.end();dit++,ii++) {
          if(forcesOnlyOnEntities.find((*dit)->get_ent())==forcesOnlyOnEntities.end()) {
            indices[ii] = -1;
          }
        }
        ierr = MatSetValues(getFEMethod()->ts_B,
        nb_row,&indices[0],
        nb_col,&col_data.getIndices()[0],
        &k(0,0),ADD_VALUES); CHKERRQ(ierr);
      } else {
        ierr = MatSetValues(getFEMethod()->ts_B,
        nb_row,&row_data.getIndices()[0],
        nb_col,&col_data.getIndices()[0],
        &k(0,0),ADD_VALUES); CHKERRQ(ierr);
      }

    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }


  ConvectiveMassElement::OpMassLhs_dM_dx::OpMassLhs_dM_dx(
    const string field_name,const string col_field,BlockData &data,CommonData &common_data
  ):
  OpMassLhs_dM_dv(field_name,col_field,data,common_data) {}

  PetscErrorCode ConvectiveMassElement::OpMassLhs_dM_dx::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
    PetscFunctionBegin;
    try {
      int nb_col = col_data.getIndices().size();
      jac.clear();
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        //h00 //h01 //h02
        jac(0,3*dd+0) += commonData.jacMass[gg](0,3+3*0+0)*diffN(dd,0);
        jac(0,3*dd+0) += commonData.jacMass[gg](0,3+3*0+1)*diffN(dd,1);
        jac(0,3*dd+0) += commonData.jacMass[gg](0,3+3*0+2)*diffN(dd,2);
        //h10 //h11 //h12
        jac(1,3*dd+1) += commonData.jacMass[gg](1,3+3*1+0)*diffN(dd,0);
        jac(1,3*dd+1) += commonData.jacMass[gg](1,3+3*1+1)*diffN(dd,1);
        jac(1,3*dd+1) += commonData.jacMass[gg](1,3+3*1+2)*diffN(dd,2);
        //h20 //h21 //h22
        jac(2,3*dd+2) += commonData.jacMass[gg](2,3+3*2+0)*diffN(dd,0);
        jac(2,3*dd+2) += commonData.jacMass[gg](2,3+3*2+1)*diffN(dd,1);
        jac(2,3*dd+2) += commonData.jacMass[gg](2,3+3*2+2)*diffN(dd,2);
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  ConvectiveMassElement::OpMassLhs_dM_dX::OpMassLhs_dM_dX(
    const string field_name,const string col_field,BlockData &data,CommonData &common_data
  ):
  OpMassLhs_dM_dv(field_name,col_field,data,common_data) {

  }

  PetscErrorCode ConvectiveMassElement::OpMassLhs_dM_dX::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
    PetscFunctionBegin;
    try {
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacVel[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        for(int nn = 0;nn<3;nn++) {
          jac(0,3*dd+nn) = commonData.jacMass[gg](0,3+9+9+nn)*N(dd)*getFEMethod()->ts_a;
          jac(1,3*dd+nn) = commonData.jacMass[gg](1,3+9+9+nn)*N(dd)*getFEMethod()->ts_a;
          jac(2,3*dd+nn) = commonData.jacMass[gg](2,3+9+9+nn)*N(dd)*getFEMethod()->ts_a;
        }
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        //h00 //h01 //h02
        jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+9+3+3*0+0)*diffN(dd,0);
        jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+9+3+3*0+1)*diffN(dd,1);
        jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+9+3+3*0+2)*diffN(dd,2);
        //h10 //h11 //h12
        jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+9+3+3*1+0)*diffN(dd,0);
        jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+9+3+3*1+1)*diffN(dd,1);
        jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+9+3+3*1+2)*diffN(dd,2);
        //h20 //h21 //h22
        jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+9+3+3*2+0)*diffN(dd,0);
        jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+9+3+3*2+1)*diffN(dd,1);
        jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+9+3+3*2+2)*diffN(dd,2);
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }


ConvectiveMassElement::OpEnergy::OpEnergy(
  const string field_name,BlockData &data,CommonData &common_data,Vec *v_ptr,bool linear
):
VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
dAta(data),commonData(common_data),Vptr(v_ptr),lInear(linear) {

}

PetscErrorCode ConvectiveMassElement::OpEnergy::doWork(
  int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    if(row_type != MBVERTEX) {
      PetscFunctionReturn(0);
    }
    if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
      PetscFunctionReturn(0);
    }

    try {

      for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
        double val = getVolume()*getGaussPts()(3,gg);
        if(getHoGaussPtsDetJac().size()>0) {
          val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
        }
        double rho0 = dAta.rho0;
        double rho;
        if(lInear) {
          rho = rho0;
        } else {
          h.resize(3,3);
          noalias(h) = (commonData.gradAtGaussPts[commonData.spatialPositions][gg]);
          if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
            H.resize(3,3);
            noalias(H) = (commonData.gradAtGaussPts[commonData.meshPositions][gg]);
            double detH;
            ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
            invH.resize(3,3);
            ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
            F.resize(3,3);
            noalias(F) = prod(h,invH);
          } else {
            F.resize(3,3);
            noalias(F) = h;
          }
          double detF;
          ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
          rho = detF*rho0;
        }
        v.resize(3);
        noalias(v) = commonData.dataAtGaussPts[commonData.spatialVelocities][gg];
        double energy = 0.5*rho*inner_prod(v,v);
        ierr = VecSetValue(*Vptr,0,val*energy,ADD_VALUES); CHKERRQ(ierr);
      }

    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

  ConvectiveMassElement::OpVelocityJacobian::OpVelocityJacobian(
    const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian
  ):
  VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
  dAta(data),
  commonData(common_data),
  tAg(tag),
  jAcobian(jacobian),
  fieldDisp(false) { }

  PetscErrorCode ConvectiveMassElement::OpVelocityJacobian::doWork(
    int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
  ) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
      PetscFunctionReturn(0);
    }

    //do it only once, no need to repeat this for edges,faces or tets
    if(row_type != MBVERTEX) PetscFunctionReturn(0);

    int nb_dofs = row_data.getIndices().size();
    if(nb_dofs==0) PetscFunctionReturn(0);

    try {

      v.resize(3);
      dot_w.resize(3);
      h.resize(3,3);
      h.clear();
      F.resize(3,3);
      dot_W.resize(3);
      dot_W.clear();
      H.resize(3,3);
      H.clear();
      invH.resize(3,3);
      invH.clear();
      dot_u.resize(3);
      for(int dd = 0;dd<3;dd++) {
        H(dd,dd) = 1;
        invH(dd,dd) = 1;
      }

      a_res.resize(3);
      int nb_gauss_pts = row_data.getN().size1();
      commonData.valVel.resize(nb_gauss_pts);
      commonData.jacVelRowPtr.resize(nb_gauss_pts);
      commonData.jacVel.resize(nb_gauss_pts);

      int nb_active_vars = 0;
      for(int gg = 0;gg<nb_gauss_pts;gg++) {

        if(gg == 0) {

          trace_on(tAg);

          for(int nn1 = 0;nn1<3;nn1++) { //0
            v[nn1] <<= commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1];
            nb_active_vars++;
          }
          for(int nn1 = 0;nn1<3;nn1++) { //3
            dot_w[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.spatialPositions][gg][nn1];
            nb_active_vars++;
          }
          if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
            for(int nn1 = 0;nn1<3;nn1++) { //3+3 = 6
              for(int nn2 = 0;nn2<3;nn2++) {
                h(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2);
                if(fieldDisp) {
                  if(nn1==nn2) {
                    h(nn1,nn2) += 1;
                  }
                }
                nb_active_vars++;
              }
            }
            for(int nn1 = 0;nn1<3;nn1++) { //3+3+9
              dot_W[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg][nn1];
              nb_active_vars++;
            }
          }
          if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
            for(int nn1 = 0;nn1<3;nn1++) { //3+3+9+3
              for(int nn2 = 0;nn2<3;nn2++) {
                H(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2);
                nb_active_vars++;
              }
            }
          }
          detH = 1;
          if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
            ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
            ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
            noalias(F) = prod(h,invH);
          } else {
            noalias(F) = h;
          }
          noalias(dot_u) = dot_w-prod(F,dot_W);
          noalias(a_res) = (v - dot_u)*detH;
          //dependant
          ublas::vector<double>& res = commonData.valVel[gg];
          res.resize(3);
          for(int rr = 0;rr<3;rr++) {
            a_res[rr] >>= res[rr];
          }
          trace_off();

        }

        active.resize(nb_active_vars);
        int aa = 0;
        for(int nn1 = 0;nn1<3;nn1++) {
          active[aa++] = commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1];
        }
        for(int nn1 = 0;nn1<3;nn1++) {
          active[aa++] = commonData.dataAtGaussPts["DOT_"+commonData.spatialPositions][gg][nn1];
        }
        if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
          for(int nn1 = 0;nn1<3;nn1++) {
            for(int nn2 = 0;nn2<3;nn2++) {
              if(fieldDisp&&nn1 == nn2) {
                active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2)+1;
              } else {
                active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2);
              }
            }
          }
          for(int nn1 = 0;nn1<3;nn1++) {
            active[aa++] = commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg][nn1];
          }
        }
        if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
          for(int nn1 = 0;nn1<3;nn1++) {
            for(int nn2 = 0;nn2<3;nn2++) {
              active[aa++] = commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2);
            }
          }
        }

        if(!jAcobian) {
          ublas::vector<double>& res = commonData.valVel[gg];
          if(gg>0) {
            res.resize(3);
            int r;
            r = ::function(tAg,3,nb_active_vars,&active[0],&res[0]);
            if(r!=3) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
            }
          }
          double val = getVolume()*getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
          }
          res *= val;
        } else {
          commonData.jacVelRowPtr[gg].resize(3);
          commonData.jacVel[gg].resize(3,nb_active_vars);
          for(int nn1 = 0;nn1<3;nn1++) {
            (commonData.jacVelRowPtr[gg])[nn1] = &(commonData.jacVel[gg](nn1,0));
          }
          int r;
          r = jacobian(
            tAg,3,nb_active_vars,
            &active[0],&(commonData.jacVelRowPtr[gg])[0]);
            if(r!=3) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
            }
            double val = getVolume()*getGaussPts()(3,gg);
            if(getHoGaussPtsDetJac().size()>0) {
              val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }
            commonData.jacVel[gg] *= val;
            //cerr << gg << " : " << commonData.jacVel[gg] << endl;
          }
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    ConvectiveMassElement::OpVelocityRhs::OpVelocityRhs(
      const string field_name,BlockData &data,CommonData &common_data
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
    dAta(data),
    commonData(common_data) {
    }

    PetscErrorCode ConvectiveMassElement::OpVelocityRhs::doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
        PetscFunctionReturn(0);
      }
      int nb_dofs = row_data.getIndices().size();
      if(nb_dofs==0) PetscFunctionReturn(0);

      try {

        nf.resize(nb_dofs);
        nf.clear();

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
          ublas::vector<double>& res = commonData.valVel[gg];
          for(int dd = 0;dd<nb_dofs/3;dd++) {
            for(int rr = 0;rr<3;rr++) {
              nf[3*dd+rr] += row_data.getN()(gg,dd)*res[rr];
            }
          }
        }

        if(row_data.getIndices().size() > 3*row_data.getN().size2()) {
          SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        }
        ierr = VecSetValues(getFEMethod()->ts_F,row_data.getIndices().size(),
        &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }


    ConvectiveMassElement::OpVelocityLhs_dV_dv::OpVelocityLhs_dV_dv(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data
    ):
    OpMassLhs_dM_dv(vel_field,field_name,data,common_data) {
    }

    PetscErrorCode ConvectiveMassElement::OpVelocityLhs_dV_dv::getJac(
      DataForcesAndSurcesCore::EntData &col_data,int gg
    ) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      //cerr << commonData.jacVel[gg] << endl;
      for(int dd = 0;dd<nb_col/3;dd++) {
        for(int nn = 0;nn<3;nn++) {
          jac(0,3*dd+nn) = commonData.jacVel[gg](0,nn)*N(dd);
          jac(1,3*dd+nn) = commonData.jacVel[gg](1,nn)*N(dd);
          jac(2,3*dd+nn) = commonData.jacVel[gg](2,nn)*N(dd);
        }
      }
      PetscFunctionReturn(0);
    }

    ConvectiveMassElement::OpVelocityLhs_dV_dx::OpVelocityLhs_dV_dx(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data
    ):
    OpVelocityLhs_dV_dv(vel_field,field_name,data,common_data) {
    }

    PetscErrorCode ConvectiveMassElement::OpVelocityLhs_dV_dx::getJac(
      DataForcesAndSurcesCore::EntData &col_data,int gg
    ) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacVel[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        for(int nn = 0;nn<3;nn++) {
          jac(0,3*dd+nn) = commonData.jacVel[gg](0,3+nn)*N(dd)*getFEMethod()->ts_a;
          jac(1,3*dd+nn) = commonData.jacVel[gg](1,3+nn)*N(dd)*getFEMethod()->ts_a;
          jac(2,3*dd+nn) = commonData.jacVel[gg](2,3+nn)*N(dd)*getFEMethod()->ts_a;
        }
      }
      if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
        ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
        for(int dd = 0;dd<nb_col/3;dd++) {
          //h00 //h01 //h02
          jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+3*0+0)*diffN(dd,0);
          jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+3*0+1)*diffN(dd,1);
          jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+3*0+2)*diffN(dd,2);
          //h10 //h11 //h12
          jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+3*1+0)*diffN(dd,0);
          jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+3*1+1)*diffN(dd,1);
          jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+3*1+2)*diffN(dd,2);
          //h20 //h21 //h22
          jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+3*2+0)*diffN(dd,0);
          jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+3*2+1)*diffN(dd,1);
          jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+3*2+2)*diffN(dd,2);
        }
      }
      //cerr << row_field_name << " " << col_field_name << endl;
      PetscFunctionReturn(0);
    }

    ConvectiveMassElement::OpVelocityLhs_dV_dX::OpVelocityLhs_dV_dX(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data
    ):
    OpVelocityLhs_dV_dv(vel_field,field_name,data,common_data) {}

    PetscErrorCode ConvectiveMassElement::OpVelocityLhs_dV_dX::getJac(
      DataForcesAndSurcesCore::EntData &col_data,int gg
    ) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacVel[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        for(int nn = 0;nn<3;nn++) {
          jac(0,3*dd+nn) = commonData.jacVel[gg](0,3+3+9+nn)*N(dd)*getFEMethod()->ts_a;
          jac(1,3*dd+nn) = commonData.jacVel[gg](1,3+3+9+nn)*N(dd)*getFEMethod()->ts_a;
          jac(2,3*dd+nn) = commonData.jacVel[gg](2,3+3+9+nn)*N(dd)*getFEMethod()->ts_a;
        }
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        //h00 //h01 //h02
        jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+9+3+3*0+0)*diffN(dd,0);
        jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+9+3+3*0+1)*diffN(dd,1);
        jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+9+3+3*0+2)*diffN(dd,2);
        //h10 //h11 //h12
        jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+9+3+3*1+0)*diffN(dd,0);
        jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+9+3+3*1+1)*diffN(dd,1);
        jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+9+3+3*1+2)*diffN(dd,2);
        //h20 //h21 //h22
        jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+9+3+3*2+0)*diffN(dd,0);
        jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+9+3+3*2+1)*diffN(dd,1);
        jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+9+3+3*2+2)*diffN(dd,2);
      }

      //cerr << row_field_name << " " << col_field_name << endl;

      PetscFunctionReturn(0);
    }

    ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumJacobian::OpEshelbyDynamicMaterialMomentumJacobian(
      const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
    dAta(data),
    commonData(common_data),
    tAg(tag),
    jAcobian(jacobian),
    fieldDisp(false) {}

    PetscErrorCode ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumJacobian::doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
        PetscFunctionReturn(0);
      }

      //do it only once, no need to repeat this for edges,faces or tets
      if(row_type != MBVERTEX) PetscFunctionReturn(0);

      int nb_dofs = row_data.getIndices().size();
      if(nb_dofs==0) PetscFunctionReturn(0);

      try {

        a.resize(3);
        v.resize(3);
        g.resize(3,3);
        G.resize(3,3);
        h.resize(3,3);
        F.resize(3,3);
        H.resize(3,3);
        H.clear();
        invH.resize(3,3);
        invH.clear();
        for(int dd = 0;dd<3;dd++) {
          H(dd,dd) = 1;
          invH(dd,dd) = 1;
        }

        int nb_gauss_pts = row_data.getN().size1();
        commonData.valT.resize(nb_gauss_pts);
        commonData.jacTRowPtr.resize(nb_gauss_pts);
        commonData.jacT.resize(nb_gauss_pts);

        int nb_active_vars = 0;
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          if(gg == 0) {

            trace_on(tAg);

            for(int nn1 = 0;nn1<3;nn1++) { //0
              a[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg][nn1]; nb_active_vars++;
            }

            for(int nn1 = 0;nn1<3;nn1++) { //3
              v[nn1] <<= commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1]; nb_active_vars++;
            }
            for(int nn1 = 0;nn1<3;nn1++) { //3+3
              for(int nn2 = 0;nn2<3;nn2++) {
                g(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialVelocities][gg](nn1,nn2); nb_active_vars++;
              }
            }
            for(int nn1 = 0;nn1<3;nn1++) { //3+3+9
              for(int nn2 = 0;nn2<3;nn2++) {
                h(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2); nb_active_vars++;
                if(fieldDisp) {
                  if(nn1==nn2) {
                    h(nn1,nn2) += 1;
                  }
                }
              }
            }
            if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
              for(int nn1 = 0;nn1<3;nn1++) { //3+3+9+9
                for(int nn2 = 0;nn2<3;nn2++) {
                  H(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2); nb_active_vars++;
                }
              }
            }
            adouble detH;
            detH = 1;
            if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
              ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
            }
            ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
            noalias(F) = prod(h,invH);
            noalias(G) = prod(g,invH);
            double rho0 = dAta.rho0;
            a_T.resize(3);
            noalias(a_T) = -rho0*(prod(trans(F),a)+prod(trans(G),v))*detH;
            commonData.valT[gg].resize(3);
            for(int nn = 0;nn<3;nn++) {
              a_T[nn] >>= (commonData.valT[gg])[nn];
            }
            trace_off();

          }

          active.resize(nb_active_vars);
          int aa = 0;
          for(int nn1 = 0;nn1<3;nn1++) { //0
            active[aa++] = commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg][nn1];
          }

          for(int nn1 = 0;nn1<3;nn1++) { //3
            active[aa++] = commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1];
          }
          for(int nn1 = 0;nn1<3;nn1++) { //3+3
            for(int nn2 = 0;nn2<3;nn2++) {
              active[aa++] = commonData.gradAtGaussPts[commonData.spatialVelocities][gg](nn1,nn2);
            }
          }
          for(int nn1 = 0;nn1<3;nn1++) { //3+3+9
            for(int nn2 = 0;nn2<3;nn2++) {
              if(fieldDisp&&nn1 == nn2) {
                active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2)+1;
              } else {
                active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2);
              }
            }
          }
          if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
            for(int nn1 = 0;nn1<3;nn1++) { //3+3+9+9
              for(int nn2 = 0;nn2<3;nn2++) {
                active[aa++] = commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2);
              }
            }
          }

          if(!jAcobian) {
            ublas::vector<double>& res = commonData.valT[gg];
            if(gg>0) {
              res.resize(3);
              int r;
              r = ::function(tAg,3,nb_active_vars,&active[0],&res[0]);
              if(r!=3) { // function is locally analytic
                SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
              }
            }
            double val = getVolume()*getGaussPts()(3,gg);
            if(getHoGaussPtsDetJac().size()>0) {
              val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
            }
            res *= val;
          } else {
            commonData.jacTRowPtr[gg].resize(3);
            commonData.jacT[gg].resize(3,nb_active_vars);
            for(int nn1 = 0;nn1<3;nn1++) {
              (commonData.jacTRowPtr[gg])[nn1] = &(commonData.jacT[gg](nn1,0));
            }
            int r;
            r = jacobian(
              tAg,3,nb_active_vars,
              &active[0],&(commonData.jacTRowPtr[gg])[0]);
              if(r!=3) {
                SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
              }
              double val = getVolume()*getGaussPts()(3,gg);
              if(getHoGaussPtsDetJac().size()>0) {
                val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
              }
              commonData.jacT[gg] *= val;
            }

          }

        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }

      PetscFunctionReturn(0);
    }

    ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumRhs::OpEshelbyDynamicMaterialMomentumRhs(
      const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
    dAta(data),
    commonData(common_data) {
      if(forcesonlyonentities_ptr!=NULL) {
        forcesOnlyOnEntities = *forcesonlyonentities_ptr;
      }
    }

    PetscErrorCode ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumRhs::doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
        PetscFunctionReturn(0);
      }
      int nb_dofs = row_data.getIndices().size();
      if(nb_dofs==0) PetscFunctionReturn(0);

      try {

        nf.resize(nb_dofs);
        nf.clear();

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
          ublas::vector<double>& res = commonData.valT[gg];
          //cerr << res << endl;
          for(int dd = 0;dd<nb_dofs/3;dd++) {
            for(int rr = 0;rr<3;rr++) {
              nf[3*dd+rr] += row_data.getN()(gg,dd)*res[rr];
            }
          }
        }

        if(row_data.getIndices().size() > 3*row_data.getN().size2()) {
          SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        }
        if(!forcesOnlyOnEntities.empty()) {
          ublas::vector<DofIdx> indices = row_data.getIndices();
          ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
          ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
          for(int ii = 0;dit!=dofs.end();dit++,ii++) {
            if(forcesOnlyOnEntities.find((*dit)->get_ent())==forcesOnlyOnEntities.end()) {
              //cerr << **dit << endl;
              indices[ii] = -1;
            }
          }
          //cerr << indices << endl;
          ierr = VecSetValues(getFEMethod()->ts_F,indices.size(),&indices[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);
        } else {
          ierr = VecSetValues(getFEMethod()->ts_F,row_data.getIndices().size(),
          &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dv::OpEshelbyDynamicMaterialMomentumLhs_dv(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,
      Range *forcesonlyonentities_ptr
    ):
    ConvectiveMassElement::OpMassLhs_dM_dv(
      vel_field,field_name,data,common_data,forcesonlyonentities_ptr
    ) {

    }


    PetscErrorCode ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dv::getJac(
      DataForcesAndSurcesCore::EntData &col_data,int gg
    ) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacT[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        for(int nn = 0;nn<3;nn++) {
          jac(0,3*dd+nn) = commonData.jacT[gg](0,nn)*N(dd)*getFEMethod()->ts_a;
          jac(1,3*dd+nn) = commonData.jacT[gg](1,nn)*N(dd)*getFEMethod()->ts_a;
          jac(2,3*dd+nn) = commonData.jacT[gg](2,nn)*N(dd)*getFEMethod()->ts_a;
        }
      }
      for(int dd = 0;dd<nb_col/3;dd++) {
        for(int nn = 0;nn<3;nn++) {
          jac(0,3*dd+nn) += commonData.jacT[gg](0,3+nn)*N(dd);
          jac(1,3*dd+nn) += commonData.jacT[gg](1,3+nn)*N(dd);
          jac(2,3*dd+nn) += commonData.jacT[gg](2,3+nn)*N(dd);
        }
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        //h00 //h01 //h02
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+3*0+0)*diffN(dd,0);
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+3*0+1)*diffN(dd,1);
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+3*0+2)*diffN(dd,2);
        //h10 //h11 //h12
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+3*1+0)*diffN(dd,0);
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+3*1+1)*diffN(dd,1);
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+3*1+2)*diffN(dd,2);
        //h20 //h21 //h22
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+3*2+0)*diffN(dd,0);
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+3*2+1)*diffN(dd,1);
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+3*2+2)*diffN(dd,2);
      }
      PetscFunctionReturn(0);
    }

    ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dx::OpEshelbyDynamicMaterialMomentumLhs_dx(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
    ):
    ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dv(
      vel_field,field_name,data,common_data,forcesonlyonentities_ptr
    ) {

    }

    PetscErrorCode ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dx::getJac(
      DataForcesAndSurcesCore::EntData &col_data,int gg
    ) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacT[gg] << endl;
      //cerr << jac << endl;
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        //h00 //h01 //h02
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+3*0+0)*diffN(dd,0);
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+3*0+1)*diffN(dd,1);
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+3*0+2)*diffN(dd,2);
        //h10 //h11 //h12
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+3*1+0)*diffN(dd,0);
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+3*1+1)*diffN(dd,1);
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+3*1+2)*diffN(dd,2);
        //h20 //h21 //h22
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+3*2+0)*diffN(dd,0);
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+3*2+1)*diffN(dd,1);
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+3*2+2)*diffN(dd,2);
      }
      PetscFunctionReturn(0);
    }


    ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dX::OpEshelbyDynamicMaterialMomentumLhs_dX(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
    ):
    ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dv(
      vel_field,field_name,data,common_data,forcesonlyonentities_ptr
    ) {

    }

    PetscErrorCode ConvectiveMassElement::OpEshelbyDynamicMaterialMomentumLhs_dX::getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacT[gg] << endl;
      //cerr << jac << endl;
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
        //h00 //h01 //h02
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+9+3*0+0)*diffN(dd,0);
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+9+3*0+1)*diffN(dd,1);
        jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+9+3*0+2)*diffN(dd,2);
        //h10 //h11 //h12
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+9+3*1+0)*diffN(dd,0);
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+9+3*1+1)*diffN(dd,1);
        jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+9+3*1+2)*diffN(dd,2);
        //h20 //h21 //h22
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+9+3*2+0)*diffN(dd,0);
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+9+3*2+1)*diffN(dd,1);
        jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+9+3*2+2)*diffN(dd,2);
      }
      PetscFunctionReturn(0);
    }


    ConvectiveMassElement::UpdateAndControl::UpdateAndControl(
      FieldInterface& _mField,TS _ts,
      const string velocity_field,
      const string spatial_position_field
    ):
    mField(_mField),tS(_ts),
    velocityField(velocity_field),
    spatialPositionField(spatial_position_field),
    jacobianLag(-1) {

    }

    PetscErrorCode ConvectiveMassElement::UpdateAndControl::preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      switch (ts_ctx) {
        case CTX_TSSETIFUNCTION: {
          snes_ctx = CTX_SNESSETFUNCTION;
          snes_f = ts_F;
          break;
        }
        case CTX_TSSETIJACOBIAN: {
          snes_ctx = CTX_SNESSETJACOBIAN;
          snes_B = ts_B;
          break;
        }
        default:
        break;
      }

      //ierr = mField.set_other_local_ghost_vector(problemPtr,velocityField,"DOT_"+velocityField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      //ierr = mField.set_other_local_ghost_vector(problemPtr,spatialPositionField,"DOT_"+spatialPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      //FIXME: This global scattering because Kuu problem and Dynamic problem
      //not share partitions. Both problem should use the same partitioning to
      //resolve this problem.
      ierr = mField.set_global_ghost_vector(problemPtr,COL,ts_u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = mField.set_other_global_ghost_vector(problemPtr,velocityField,"DOT_"+velocityField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = mField.set_other_global_ghost_vector(problemPtr,spatialPositionField,"DOT_"+spatialPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    PetscErrorCode ConvectiveMassElement::UpdateAndControl::postProcess() {
      PetscFunctionBegin;
      //PetscErrorCode ierr;
      //SNES snes;
      //ierr = TSGetSNES(tS,&snes); CHKERRQ(ierr);
      //ierr = SNESSetLagJacobian(snes,jacobianLag); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }



    PetscErrorCode ConvectiveMassElement::setBlocks() {
      PetscFunctionBegin;
      ErrorCode rval;
      PetscErrorCode ierr;

      Range added_tets;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|BODYFORCESSET,it)) {
        int id = it->get_msId();
        EntityHandle meshset = it->get_meshset();
        rval = mField.get_moab().get_entities_by_type(meshset,MBTET,setOfBlocks[id].tEts,true); CHKERR_PETSC(rval);
        added_tets.merge(setOfBlocks[id].tEts);
        Block_BodyForces mydata;
        ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
        setOfBlocks[id].rho0 = mydata.data.density;
        setOfBlocks[id].a0.resize(3);
        setOfBlocks[id].a0[0] = mydata.data.acceleration_x;
        setOfBlocks[id].a0[1] = mydata.data.acceleration_y;
        setOfBlocks[id].a0[2] = mydata.data.acceleration_z;
        //cerr << setOfBlocks[id].tEts << endl;
      }

      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        Mat_Elastic mydata;
        ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
        if(mydata.data.User1 == 0) continue;
        Range tets;
        EntityHandle meshset = it->get_meshset();
        rval = mField.get_moab().get_entities_by_type(meshset,MBTET,tets,true); CHKERR_PETSC(rval);
        tets = subtract(tets,added_tets);
        if(tets.empty()) continue;
        int id = it->get_msId();
        setOfBlocks[-id].tEts = tets;
        setOfBlocks[-id].rho0 = mydata.data.User1;
        setOfBlocks[-id].a0.resize(3);
        setOfBlocks[-id].a0[0] = mydata.data.User2;
        setOfBlocks[-id].a0[1] = mydata.data.User3;
        setOfBlocks[-id].a0[2] = mydata.data.User4;
        //cerr << setOfBlocks[id].tEts << endl;
      }

      PetscFunctionReturn(0);
    }

  PetscErrorCode ConvectiveMassElement::addConvectiveMassElement(
    string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,
    bool ale,BitRefLevel bit) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
    if(mField.check_field(material_position_field_name)) {
      if(ale) {
        ierr = mField.modify_finite_element_add_field_row(element_name,material_position_field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+material_position_field_name); CHKERRQ(ierr);
      }
      ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
    }
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+spatial_position_field_name); CHKERRQ(ierr);

    Range tets;
    if(bit.any()) {
      ierr = mField.get_entities_by_type_and_ref_level(bit,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
    }

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      Range add_tets = sit->second.tEts;
      if(!tets.empty()) {
        add_tets = intersect(add_tets,tets);
      }
      ierr = mField.add_ents_to_finite_element_by_TETs(add_tets,element_name); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::addVelocityElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,
    bool ale,BitRefLevel bit) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      //ErrorCode rval;

      ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(element_name,velocity_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
      if(mField.check_field(material_position_field_name)) {
        if(ale) {
          ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
          ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+material_position_field_name); CHKERRQ(ierr);
        }
        ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
      }
      ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+velocity_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+spatial_position_field_name); CHKERRQ(ierr);

      Range tets;
      if(bit.any()) {
        ierr = mField.get_entities_by_type_and_ref_level(bit,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
      }

      map<int,BlockData>::iterator sit = setOfBlocks.begin();
      for(;sit!=setOfBlocks.end();sit++) {
        Range add_tets = sit->second.tEts;
        if(!tets.empty()) {
          add_tets = intersect(add_tets,tets);
        }
        ierr = mField.add_ents_to_finite_element_by_TETs(add_tets,element_name); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

  PetscErrorCode ConvectiveMassElement::addEshelbyDynamicMaterialMomentum(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,
    bool ale,
    BitRefLevel bit,
    Range *intersected
  ) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
    if(mField.check_field(material_position_field_name)) {
      if(ale) {
        ierr = mField.modify_finite_element_add_field_row(element_name,material_position_field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+material_position_field_name); CHKERRQ(ierr);
      }
      ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
    }
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+spatial_position_field_name); CHKERRQ(ierr);

    Range tets;
    if(bit.any()) {
      ierr = mField.get_entities_by_type_and_ref_level(bit,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
    }
    if(intersected!=NULL) {
      if(tets.empty()) {
        tets = *intersected;
      } else {
        tets = intersect(*intersected,tets);
      }
    }

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      Range add_tets = sit->second.tEts;
      if(!tets.empty()) {
        add_tets = intersect(add_tets,tets);
      }
      ierr = mField.add_ents_to_finite_element_by_TETs(add_tets,element_name); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::setConvectiveMassOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,
    bool ale,
    bool linear
  ) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
        feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      } else {
        feMassRhs.meshPositionsFieldName = material_position_field_name;
      }
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassRhs.getOpPtrVector().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,false,linear));
      feMassRhs.getOpPtrVector().push_back(new OpMassRhs(spatial_position_field_name,sit->second,commonData));
    }

    //Lhs
    feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
        feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      } else {
        feMassLhs.meshPositionsFieldName = material_position_field_name;
      }
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassLhs.getOpPtrVector().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,true,linear));
      feMassLhs.getOpPtrVector().push_back(new OpMassLhs_dM_dv(spatial_position_field_name,velocity_field_name,sit->second,commonData));
      feMassLhs.getOpPtrVector().push_back(new OpMassLhs_dM_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
        if(ale) {
          feMassLhs.getOpPtrVector().push_back(new OpMassLhs_dM_dX(spatial_position_field_name,material_position_field_name,sit->second,commonData));
        } else {
          feMassLhs.meshPositionsFieldName = material_position_field_name;
        }
      }
    }

    //Energy
    feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feEnergy.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feEnergy.getOpPtrVector().push_back(new OpEnergy(spatial_position_field_name,sit->second,commonData,&feEnergy.V,linear));
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::setVelocityOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,
    bool ale
  ) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feVelRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feVelRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feVelRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feVelRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
      feVelRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
        feVelRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      } else {
        feVelRhs.meshPositionsFieldName = material_position_field_name;
      }
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feVelRhs.getOpPtrVector().push_back(new OpVelocityJacobian(velocity_field_name,sit->second,commonData,tAg,false));
      feVelRhs.getOpPtrVector().push_back(new OpVelocityRhs(velocity_field_name,sit->second,commonData));
    }

    //Lhs
    feVelLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feVelLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feVelLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feVelLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
      feVelLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
        feVelLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      } else {
        feVelLhs.meshPositionsFieldName = material_position_field_name;
      }
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feVelLhs.getOpPtrVector().push_back(new OpVelocityJacobian(velocity_field_name,sit->second,commonData,tAg));
      feVelLhs.getOpPtrVector().push_back(new OpVelocityLhs_dV_dv(velocity_field_name,velocity_field_name,sit->second,commonData));
      feVelLhs.getOpPtrVector().push_back(new OpVelocityLhs_dV_dx(velocity_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
        if(ale) {
          feVelLhs.getOpPtrVector().push_back(new OpVelocityLhs_dV_dX(velocity_field_name,material_position_field_name,sit->second,commonData));
        } else {
          feVelLhs.meshPositionsFieldName = material_position_field_name;
        }
      }
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode ConvectiveMassElement::setKinematicEshelbyOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,
    Range *forces_on_entities_ptr
  ) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feTRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feTRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feTRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    feTRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feTRhs.getOpPtrVector().push_back(
        new OpEshelbyDynamicMaterialMomentumJacobian(material_position_field_name,sit->second,commonData,tAg,false)
      );
      feTRhs.getOpPtrVector().push_back(
        new OpEshelbyDynamicMaterialMomentumRhs(material_position_field_name,sit->second,commonData,forces_on_entities_ptr)
      );
    }

    //Lhs
    feTLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feTLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feTLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feTLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feTLhs.getOpPtrVector().push_back(new OpEshelbyDynamicMaterialMomentumJacobian(material_position_field_name,sit->second,commonData,tAg));
      feTLhs.getOpPtrVector().push_back(
        new OpEshelbyDynamicMaterialMomentumLhs_dv(material_position_field_name,velocity_field_name,sit->second,commonData,forces_on_entities_ptr));
      feTLhs.getOpPtrVector().push_back(
        new OpEshelbyDynamicMaterialMomentumLhs_dx(material_position_field_name,spatial_position_field_name,sit->second,commonData,forces_on_entities_ptr));
      feTLhs.getOpPtrVector().push_back(
        new OpEshelbyDynamicMaterialMomentumLhs_dX(material_position_field_name,material_position_field_name,sit->second,commonData,forces_on_entities_ptr));
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::setShellMatrixMassOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,
    bool linear) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassRhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feMassRhs.meshPositionsFieldName = material_position_field_name;
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassRhs.getOpPtrVector().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,false,linear));
      feMassRhs.getOpPtrVector().push_back(new OpMassRhs(spatial_position_field_name,sit->second,commonData));
    }

    //Lhs
    feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feMassLhs.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassLhs.getOpPtrVector().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,true,linear));
      feMassLhs.getOpPtrVector().push_back(new OpMassLhs_dM_dv(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
        feMassLhs.meshPositionsFieldName = material_position_field_name;
      }
    }

    //Aux Lhs
    feMassAuxLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassAuxLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassAuxLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassAuxLhs.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feMassAuxLhs.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassAuxLhs.getOpPtrVector().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,true,linear));
      feMassAuxLhs.getOpPtrVector().push_back(new OpMassLhs_dM_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
        feMassAuxLhs.meshPositionsFieldName = material_position_field_name;
      }
    }

    //Energy E=0.5*rho*v*v
    feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feEnergy.getOpPtrVector().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feEnergy.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feEnergy.getOpPtrVector().push_back(new OpEnergy(spatial_position_field_name,sit->second,commonData,&feEnergy.V,linear));
    }

    PetscFunctionReturn(0);
  }

  ConvectiveMassElement::MatShellCtx::MatShellCtx(): iNitialized(false) {}
  ConvectiveMassElement::MatShellCtx::~MatShellCtx() {
    if(iNitialized) {
      PetscErrorCode ierr;
      ierr = dEstroy(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
  }

  PetscErrorCode ConvectiveMassElement::MatShellCtx::iNit() {
    PetscFunctionBegin;
    if(!iNitialized) {
      PetscErrorCode ierr;
      #if PETSC_VERSION_GE(3,5,3)
      ierr = MatCreateVecs(K,&u,&Ku); CHKERRQ(ierr);
      ierr = MatCreateVecs(M,&v,&Mv); CHKERRQ(ierr);
      #else
      ierr = MatGetVecs(K,&u,&Ku); CHKERRQ(ierr);
      ierr = MatGetVecs(M,&v,&Mv); CHKERRQ(ierr);
      #endif
      ierr = MatDuplicate(K,MAT_SHARE_NONZERO_PATTERN,&barK); CHKERRQ(ierr);
      iNitialized = true;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::MatShellCtx::dEstroy() {
    PetscFunctionBegin;
    if(iNitialized) {
      PetscErrorCode ierr;
      ierr = VecDestroy(&u); CHKERRQ(ierr);
      ierr = VecDestroy(&Ku); CHKERRQ(ierr);
      ierr = VecDestroy(&v); CHKERRQ(ierr);
      ierr = VecDestroy(&Mv); CHKERRQ(ierr);
      ierr = MatDestroy(&barK); CHKERRQ(ierr);
      iNitialized = false;
    }
    PetscFunctionReturn(0);
  }


  PetscErrorCode ConvectiveMassElement::PCShellCtx::iNit() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if(!initPC) {
      MPI_Comm comm;
      ierr = PetscObjectGetComm((PetscObject)shellMat,&comm); CHKERRQ(ierr);
      ierr = PCCreate(comm,&pC); CHKERRQ(ierr);
      initPC = true;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::PCShellCtx::dEstroy() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if(initPC) {
      ierr = PCDestroy(&pC); CHKERRQ(ierr);
      initPC = false;
    }
    PetscFunctionReturn(0);
  }

  ConvectiveMassElement::ShellResidualElement::ShellResidualElement(
    FieldInterface &m_field
  ): mField(m_field) {

  }

  PetscErrorCode ConvectiveMassElement::ShellResidualElement::preProcess() {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    if(ts_ctx != CTX_TSSETIFUNCTION) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"It is used to residual of velocities");
    }
    if(!shellMatCtx->iNitialized) {
      ierr = shellMatCtx->iNit(); CHKERRQ(ierr);
    }
    ierr = VecScatterBegin(shellMatCtx->scatterU,ts_u_t,shellMatCtx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shellMatCtx->scatterU,ts_u_t,shellMatCtx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterBegin(shellMatCtx->scatterV,ts_u,shellMatCtx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shellMatCtx->scatterV,ts_u,shellMatCtx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecAXPY(shellMatCtx->v,-1,shellMatCtx->u); CHKERRQ(ierr);
    ierr = VecScatterBegin(shellMatCtx->scatterV,shellMatCtx->v,ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(shellMatCtx->scatterV,shellMatCtx->v,ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //VecView(shellMatCtx->v,PETSC_VIEWER_STDOUT_WORLD);

    PetscFunctionReturn(0);
  }

  PetscErrorCode ConvectiveMassElement::ShellResidualElement::postProcess() {
    PetscFunctionBegin;

    /*PetscErrorCode ierr;
    if(ts_ctx != CTX_TSSETIFUNCTION) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"It is used to residual of velocities");
    }
    if(!shellMatCtx->iNitialized) {
      ierr = shellMatCtx->iNit(); CHKERRQ(ierr);
    }
    ierr = VecScatterBegin(shellMatCtx->scatterU,ts_F,shellMatCtx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shellMatCtx->scatterU,ts_F,shellMatCtx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterBegin(shellMatCtx->scatterV,ts_F,shellMatCtx->Mv,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shellMatCtx->scatterV,ts_F,shellMatCtx->Mv,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    double nrm2_ku,nrm2_mv;
    ierr = VecNorm(shellMatCtx->Ku,NORM_INFINITY,&nrm2_ku); CHKERRQ(ierr);
    ierr = VecNorm(shellMatCtx->Mv,NORM_INFINITY,&nrm2_mv); CHKERRQ(ierr);
    PetscPrintf(mField.get_comm(),"nrm2 U = %6.4e nrm2 V = %6.4e scale = %6.4e\n",nrm2_ku,nrm2_mv,nrm2_ku/fmax(nrm2_mv,nrm2_ku));
    //shellMatCtx->scale = nrm2_ku/fmax(nrm2_mv,nrm2_ku);
    //ierr = VecScale(shellMatCtx->Mv,shellMatCtx->scale); CHKERRQ(ierr);
    ierr = VecScatterBegin(shellMatCtx->scatterV,shellMatCtx->Mv,ts_F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(shellMatCtx->scatterV,shellMatCtx->Mv,ts_F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);*/

    PetscFunctionReturn(0);

  }

  #ifdef __DIRICHLETBC_HPP__

  ConvectiveMassElement::ShellMatrixElement::ShellMatrixElement(
    FieldInterface &m_field
  ): mField(m_field) {

  }

  PetscErrorCode ConvectiveMassElement::ShellMatrixElement::preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    if(ts_ctx != CTX_TSSETIJACOBIAN) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"It is used to calculate shell matrix only");
    }

    shellMatCtx->ts_a = ts_a;
    DirichletBcPtr->copy_ts(*((TSMethod*)this)); //copy context for TSMethod

    DirichletBcPtr->dIag = 1;
    DirichletBcPtr->ts_B = shellMatCtx->K;
    ierr = MatZeroEntries(shellMatCtx->K); CHKERRQ(ierr);
    ierr = mField.problem_basic_method_preProcess(problemName,*DirichletBcPtr); CHKERRQ(ierr);
    LoopsToDoType::iterator itk = loopK.begin();
    for(;itk!=loopK.end();itk++) {
      itk->second->copy_ts(*((TSMethod*)this));
      itk->second->ts_B = shellMatCtx->K;
      ierr = mField.loop_finite_elements(problemName,itk->first,*itk->second); CHKERRQ(ierr);
    }
    LoopsToDoType::iterator itam = loopAuxM.begin();
    for(;itam!=loopAuxM.end();itam++) {
      itam->second->copy_ts(*((TSMethod*)this));
      itam->second->ts_B = shellMatCtx->K;
      ierr = mField.loop_finite_elements(problemName,itam->first,*itam->second); CHKERRQ(ierr);
    }
    ierr = mField.problem_basic_method_postProcess(problemName,*DirichletBcPtr); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(shellMatCtx->K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(shellMatCtx->K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    DirichletBcPtr->dIag = 0;
    DirichletBcPtr->ts_B = shellMatCtx->M;
    ierr = MatZeroEntries(shellMatCtx->M); CHKERRQ(ierr);
    //ierr = mField.problem_basic_method_preProcess(problemName,*DirichletBcPtr); CHKERRQ(ierr);
    LoopsToDoType::iterator itm = loopM.begin();
    for(;itm!=loopM.end();itm++) {
      itm->second->copy_ts(*((TSMethod*)this));
      itm->second->ts_B = shellMatCtx->M;
      ierr = mField.loop_finite_elements(problemName,itm->first,*itm->second); CHKERRQ(ierr);
    }
    ierr = mField.problem_basic_method_postProcess(problemName,*DirichletBcPtr); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(shellMatCtx->M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(shellMatCtx->M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    //barK
    ierr = MatZeroEntries(shellMatCtx->barK); CHKERRQ(ierr);
    ierr = MatCopy(shellMatCtx->K,shellMatCtx->barK,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = MatAXPY(shellMatCtx->barK,ts_a,shellMatCtx->M,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(shellMatCtx->barK,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(shellMatCtx->barK,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    //Matrix View
    //MatView(shellMatCtx->barK,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    //std::string wait;
    //std::cin >> wait;

    PetscFunctionReturn(0);
  }

  #endif //__DIRICHLETBC_HPP__
