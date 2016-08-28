/** \file AnalyticalDirichlet.cpp

  Enforce Dirichlet boundary condition for given analytical function,

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
#include <MethodForForceScaling.hpp>
#include <DirichletBC.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <AnalyticalDirichlet.hpp>

AnalyticalDirichletBC::ApproxField::OpHoCoord::OpHoCoord(const std::string field_name,ublas::matrix<double> &ho_coords):
FaceElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
hoCoords(ho_coords) {}

PetscErrorCode AnalyticalDirichletBC::ApproxField::OpHoCoord::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  try {

    if(data.getFieldData().size()==0) PetscFunctionReturn(0);

    hoCoords.resize(data.getN().size1(),3);
    if(type == MBVERTEX) {
      hoCoords.clear();
    }

    int nb_dofs = data.getFieldData().size();
    for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
      for(int dd = 0;dd<3;dd++) {
        hoCoords(gg,dd) += cblas_ddot(nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
      }
    }

  } catch (const std::exception& ex) {
    std::ostringstream ss;
    ss << "throw in method: " << ex.what() << std::endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

AnalyticalDirichletBC::ApproxField::OpLhs::OpLhs(const std::string field_name,ublas::matrix<double> &ho_coords):
FaceElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROWCOL),
hoCoords(ho_coords)
{

}

PetscErrorCode AnalyticalDirichletBC::ApproxField::OpLhs::doWork(
  int row_side,int col_side,
  EntityType row_type,EntityType col_type,
  DataForcesAndSurcesCore::EntData &row_data,
  DataForcesAndSurcesCore::EntData &col_data
) {
  PetscFunctionBegin;

  if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
  if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  const FENumeredDofEntity *dof_ptr;
  ierr = getNumeredEntFiniteElementPtr()->getRowDofsByPetscGlobalDofIdx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
  int rank = dof_ptr->getNbOfCoeffs();

  int nb_row_dofs = row_data.getIndices().size()/rank;
  int nb_col_dofs = col_data.getIndices().size()/rank;

  NN.resize(nb_row_dofs,nb_col_dofs);
  NN.clear();

  unsigned int nb_gauss_pts = row_data.getN().size1();
  for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {

    double w = getGaussPts()(2,gg);
    if(hoCoords.size1() == row_data.getN().size1()) {

      // higher order element
      double area = norm_2(getNormalsAtGaussPt(gg))*0.5;
      w *= area;

    } else {

      //linear element
      w *= getArea();

    }

    cblas_dger(CblasRowMajor,
      nb_row_dofs,nb_col_dofs,
      w,&row_data.getN()(gg,0),1,&col_data.getN()(gg,0),1,
      &*NN.data().begin(),nb_col_dofs);

    }

    if( (row_type != col_type) || (row_side != col_side) ) {
      transNN.resize(nb_col_dofs,nb_row_dofs);
      ublas::noalias(transNN) = trans(NN);
    }

    double *data = &*NN.data().begin();
    double *trans_data = &*transNN.data().begin();
    ublas::vector<DofIdx> row_indices,col_indices;
    row_indices.resize(nb_row_dofs);
    col_indices.resize(nb_col_dofs);

    for(int rr = 0;rr < rank; rr++) {

      if((row_data.getIndices().size()%rank)!=0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }

      if((col_data.getIndices().size()%rank)!=0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }

      unsigned int nb_rows;
      unsigned int nb_cols;
      int *rows;
      int *cols;


      if(rank > 1) {

        ublas::noalias(row_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
        (row_data.getIndices(), ublas::slice(rr, rank, row_data.getIndices().size()/rank));
        ublas::noalias(col_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
        (col_data.getIndices(), ublas::slice(rr, rank, col_data.getIndices().size()/rank));

        nb_rows = row_indices.size();
        nb_cols = col_indices.size();
        rows = &*row_indices.data().begin();
        cols = &*col_indices.data().begin();

      } else {

        nb_rows = row_data.getIndices().size();
        nb_cols = col_data.getIndices().size();
        rows = &*row_data.getIndices().data().begin();
        cols = &*col_data.getIndices().data().begin();

      }

      if(nb_rows != NN.size1()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      if(nb_cols != NN.size2()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }

      ierr = MatSetValues(getFEMethod()->snes_B,nb_rows,rows,nb_cols,cols,data,ADD_VALUES); CHKERRQ(ierr);
      if( (row_type != col_type) || (row_side != col_side) ) {
        if(nb_rows != transNN.size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        if(nb_cols != transNN.size1()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        ierr = MatSetValues(getFEMethod()->snes_B,nb_cols,cols,nb_rows,rows,trans_data,ADD_VALUES); CHKERRQ(ierr);
      }

    }

    PetscFunctionReturn(0);
  }

  AnalyticalDirichletBC::DirichletBC::DirichletBC(
    MoFEM::Interface& m_field,const std::string &field,Mat A,Vec X,Vec F
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field,A,X,F),
  tRis_ptr(NULL)
  {

  }

  AnalyticalDirichletBC::DirichletBC::DirichletBC(
    MoFEM::Interface& m_field,const std::string &field
  ):
  DisplacementBCFEMethodPreAndPostProc(m_field,field),
  tRis_ptr(NULL) {

  }


  PetscErrorCode AnalyticalDirichletBC::DirichletBC::iNitalize() {
    PetscFunctionBegin;
    if(mapZeroRows.empty()) {
      if(tRis_ptr == NULL) {
        SETERRQ(PETSC_COMM_SELF,1,"need to initialised from AnalyticalDirichletBC::solveProblem");
      }
      ierr = iNitalize(*tRis_ptr); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode AnalyticalDirichletBC::DirichletBC::iNitalize(Range &tris) {
    PetscFunctionBegin;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    Range ents;
    rval = mField.get_moab().get_connectivity(tris,ents,true); CHKERRQ_MOAB(rval);
    ierr = mField.get_moab().get_adjacencies(tris,1,false,ents,moab::Interface::UNION); CHKERRQ(ierr);
    ents.merge(tris);
    for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof)) {
        mapZeroRows[dof->get()->getPetscGlobalDofIdx()] = dof->get()->getFieldData();
      }
    }
    dofsIndices.resize(mapZeroRows.size());
    dofsValues.resize(mapZeroRows.size());
    int ii = 0;
    std::map<DofIdx,FieldData>::iterator mit = mapZeroRows.begin();
    for(;mit!=mapZeroRows.end();mit++,ii++) {
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
    PetscFunctionReturn(0);
  }

  AnalyticalDirichletBC::AnalyticalDirichletBC(MoFEM::Interface& m_field): approxField(m_field) {};

  PetscErrorCode AnalyticalDirichletBC::initializeProblem(
    MoFEM::Interface &m_field,string fe,string field,Range& tris,string nodals_positions
  ) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = m_field.add_finite_element(fe,MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row(fe,field); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(fe,field); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(fe,field); CHKERRQ(ierr);
    if(m_field.check_field(nodals_positions)) {
      ierr = m_field.modify_finite_element_add_field_data(fe,nodals_positions); CHKERRQ(ierr);
    }
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,fe); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode AnalyticalDirichletBC::setProblem(
    MoFEM::Interface &m_field,string problem
  ) {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = m_field.VecCreateGhost(problem,ROW,&F); CHKERRQ(ierr);
    ierr = m_field.VecCreateGhost(problem,COL,&D); CHKERRQ(ierr);
    ierr = m_field.MatCreateMPIAIJWithArrays(problem,&A); CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&kspSolver); CHKERRQ(ierr);
    ierr = KSPSetOperators(kspSolver,A,A); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(kspSolver); CHKERRQ(ierr);

    //PC pc;
    //ierr = KSPGetPC(kspSolver,&pc); CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCLU); CHKERRQ(ierr);
    //ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS); CHKERRQ(ierr);
    //ierr = PCFactorSetUpMatSolverPackage(pc);  CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode AnalyticalDirichletBC::solveProblem(
    MoFEM::Interface &m_field,string problem,string fe,DirichletBC &bc,Range &tris
  ) {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = MatZeroEntries(A); CHKERRQ(ierr);

    approxField.getLoopFeApprox().snes_B = A;
    approxField.getLoopFeApprox().snes_f = F;
    ierr = m_field.loop_finite_elements(problem,fe,approxField.getLoopFeApprox()); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = KSPSolve(kspSolver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ierr = m_field.set_global_ghost_vector(problem,ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    bc.tRis_ptr = &tris;
    bc.mapZeroRows.clear();
    bc.dofsIndices.clear();
    bc.dofsValues.clear();

    PetscFunctionReturn(0);
  }

  PetscErrorCode AnalyticalDirichletBC::destroyProblem() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = KSPDestroy(&kspSolver); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
