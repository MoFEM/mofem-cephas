/** \file AnalyticalDirichlet.hpp

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

#ifndef __ANALYTICALDIRICHLETBC_HPP__
#define __ANALYTICALDIRICHLETBC_HPP__

using namespace boost::numeric;
using namespace MoFEM;

/** \brief Analytical Dirichlet boundary conditions
  \ingroup mofem_forces_and_sources
  */
struct AnalyticalDirichletBC {

  /** \brief finite element to approximate analytical solution on surface
    */
  struct ApproxField {


    struct MyTriFE: public FaceElementForcesAndSourcesCore {

      int addToRule; ///< this is add to integration rule if 2nd order geometry approximation
      MyTriFE(FieldInterface &m_field): FaceElementForcesAndSourcesCore(m_field),addToRule(1) {}
      int getRule(int order) { return order+addToRule; };

    };

    ApproxField(FieldInterface &m_field): feApprox(m_field) {}
    virtual ~ApproxField() {}

    MyTriFE feApprox;
    MyTriFE& getLoopFeApprox() { return feApprox; }

    ublas::matrix<double> hoCoords;
    struct OpHoCoord: public FaceElementForcesAndSourcesCore::UserDataOperator {

      ublas::matrix<double> &hoCoords;
      OpHoCoord(const string field_name,ublas::matrix<double> &ho_coords);

      PetscErrorCode doWork(
        int side,EntityType type,DataForcesAndSurcesCore::EntData &data
      );

    };


    /** \brief Lhs operator used to build matrix
      */
    struct OpLhs:public FaceElementForcesAndSourcesCore::UserDataOperator {

      ublas::matrix<double> &hoCoords;
      OpLhs(const string field_name,ublas::matrix<double> &ho_coords);

      ublas::matrix<FieldData> NN,transNN;
      PetscErrorCode doWork(
        int row_side,int col_side,
        EntityType row_type,EntityType col_type,
        DataForcesAndSurcesCore::EntData &row_data,
        DataForcesAndSurcesCore::EntData &col_data
      );

    };

    /** \brief Rhs operator used to build matrix
      */
    template<typename FUNEVAL>
    struct OpRhs:public FaceElementForcesAndSourcesCore::UserDataOperator {

      Range tRis;
      ublas::matrix<double> &hoCoords;
      boost::shared_ptr<FUNEVAL> functionEvaluator;
      int fieldNumber;

      OpRhs(const string field_name,Range tris,
        ublas::matrix<double> &ho_coords,
        boost::shared_ptr<FUNEVAL> function_evaluator,int field_number
      ):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
      tRis(tris),
      hoCoords(ho_coords),functionEvaluator(function_evaluator),
      fieldNumber(field_number)
      {

      }

      ublas::vector<FieldData> NTf;
      ublas::vector<DofIdx> iNdices;

      PetscErrorCode doWork(
        int side,EntityType type,DataForcesAndSurcesCore::EntData &data
      ) {
        PetscFunctionBegin;
        PetscErrorCode ierr;

        try {

          unsigned int nb_row = data.getIndices().size();
          if(nb_row==0) PetscFunctionReturn(0);
          if(tRis.find(getMoFEMFEPtr()->get_ent()) == tRis.end()) {
            PetscFunctionReturn(0);
          }

          const FENumeredDofMoFEMEntity *dof_ptr;
          ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
          unsigned int rank = dof_ptr->get_max_rank();

          NTf.resize(nb_row/rank);
          iNdices.resize(nb_row/rank);

          for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

            double x,y,z;
            double val = getGaussPts()(2,gg);
            if(hoCoords.size1() == data.getN().size1()) {
              double area = norm_2(getNormals_at_GaussPt(gg))*0.5;
              val *= area;
              x = hoCoords(gg,0);
              y = hoCoords(gg,1);
              z = hoCoords(gg,2);
            } else {
              val *= getArea();
              x = getCoordsAtGaussPts()(gg,0);
              y = getCoordsAtGaussPts()(gg,1);
              z = getCoordsAtGaussPts()(gg,2);
            }

            ublas::vector<double> a;
            try {

              a = (*functionEvaluator)(x,y,z)[fieldNumber];

            } catch (exception& ex) {
              ostringstream ss;
              ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
              SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
            }
            if(a.size()!=rank) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
            }


            for(unsigned int rr = 0;rr<rank;rr++) {

              ublas::noalias(iNdices) = ublas::vector_slice<ublas::vector<int> >
              (data.getIndices(), ublas::slice(rr, rank, data.getIndices().size()/rank));

              noalias(NTf) = data.getN(gg,nb_row/rank)*a[rr]*val;
              ierr = VecSetValues(getFEMethod()->snes_f,iNdices.size(),
              &iNdices[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);

            }

          }


        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }

        PetscFunctionReturn(0);
      }



    };

  };

  struct DirichletBC : public DisplacementBCFEMethodPreAndPostProc {

    DirichletBC(
      FieldInterface& m_field,const string &field,Mat A,Vec X,Vec F
    );

    DirichletBC(
      FieldInterface& m_field,const string &field
    );

    Range *tRis_ptr;

    PetscErrorCode iNitalize();
    PetscErrorCode iNitalize(Range &tris);

  };

  ApproxField approxField;
  AnalyticalDirichletBC(FieldInterface& m_field);

  template<typename FUNEVAL>
  PetscErrorCode setApproxOps(
    FieldInterface &m_field,
    string field_name,Range& tris,
    boost::shared_ptr<FUNEVAL> funtcion_evaluator,int field_number = 0,
    string nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      if(approxField.getLoopFeApprox().getOpPtrVector().empty()) {
        if(m_field.check_field(nodals_positions)) {
          approxField.getLoopFeApprox().getOpPtrVector().push_back(new ApproxField::OpHoCoord(nodals_positions,approxField.hoCoords));
        }
        approxField.getLoopFeApprox().getOpPtrVector().push_back(new ApproxField::OpLhs(field_name,approxField.hoCoords));
      }
      approxField.getLoopFeApprox().getOpPtrVector().push_back(new ApproxField::OpRhs<FUNEVAL>(field_name,tris,approxField.hoCoords,funtcion_evaluator,field_number));
      PetscFunctionReturn(0);
    }

  PetscErrorCode initializeProblem(
    FieldInterface &m_field,string fe,string field,Range& tris,
    string nodals_positions = "MESH_NODE_POSITIONS"
  );

  Mat A;
  Vec D,F;
  KSP kspSolver;
  PetscErrorCode setProblem(
    FieldInterface &m_field,string problem
  );

  PetscErrorCode solveProblem(
    FieldInterface &m_field,string problem,string fe,DirichletBC &bc,Range &tris
  );

  PetscErrorCode destroyProblem();


};

#endif //__ANALYTICALDIRICHLETBC_HPP__
