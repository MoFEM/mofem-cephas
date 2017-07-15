/**
 * \file PoissonOperators.hpp
 * \example PoissonOperators.hpp
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

#ifndef __POISSONOPERATORS_HPP__
#define __POISSONOPERATORS_HPP__

namespace PoissonExample {

  /**
   * \brief Calculate the grad-grad operator and assemble matrix
   *
   * Calculate
   * \f[
   * \mathbf{K}=\int_\Omega \nabla \boldsymbol\phi \cdot \nabla \boldsymbol\phi \textrm{d}\Omega
   * \f]
   * and assemble to global matrix.
   *
   * This operator is executed on element for each unique combination of entities.
   *
   */
  struct OpGradGrad: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    OpGradGrad():
    VolumeElementForcesAndSourcesCore::UserDataOperator("U","U",OPROWCOL,true) {
    }

    /**
     * \brief Do calculations for give operator
     * @param  row_side row side number (local number) of entity on element
     * @param  col_side column side number (local number) of entity on element
     * @param  row_type type of row entity MBVERTEX, MBEDGE, MBTRI or MBTET
     * @param  col_type type of column entity MBVERTEX, MBEDGE, MBTRI or MBTET
     * @param  row_data data for row
     * @param  col_data data for column
     * @return          error code
     */
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // get number of dofs on row
      nbRows = row_data.getIndices().size();
      // if no dofs on row, exit that work, nothing to do here
      if(!nbRows) PetscFunctionReturn(0);
      // get number of dofs on column
      nbCols = col_data.getIndices().size();
      // if no dofs on Columbia, exit nothing to do here
      if(!nbCols) PetscFunctionReturn(0);
      // get number of integration points
      nbIntegrationPts = getGaussPts().size2();
      // chekk if entity block is on matrix diagonal
      if(
        row_side==col_side&&
        row_type==col_type
      ) {
        isDiag = true; // yes, it is
      } else {
        isDiag = false;
      }
      // integrate local matrix for entity block
      ierr = iNtegrte(row_data,col_data); CHKERRQ(ierr);
      // asseble local matrix
      ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  private:

    PetscErrorCode ierr;  ///< error code

    int nbRows;           ///< number of dofs on rows
    int nbCols;           ///< number if dof on column
    int nbIntegrationPts; ///< number of integration points
    bool isDiag;          ///< true if this block is on diagonal

    FTensor::Index<'i',3> i;  ///< summit Index
    MatrixDouble locMat;      ///< local entity block matrix

    /**
     * \brief Integrate grad-grad operator
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
    inline PetscErrorCode iNtegrte(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // set size of local entity bock
      locMat.resize(nbRows,nbCols,false);
      // clear matrux
      locMat.clear();
      // get element volume
      double vol = getVolume();
      // get integration weigths
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      // get base function gradient on rows
      FTensor::Tensor1<double*,3> t_row_grad = row_data.getFTensor1DiffN<3>();
      // loop over integration points
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        // take into account Jacobean
        const double alpha = t_w*vol;
        // noalias(locMat) += alpha*prod(row_data.getDiffN(gg),trans(col_data.getDiffN(gg)));
        // take fist element to local matrix
        FTensor::Tensor0<double*> a(&*locMat.data().begin());
        // loop over rows base functions
        for(int rr = 0;rr!=nbRows;rr++) {
          // get column base functions gradient at gauss point gg
          FTensor::Tensor1<double*,3> t_col_grad = col_data.getFTensor1DiffN<3>(gg,0);
          // loop over columbs
          for(int cc = 0;cc!=nbCols;cc++) {
            // calculate element of loacl matrix
            a += alpha*(t_row_grad(i)*t_col_grad(i));
            ++t_col_grad; // move to another gradient of base function on column
            ++a;  // move to another element of local matrix in column
          }
          ++t_row_grad; // move to another element of gradient of base function on row
        }
        ++t_w; // move to another integration weight
      }
      PetscFunctionReturn(0);
    }

    /**
     * \brief Assemble local entity block matrix
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
    inline PetscErrorCode aSsemble(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // get pointer to first global index on row
      const int* row_indices = &*row_data.getIndices().data().begin();
      // get pointer to first global index on column
      const int* col_indices = &*col_data.getIndices().data().begin();
      // assemble local matrix
      ierr = MatSetValues(
        getFEMethod()->ksp_B,
        nbRows,row_indices,
        nbCols,col_indices,
        &*locMat.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
      if(!isDiag) {
        // if not diagonal term and since global matrix is symmetric assemble
        // transpose term.
        locMat = trans(locMat);
        ierr = MatSetValues(
          getFEMethod()->ksp_B,
          nbCols,col_indices,
          nbRows,row_indices,
          &*locMat.data().begin(),ADD_VALUES
        ); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief template class for integration oh the right hand side
   */
  template<typename OPBASE>
  struct OpBaseRhs: public OPBASE {

    OpBaseRhs(const std::string field_name):
    OPBASE(field_name,OPBASE::OPROW) {
    }

    /**
     * \brief This function is called by finite element
     *
     * Do work is composed from two operations, integrate and assembly. Also,
     * it set values nbRows, and nbIntegrationPts.
     *
     */
    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;
      // get number of dofs on row
      nbRows = row_data.getIndices().size();
      if(!nbRows) PetscFunctionReturn(0);
      // get number of integration points
      nbIntegrationPts = OPBASE::getGaussPts().size2();
      // integrate local vector
      ierr = iNtegrte(row_data); CHKERRQ(ierr);
      // assemble local vector
      ierr = aSsemble(row_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    /**
     * \brief Class dedicated to integrate operator
     * @param  data entity data on element row
     * @return      error code
     */
    virtual PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) = 0;

    /**
     * \brief Class dedicated to assemble operator to global system vector
     * @param  data entity data (indices, base functions, etc. ) on element row
     * @return      error code
     */
    virtual PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) = 0;

  protected:

    PetscErrorCode ierr;    ///< error code
    int nbRows;             ///< number of dofs on row
    int nbIntegrationPts;   ///< number of integration points

  };

  /**
   * \brief Operator calculate source term,
   *
   * \f[
   * \mathbf{F} = \int_\Omega \boldsymbol\phi f \textrm{d}\Omega
   * \f]
   *
   */
  struct OpVF: public OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator> {

    typedef boost::function<double (const double,const double,const double)> FSource;

    OpVF(FSource f_source):
    OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator>("U"),
    fSource(f_source) {
    }

  private:

    PetscErrorCode ierr;
    FTensor::Number<0> NX;
    FTensor::Number<1> NY;
    FTensor::Number<2> NZ;
    FSource fSource;

    VectorDouble locVec;

    /**
     * \brief Integrate local entity vector
     * @param  data entity data on element row
     * @return      error code
     */
    PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      // set size of local vector
      locVec.resize(nbRows,false);
      // clear local entity vector
      locVec.clear();
      // get finite element volume
      double vol = getVolume();
      // get integration weights
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      // get base functions on entity
      FTensor::Tensor0<double*> t_v = data.getFTensor0N();
      // get coordinates at integration points
      FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
      // loop over all integration points
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        // evaluate constant term
        const double alpha = vol*t_w*fSource(t_coords(NX),t_coords(NY),t_coords(NZ));
        // get element of local vector
        FTensor::Tensor0<double*> t_a(&*locVec.data().begin());
        // loop over base functions
        for(int rr = 0;rr!=nbRows;rr++) {
          // add to local vector source term
          t_a -= alpha*t_v;
          ++t_a;  // move to next element of local vector
          ++t_v;  // move to next base function
        }
        ++t_w;  // move to next integration weight
        ++t_coords; // move to next physical coordinates at integration point
      }
      PetscFunctionReturn(0);
    }

    /**
     * \brief assemble local entity vector to the global right hand side
     * @param  data entity data, i.e. global indices of local vector
     * @return      error code
     */
    PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      // get global indices of local vector
      const int* indices = &*data.getIndices().data().begin();
      // get values from local vector
      const double* vals = &*locVec.data().begin();
      // assemble vector
      ierr = VecSetValues(
        getFEMethod()->ksp_f,nbRows,indices,vals,ADD_VALUES
      ); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief Calculate constrains matrix
   *
   * \f[
   * \mathbf{C} = \int_{\partial\Omega} \boldsymbol\psi \boldsymbol\phi \textrm{d}\partial\Omega
   * \f]
   * where \f$\lambda \f$ is base function on boundary
   *
   */
  struct OpLU: public FaceElementForcesAndSourcesCore::UserDataOperator {

    OpLU(const bool assemble_transpose):
    FaceElementForcesAndSourcesCore::UserDataOperator("L","U",OPROWCOL,false),
    assembleTraspose(assemble_transpose) {
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // get number of dofs on row
      nbRows = row_data.getIndices().size();
      // exit here if no dofs on row, nothing to do
      if(!nbRows) PetscFunctionReturn(0);
      // get number of dofs on column,
      nbCols = col_data.getIndices().size();
      // exit here if no dofs on roe, nothing to do
      if(!nbCols) PetscFunctionReturn(0);
      // get number of integration points
      nbIntegrationPts = getGaussPts().size2();
      // integrate local constrains matrix
      ierr = iNtegrte(row_data,col_data); CHKERRQ(ierr);
      // assemble local constrains matrix
      ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  private:

    PetscErrorCode ierr; ///< error code

    int nbRows;            ///< number of dofs on row
    int nbCols;            ///< number of dofs on column
    int nbIntegrationPts;  ///< number of integration points
    const bool assembleTraspose;  ///< assemble transpose, i.e. CT if set to true

    MatrixDouble locMat;   ///< local constrains matrxi

    /** \brief Integrate local constrains matrix
     */
    inline PetscErrorCode iNtegrte(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // set size of local constrains matrix
      locMat.resize(nbRows,nbCols,false);
      // clear matrix
      locMat.clear();
      // get area of element
      const double area = getArea();
      // get integration weights
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      // get base functions on entity
      FTensor::Tensor0<double*> t_row = row_data.getFTensor0N();
      // run over integration points
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        const double alpha = area*t_w;
        // get element of local matrix
        FTensor::Tensor0<double*> c(&*locMat.data().begin());
        // run over base functions on rows
        for(int rr = 0;rr!=nbRows;rr++) {
          // get first base functions on column for integration point gg
          FTensor::Tensor0<double*> t_col = col_data.getFTensor0N(gg,0);
          // run over base function on column
          for(int cc = 0;cc!=nbCols;cc++) {
            // integrate element of constrains matrix
            c += alpha*t_row*t_col;
            ++t_col; // move to next base function on column
            ++c;  // move to next element of local matrix
          }
          ++t_row; // move to next base function on row
        }
        ++t_w; // move to next integrate weight
      }
      PetscFunctionReturn(0);
    }

    /**
     * \brief integrate local constrains matrix
     */
    inline PetscErrorCode aSsemble(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // get indices on row
      const int* row_indices = &*row_data.getIndices().data().begin();
      // get indices on column
      const int* col_indices = &*col_data.getIndices().data().begin();
      // assemble local matrix
      ierr = MatSetValues(
        getFEMethod()->ksp_B,
        nbRows,row_indices,
        nbCols,col_indices,
        &*locMat.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
      // cerr << locMat << endl;
      if(assembleTraspose) {
        // assmble transpose of local matrix
        locMat = trans(locMat);
        ierr = MatSetValues(
          getFEMethod()->ksp_B,
          nbCols,col_indices,
          nbRows,row_indices,
          &*locMat.data().begin(),ADD_VALUES
        ); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief Assemble constrains vector
   *
   * \f[
   * \mathbf{g} = \int_{\partial\Omega} \boldsymbol\psi \overline{u} \textrm{d}\partial\Omega
   * \f]
   *
   */
  struct OpLU_exact: public OpBaseRhs<FaceElementForcesAndSourcesCore::UserDataOperator> {

    typedef boost::function<double (const double,const double,const double)> FVal;

    OpLU_exact(FVal f_value,const string field_name = "L",const double beta = 1):
    OpBaseRhs<FaceElementForcesAndSourcesCore::UserDataOperator>(field_name),
    fValue(f_value),
    bEta(beta) {
    }

  private:

    FTensor::Number<0> NX; ///< x-direction index
    FTensor::Number<1> NY; ///< y-direction index
    FTensor::Number<2> NZ; ///< z-direction index
    FVal fValue;           ///< Function pointer evaluating values of "U" at the boundary

    VectorDouble locVec;
    const double bEta;

    /**
     * \brief Integrate local constrains vector
     */
    PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      // set size to local vector
      locVec.resize(nbRows,false);
      // clear loacl vector
      locVec.clear();
      // get face area
      const double area = getArea()*bEta;
      // get integration wiegth
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      // get base function
      FTensor::Tensor0<double*> t_l = data.getFTensor0N();
      // get coordinates at integration point
      FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
      // make loop over integration points
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        // evalue function on boundary and scale it by area and integration weight
        double alpha = area*t_w*fValue(t_coords(NX),t_coords(NY),t_coords(NZ));
        // get element of vector
        FTensor::Tensor0<double*> t_a(&*locVec.data().begin());
        //
        for(int rr = 0;rr!=nbRows;rr++) {
          t_a += alpha*t_l;
          ++t_a;
          ++t_l;
        }
        ++t_w;
        ++t_coords;
      }
      PetscFunctionReturn(0);
    }

    /**
     * \brief assemble constrains vectors
     */
    PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      const int* indices = &*data.getIndices().data().begin();
      const double* vals = &*locVec.data().begin();
      ierr = VecSetValues(
        getFEMethod()->ksp_f,nbRows,indices,&*vals,ADD_VALUES
      ); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief Evaluate error
   */
  struct OpError: public OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator> {

    typedef boost::function<double (const double,const double,const double)> UVal;
    typedef boost::function<FTensor::Tensor1<double,3> (const double,const double,const double)> GVal;

    OpError(
      UVal u_value,
      GVal g_value,
      boost::shared_ptr<VectorDouble>& field_vals,
      boost::shared_ptr<MatrixDouble>& grad_vals,
      Vec global_error
    ):
    OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator>("ERROR"),
    globalError(global_error),
    uValue(u_value),
    gValue(g_value),
    fieldVals(field_vals),
    gradVals(grad_vals) {
    }

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;
      nbRows = row_data.getFieldData().size();
      if(!nbRows) PetscFunctionReturn(0);
      nbIntegrationPts = getGaussPts().size2();
      ierr = iNtegrte(row_data); CHKERRQ(ierr);
      ierr = aSsemble(row_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  private:

    Vec globalError;  ///< ghost vector with global (integrated over volume) error

    FTensor::Number<0> NX;
    FTensor::Number<1> NY;
    FTensor::Number<2> NZ;
    FTensor::Index<'i',3> i;
    UVal uValue;  ///< function with exact solution
    GVal gValue;  ///< function with exact solution for gradient

    boost::shared_ptr<VectorDouble> fieldVals;
    boost::shared_ptr<MatrixDouble> gradVals;

    /**
     * \brief Integrate error
     */
    PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      // clear field dofs
      data.getFieldData().clear();
      // get volume of element
      const double vol = getVolume();
      // get integration weight
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      // get solution at integration point
      FTensor::Tensor0<double*> t_u = getTensor0FormData(*fieldVals);
      // get solution at integration point
      FTensor::Tensor1<double*,3> t_grad = getTensor1FormData<3>(*gradVals);
      // get coordinates at integration point
      FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
      // keep exact gradient and error or gradient
      FTensor::Tensor1<double,3> t_exact_grad,t_error_grad;
      // integrate over
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        double alpha = vol*t_w;
        // evalue exact value
        double exact_u = uValue(t_coords(NX),t_coords(NY),t_coords(NZ));
        // evalue exact hradient
        t_exact_grad = gValue(t_coords(NX),t_coords(NY),t_coords(NZ));
        // calculate gradient errro
        t_error_grad(i) = t_grad(i)-t_exact_grad(i);
        // error
        double error = pow(t_u-exact_u,2)+t_error_grad(i)*t_error_grad(i);
        // iterate over base functions
        data.getFieldData()[0] += alpha*error;
        ++t_w;      // move to next integration point
        ++t_u;      // next value of function at integration point
        ++t_grad;   // next gradient at integration point
        ++t_coords; // next coordinate at integration point
      }
      PetscFunctionReturn(0);
    }

    /**
     * \brief Assemble error
     */
    PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      // set error on mesh
      data.getFieldDofs()[0]->getFieldData() = sqrt(data.getFieldData()[0]);
      // assemble vector to global error
      ierr = VecSetValue(globalError,0,data.getFieldData()[0],ADD_VALUES); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  /**
   * \brief Set integration rule to volume elements
   *
   * This rule is used to integrate \f$\nabla v \cdot \nabla u\f$, thus
   * if approximation field and testing field are polynomial order "p",
   * then rule for exact integration is 2*(p-1).
   *
   * Integration rule is order of polynomial which is calculated exactly. Finite element
   * selects integration method based on return of this function.
   *
   */
  struct VolRule {
    int operator()(int,int,int p) const {
      return 2*(p-1);
    }
  };

  /**
   * \brief Set integration rule to boundary elements
   *
   * This is uses to integrate values on the face. Is used to integrate
   * \f$(\mathbf{n} \cdot \lambda) u\f$, where Lagrange multiplayer
   * is order "p_row" and approximate function is order "p_col".
   *
   * Integration rule is order of polynomial which is calculated exactly. Finite element
   * selects integration method based on return of this function.
   *
   */
  struct FaceRule {
    int operator()(int p_row,int p_col,int p_data) const {
      return 2*p_data+1;
    }
  };

  struct OpUU: public FaceElementForcesAndSourcesCore::UserDataOperator {

    OpUU(const bool beta = 1):
    FaceElementForcesAndSourcesCore::UserDataOperator("U","U",OPROWCOL,true),
    bEta(beta) {
    }

    /**
     * \brief Do calculations for give operator
     * @param  row_side row side number (local number) of entity on element
     * @param  col_side column side number (local number) of entity on element
     * @param  row_type type of row entity MBVERTEX, MBEDGE, MBTRI or MBTET
     * @param  col_type type of column entity MBVERTEX, MBEDGE, MBTRI or MBTET
     * @param  row_data data for row
     * @param  col_data data for column
     * @return          error code
     */
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // get number of dofs on row
      nbRows = row_data.getIndices().size();
      // if no dofs on row, exit that work, nothing to do here
      if(!nbRows) PetscFunctionReturn(0);
      // get number of dofs on column
      nbCols = col_data.getIndices().size();
      // if no dofs on Columbia, exit nothing to do here
      if(!nbCols) PetscFunctionReturn(0);
      // get number of integration points
      nbIntegrationPts = getGaussPts().size2();
      // chekk if entity block is on matrix diagonal
      if(
        row_side==col_side&&
        row_type==col_type
      ) {
        isDiag = true; // yes, it is
      } else {
        isDiag = false;
      }
      // integrate local matrix for entity block
      ierr = iNtegrte(row_data,col_data); CHKERRQ(ierr);
      // asseble local matrix
      ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  private:

    const double bEta;

    PetscErrorCode ierr;  ///< error code

    int nbRows;           ///< number of dofs on rows
    int nbCols;           ///< number if dof on column
    int nbIntegrationPts; ///< number of integration points
    bool isDiag;          ///< true if this block is on diagonal

    FTensor::Index<'i',3> i;  ///< summit Index
    MatrixDouble locMat;      ///< local entity block matrix

    /**
     * \brief Integrate grad-grad operator
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
    inline PetscErrorCode iNtegrte(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // set size of local entity bock
      locMat.resize(nbRows,nbCols,false);
      // clear matrux
      locMat.clear();
      // get element area
      double area = getArea()*bEta;
      // get integration weigths
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      // get base function gradient on rows
      FTensor::Tensor0<double*> t_row_base = row_data.getFTensor0N();
      // loop over integration points
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        // take into account Jacobean
        const double alpha = t_w*area;
        // noalias(locMat) += alpha*prod(row_data.getDiffN(gg),trans(col_data.getDiffN(gg)));
        // take fist element to local matrix
        FTensor::Tensor0<double*> a(&*locMat.data().begin());
        // loop over rows base functions
        for(int rr = 0;rr!=nbRows;rr++) {
          // get column base functions gradient at gauss point gg
          FTensor::Tensor0<double*> t_col_base = col_data.getFTensor0N(gg,0);
          // loop over columbs
          for(int cc = 0;cc!=nbCols;cc++) {
            // calculate element of loacl matrix
            a += alpha*t_row_base*t_col_base;
            ++t_col_base; // move to another gradient of base function on column
            ++a;  // move to another element of local matrix in column
          }
          ++t_row_base; // move to another element of gradient of base function on row
        }
        ++t_w; // move to another integration weight
      }
      PetscFunctionReturn(0);
    }

    /**
     * \brief Assemble local entity block matrix
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
    inline PetscErrorCode aSsemble(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      // get pointer to first global index on row
      const int* row_indices = &*row_data.getIndices().data().begin();
      // get pointer to first global index on column
      const int* col_indices = &*col_data.getIndices().data().begin();
      // assemble local matrix
      ierr = MatSetValues(
        getFEMethod()->ksp_B,
        nbRows,row_indices,
        nbCols,col_indices,
        &*locMat.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
      if(!isDiag) {
        // if not diagonal term and since global matrix is symmetric assemble
        // transpose term.
        locMat = trans(locMat);
        ierr = MatSetValues(
          getFEMethod()->ksp_B,
          nbCols,col_indices,
          nbRows,row_indices,
          &*locMat.data().begin(),ADD_VALUES
        ); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };


  /**
   * \brief Create finite elements instances
   *
   * Create finite element instances and add operators to finite elements.
   *
   */
  struct CreateFiniteElementes {

    CreateFiniteElementes(MoFEM::Interface &m_field):
    mField(m_field) {
    }

    /**
     * \brief Create finite element to calculate matrix and vectors
     */
    PetscErrorCode createFEToAssmbleMatrceAndVector(
      boost::function<double (const double,const double,const double)> f_u,
      boost::function<double (const double,const double,const double)> f_source,
      boost::shared_ptr<ForcesAndSurcesCore>& domain_lhs_fe,
      boost::shared_ptr<ForcesAndSurcesCore>& boundary_lhs_fe,
      boost::shared_ptr<ForcesAndSurcesCore>& domain_rhs_fe,
      boost::shared_ptr<ForcesAndSurcesCore>& boundary_rhs_fe,
      bool trans = true
    ) const {
      PetscFunctionBegin;

      // Create elements element instances
      domain_lhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(mField));
      boundary_lhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(mField));
      domain_rhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new VolumeElementForcesAndSourcesCore(mField));
      boundary_rhs_fe = boost::shared_ptr<ForcesAndSurcesCore>(new FaceElementForcesAndSourcesCore(mField));

      // Set integration rule to elements instances
      domain_lhs_fe->getRuleHook = VolRule();
      domain_rhs_fe->getRuleHook = VolRule();
      boundary_lhs_fe->getRuleHook = FaceRule();
      boundary_rhs_fe->getRuleHook = FaceRule();

      // Ass operators to element instances
      // Add operator grad-grad for calualte matrix
      domain_lhs_fe->getOpPtrVector().push_back(new OpGradGrad());
      // Add operator to calculate source terms
      domain_rhs_fe->getOpPtrVector().push_back(new OpVF(f_source));
      // Add operator calculating constrains matrix
      boundary_lhs_fe->getOpPtrVector().push_back(new OpLU(trans));
      // Add operator calculating constrains vector
      boundary_rhs_fe->getOpPtrVector().push_back(new OpLU_exact(f_u));

      PetscFunctionReturn(0);
    }

    /**
     * \brief Create finite element to calculate error
     */
    PetscErrorCode createFEToEvaluateError(
      boost::function<double (const double,const double,const double)> f_u,
      boost::function<FTensor::Tensor1<double,3> (const double,const double,const double)> g_u,
      Vec global_error,
      boost::shared_ptr<ForcesAndSurcesCore>& domain_error
    ) const {
      PetscFunctionBegin;
      // Create finite element instance to calualte error
      domain_error = boost::shared_ptr<ForcesAndSurcesCore>(
        new VolumeElementForcesAndSourcesCore(mField)
      );
      domain_error->getRuleHook = VolRule();
      // Set integration rule
      // Crate shared vector storing values of field "u" on integration points on element. element
      // is local and is used to exchange data between operators.
      boost::shared_ptr<VectorDouble> values_at_integation_ptr = boost::make_shared<VectorDouble>();
      // Storing gradients of field
      boost::shared_ptr<MatrixDouble> grad_at_integation_ptr = boost::make_shared<MatrixDouble>();
      // Add default operator to calculate field values at integration points
      domain_error->getOpPtrVector().push_back(new OpCalculateScalarFieldValues("U",values_at_integation_ptr));
      // Add default operator to calculate field gradient at integration points
      domain_error->getOpPtrVector().push_back(new OpCalculateScalarFieldGradient<3>("U",grad_at_integation_ptr));
      // Add operator to integrate error element by element.
      domain_error->getOpPtrVector().push_back(
        new OpError(f_u,g_u,values_at_integation_ptr,grad_at_integation_ptr,global_error)
      );
      PetscFunctionReturn(0);
    }

    /**
     * \brief Create finite element to post-process results
     */
    PetscErrorCode creatFEToPostProcessResults(
      boost::shared_ptr<ForcesAndSurcesCore>& post_proc_volume
    ) const {
      PetscErrorCode ierr;
      PetscFunctionBegin;

      // Note that user can stack together arbitrary number of operators to compose
      // complex PDEs.

      // Post-process results. This is standard element, with functionality
      // enabling refining mesh for post-processing. In addition in PostProcOnRefMesh.hpp
      // are implanted set of  users operators to post-processing fields. Here
      // using simplified mechanism for post-processing finite element, we
      // add operators to save data from field on mesh tags for ParaView
      // visualization.
      post_proc_volume = boost::shared_ptr<ForcesAndSurcesCore>(new PostProcVolumeOnRefinedMesh(mField));
      // Add operators to the elements, starting with some generic
      ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
      generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
      addFieldValuesPostProc("U"); CHKERRQ(ierr);
      ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
      addFieldValuesPostProc("ERROR"); CHKERRQ(ierr);
      ierr = boost::static_pointer_cast<PostProcVolumeOnRefinedMesh>(post_proc_volume)->
      addFieldValuesGradientPostProc("U"); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  private:

    MoFEM::Interface &mField;

  };


}

#endif //__POISSONOPERATORS_HPP__
