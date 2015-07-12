/** \file ConvectiveMassElement.hpp
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

#ifndef __CONVECTIVE_MASS_ELEMENT_HPP
#define __CONVECTIVE_MASS_ELEMENT_HPP

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

/** \brief structure grouping operators and data used for calculation of mass (convective) element
  * \ingroup convective_mass_elem
  * \ingroup nonlinear_elastic_elem
  *
  * In order to assemble matrices and right hand vectors, the loops over
  * elements, entities over that elements and finally loop over integration
  * points are executed.
  *
  * Following implementation separate those three celeries of loops and to each
  * loop attach operator.
  *
  */
struct ConvectiveMassElement {

  /// \brief  definition of volume element
  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {

    Mat A;
    Vec F;
    bool initV; ///< check if ghost vector used to accumalte Kinetin energy is created

    MyVolumeFE(FieldInterface &m_field);

    /** \brief it is used to calculate nb. of Gauss integration points
     *
     * for more details pleas look
     *   Reference:
     *
     * Albert Nijenhuis, Herbert Wilf,
     * Combinatorial Algorithms for Computers and Calculators,
     * Second Edition,
     * Academic Press, 1978,
     * ISBN: 0-12-519260-6,
     * LC: QA164.N54.
     *
     * More details about algorithm
     * http://people.sc.fsu.edu/~jburkardt/cpp_src/gm_rule/gm_rule.html
    **/
      int getRule(int order);

    Vec V;
    double eNergy;

    PetscErrorCode preProcess();
    PetscErrorCode postProcess();


  };

  MyVolumeFE feMassRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeMassRhs() { return feMassRhs; } ///< get rhs volume element
  MyVolumeFE feMassLhs; ///< calculate left hand side for tetrahedral elements,i.e. mass element
  MyVolumeFE& getLoopFeMassLhs() { return feMassLhs; } ///< get lhs volume element
  MyVolumeFE feMassAuxLhs; ///< calculate left hand side for tetrahedral elements for Kuu shell matrix
  MyVolumeFE& getLoopFeMassAuxLhs() { return feMassAuxLhs; } ///< get lhs volume element for Kuu shell matrix

  MyVolumeFE feVelRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeVelRhs() { return feVelRhs; } ///< get rhs volume element
  MyVolumeFE feVelLhs; ///< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeVelLhs() { return feVelLhs; } ///< get lhs volume element

  MyVolumeFE feTRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeTRhs() { return feTRhs; } ///< get rhs volume element
  MyVolumeFE feTLhs; ///< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeTLhs() { return feTLhs; } ///< get lhs volume element

  MyVolumeFE feEnergy; ///< calculate kinetic energy
  MyVolumeFE& getLoopFeEnergy() { return feEnergy; } ///< get kinetic energy element

  FieldInterface &mField;
  short int tAg;

  ConvectiveMassElement(FieldInterface &m_field,short int tag);

  /** \brief data for calculation inertia forces
    * \ingroup user_modules
    */
  struct BlockData {
    double rho0; ///< reference density
    ublas::vector<double> a0; ///< constant acceleration
    Range tEts; ///< elements in block set
  };
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropriate BlockData

  /** \brief common data used by volume elements
    * \ingroup user_modules
    */
  struct CommonData {
    map<string,vector<ublas::vector<double> > > dataAtGaussPts;
    map<string,vector<ublas::matrix<double> > > gradAtGaussPts;
    string spatialPositions;
    string meshPositions;
    string spatialVelocities;
    vector<ublas::vector<double> > valVel;
    vector<vector<double*> > jacVelRowPtr;
    vector<ublas::matrix<double> > jacVel;
    vector<ublas::vector<double> > valMass;
    vector<vector<double*> > jacMassRowPtr;
    vector<ublas::matrix<double> > jacMass;
    vector<ublas::vector<double> > valT;
    vector<vector<double*> > jacTRowPtr;
    vector<ublas::matrix<double> > jacT;

  };
  CommonData commonData;

  struct OpGetDataAtGaussPts: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    vector<ublas::vector<double> > &valuesAtGaussPts;
    vector<ublas::matrix<double> > &gradientAtGaussPts;
    const EntityType zeroAtType;

    OpGetDataAtGaussPts(const string field_name,
      vector<ublas::vector<double> > &values_at_gauss_pts,
      vector<ublas::matrix<double> > &gardient_at_gauss_pts
    );

    /** \brief operator calculating deformation gradient
      *
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data
    );

  };

  struct OpGetCommonDataAtGaussPts: public OpGetDataAtGaussPts {
    OpGetCommonDataAtGaussPts(const string field_name,CommonData &common_data);
  };

  struct CommonFunctions {

    template<typename TYPE>
    PetscErrorCode dEterminatnt(ublas::matrix<TYPE> a,TYPE &det) {
      PetscFunctionBegin;
      //a11a22a33
      //+a21a32a13
      //+a31a12a23
      //-a11a32a23
      //-a31a22a13
      //-a21a12a33
      //http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
      //http://mathworld.wolfram.com/MatrixInverse.html
      det = a(0,0)*a(1,1)*a(2,2)
        +a(1,0)*a(2,1)*a(0,2)
        +a(2,0)*a(0,1)*a(1,2)
        -a(0,0)*a(2,1)*a(1,2)
        -a(2,0)*a(1,1)*a(0,2)
        -a(1,0)*a(0,1)*a(2,2);
      PetscFunctionReturn(0);
    }

    template<typename TYPE>
    PetscErrorCode iNvert(TYPE det,ublas::matrix<TYPE> a,ublas::matrix<TYPE> &inv_a) {
      PetscFunctionBegin;
      //PetscErrorCode ierr;
      inv_a.resize(3,3);
      //http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
      //http://mathworld.wolfram.com/MatrixInverse.html
      inv_a(0,0) = a(1,1)*a(2,2)-a(1,2)*a(2,1);
      inv_a(0,1) = a(0,2)*a(2,1)-a(0,1)*a(2,2);
      inv_a(0,2) = a(0,1)*a(1,2)-a(0,2)*a(1,1);
      inv_a(1,0) = a(1,2)*a(2,0)-a(1,0)*a(2,2);
      inv_a(1,1) = a(0,0)*a(2,2)-a(0,2)*a(2,0);
      inv_a(1,2) = a(0,2)*a(1,0)-a(0,0)*a(1,2);
      inv_a(2,0) = a(1,0)*a(2,1)-a(1,1)*a(2,0);
      inv_a(2,1) = a(0,1)*a(2,0)-a(0,0)*a(2,1);
      inv_a(2,2) = a(0,0)*a(1,1)-a(0,1)*a(1,0);
      inv_a /= det;
      PetscFunctionReturn(0);
    }

  };

  struct OpMassJacobian: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian;
    bool lInear;
    bool fieldDisp;

    ublas::vector<adouble> a,dot_W,dp_dt,a_res;
    ublas::matrix<adouble> h,H,invH,F,g,G;
    vector<double> active;
    OpMassJacobian(
      const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true,bool linear = false
    );

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpMassRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;

    OpMassRhs(const string field_name,BlockData &data,CommonData &common_data);

    ublas::vector<double> nf;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpMassLhs_dM_dv: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    Range forcesOnlyOnEntities;

    OpMassLhs_dM_dv(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr = NULL
    );

    ublas::matrix<double> k,jac;

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    );


  };

  struct OpMassLhs_dM_dx: public OpMassLhs_dM_dv {

    OpMassLhs_dM_dx(const string field_name,const string col_field,BlockData &data,CommonData &common_data);

    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct OpMassLhs_dM_dX: public OpMassLhs_dM_dv  {

    OpMassLhs_dM_dX(const string field_name,const string col_field,BlockData &data,CommonData &common_data);

    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct OpEnergy: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    Vec *Vptr;
    bool lInear;

    OpEnergy(const string field_name,BlockData &data,CommonData &common_data,Vec *v_ptr,bool linear = false);

    ublas::matrix<double> h,H,invH,F;
    ublas::vector<double> v;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpVelocityJacobian: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian,fieldDisp;

    OpVelocityJacobian(const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true);

    ublas::vector<adouble> a_res;
    ublas::vector<adouble> v,dot_w,dot_W;
    ublas::matrix<adouble> h,H,invH,F;
    ublas::vector<adouble> dot_u;
    adouble detH;

    vector<double> active;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpVelocityRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;

    OpVelocityRhs(const string field_name,BlockData &data,CommonData &common_data);

    ublas::vector<double> nf;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpVelocityLhs_dV_dv: public OpMassLhs_dM_dv {

    OpVelocityLhs_dV_dv(const string vel_field,const string field_name,BlockData &data,CommonData &common_data);

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct OpVelocityLhs_dV_dx: public OpVelocityLhs_dV_dv {

    OpVelocityLhs_dV_dx(const string vel_field,const string field_name,BlockData &data,CommonData &common_data);

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };


  struct OpVelocityLhs_dV_dX: public OpVelocityLhs_dV_dv {

    OpVelocityLhs_dV_dX(const string vel_field,const string field_name,BlockData &data,CommonData &common_data);

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct OpEshelbyDynamicMaterialMomentumJacobian: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian;
    bool fieldDisp;

    OpEshelbyDynamicMaterialMomentumJacobian(
      const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true
    );

    ublas::vector<adouble> a,v,a_T;
    ublas::matrix<adouble> g,H,invH,h,F,G;
    ublas::vector<double> active;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpEshelbyDynamicMaterialMomentumRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    Range forcesOnlyOnEntities;

    OpEshelbyDynamicMaterialMomentumRhs(
      const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
    );

    ublas::vector<double> nf;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    );

  };

  struct OpEshelbyDynamicMaterialMomentumLhs_dv: public OpMassLhs_dM_dv {

    OpEshelbyDynamicMaterialMomentumLhs_dv(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
    );

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct OpEshelbyDynamicMaterialMomentumLhs_dx: public OpEshelbyDynamicMaterialMomentumLhs_dv {

    OpEshelbyDynamicMaterialMomentumLhs_dx(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
    );

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct OpEshelbyDynamicMaterialMomentumLhs_dX: public OpEshelbyDynamicMaterialMomentumLhs_dv {

    OpEshelbyDynamicMaterialMomentumLhs_dX(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr
    );

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg);

  };

  struct UpdateAndControl: public FEMethod {

    FieldInterface& mField;
    TS tS;
    const string velocityField;
    const string spatialPositionField;

    int jacobianLag;
    UpdateAndControl(FieldInterface& m_field,TS _ts,
      const string velocity_field,
      const string spatial_position_field
    );

    PetscErrorCode preProcess();
    PetscErrorCode postProcess();

  };

  PetscErrorCode setBlocks();

  PetscErrorCode addConvectiveMassElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,BitRefLevel bit = BitRefLevel()
  );

  PetscErrorCode addVelocityElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,
    BitRefLevel bit = BitRefLevel());

  PetscErrorCode addEshelbyDynamicMaterialMomentum(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,
    BitRefLevel bit = BitRefLevel(),
    Range *intersected = NULL);

  PetscErrorCode setConvectiveMassOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,
    bool linear = false
  );

  PetscErrorCode setVelocityOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false
  );

  PetscErrorCode setKinematicEshelbyOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    Range *forces_on_entities_ptr = NULL
  );

  PetscErrorCode setShellMatrixMassOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool linear = false
  );

  struct MatShellCtx {

    Mat K,M;
    VecScatter scatterU,scatterV;
    double ts_a;//,scale;

    bool iNitialized;
    MatShellCtx();
    virtual ~MatShellCtx();

    Mat barK;
    Vec u,v,Ku,Mv;
    PetscErrorCode iNit();

    PetscErrorCode dEstroy();

    friend PetscErrorCode MultOpA(Mat A,Vec x,Vec f);
    friend PetscErrorCode ZeroEntriesOp(Mat A);

  };

  /** \brief Mult operator for shell matrix
    *
    * \f[
    \left[
    \begin{array}{cc}
    \mathbf{M} & \mathbf{K} \\
    \mathbf{I} & -\mathbf{I}a
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \mathbf{v} \\
    \mathbf{u}
    \end{array}
    \right] =
    \left[
    \begin{array}{c}
    \mathbf{r}_u \\
    \mathbf{r}_v
    \end{array}
    \right]
    * \f]
    *
    */
  static PetscErrorCode MultOpA(Mat A,Vec x,Vec f) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
    MatShellCtx *ctx = (MatShellCtx*)void_ctx;
    if(!ctx->iNitialized) {
      ierr = ctx->iNit(); CHKERRQ(ierr);
    }
    ierr = VecZeroEntries(f); CHKERRQ(ierr);
    //Mult Ku
    ierr = VecScatterBegin(ctx->scatterU,x,ctx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterU,x,ctx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatMult(ctx->K,ctx->u,ctx->Ku); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx->scatterU,ctx->Ku,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterU,ctx->Ku,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //Mult Mv
    ierr = VecScatterBegin(ctx->scatterV,x,ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterV,x,ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatMult(ctx->M,ctx->v,ctx->Mv); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx->scatterU,ctx->Mv,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterU,ctx->Mv,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //Velocities
    ierr = VecAXPY(ctx->v,-ctx->ts_a,ctx->u); CHKERRQ(ierr);
    //ierr = VecScale(ctx->v,ctx->scale); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx->scatterV,ctx->v,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterV,ctx->v,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //Assemble
    ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  static PetscErrorCode ZeroEntriesOp(Mat A) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
    MatShellCtx *ctx = (MatShellCtx*)void_ctx;
    ierr = MatZeroEntries(ctx->K); CHKERRQ(ierr);
    ierr = MatZeroEntries(ctx->M); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  struct PCShellCtx {

    Mat shellMat;
    bool initPC; ///< check if PC is initialized

    PCShellCtx(Mat shell_mat):
      shellMat(shell_mat),initPC(false) {
    }

    PC pC;

    PetscErrorCode iNit();

    PetscErrorCode dEstroy();

    friend PetscErrorCode PCShellSetUpOp(PC pc);
    friend PetscErrorCode PCShellDestroy(PC pc);
    friend PetscErrorCode PCShellApplyOp(PC pc,Vec f,Vec x);

  };

  static PetscErrorCode PCShellSetUpOp(PC pc) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
    PCShellCtx *ctx = (PCShellCtx*)void_ctx;
    ierr = ctx->iNit(); CHKERRQ(ierr);
    MatShellCtx *shell_mat_ctx;
    ierr = MatShellGetContext(ctx->shellMat,&shell_mat_ctx); CHKERRQ(ierr);
    ierr = PCSetFromOptions(ctx->pC); CHKERRQ(ierr);
    ierr = PCSetOperators(ctx->pC,shell_mat_ctx->barK,shell_mat_ctx->barK); CHKERRQ(ierr);
    ierr = PCSetUp(ctx->pC); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  static PetscErrorCode PCShellDestroy(PC pc) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
    PCShellCtx *ctx = (PCShellCtx*)void_ctx;
    ierr = ctx->dEstroy(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** \brief apply pre-conditioner for shell matrix
    *
    * \f[
    \left[
    \begin{array}{cc}
    \mathbf{M} & \mathbf{K} \\
    \mathbf{I} & -\mathbf{I}a
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \mathbf{v} \\
    \mathbf{u}
    \end{array}
    \right] =
    \left[
    \begin{array}{c}
    \mathbf{r}_u \\
    \mathbf{r}_v
    \end{array}
    \right]
    * \f]
    *
    * where \f$\mathbf{v} = \mathbf{r}_v + a\mathbf{u}\f$ and \f$\mathbf{u}=(a\mathbf{M}+\mathbf{K})^{-1}(\mathbf{r}_u - \mathbf{M}\mathbf{r}_v\f$.
    *
    */
  static PetscErrorCode PCShellApplyOp(PC pc,Vec f,Vec x) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
    PCShellCtx *ctx = (PCShellCtx*)void_ctx;
    MatShellCtx *shell_mat_ctx;
    ierr = MatShellGetContext(ctx->shellMat,&shell_mat_ctx); CHKERRQ(ierr);
    //forward
    ierr = VecScatterBegin(shell_mat_ctx->scatterU,f,shell_mat_ctx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterU,f,shell_mat_ctx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterBegin(shell_mat_ctx->scatterV,f,shell_mat_ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterV,f,shell_mat_ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //ierr = VecScale(shell_mat_ctx->v,1/shell_mat_ctx->scale); CHKERRQ(ierr);
    //apply pre-conditioner and calculate u
    ierr = MatMult(shell_mat_ctx->M,shell_mat_ctx->v,shell_mat_ctx->Mv); CHKERRQ(ierr); // Mrv
    ierr = VecAXPY(shell_mat_ctx->Ku,-1,shell_mat_ctx->Mv); CHKERRQ(ierr); // f-Mrv
    ierr = PCApply(ctx->pC,shell_mat_ctx->Ku,shell_mat_ctx->u); CHKERRQ(ierr); //u = (aM+K)^(-1)(ru-Mrv)
    //VecView(shell_mat_ctx->u,PETSC_VIEWER_STDOUT_WORLD);
    //calculate velocities
    ierr = VecAXPY(shell_mat_ctx->v,shell_mat_ctx->ts_a,shell_mat_ctx->u); CHKERRQ(ierr); // v = v + a*u
    //VecView(shell_mat_ctx->v,PETSC_VIEWER_STDOUT_WORLD);
    //reverse
    ierr = VecZeroEntries(x); CHKERRQ(ierr);
    ierr = VecScatterBegin(shell_mat_ctx->scatterU,shell_mat_ctx->u,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterU,shell_mat_ctx->u,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterBegin(shell_mat_ctx->scatterV,shell_mat_ctx->v,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterV,shell_mat_ctx->v,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  struct ShellResidualElement: public FEMethod {
    FieldInterface &mField;
    ShellResidualElement(FieldInterface &m_field);

    //variables bellow need to be set by user
    MatShellCtx *shellMatCtx; 					///< pointer to shell matrix

    PetscErrorCode preProcess();

    PetscErrorCode postProcess();

  };

  #ifdef __DIRICHLETBC_HPP__

  /** \brief blocked element/problem
    *
    * Blocked element run loops for different problem than TS problem. It is
    * used to calculate matrices of shell matrix.
    *
    */
  struct ShellMatrixElement: public FEMethod {

    FieldInterface &mField;
    ShellMatrixElement(FieldInterface &m_field);

    typedef pair<string,FEMethod*> LoopPairType;
    typedef vector<LoopPairType > LoopsToDoType;
    LoopsToDoType loopK; 	///< methods to calculate K shell matrix
    LoopsToDoType loopM; 	///< methods to calculate M shell matrix
    LoopsToDoType loopAuxM; 	///< methods to calculate derivatives of inertia forces over displacements shell matrix

    //variables bellow need to be set by user
    string problemName; 					///< name of shell problem
    MatShellCtx *shellMatCtx; 					///< pointer to shell matrix
    SpatialPositionsBCFEMethodPreAndPostProc *DirichletBcPtr; 	///< boundary conditions

    PetscErrorCode preProcess();

  };

  #endif //__DIRICHLETBC_HPP__

};


#endif //__CONVECTIVE_MASS_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup convective_mass_elem Mass Element
 * \ingroup user_modules
 ******************************************************************************/
