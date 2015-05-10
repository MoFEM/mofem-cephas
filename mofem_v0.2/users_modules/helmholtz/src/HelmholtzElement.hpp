/** \file HelmholtzElement.hpp
 \ingroup mofem_helmholtz_elem

 \brief Operators and data structures for wave propagation analyze (Galerkin Element)

 Implementation of Helmholtz element for wave propagation problem

 */

/*
  This work is part of PhD thesis by on Micro-fluids: Thomas Felix Xuan Meng
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>.
 * The header file should contains as least #include as possible for speed
 */

#ifndef __HELMHOLTZ_ELEMENT_HPP
#define __HELMHOLTZ_ELEMENT_HPP

/** \brief Structure grouping operators and data used for wave propagation problem
  \ingroup mofem_helmholtz_elem

  Complex-real transformation, See Ercegovac, Milos, and Jean-Michel Muller.
  "Solving Systems of Linear Equations in Complex Domain: Complex E-Method."
  (2007)"

  Three dimensional homogeneous isotropic medium and time harmonic Helmholtz operator
*/
struct HelmholtzElement {


  /// \brief  Volume element
  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
    int addToRank; ///< default value 1, i.e. assumes that geometry is approx. by quadratic functions.
    MyVolumeFE(FieldInterface &_mField,int add_to_rank): VolumeElementForcesAndSourcesCore(_mField),addToRank(add_to_rank) {}
    int getRule(int order) { return order+addToRank; };
  };

  /// \brief Surface element
  struct MySurfaceFE: public FaceElementForcesAndSourcesCore {
    int addToRank; ///< default value 1, i.e. assumes that geometry is approx. by quadratic functions.
    MySurfaceFE(FieldInterface &_mField,int add_to_rank): FaceElementForcesAndSourcesCore(_mField),addToRank(add_to_rank) {}
    int getRule(int order) { return order+addToRank; };
  };

  boost::ptr_map<string,ForcesAndSurcesCore> feRhs; // surface element for LHS
  boost::ptr_map<string,ForcesAndSurcesCore> feLhs; // surface element for RHS

  /** \brief Volume element data
  * \ingroup mofem_helmholtz_elem
  */
  struct VolumeData {
    double waveNumber;
    Range tEts; ///< contains elements in block set
  };
  map<int,VolumeData> volumeData;

  /** \brief Surface element data
  * \ingroup mofem_helmholtz_elem
  */
  struct SurfaceData {

    double aDmittance_real;
    double aDmittance_imag;
    Range tRis; ///< surface triangles where hate flux is applied

  };
  map<int,SurfaceData> surfaceIncidentWaveBcData;
  map<int,SurfaceData> sommerfeldBcData;
  map<int,SurfaceData> baylissTurkelBcData;

  struct GlobalParameters {
    pair<double,PetscBool> waveNumber;
    pair<double,PetscBool> surfaceAdmittance;
    pair<double,PetscBool> powerOfIncidentWaveReal;
    pair<double,PetscBool> powerOfIncidentWaveImag;

    pair<ublas::vector<double>,PetscBool> waveDirection;
  };
  GlobalParameters globalParameters;

  PetscErrorCode getGlobalParametersFromLineCommandOptions() {
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = PetscOptionsBegin(mField.get_comm(),NULL,"Helmholtz problem options","none"); CHKERRQ(ierr);

    globalParameters.waveNumber.first = 1;
    ierr = PetscOptionsReal("-wave_number","wave number","",
      globalParameters.waveNumber.first,
      &globalParameters.waveNumber.first,&globalParameters.waveNumber.second); CHKERRQ(ierr);
    if(!globalParameters.waveNumber.second) {

      SETERRQ(PETSC_COMM_SELF,1,"wave number not given, set in line command -wave_number to fix problem");

    }

    globalParameters.surfaceAdmittance.first = 0;
    ierr = PetscOptionsReal("-surface_admittance","surface admitance applied to all surface elements on MIX_INCIDENT_WAVE_BC","",
      globalParameters.surfaceAdmittance.first,
      &globalParameters.surfaceAdmittance.first,&globalParameters.surfaceAdmittance.second); CHKERRQ(ierr);

    globalParameters.powerOfIncidentWaveReal.first = 1;
    ierr = PetscOptionsReal("-power_of_incident_wave",
      "power of incident wave applied to all surface elements on MIX_INCIDENT_WAVE_BC and HARD_INCIDENT_WAVE_BC","",
      globalParameters.powerOfIncidentWaveReal.first,
      &globalParameters.powerOfIncidentWaveReal.first,&globalParameters.powerOfIncidentWaveReal.second); CHKERRQ(ierr);

    globalParameters.powerOfIncidentWaveImag.first = 0;
    globalParameters.powerOfIncidentWaveImag.second = PETSC_FALSE;

    globalParameters.waveDirection.first.resize(3);
    globalParameters.waveDirection.first.clear();
    globalParameters.waveDirection.first[2] = 1;
    int nmax = 3;
    ierr = PetscOptionsRealArray(
      "-wave_direction","direction of incident wave","",&globalParameters.waveDirection.first[0],&nmax,&globalParameters.waveDirection.second
    ); CHKERRQ(ierr);
    if(globalParameters.waveDirection.second) {
      if(nmax > 0 && nmax != 3) {

        SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"*** ERROR -wave_direction [3*1 vector] default:X direction [0,0,1]");

      }
    }

    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** \brief Common data used by volume and surface elements
  * \ingroup mofem_helmholtz_elem
  */
  struct CommonData {

    map<string,ublas::vector<double> > pressureAtGaussPts;
    map<string,ublas::matrix<double> > gradPressureAtGaussPts;
    ublas::matrix<double> hoCoords;

    map<EntityType, vector< ublas::vector<int> > > imIndices;

  };
  CommonData commonData;

  FieldInterface &mField;
  int addToRank; ///< default value 1, i.e. assumes that geometry is approx. by quadratic functions.

  HelmholtzElement(
    FieldInterface &m_field):
    mField(m_field),addToRank(1) {}

  struct OpGetImIndices: public ForcesAndSurcesCore::UserDataOperator  {

    CommonData &commonData;
    const string reFieldName,imFieldName;
    bool takeIndicesFromElementRowIndices;

    OpGetImIndices(
      const string re_field_name,const string im_field_name,CommonData &common_data):
      ForcesAndSurcesCore::UserDataOperator(re_field_name),
      commonData(common_data),
      reFieldName(re_field_name),imFieldName(im_field_name) {

      if(reFieldName!=imFieldName) {

        takeIndicesFromElementRowIndices = false;

      } else {

        takeIndicesFromElementRowIndices = true;

      }

    }

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      int nb_row_dofs = data.getIndices().size();
      if(nb_row_dofs==0) PetscFunctionReturn(0);

      if(!side) {
        commonData.imIndices[type].resize(6);
      }

      if(takeIndicesFromElementRowIndices) {

        (commonData.imIndices[type])[side] = data.getIndices();

      } else {

      	// Get rows and cols indices of imaginary part and assemble matrix.
      	// Note: However HELMHOLTZ_IMIM_FE element is not calculated, since
      	// matrix A on real elements is equal to matrix in imaginary elements,
      	// it need be to declared. Declaration indicate that on imaginary part,
      	// assembled matrix has non-zero values;
      	PetscErrorCode ierr;
      	ierr = getPorblemRowIndices(imFieldName,type,side,(commonData.imIndices[type])[side]); CHKERRQ(ierr);

      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Calculate pressure and gradient of pressure in volume
    */
  struct OpGetValueAndGradAtGaussPts: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    const string fieldName;
    OpGetValueAndGradAtGaussPts(const string field_name,
      CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),fieldName(field_name) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        int nb_dofs = data.getFieldData().size();
        if(nb_dofs==0) PetscFunctionReturn(0);
        int nb_gauss_pts = data.getN().size1();

        ublas::vector<double> &value = commonData.pressureAtGaussPts[fieldName];
        ublas::matrix<double> &gradient = commonData.gradPressureAtGaussPts[fieldName];

        // initialize
        value.resize(nb_gauss_pts);
        gradient.resize(nb_gauss_pts,3);
        if(type == MBVERTEX) {
          gradient.clear();
          value.clear();
        }

        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          value[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
          //ublas::noalias(ublas::matrix_row<ublas::matrix<double> >(gradient,gg)) +=
          //prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
          cblas_dgemv(CblasRowMajor,CblasTrans,
            nb_dofs,3,1,
            &data.getDiffN()(gg,0),3,
            &data.getFieldData()[0],1,
            1,&gradient(gg,0),1
          );

        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Calculate pressure on surface
    */
  struct OpGetValueAtGaussPts: public FaceElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    const string fieldName;
    OpGetValueAtGaussPts(const string field_name,
      CommonData &common_data):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),fieldName(field_name) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        int nb_dofs = data.getFieldData().size();
        if(nb_dofs==0) PetscFunctionReturn(0);
        int nb_gauss_pts = data.getN().size1();

        ublas::vector<double> &value = commonData.pressureAtGaussPts[fieldName];

        // Initialize
        value.resize(nb_gauss_pts);
        if(type == MBVERTEX) {
          value.clear();
        }

        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          value[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };


  /** \brief Coordinates at integration points, on surface.

    This takes into account HO approximation for geometry

    */
  struct OpHoCoordTri: public FaceElementForcesAndSourcesCore::UserDataOperator {

    ublas::matrix<double> &hoCoordsTri;
    OpHoCoordTri(const string field_name,ublas::matrix<double> &ho_coords):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      hoCoordsTri(ho_coords) {}

    /*
    Cartesian coordinates for integration points inside elements
    X^coordinates = DOF dot* N
    */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

        int nb_dofs = data.getFieldData().size();
        if(nb_dofs==0) PetscFunctionReturn(0);

        hoCoordsTri.resize(data.getN().size1(),3);
        if(type == MBVERTEX) {
          hoCoordsTri.clear();
        }

        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	        // Calculate x,y,z at each GaussPts
          for(int dd = 0;dd<3;dd++) {
            hoCoordsTri(gg,dd) += cblas_ddot(nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
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

  /** \brief Rhs vector for Helmholtz operator
    \ingroup mofem_helmholtz_elem

    \f[
    F_i = \int_{\Omega^e} \frac{\partial N_i}{\partial X_j} \frac{\partial p}{\partial X_j} - k^2 N_i p \textrm{d}V
    \f]

  */
  struct OpHelmholtzRhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    VolumeData &dAta;
    CommonData &commonData;
    const string fieldName;
    Vec F;

    OpHelmholtzRhs(
      const string field_name,Vec _F,VolumeData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),fieldName(field_name),F(_F) { }

    ublas::vector<double> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      int nb_row_dofs = data.getIndices().size();
      if(nb_row_dofs==0) PetscFunctionReturn(0);
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
        PetscFunctionReturn(0);
      }

      try {

        PetscErrorCode ierr;

        ublas::vector<double> &pressure = commonData.pressureAtGaussPts[fieldName];
        ublas::matrix<double> &grad_p = commonData.gradPressureAtGaussPts[fieldName];

        Nf.resize(nb_row_dofs);
        Nf.clear();

        // wave number "k" is the proportional to the frequency of incident wave
        // and represents number of waves per wave length 2Pi - 2Pi/K
        double k_pow2 = dAta.waveNumber*dAta.waveNumber;

        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

          double val = getVolume()*getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; // higher order geometry
          }

          /// Integrate diffN^T grad_p - k^2 N^T p dV

          //const ublas::matrix_row<ublas::matrix<double> > gard_p_at_gauss_pt(grad_p,gg);
          //ublas::noalias(Nf) += val*prod(data.getDiffN(gg,nb_row_dofs),gard_p_at_gauss_pt);
          cblas_dgemv(CblasRowMajor,CblasNoTrans,nb_row_dofs,3,val,
            &data.getDiffN()(gg,0),3,
            &grad_p(gg,0),1,
            1.,&Nf[0],1
          );
          ublas::noalias(Nf) -= val*k_pow2*data.getN(gg,nb_row_dofs)*pressure[gg];

        }

        ierr = VecSetValues(F,data.getIndices().size(),
        &data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Lhs for helmholtz operator

    \ingroup mofem_helmholtz_elem

    \f[
    A_{ik} = \int_{\Omega^e} \frac{\partial N_i}{\partial X_j} \frac{\partial N_k}{\partial X_j} - k^2 N_i N_k \textrm{d}V
    \f]

    */
  struct OpHelmholtzLhs: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    VolumeData &dAta;
    CommonData &commonData;
    const string imFieldName;
    Mat A;

    OpHelmholtzLhs(
      const string &re_field_name,const string &im_field_name,
      Mat _A,VolumeData &data,CommonData &common_data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(re_field_name,re_field_name),
      dAta(data),commonData(common_data),imFieldName(im_field_name),A(_A) {}

    ublas::matrix<double> K,transK;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      int nb_rows = row_data.getIndices().size();
      if(nb_rows==0) PetscFunctionReturn(0);
      int nb_cols = col_data.getIndices().size();
      if(nb_cols==0) PetscFunctionReturn(0);
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
        PetscFunctionReturn(0);
      }

      try {

        K.resize(nb_rows,nb_cols);
        K.clear();

        double k_pow2 = dAta.waveNumber*dAta.waveNumber;

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

          double val = getVolume()*getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
          }

          const double *diff_row_mat_ptr = &row_data.getDiffN()(gg,0);
          const double *diff_col_mat_ptr = &col_data.getDiffN()(gg,0);
          cblas_dgemm(
            CblasRowMajor,CblasNoTrans,CblasTrans,
            nb_rows,nb_cols,3,
            +val,diff_row_mat_ptr,3,
            diff_col_mat_ptr,3,
            1,&K(0,0),nb_cols
          );
          const double *row_mat_ptr = &row_data.getN()(gg,0);
          const double *col_mat_ptr = &col_data.getN()(gg,0);
          cblas_dger(CblasRowMajor,
            nb_rows,nb_cols,
            -val*k_pow2,
            row_mat_ptr,1,
            col_mat_ptr,1,
            &K(0,0),nb_cols
          );
          //noalias(K) -= val*k_pow2*outer_prod( row_data.getN(gg,nb_rows),col_data.getN(gg,nb_cols) );
          //noalias(K) += val*prod(row_data.getDiffN(gg,nb_rows),trans(col_data.getDiffN(gg,nb_cols)));
        }

        PetscErrorCode ierr;
        ierr = MatSetValues(
          A,
          nb_rows,&row_data.getIndices()[0],
          nb_cols,&col_data.getIndices()[0],
          &K(0,0),ADD_VALUES
        ); CHKERRQ(ierr);

        if(row_side != col_side || row_type != col_type) {
          transK.resize(nb_cols,nb_rows);
          noalias(transK) = trans(K);
          ierr = MatSetValues(
            A,
            nb_cols,&col_data.getIndices()[0],
            nb_rows,&row_data.getIndices()[0],
            &transK(0,0),ADD_VALUES
          ); CHKERRQ(ierr);
        }

        ierr = MatSetValues(
          A,
          nb_rows,&((commonData.imIndices[row_type])[row_side])[0],
          nb_cols,&((commonData.imIndices[col_type])[col_side])[0],
          &K(0,0),ADD_VALUES
        ); CHKERRQ(ierr);

        if(row_side != col_side || row_type != col_type) {
          ierr = MatSetValues(
            A,
            nb_cols,&((commonData.imIndices[col_type])[col_side])[0],
            nb_rows,&((commonData.imIndices[row_type])[row_side])[0],
            &transK(0,0),ADD_VALUES
          ); CHKERRQ(ierr);
        }


      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct ZeroFunVal {

    ublas::vector<double> vAl;

    ublas::vector<double>& operator()(double x,double y,double z,ublas::vector<double> &normal) {
      vAl.resize(2);
      vAl[0] = 0;
      vAl[1] = 0;
      return vAl;
    }

  };
  ZeroFunVal zeroFunVal;

  struct BaylissTurkel {

    BaylissTurkel() {}

    ublas::vector<double> vAl;
    ublas::vector<double>& operator()(double x,double y,double z,ublas::vector<double> &normal) {

      const complex< double > i( 0.0, 1.0 );
      double x2=x*x,y2=y*y,z2=z*z;
      double R = sqrt(x2+y2+z2);

      complex< double > result = 1.0/(2.0*R);

      vAl.resize(2);
      vAl[0] = std::real(result);
      vAl[1] = std::imag(result);

      return vAl;
    }

  };




  /** \brief Calculate incident wave scattered on hard surface

    \bug Assumes that normal sf surface pointing outward.
    */
  struct IncidentWaveNeumannF2 {

    GlobalParameters &globalParameters;
    IncidentWaveNeumannF2(GlobalParameters &global_parameters): globalParameters(global_parameters) {}

    ublas::vector<double> cOordinate;
    ublas::vector<double> vAl;

    ublas::vector<double>& operator()(double x,double y,double z,ublas::vector<double> &normal) {

      const complex< double > i( 0.0, 1.0 );

      cOordinate.resize(3);
      cOordinate[0] = x;
      cOordinate[1] = y;
      cOordinate[2] = z;

      double x1d = inner_prod(globalParameters.waveDirection.first,cOordinate);
      complex<double> power = globalParameters.powerOfIncidentWaveReal.first+i*globalParameters.powerOfIncidentWaveImag.first;
      complex<double> angle = globalParameters.waveNumber.first*(x1d);
      complex<double> p_inc = power*exp(i*angle);

      ublas::vector<complex<double > > grad(3);
      for(int ii = 0;ii<3;ii++) {
        grad[ii] = i*power*globalParameters.waveNumber.first*globalParameters.waveDirection.first[ii]*p_inc;
      }
      complex<double > grad_n = inner_prod(grad,normal);

      //// check if normal pointing to ceneter;
      //double dot = -inner_prod(normal,cOordinate);
      //if(dot < 0) grad_n *= -1;

      vAl.resize(2);
      vAl[0] = std::real(grad_n);
      vAl[1] = std::imag(grad_n);

      return vAl;
    }

  };

  /** \brief Rhs vector for Helmholtz operator
    \ingroup mofem_helmholtz_elem

    Operator is build using two template functions, see equations below.
    Depending on returning values of those functions user can apply, Nuemman, Mix or
    any variant of above conditions.

    \f[
    \left.\left\{ \frac{\partial p}{\partial n} +
      (\textrm{re}[\sigma]+i*\textrm{im}[\sigma] + \textrm{re}[F1]+i\textrm{im}[F1]) p
      + \textrm{re}[F2]+i\textrm{im}[F2] \right\}\right|_\Gamma = 0
    \f]
    where \f$F1(x,y,z,\mathbf{n})\f$ and \f$F2(x,y,z,\mathbf{n})\f$ are template function evaluators given by user.

    \f[
    F_i =
      \int_{\Omega^e} \frac{\partial N_i}{\partial x_j}\frac{\partial p}{\partial x_j} \textrm{d}V +
      \int_{\Gamma^e} N_i \left\{
	(\textrm{re}[\sigma]+i*\textrm{im}[\sigma] + \textrm{re}[F1]+i\textrm{im}[F1]) p + \textrm{re}[F2]+i\textrm{im}[F2]
      \right\} \textrm{d}V
    \f]

    Note: that first term is volume interval and is taken into account in
    volume operator. Expanding surface integral and taking into account that
    admittance and shape function are complex numbers, we get:
    \f[
    \textrm{re}[F^\Gamma_i] = \int_{\Gamma^e} N_i
      \left\{
      (\textrm{re}[\sigma]+\textrm{re}[F1])\textrm{re}[p] - (\textrm{im}[\sigma]+\textrm{im}[F1])\textrm{im}[p] +  \textrm{re}[F2])
      \right\}
    \textrm{d}V
    \f]
    and
    \f[
    \textrm{im}[F^\Gamma_i] = \int_{\Gamma^e} N_i
      \left\{
      (\textrm{re}[\sigma]+\textrm{re}[F1])\textrm{im}[p] + (\textrm{im}[\sigma]+\textrm{im}[F1])\textrm{re}[p] +  \textrm{im}[F2])
      \right\}
    \textrm{d}V
    \f]


  */
  template<typename FUNEVAL1,typename FUNEVAL2>
  struct OpHelmholtzMixBCRhs: public FaceElementForcesAndSourcesCore::UserDataOperator {

    SurfaceData &dAta;
    CommonData &commonData;
    Vec F;

    const string rePressure;
    const string imPressure;

    boost::shared_ptr<FUNEVAL1> functionEvaluator1;
    boost::shared_ptr<FUNEVAL2> functionEvaluator2;

    OpHelmholtzMixBCRhs(
      const string re_field_name,const string im_field_name,
      Vec _F,SurfaceData &data,CommonData &common_data,
      boost::shared_ptr<FUNEVAL1> function_evaluator1,boost::shared_ptr<FUNEVAL2> function_evaluator2):
      FaceElementForcesAndSourcesCore::UserDataOperator(re_field_name),
      dAta(data),commonData(common_data),F(_F),
      rePressure(re_field_name),imPressure(im_field_name),
      functionEvaluator1(function_evaluator1),functionEvaluator2(function_evaluator2) {}

    ublas::vector<double> reNf,imNf;

    double reResidual,imResidual;
    virtual PetscErrorCode calculateResidualRe(int gg,ublas::vector<double> &f1,ublas::vector<double> &f2) {
      PetscFunctionBegin;

      double re_pressure = commonData.pressureAtGaussPts[rePressure](gg);
      double im_pressure = commonData.pressureAtGaussPts[imPressure](gg);

      reResidual = (dAta.aDmittance_real+f1[0])*re_pressure - (dAta.aDmittance_imag+f1[1])*im_pressure + f2[0];
      imResidual = (dAta.aDmittance_real+f1[0])*im_pressure + (dAta.aDmittance_imag+f1[1])*re_pressure + f2[1];

      PetscFunctionReturn(0);
    }

    ublas::vector<double> nOrmal;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {

        PetscErrorCode ierr;

        int nb_row_dofs = data.getIndices().size();
        if(!nb_row_dofs) {
          PetscFunctionReturn(0);
        }

        if(dAta.tRis.find(getMoFEMFEPtr()->get_ent()) == dAta.tRis.end()) {
          PetscFunctionReturn(0);
        }

        reNf.resize(nb_row_dofs);
        reNf.clear();
        imNf.resize(nb_row_dofs);
        imNf.clear();
        nOrmal.resize(3);

        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

          double area = getArea();
          if(getNormals_at_GaussPt().size1()) {

            noalias(nOrmal) = getNormals_at_GaussPt(gg);
            area = ublas::norm_2(nOrmal)*0.5;
            nOrmal /= 2*area;

          }
          double val = area*getGaussPts()(2,gg);

          double x,y,z;
          if(commonData.hoCoords.size1()) {
            x = commonData.hoCoords(gg,0);
            y = commonData.hoCoords(gg,1);
            z = commonData.hoCoords(gg,2);
          } else {
            x = getCoordsAtGaussPts()(gg,0);
            y = getCoordsAtGaussPts()(gg,1);
            z = getCoordsAtGaussPts()(gg,2);
            if(gg == 0) {
              noalias(nOrmal) = getNormal();
              nOrmal /= 2*area;
            }
          }

          try {
            ublas::vector<double>& f1 = (*functionEvaluator1)(x,y,z,nOrmal);
            ublas::vector<double>& f2 = (*functionEvaluator2)(x,y,z,nOrmal);

            ierr = calculateResidualRe(gg,f1,f2); CHKERRQ(ierr);

            noalias(reNf) += val*(reResidual*data.getN(gg,nb_row_dofs));
            noalias(imNf) += val*(imResidual*data.getN(gg,nb_row_dofs));

          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "throw in method: " << ex.what() << endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
          }

        }

        ierr = VecSetValues(F,
          data.getIndices().size(),&
          data.getIndices()[0],
          &reNf[0],ADD_VALUES
        ); CHKERRQ(ierr);
        ierr = VecSetValues(F,
          (commonData.imIndices[type][side]).size(),
          &(commonData.imIndices[type][side])[0],
          &imNf[0],ADD_VALUES
        ); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Mix boundary conditions
    * \ingroup mofem_helmholtz_elem

  \f[
  C_{ij}= \int_{\Gamma^e} (\sigma+F1) N_i N_j \textrm{d}S
  \f]

  Note that off diagonal imaginary part is assembled on the fly. C matrix is
  transposed and assembled with negative sign.

    */
  template<typename FUNEVAL1>
  struct OpHelmholtzMixBCLhs:public FaceElementForcesAndSourcesCore::UserDataOperator {

    SurfaceData &dAta;
    CommonData &commonData;
    const string rePressure;
    const string imPressure;

    boost::shared_ptr<FUNEVAL1> functionEvaluator1;

    Mat A;
    OpHelmholtzMixBCLhs(
      const string re_field_name,const string im_field_name,
      Mat _A,SurfaceData &data,CommonData &common_data,
      boost::shared_ptr<FUNEVAL1> function_evaluator1):
      FaceElementForcesAndSourcesCore::UserDataOperator(re_field_name,re_field_name), // Note: operator is real-real
      dAta(data),commonData(common_data),rePressure(re_field_name),imPressure(im_field_name),
      functionEvaluator1(function_evaluator1),A(_A) {

      sYmm = false; /// This operator is not symmetric

    }

    ublas::matrix<double> K,K0,K1,reF1K,imF1K;
    ublas::vector<double> nOrmal;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      try {

        int nb_rows = row_data.getIndices().size();
        int nb_cols = col_data.getIndices().size();
        if(nb_rows==0) PetscFunctionReturn(0);
        if(nb_cols==0) PetscFunctionReturn(0);

        if(dAta.tRis.find(getMoFEMFEPtr()->get_ent()) == dAta.tRis.end()) {
          PetscFunctionReturn(0);
        }

        K.resize(nb_rows,nb_cols);
        K.clear();
        reF1K.resize(nb_rows,nb_cols);
        reF1K.clear();
        imF1K.resize(nb_rows,nb_cols);
        imF1K.clear();

        K0.resize(nb_rows,nb_cols);
        nOrmal.resize(3);

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

          double area = getArea();
          if(getNormals_at_GaussPt().size1()) {
            noalias(nOrmal) = getNormals_at_GaussPt(gg);
            area = ublas::norm_2(nOrmal)*0.5;
          }
          double val = area*getGaussPts()(2,gg);

          double x,y,z;
          if(commonData.hoCoords.size1()) {
            x = commonData.hoCoords(gg,0);
            y = commonData.hoCoords(gg,1);
            z = commonData.hoCoords(gg,2);
          } else {
            x = getCoordsAtGaussPts()(gg,0);
            y = getCoordsAtGaussPts()(gg,1);
            z = getCoordsAtGaussPts()(gg,2);
            if(gg == 0) {
              noalias(nOrmal) = getNormal();
            }
          }

          noalias(K0) = outer_prod(row_data.getN(gg,nb_rows),col_data.getN(gg,nb_cols));
          noalias(K) += val*K0;

          ublas::vector<double>& f1 = (*functionEvaluator1)(x,y,z,nOrmal);
          if(f1[0]!=0) {
            noalias(reF1K) += val*f1[0]*K0;
          }
          if(f1[1]!=0) {
            noalias(imF1K) += val*f1[1]*K0;
          }


        }

        if(row_data.getIndices().size()!=(commonData.imIndices[row_type][row_side]).size()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        }

        if(col_data.getIndices().size()!=(commonData.imIndices[col_type][col_side]).size()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        }

        K1.resize(nb_rows,nb_cols);

        PetscErrorCode ierr;

        // real-real
        noalias(K1) = dAta.aDmittance_real*K+reF1K;
        ierr = MatSetValues(
          A,
          nb_rows,&row_data.getIndices()[0],
          nb_cols,&col_data.getIndices()[0],
          &K1(0,0),ADD_VALUES); CHKERRQ(ierr);

        // imag-imag
        noalias(K1) = dAta.aDmittance_real*K+reF1K;
        ierr = MatSetValues(
          A,
          nb_rows,&(commonData.imIndices[row_type][row_side])[0],
          nb_cols,&(commonData.imIndices[col_type][col_side])[0],
          &K1(0,0),ADD_VALUES); CHKERRQ(ierr);

        // real-imag
        noalias(K1) = -(dAta.aDmittance_imag*K+imF1K);
        ierr = MatSetValues(
          A,
          nb_rows,&row_data.getIndices()[0],
          nb_cols,&(commonData.imIndices[col_type][col_side])[0],
          &K1(0,0),ADD_VALUES); CHKERRQ(ierr);

        // imag-real
        noalias(K1) = dAta.aDmittance_imag*K+imF1K;
        ierr = MatSetValues(
          A,
          nb_rows,&(commonData.imIndices[row_type][row_side])[0],
          nb_cols,&col_data.getIndices()[0],
          &K1(0,0),ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /* \brief Infinite Helmholtz Element

  Implementation based on Demkowicz Book:
  Computing with Adaptive HP-Elements Volume 2
  Page: 225

  Element is integrated on outer mesh surface \f$S_a\f$ up to the infinity where
  Sommerfeld radiation condition is applied. Integration is done which-in
  element.

  Moreover approximations functions are scaled in such a way, that despite they
  are different outside discretized volume, on outer mesh surface \f$S_a\f$,
  those functions looks as a shape function which we use to approximate filed in
  the interior.

  The trick is here that this is Bubnov-Galerkin formulation, so we using the
  same functions to approximate pressures on surface \f$S_a\f$ like we use in
  standard FE formulation.

  The Cartesian position in space is parametrised by surface \f$S_a\f$
  parametrisation \f$(\xi_1,\xi_2\f$), which in practice is local approximation
  on triangle on surface and radial parametrisation, as follows
  \f[
  x(\xi_1,\xi_2,\xi_3) = \xi_3^{-1}x_a(\xi_1,\xi_2),\quad a/R < \xi_3 < 1
  \f]

  Functions in the infinite element domain are constructed by tensor product of
  one-divisional radial shape functions of coordinate \f$\xi_3\f$ and standard FE
  hierarchical shape functions.
  \f[
  P(\xi_1,\xi_2,\xi_3) = \sum_{kl} P_{kl} \psi_k(\xi_3) e_l(\xi_1,\xi_2)
  \f]
  where \f$e(\xi_1,\xi_2)\f$ are standard shape functions on element face and
  \f$\psi(\xi_3)\f$ is a radial shape function. The radial shape functions are
  generated as follows
  \f[
  \psi_j =
  \left\{
    \begin{array}{cc}
      1, & j = 0 \\
      P_j(2\xi_3-1), & j > 0
    \end{array}
  \right.
  \f]
  Note that function for \f$j>0\f$ are bubble functions so there are not
  assembled into global coordinate system but are condensed in the finite element
  procedure and then element is assembled. This makes a trick mentioned above,
  that to global system are assembled additional terms but size remains
  unchanged, moreover is symmetric.

  To make all above true, pleas note that
  \f[
  p=\xi_3 e^{-ik(\xi_3^{-1}-1)}P
  \f]
  where for \f$\xi = 1\f$, \f$p(\xi_1,\xi_2) = P(\xi_1,\xi_2,1)\f$

  Infinite part,
  \f[
  b_\infty(q,p) =
    \frac{1}{a} \int_{S_a}
      \left\{
      \int_0^1
	\frac{\partial P}{\partial \xi_3}
	\left(
	\xi_3^2\frac{\partial \overline{Q}}{\partial \xi_3} - i2ka\overline{Q}
	\right)
      \textrm{d}\xi_3
      +(1+ika)P\overline{Q}|_{\xi_3=1}
    \right\}
  \textrm{d}S_a
  +
  a\int_{S_a}
    \int_0^1
      \nabla_{S_a} P\cdot
      \nabla_{S_a} \overline{Q}
    \textrm{d}xi_3
  \textrm{d}S_a
  \f]


  */

  struct InfiniteHelmholtz {

    /*ublas::matrix<double> gaussPointsInRadialDirection;
    ublas::matrix<double> radialShapeFunctions;
    ublas::matrix<double> direvativesOfRadialShapeFunctions;

    ublas::matrix<double> elementShapeFunctions;
    ublas::matrxi<double> direvativeOfElementShapeFunctions;*/

    /** \brief Generate approximation space for infinite element

    It is ordered in such a way, that first indices are loop over surface
    degrees of freedom then radial degrees of freedom. So that
    n_surface_dofs*n_radial_dofs, wher first n_surface_dofs are for surface dofs
    and constant radial shape function. Remaining dofs then could be statically
    condensed.

    */
    //PetscErrorCode generateBase(
      //VectorAdaptor &N,MatrixAdaptor &diffN) {
      //PetscFunctionBegin;

      /*int nb_gauss_points_on_surface = N.size1();
      int nb_gauss_pts_radial = radialShapeFunctions.size1();
      int nb_surface_dofs = N.size2();
      int nb_radial_dofs = radialShapeFunctions.size2();

      elementShapeFunctions.resize(
	nb_gauss_points_on_surface*nb_gauss_pts_radial,
	nb_surface_dofs*nb_radial_dofs);

      for(int ii = 0;ii<nb_radial_dofs;ii++) {
	for(int jj = 0;jj<nb_surface_dofs;jj++) {

	  for(int gg_s = 0;gg_s<nb_gauss_points_on_surface;gg++) {
	    for(int gg_r = 0;gg_r<nb_gauss_pts_radial;gg_r++) {

	    double direvative_shf_in_radial_direction =
	      direvativesOfRadialShapeFunctions(gg_s,ii)*N(gg_s,jj);

	      elementShapeFunctions(
		gg_r*nb_gauss_points_on_surface+gg_s,
		ii*nb_surface_dofs+jj) =
		radialShapeFunctions(gg_s,ii)*N(gg_r,jj);

	      direvativeOfElementShapeFunctions(
		gg_r*nb_gauss_points_on_surface+gg_s,
		3*ii*nb_surface_dofs+3*jj+0) =
		radialShapeFunctions(gg_s,ii)*diffN(gg_r,3*jj+0);
	      direvativeOfElementShapeFunctions(
		gg_r*nb_gauss_points_on_surface+gg_s,
		3*ii*nb_surface_dofs+3*jj+1) =
		radialShapeFunctions(gg_s,ii)*diffN(gg_r,3*jj+1);
	      direvativeOfElementShapeFunctions(
		gg_r*nb_gauss_points_on_surface+gg_s,
		3*ii*nb_surface_dofs+3*jj+3) =
		direvativesOfRadialShapeFunctions(gg_s,ii)*N(gg_r,jj);

	    }
	  }

	}
      }*/

      //PetscFunctionReturn(0);
    //}

  };

  /** \brief Add Helmholtz elements to problem
    * \ingroup mofem_helmholtz_elem
    *
    * It get data from block set and define element in moab
    *w
    * \param problem name
    * \param field name
    * \param name of mesh nodal positions (if not defined nodal coordinates are used)
    */
  PetscErrorCode addHelmholtzElements(
    const string re_field_name,const string im_field_name,
    const string mesh_nodals_positions = "MESH_NODE_POSITIONS",
    const string pressure_field = "P") {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    if(mField.check_field(pressure_field)) {

      ierr = mField.add_finite_element("PRESSURE_FE",MF_ZERO); CHKERRQ(ierr );
      ierr = mField.modify_finite_element_add_field_row("PRESSURE_FE",pressure_field); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("PRESSURE_FE",pressure_field); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("PRESSURE_FE",pressure_field); CHKERRQ(ierr);

    }

    ierr = mField.add_finite_element("HELMHOLTZ_RERE_FE",MF_ZERO); CHKERRQ(ierr );
    ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_RERE_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_RERE_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE",im_field_name); CHKERRQ(ierr);

    ierr = mField.add_finite_element("HELMHOLTZ_IMIM_FE",MF_ZERO); CHKERRQ(ierr );
    ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_IMIM_FE",im_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_IMIM_FE",im_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_IMIM_FE",im_field_name); CHKERRQ(ierr);

    ierr = mField.add_finite_element("HELMHOLTZ_REIM_FE",MF_ZERO); CHKERRQ(ierr );
    ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_REIM_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_REIM_FE",im_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_REIM_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_REIM_FE",im_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_REIM_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_REIM_FE",im_field_name); CHKERRQ(ierr);

    if(mField.check_field(mesh_nodals_positions)) {

      ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE",mesh_nodals_positions); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_IMIM_FE",mesh_nodals_positions); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_REIM_FE",mesh_nodals_positions); CHKERRQ(ierr);

      if(mField.check_field(pressure_field)) {

        ierr = mField.modify_finite_element_add_field_data("PRESSURE_FE",mesh_nodals_positions); CHKERRQ(ierr);

      }

    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {

      if(it->get_name().compare(0,13,"MAT_HELMHOLTZ") == 0) {

        volumeData[it->get_msId()].waveNumber = globalParameters.waveNumber.first;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,volumeData[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TETs(volumeData[it->get_msId()].tEts,"HELMHOLTZ_RERE_FE"); CHKERRQ(ierr);
        ierr = mField.add_ents_to_finite_element_by_TETs(volumeData[it->get_msId()].tEts,"HELMHOLTZ_IMIM_FE"); CHKERRQ(ierr);

        if(mField.check_field(pressure_field)) {

          ierr = mField.add_ents_to_finite_element_by_TETs(volumeData[it->get_msId()].tEts,"PRESSURE_FE"); CHKERRQ(ierr);

        }

      }

    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {

      if(it->get_name().compare(0,23,"HARD_INCIDENT_WAVE_BC") == 0) {

        surfaceIncidentWaveBcData[it->get_msId()].aDmittance_real = 0;
        surfaceIncidentWaveBcData[it->get_msId()].aDmittance_imag = 0;

        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,surfaceIncidentWaveBcData[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TRIs(surfaceIncidentWaveBcData[it->get_msId()].tRis,"HELMHOLTZ_REIM_FE"); CHKERRQ(ierr);

      }

      if(it->get_name().compare(0,22,"MIX_INCIDENT_WAVE_BC") == 0) {

        //get block attributes
        vector<double> attributes;
        ierr = it->get_attributes(attributes); CHKERRQ(ierr);
        if(attributes.size()<1) {
          SETERRQ1(PETSC_COMM_SELF,1,"first block attribute should define surface admitance",attributes.size());
        }

        surfaceIncidentWaveBcData[it->get_msId()].aDmittance_real = 0;
        surfaceIncidentWaveBcData[it->get_msId()].aDmittance_imag = attributes[0];
        if(globalParameters.surfaceAdmittance.second) {

          surfaceIncidentWaveBcData[it->get_msId()].aDmittance_imag = globalParameters.surfaceAdmittance.first;

        }

        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,surfaceIncidentWaveBcData[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TRIs(surfaceIncidentWaveBcData[it->get_msId()].tRis,"HELMHOLTZ_REIM_FE"); CHKERRQ(ierr);

      }

      if(it->get_name().compare(0,13,"SOMMERFELD_BC") == 0) {

        sommerfeldBcData[it->get_msId()].aDmittance_real = 0;
        sommerfeldBcData[it->get_msId()].aDmittance_imag = -globalParameters.waveNumber.first;

        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,sommerfeldBcData[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TRIs(sommerfeldBcData[it->get_msId()].tRis,"HELMHOLTZ_REIM_FE"); CHKERRQ(ierr);

      }

      if(it->get_name().compare(0,17,"BAYLISS_TURKEL_BC") == 0) {

        baylissTurkelBcData[it->get_msId()].aDmittance_real = 0;
        baylissTurkelBcData[it->get_msId()].aDmittance_imag = -globalParameters.waveNumber.first;


        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,baylissTurkelBcData[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TRIs(baylissTurkelBcData[it->get_msId()].tRis,"HELMHOLTZ_REIM_FE"); CHKERRQ(ierr);

      }

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setOperators(
    Mat A,Vec F,
    const string re_field_name,const string im_field_name,
    const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    string fe_name;

    fe_name = "HELMHOLTZ_RERE_FE"; feLhs.insert(fe_name,new MyVolumeFE(mField,addToRank));
    fe_name = "HELMHOLTZ_RERE_FE"; feRhs.insert(fe_name,new MyVolumeFE(mField,addToRank));
    fe_name = "HELMHOLTZ_IMIM_FE"; feRhs.insert(fe_name,new MyVolumeFE(mField,addToRank));

    feLhs.at("HELMHOLTZ_RERE_FE").getRowOpPtrVector().push_back(
      new OpGetImIndices(re_field_name,im_field_name,commonData));

    /* real field and imag field */
    feRhs.at("HELMHOLTZ_RERE_FE").getRowOpPtrVector().push_back(
      new OpGetValueAndGradAtGaussPts(re_field_name,commonData));
    feRhs.at("HELMHOLTZ_IMIM_FE").getRowOpPtrVector().push_back(
     new OpGetValueAndGradAtGaussPts(im_field_name,commonData));

    map<int,VolumeData>::iterator mit = volumeData.begin();
    for(;mit!=volumeData.end();mit++) {

      feLhs.at("HELMHOLTZ_RERE_FE").getRowColOpPtrVector().push_back(
	       new OpHelmholtzLhs(re_field_name,im_field_name,A,mit->second,commonData));

      feRhs.at("HELMHOLTZ_RERE_FE").getRowOpPtrVector().push_back(
	       new OpHelmholtzRhs(re_field_name,F,mit->second,commonData));
      feRhs.at("HELMHOLTZ_IMIM_FE").getRowOpPtrVector().push_back(
	       new OpHelmholtzRhs(im_field_name,F,mit->second,commonData));

    }

    fe_name = "HELMHOLTZ_REIM_FE";
    feLhs.insert(fe_name,new MySurfaceFE(mField,addToRank));
    feRhs.insert(fe_name,new MySurfaceFE(mField,addToRank));

    if(mField.check_field(mesh_nodals_positions)) {

      feLhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(
        new OpHoCoordTri(mesh_nodals_positions,commonData.hoCoords)
      );

      feRhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(
        new OpHoCoordTri(mesh_nodals_positions,commonData.hoCoords)
      );

    }

    // Get im indices
    feLhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(new OpGetImIndices(im_field_name,im_field_name,commonData));
    feRhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(new OpGetImIndices(im_field_name,im_field_name,commonData));

    // Get field values
    feRhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(new OpGetValueAtGaussPts(re_field_name,commonData));
    feRhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(new OpGetValueAtGaussPts(im_field_name,commonData));

    boost::shared_ptr<ZeroFunVal> zero_function = boost::shared_ptr<ZeroFunVal>(new ZeroFunVal());

    map<int,SurfaceData>::iterator miit = surfaceIncidentWaveBcData.begin();
    for(;miit!=surfaceIncidentWaveBcData.end();miit++) {

      boost::shared_ptr<IncidentWaveNeumannF2> incident_wave_neumann_bc =
        boost::shared_ptr<IncidentWaveNeumannF2>(new IncidentWaveNeumannF2(globalParameters));

      if(miit->second.aDmittance_imag!=0) {

        feRhs.at("HELMHOLTZ_REIM_FE").getRowColOpPtrVector().push_back(
          new OpHelmholtzMixBCLhs<ZeroFunVal>(re_field_name,im_field_name,A,miit->second,commonData, zero_function)
        );

      }

      // assembled to the right hand vector
      feRhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(
        new OpHelmholtzMixBCRhs<ZeroFunVal,IncidentWaveNeumannF2>(
          re_field_name,im_field_name,F,miit->second,commonData,
          zero_function,incident_wave_neumann_bc
        )
      );

    }

    miit = sommerfeldBcData.begin();
    for(;miit!=sommerfeldBcData.end();miit++) {

      feLhs.at("HELMHOLTZ_REIM_FE").getRowColOpPtrVector().push_back(
        new OpHelmholtzMixBCLhs<ZeroFunVal>(
          re_field_name,im_field_name,A,miit->second,commonData,
          zero_function
        )
      );

      // FIXME: need to add second functions so that residual is calculated properly
      feRhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(
        new OpHelmholtzMixBCRhs<ZeroFunVal,ZeroFunVal>(
          re_field_name,im_field_name,F,miit->second,commonData,
          zero_function,zero_function
        )
      );

    }

    boost::shared_ptr<BaylissTurkel> bayliss_turkel_bc =
    boost::shared_ptr<BaylissTurkel>(new BaylissTurkel());

    miit = baylissTurkelBcData.begin();
    for(;miit!=baylissTurkelBcData.end();miit++) {

      feLhs.at("HELMHOLTZ_REIM_FE").getRowColOpPtrVector().push_back(
        new OpHelmholtzMixBCLhs<BaylissTurkel>(
          re_field_name,im_field_name,A,miit->second,commonData,
          bayliss_turkel_bc
        )
      );

      // FIXME: need to add second functions so that residual is calculated properly
      feRhs.at("HELMHOLTZ_REIM_FE").getRowOpPtrVector().push_back(
        new OpHelmholtzMixBCRhs<BaylissTurkel,ZeroFunVal>(
          re_field_name,im_field_name,F,miit->second,commonData,
          bayliss_turkel_bc,zero_function
        )
      );

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode calculateA(const string problem_name) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    boost::ptr_map<string,ForcesAndSurcesCore>::iterator mit = feLhs.begin();
    for(;mit!=feLhs.end();mit++) {

      ierr = mField.loop_finite_elements(problem_name,mit->first,*(mit->second)); CHKERRQ(ierr); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode calculateF(const string problem_name) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    boost::ptr_map<string,ForcesAndSurcesCore>::iterator mit = feRhs.begin();
    for(;mit!=feRhs.end();mit++) {

      ierr = mField.loop_finite_elements(problem_name,mit->first,*(mit->second)); CHKERRQ(ierr); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

};

#endif //__HELMHOLTZ_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup mofem_helmholtz_elem Helmholtz element
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
