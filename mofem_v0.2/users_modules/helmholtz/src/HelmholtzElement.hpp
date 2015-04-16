/** \file HelmholtzElement.hpp
\ingroup mofem_helmholtz_elem

 \brief Operators and data structures for wave propagation analyse (Galerkin Element)

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
  struct MySurfaceFE: public TriElementForcesAndSurcesCore {
    int addToRank; ///< default value 1, i.e. assumes that geometry is approx. by quadratic functions.
    MySurfaceFE(FieldInterface &_mField,int add_to_rank): TriElementForcesAndSurcesCore(_mField),addToRank(add_to_rank) {}
    int getRule(int order) { return order+addToRank; };
  };

  boost::ptr_map<string,ForcesAndSurcesCore> feRhs; // surface element for LHS 
  boost::ptr_map<string,ForcesAndSurcesCore> feLhs; // surface element for RHS

  /** \brief Volume element data
  * \ingroup mofem_helmholtz_elem
  */
  struct VolumeData {
    double wAvenumber;
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
  map<int,SurfaceData> surfaceData;

  /** \brief Common data used by volume and surface elements
  * \ingroup mofem_helmholtz_elem
  */
  struct CommonData {
    
    map<string,ublas::vector<double> > pressureAtGaussPts; 
    map<string,ublas::matrix<double> > gradPressureAtGaussPts;
    ublas::matrix<double> hoCoords;

  };
  CommonData commonData;
  
  FieldInterface &mField;
  int addToRank; ///< default value 1, i.e. assumes that geometry is approx. by quadratic functions.

  HelmholtzElement(
    FieldInterface &m_field):
    mField(m_field),addToRank(1) {}
  
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
  
        // Initialize
	value.resize(nb_gauss_pts);
        gradient.resize(nb_gauss_pts,3);
        if(type == MBVERTEX) {
	  gradient.clear();
	  value.clear();
        }
  
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          value[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
          ublas::noalias(ublas::matrix_row<ublas::matrix<double> >(gradient,gg))
	    += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
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
  struct OpGetValueAtGaussPts: public TriElementForcesAndSurcesCore::UserDataOperator {
  
    CommonData &commonData;
    const string fieldName;
    OpGetValueAtGaussPts(const string field_name,
      CommonData &common_data):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
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
  struct OpHoCoordTri: public TriElementForcesAndSurcesCore::UserDataOperator {
  
    ublas::matrix<double> &hoCoordsTri;
    OpHoCoordTri(const string field_name,ublas::matrix<double> &ho_coords): 
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
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
        double k_pow2 = dAta.wAvenumber*dAta.wAvenumber;

        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
  
          double val = getVolume()*getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; // higher order geometry
          }

					const ublas::matrix_row<ublas::matrix<double> > gard_p_at_gauss_pt(grad_p,gg);

					/// Integrate diffN^T grad_p - k^2 N^T p dV
          ublas::noalias(Nf) += val*prod(data.getDiffN(gg,nb_row_dofs),gard_p_at_gauss_pt);
					ublas::noalias(Nf) -= val*k_pow2*data.getN(gg)*pressure[gg];

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
  
  
  /** \brief Mix boundary conditions on surface element

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
    ublas::vector<int> imRowIndices,imColIndices;
  
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

        double k_pow2 = dAta.wAvenumber*dAta.wAvenumber;

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
          

          double val = getVolume()*getGaussPts()(3,gg);
          if(getHoGaussPtsDetJac().size()>0) {
            val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
          }

          noalias(K) += val*prod(row_data.getDiffN(gg,nb_rows),trans(col_data.getDiffN(gg,nb_cols)));
	  noalias(K) -= val*k_pow2*outer_prod( row_data.getN(gg,nb_rows),col_data.getN(gg,nb_cols) );

        }

        PetscErrorCode ierr;
        ierr = MatSetValues(
          A,  
          nb_rows,&row_data.getIndices()[0],
          nb_cols,&col_data.getIndices()[0],
	  &K(0,0),ADD_VALUES); CHKERRQ(ierr);
        
        if(row_side != col_side || row_type != col_type) {
          transK.resize(nb_cols,nb_rows);
          noalias(transK) = trans(K);
          ierr = MatSetValues(
            A,
            nb_cols,&col_data.getIndices()[0],
            nb_rows,&row_data.getIndices()[0],
            &transK(0,0),ADD_VALUES); CHKERRQ(ierr);
        }  

	// Get rows and cols indices of imaginary part and assemble matrix.
	// Note: However HELMHOLTZ_IMIM_FE element is not calculated, since
	// matrix A on real elements is equal to matrix in imaginary elements,
	// it need be to declared. Declaration indicate that on imaginary part,
	// assembled matrix has non-zero values;
	ierr = getPorblemRowIndices(imFieldName,row_type,row_side,imRowIndices); CHKERRQ(ierr);
	ierr = getPorblemColIndices(imFieldName,col_type,col_side,imColIndices); CHKERRQ(ierr);

        ierr = MatSetValues(
          A,  
          nb_rows,&imRowIndices[0],
          nb_cols,&imColIndices[0],
	  &K(0,0),ADD_VALUES); CHKERRQ(ierr);
        
        if(row_side != col_side || row_type != col_type) {
          ierr = MatSetValues(
            A,
            nb_cols,&imColIndices[0],
            nb_rows,&imRowIndices[0],
            &transK(0,0),ADD_VALUES); CHKERRQ(ierr);
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

  ZeroFunVal zero_fun_val;
	
	struct BT2FunVal {
    
    ublas::vector<double> vAl;    
	
    ublas::vector<double>& operator()(double x,double y,double z,ublas::vector<double> &normal) {
      vAl.resize(2);
			const complex< double > i( 0.0, 1.0 );
      vAl[0] = 0;
      vAl[1] = 0;
      return vAl;
    }
    
  };
	
	BT2FunVal bt2_fun_val;

  /** \brief Rhs vector for Helmholtz operator
    \ingroup mofem_helmholtz_elem

    \f[
    \left.\left( \frac{\partial p}{\partial n} + i\sigma p \right)\right|_\Gamma = 0
    \f]

    \f[
    F_i =  \int_{\Omega^e} \frac{\partial N_i}{\partial x_j}\frac{\partial p}{\partial x_j} \textrm{d}V + \int_{\Gamma^e} N_i i\sigma p \textrm{d}V
    \f]

    Note: that first term is volume interval and is taken into account in
    volume operator. Expanding surface integral and taking into account that
    admittance and shape function are complex numbers, we get:
    \f[
    \textrm{re}[F^\Gamma_i] = \int_{\Gamma^e} N_i (\textrm{re}[\sigma]\textrm{re}[p] - \textrm{im}[\sigma]\textrm{im}[p]) \textrm{d}V
    \f]
    and
    \f[
    \textrm{im}[F^\Gamma_i] = \int_{\Gamma^e} N_i (\textrm{im}[\sigma]\textrm{re}[p] + \textrm{re}[\sigma]\textrm{im}[p])   \textrm{d}V
    \f]

  */
  template<typename FUNVAL>
  struct OpHelmholtzMixBCRhs: public TriElementForcesAndSurcesCore::UserDataOperator {
  
    SurfaceData &dAta;
    CommonData &commonData;
    Vec F;

    const string rePressure;
    const string imPressure;

    FUNVAL &functionEvaluator;

    OpHelmholtzMixBCRhs(
      const string re_field_name,const string im_field_name,
      Vec _F,SurfaceData &data,CommonData &common_data,
      FUNVAL &function_evaluator):
      TriElementForcesAndSurcesCore::UserDataOperator(re_field_name),
      dAta(data),commonData(common_data),F(_F),
      rePressure(re_field_name),imPressure(im_field_name),
      functionEvaluator(function_evaluator) { }

    ublas::vector<double> reNf,imNf;

    double reResidual,imResidual;
    virtual PetscErrorCode calculateResidualRe(int gg) {
      PetscFunctionBegin;

      double re_pressure = commonData.pressureAtGaussPts[rePressure](gg);
      double im_pressure = commonData.pressureAtGaussPts[imPressure](gg);
    
      reResidual = dAta.aDmittance_real*re_pressure - dAta.aDmittance_imag*im_pressure;
      imResidual = dAta.aDmittance_real*im_pressure + dAta.aDmittance_imag*re_pressure;

      PetscFunctionReturn(0);
    }

    ublas::vector<int> imRowIndices; 
    ublas::vector<double> nOrmal;
  
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
  
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent()) == dAta.tRis.end()) {
        PetscFunctionReturn(0);
      }
  
      try {
  
        PetscErrorCode ierr;
       
        int nb_row_dofs = data.getIndices().size();
	if(!nb_row_dofs) {
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
	    area = ublas::norm_2(getNormals_at_GaussPt(gg))*0.5;
	  }
          double val = area*getGaussPts()(2,gg);

	  double x,y,z;
	  if(commonData.hoCoords.size1() == data.getN().size1()) {
	    x = commonData.hoCoords(gg,0);
	    y = commonData.hoCoords(gg,1);
	    z = commonData.hoCoords(gg,2);	
	    if(gg == 0) {
	      noalias(nOrmal) = getNormal();
	    }
	  } else {
	    x = getCoordsAtGaussPts()(gg,0);
	    y = getCoordsAtGaussPts()(gg,1);
	    z = getCoordsAtGaussPts()(gg,2);
	    noalias(nOrmal) = getNormals_at_GaussPt(gg);
	    
	  }
	  ublas::vector<double>& f = functionEvaluator(x,y,z,nOrmal);
    
	  ierr = calculateResidualRe(gg); CHKERRQ(ierr);
  
	  noalias(reNf) += val*((reResidual-f[0])*data.getN(gg));
	  noalias(imNf) += val*((imResidual-f[0])*data.getN(gg));

        }
  
        ierr = VecSetValues(F,data.getIndices().size(),
          &data.getIndices()[0],&reNf[0],ADD_VALUES); CHKERRQ(ierr);

	// get indices of imaginary degrees of freedoms and assemble imaginary part
	ierr = getPorblemRowIndices(imPressure,type,side,imRowIndices); CHKERRQ(ierr);
        ierr = VecSetValues(F,imRowIndices.size(),&imRowIndices[0],&imNf[0],ADD_VALUES); CHKERRQ(ierr);
  
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
  C_{ij}= \int_{\Gamma^e} N_i N_j \textrm{d}S
  \f]
  
  Note that off diagonal imaginary part is assembled on the fly. C matrix is
  transposed and assembled with negative sign.

    */
	template<typename FUNVAL>
  struct OpHelmholtzMixBCLhs:public TriElementForcesAndSurcesCore::UserDataOperator {

    SurfaceData &dAta;
    CommonData &commonData;
    const string rePressure;
    const string imPressure; 
		FUNVAL &functionEvaluator;
			
    Mat A;
    OpHelmholtzMixBCLhs(
      const string re_field_name,const string im_field_name,
      Mat _A,SurfaceData &data,CommonData &common_data,
			FUNVAL &function_evaluator):
      TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
      dAta(data),commonData(common_data),rePressure(re_field_name),imPressure(im_field_name),A(_A),
			functionEvaluator(function_evaluator) {
      symm = false;
		}
  
    ublas::matrix<double> K,K1;
    ublas::vector<int> imRowIndices,imColIndices;
    ublas::vector<double> nOrmal;
    
		PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;
  
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent()) == dAta.tRis.end()) {
        PetscFunctionReturn(0);
      }
  
      try {
  
				int nb_rows = row_data.getIndices().size();
				int nb_cols = col_data.getIndices().size();
        if(nb_rows==0) PetscFunctionReturn(0);
        if(nb_cols==0) PetscFunctionReturn(0);
        
        K.resize(nb_rows,nb_cols);
				K.clear();

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
            
					double area = getArea();
					if(getNormals_at_GaussPt().size1()) {
					area = ublas::norm_2(getNormals_at_GaussPt(gg))*0.5;
					}
					double val = area*getGaussPts()(2,gg);

					/*get cartesian coordinates */
					double x,y,z;
					if(commonData.hoCoords.size1() == row_data.getN().size1()) {
						x = commonData.hoCoords(gg,0);
						y = commonData.hoCoords(gg,1);
						z = commonData.hoCoords(gg,2);	
						if(gg == 0) {
							noalias(nOrmal) = getNormal();
						}
					} else {
						x = getCoordsAtGaussPts()(gg,0);
						y = getCoordsAtGaussPts()(gg,1);
						z = getCoordsAtGaussPts()(gg,2);
						noalias(nOrmal) = getNormals_at_GaussPt(gg);
						
					}
					ublas::vector<double>& f = functionEvaluator(x,y,z,nOrmal);
					
					noalias(K) += val*f[0]*outer_prod(row_data.getN(gg,nb_rows),col_data.getN(gg,nb_cols));

        }

        PetscErrorCode ierr;

				ierr = getPorblemRowIndices(imPressure,row_type,row_side,imRowIndices); CHKERRQ(ierr);
				ierr = getPorblemRowIndices(imPressure,col_type,col_side,imColIndices); CHKERRQ(ierr);

				K1.resize(nb_rows,nb_rows);
	
				// real-real
				if(dAta.aDmittance_real!=0) {
				noalias(K1) = dAta.aDmittance_real*K;
				ierr = MatSetValues(
					A,  
					nb_rows,&row_data.getIndices()[0],
					nb_cols,&col_data.getIndices()[0],
					&K1(0,0),ADD_VALUES); CHKERRQ(ierr);
				// imag-imag
				noalias(K1) = dAta.aDmittance_real*K;
				ierr = MatSetValues(
					A,  
					nb_rows,&imRowIndices[0],
					nb_cols,&imRowIndices[0],
					&K1(0,0),ADD_VALUES); CHKERRQ(ierr);   
				}

				// real-imag
				if(dAta.aDmittance_imag!=0) {
					noalias(K1) =  -dAta.aDmittance_imag*K;
					ierr = MatSetValues(
						A,  
						nb_rows,&row_data.getIndices()[0],
						nb_cols,&imColIndices[0],
						&K1(0,0),ADD_VALUES); CHKERRQ(ierr);
				// imag-real
				noalias(K1) = dAta.aDmittance_imag*K;
				ierr = MatSetValues(
					A,  
					nb_rows,&imRowIndices[0],
					nb_cols,&col_data.getIndices()[0],
					&K1(0,0),ADD_VALUES); CHKERRQ(ierr);   
				}
  
			} catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
  
      PetscFunctionReturn(0);
    }

      
   
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
    const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
  
    PetscErrorCode ierr;
    ErrorCode rval;
    
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
    ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_REIM_FE",im_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_REIM_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_REIM_FE",im_field_name); CHKERRQ(ierr);

    ierr = mField.add_finite_element("HELMHOLTZ_IMRE_FE",MF_ZERO); CHKERRQ(ierr ); 
    ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_IMRE_FE",re_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_IMRE_FE",im_field_name); CHKERRQ(ierr);
    
    if(mField.check_field(mesh_nodals_positions)) {

      ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_RERE_FE",mesh_nodals_positions); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_IMIM_FE",mesh_nodals_positions); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_REIM_FE",mesh_nodals_positions); CHKERRQ(ierr);

    }
   
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
        
      if(it->get_Cubit_name().compare(0,13,"MAT_HELMHOLTZ") == 0) {
        
        //get block attributes
        vector<double> attributes;
        ierr = it->get_Cubit_attributes(attributes); CHKERRQ(ierr);
        if(attributes.size()<2) {
          SETERRQ1(PETSC_COMM_SELF,1,"not enough block attributes to deffine fluid pressure element, attributes.size() = %d ",attributes.size());
        }
        
        double aNgularfreq = attributes[0];
        double sPeed = attributes[1];
        volumeData[it->get_msId()].wAvenumber = aNgularfreq/sPeed;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,volumeData[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TETs(volumeData[it->get_msId()].tEts,"HELMHOLTZ_RERE_FE"); CHKERRQ(ierr);
        ierr = mField.add_ents_to_finite_element_by_TETs(volumeData[it->get_msId()].tEts,"HELMHOLTZ_IMIM_FE"); CHKERRQ(ierr);

				ierr = PetscOptionsGetScalar(NULL,"-wave_number",&volumeData[it->get_msId()].wAvenumber,NULL); CHKERRQ(ierr);

      }
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->get_Cubit_name().compare(0,9,"IMPEDANCE") == 0) {
        
        vector<double> data;
        ierr = it->get_Cubit_attributes(data); CHKERRQ(ierr);
        if(data.size()!=2) {
          SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
        }
        surfaceData[it->get_msId()].aDmittance_real = data[0];
        surfaceData[it->get_msId()].aDmittance_imag = data[1];

        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,surfaceData[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TRIs(surfaceData[it->get_msId()].tRis,"HELMHOLTZ_REIM_FE"); CHKERRQ(ierr);
        ierr = mField.add_ents_to_finite_element_by_TRIs(surfaceData[it->get_msId()].tRis,"HELMHOLTZ_IMRE_FE"); CHKERRQ(ierr);

	ierr = PetscOptionsGetScalar(NULL,"-sigma_real",&surfaceData[it->get_msId()].aDmittance_real,NULL); CHKERRQ(ierr);
	ierr = PetscOptionsGetScalar(NULL,"-sigma_imag",&surfaceData[it->get_msId()].aDmittance_imag,NULL); CHKERRQ(ierr);

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

		/* volume elements */
    fe_name = "HELMHOLTZ_RERE_FE"; feLhs.insert(fe_name,new MyVolumeFE(mField,addToRank));
    fe_name = "HELMHOLTZ_RERE_FE"; feRhs.insert(fe_name,new MyVolumeFE(mField,addToRank));
    fe_name = "HELMHOLTZ_IMIM_FE"; feRhs.insert(fe_name,new MyVolumeFE(mField,addToRank));

    feRhs.at("HELMHOLTZ_RERE_FE").get_op_to_do_Rhs().push_back(
     new OpGetValueAndGradAtGaussPts(re_field_name,commonData));
    feRhs.at("HELMHOLTZ_IMIM_FE").get_op_to_do_Rhs().push_back(
     new OpGetValueAndGradAtGaussPts(im_field_name,commonData));

    map<int,VolumeData>::iterator mit = volumeData.begin();
    for(;mit!=volumeData.end();mit++) {

      feLhs.at("HELMHOLTZ_RERE_FE").get_op_to_do_Lhs().push_back(
			 new OpHelmholtzLhs(re_field_name,im_field_name,A,mit->second,commonData));

      feRhs.at("HELMHOLTZ_RERE_FE").get_op_to_do_Rhs().push_back(
			 new OpHelmholtzRhs(re_field_name,F,mit->second,commonData));
			
      feRhs.at("HELMHOLTZ_IMIM_FE").get_op_to_do_Rhs().push_back(
			 new OpHelmholtzRhs(im_field_name,F,mit->second,commonData));

    }

		/* surface elements */
    fe_name = "HELMHOLTZ_REIM_FE"; feLhs.insert(fe_name,new MySurfaceFE(mField,addToRank));
    fe_name = "HELMHOLTZ_REIM_FE"; feRhs.insert(fe_name,new MySurfaceFE(mField,addToRank));

    if(mField.check_field(mesh_nodals_positions)) {

      feRhs.at("HELMHOLTZ_REIM_FE").get_op_to_do_Rhs().push_back(
				new OpHoCoordTri(mesh_nodals_positions,commonData.hoCoords));

    }

    feRhs.at("HELMHOLTZ_REIM_FE").get_op_to_do_Rhs().push_back(
      new OpGetValueAtGaussPts(re_field_name,commonData));
    feRhs.at("HELMHOLTZ_REIM_FE").get_op_to_do_Rhs().push_back(
      new OpGetValueAtGaussPts(im_field_name,commonData));

    map<int,SurfaceData>::iterator miit = surfaceData.begin();
    for(;miit!=surfaceData.end();miit++) {
      feLhs.at("HELMHOLTZ_REIM_FE").get_op_to_do_Lhs().push_back(
        new OpHelmholtzMixBCLhs<ZeroFunVal>(re_field_name,im_field_name,A,
	miit->second,commonData,zero_fun_val));
      feRhs.at("HELMHOLTZ_REIM_FE").get_op_to_do_Rhs().push_back(
        new OpHelmholtzMixBCRhs<ZeroFunVal>(re_field_name,im_field_name,F,
	  miit->second,commonData,zero_fun_val));
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
