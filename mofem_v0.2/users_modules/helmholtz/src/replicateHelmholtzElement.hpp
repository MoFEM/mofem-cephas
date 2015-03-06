/** \file HelmholtzElement.hpp
 * \brief Operators and data structures for helmholtz analysis
 *
 * Implementation of Helmholtz element for unsteady and steady case.
 *
 */

/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * PhD student: Thomas Felix Xuan Meng
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. 
 * The header file should contains as least #include as possible for speed */

#ifndef __HELMHOLTZ_ELEMENT_HPP
#define __HELMHOLTZ_ELEMENT_HPP

#include<moab/Skinner.hpp>
using namespace std;
using namespace boost::math;



namespace MoFEM {

/** \brief struture grouping operators and data used for helmholtz problems
* \ingroup mofem_helmholtz_elem
*
* In order to assemble matrices and right hand side vectors, the loops over
* elements, enetities over that elememnts and finally loop over intergration
* points are executed.
*
* complex-real transformation are implemented for computational efficenicy
* See Ercegovac, Milos, and Jean-Michel Muller. "Solving Systems of Linear 
* Equations in Complex Domain: Complex E-Method." (2007). for details. 
* 
*
* Following implementation separte those three cegories of loops and to each
* loop attach operator.
*
*/
// Three dimensional homogeneous isotropic medium and time harmonic Helmholtz wave operator
struct HelmholtzElement {
	/// \brief  definition of volume element
	struct MyVolumeFE: public TetElementForcesAndSourcesCore {
		MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSourcesCore(_mField) {}
	
		/** \brief it is used to calculate nb. of Gauss integartion points
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
		int getRule(int order) { return order-1; };
	};
	MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
	MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element
	MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
	MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element
	
	/** \brief define surface element
	*
	* 
	*/
	struct MyTriFE: public TriElementForcesAndSurcesCore {
		MyTriFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
		int getRule(int order) { return ceil(order/2); };
	};
	MyTriFE feFlux; //< scalar velocity element
    MyTriFE& getLoopFeFlux() { return feFlux; } //< get scalar velocity element
	
	MyTriFE feIncidentWave; //< Incident wave flux element
    MyTriFE& getLoopfeIncidentWave() { return feIncidentWave; } //< get Incident Waveflux element
	
	//MyTriFE feImpedanceRhs; //< Robin Boundary element
    //MyTriFE& getLoopFeImpedanceRhs() { return feImpedanceRhs; } //< get Robin element
    MyTriFE feImpedanceLhs;
    MyTriFE& getLoopFeImpedanceLhs() { return feImpedanceLhs; }

	
	
	FieldInterface &mField;
    HelmholtzElement(
        FieldInterface &m_field):
        feRhs(m_field),feLhs(m_field),feFlux(m_field),feIncidentWave(m_field),feImpedanceLhs(m_field),mField(m_field) {}
	
	/** \brief data for calulation Angular Frequency and wave speed elements
	* \infroup mofem_helmholtz_elem
	*/
	struct BlockData {
		double aNgularfreq; // Angular frequency
		double sPeed;   // Wave Speed
		
		Range tEts; ///< constatins elements in block set
	};
	map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData
	BlockData blockData;

	/** \brief data for calulation heat flux
	* \infroup mofem_helmholtz_elem
	*/
	struct FluxData {
		HeatfluxCubitBcData dAta; ///< for more details look to BCMultiIndices.hpp to see details of HeatfluxCubitBcData
		Range tRis; ///< suraface triangles where hate flux is applied
	};
	map<int,FluxData> setOfFluxes; ///< maps side set id with appropiate FluxData
	
	/** \brief data for Robin Boundary condition
	* \infroup mofem_helmholtz_elem
	*/
	struct ImpedanceData {
		double g; /*the zero velocity term of the Robin BC*/
		double sIgma;  //admittance of the surface, if g and sIgma both zero the surface is called sound hard.
		Range tRis; ///< those will be on body skin, except thos with contact whith other body where pressure is applied
	};
	map<int,ImpedanceData> setOfImpedance; //< maps block set id with appropiate data
	
	
	
	/** \brief common data used by volume elements
	* \infroup mofem_helmholtz_elem
	*/
	struct CommonData {
		ublas::vector<double> pressureReAtGaussPts;   //Helmholtz Pressure
		ublas::vector<double> pressureImAtGaussPts;   //Helmholtz Pressure
		ublas::vector<double> pressureRateAtGaussPts;
		ublas::matrix<double> gradReAtGaussPts;
		ublas::matrix<double> gradImAtGaussPts;
		inline ublas::matrix_row<ublas::matrix<double> > getGradReAtGaussPts(const int gg) {
			return ublas::matrix_row<ublas::matrix<double> >(gradReAtGaussPts,gg);
		}
		inline ublas::matrix_row<ublas::matrix<double> > getGradImAtGaussPts(const int gg) {
			return ublas::matrix_row<ublas::matrix<double> >(gradImAtGaussPts,gg);
		}
	};
	CommonData commonData;
	
	/// \brief operator to calculate pressure gradient at Gauss points
	struct OpGetGradReAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
	
		CommonData &commonData;
		OpGetGradReAtGaussPts(const string field_name,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(field_name),
			commonData(common_data) {}
	
		/** \brief operator calculating pressure gradients
		  *
		  * pressure gradient is calculated multiplying derivatives of shape functions by degrees of freedom.
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			try {
				
				//if(data.getIndices().size()==0) PetscFunctionReturn(0);
				if(data.getFieldData().size()==0) PetscFunctionReturn(0);
				int nb_dofs = data.getFieldData().size();
				int nb_gauss_pts = data.getN().size1();
	
				//initialize
				commonData.gradReAtGaussPts.resize(nb_gauss_pts,3);
				if(type == MBVERTEX) {
					fill(commonData.gradReAtGaussPts.data().begin(),commonData.gradReAtGaussPts.data().end(),0);
				}
	
				for(int gg = 0;gg<nb_gauss_pts;gg++) {
					std::cout << "\n data.getFieldData = \n" << data.getFieldData() << std::endl;
					ublas::noalias(commonData.getGradReAtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
				}
	
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
	
			PetscFunctionReturn(0);
		}
	
	};
	
	/// \brief operator to calculate pressure gradient at Gauss points
	struct OpGetGradImAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
		
		CommonData &commonData;
		OpGetGradImAtGaussPts(const string field_name,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(field_name),
			commonData(common_data) {}
		
		/** \brief operator calculating pressure gradients
		  *
		  * pressure gradient is calculated multiplying derivatives of shape functions by degrees of freedom.
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			try {
				
				//if(data.getIndices().size()==0) PetscFunctionReturn(0);
				if(data.getFieldData().size()==0) PetscFunctionReturn(0);
				int nb_dofs = data.getFieldData().size();
				int nb_gauss_pts = data.getN().size1();
				
				//initialize
				commonData.gradReAtGaussPts.resize(nb_gauss_pts,3);
				if(type == MBVERTEX) {
					fill(commonData.gradImAtGaussPts.data().begin(),commonData.gradImAtGaussPts.data().end(),0);
				}
				
				for(int gg = 0;gg<nb_gauss_pts;gg++) {
					ublas::noalias(commonData.getGradImAtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
				}
				
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
			
			PetscFunctionReturn(0);
		}
		
	};
	
	/** \brief opearator to caulate pressure  and rate of pressure at Gauss points
	  * \infroup mofem_helmholtz_elem
	  */
	template<typename OP>
	struct OpGetFieldAtGaussPts: public OP::UserDataOperator {
				
		ublas::vector<double> &fieldAtGaussPts;
		OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
			OP::UserDataOperator(field_name),
			fieldAtGaussPts(field_at_gauss_pts) {}
		/** \brief operator calculating pressure and rate of pressure
		  *
		  * pressure pressure or rate of pressure is calculated multiplyingshape functions by degrees of freedom
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			try {
	
				//if(data.getIndices().size()==0) PetscFunctionReturn(0);
				if(data.getFieldData().size()==0) PetscFunctionReturn(0);
				int nb_dofs = data.getFieldData().size();
				int nb_gauss_pts = data.getN().size1();
	
				//initialize
				fieldAtGaussPts.resize(nb_gauss_pts);
				if(type == MBVERTEX) {
					//loop over shape functions on entities allways start from
					//vertices, so if nodal shape functions are processed, vector of
					//field values is zeroad at initialization
					fill(fieldAtGaussPts.begin(),fieldAtGaussPts.end(),0);
				}
	
				for(int gg = 0;gg<nb_gauss_pts;gg++) {
					std::cout << "\n data.getFieldData = \n" << data.getFieldData() << std::endl;
					fieldAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
					
				}
	
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
			}
	
			PetscFunctionReturn(0);
		}
	
	};

	/** \brief operator to calculate pressure at Gauss pts
	* \infroup mofem_helmholtz_elem
	*/
	struct OpGetTetPressureReAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetPressureReAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field_name,common_data.pressureReAtGaussPts) {}
	};
	struct OpGetTetPressureImAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetPressureImAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field_name,common_data.pressureImAtGaussPts) {}
	};
	
	/** \brief operator to calculate pressure at Gauss pts
	* \infroup mofem_helmholtz_elem
	*/
	struct OpGetTriPressureReAtGaussPts: public OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore> {
		OpGetTriPressureReAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore>(field_name,common_data.pressureReAtGaussPts) {}
	};
	
	struct OpGetTriPressureImAtGaussPts: public OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore> {
		OpGetTriPressureImAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore>(field_name,common_data.pressureImAtGaussPts) {}
	};
	
	/** \brief operator to calculate Pressure rate at Gauss pts
	  * \infroup mofem_helmholtz_elem
	  */
	struct OpGetTetRateAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetRateAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field_name,common_data.pressureRateAtGaussPts) {}
	};
	
    ublas::matrix<double> hoCoordsTri;
	
	//get the higher order approximation coordinates.
    struct OpHoCoordTri: public TriElementForcesAndSurcesCore::UserDataOperator {
	
		ublas::matrix<double> &hoCoordsTri;
		OpHoCoordTri(const string field_name,ublas::matrix<double> &ho_coords): 
			TriElementForcesAndSurcesCore::UserDataOperator(field_name),
			hoCoordsTri(ho_coords) {}
	
		/*	
		Cartesian coordinates for gaussian points inside elements
		X^coordinates = DOF dot* N
		*/
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
	
			try {
	
				if(data.getFieldData().size()==0) PetscFunctionReturn(0);
				//if(data.getIndices().size()==0) PetscFunctionReturn(0); //this will return zero
	
				hoCoordsTri.resize(data.getN().size1(),3);
				if(type == MBVERTEX) {
					hoCoordsTri.clear();
				}
	
				int nb_dofs = data.getN().size2();
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
					for(int dd = 0;dd<3;dd++) {
						hoCoordsTri(gg,dd) += cblas_ddot(nb_dofs,&data.getN(gg)[0],1,&data.getFieldData()[dd],3); //calculate x,y,z in each GaussPts
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
	
	
	/** \biref operator to calculate right hand side of stiffness terms
	* \infroup mofem_helmholtz_elem
	*/
	struct OpHelmholtzRhs: public TetElementForcesAndSourcesCore::UserDataOperator {
	
		BlockData &dAta;
		CommonData &commonData;
		bool useTsF;
		bool useImim;
		string field_name;
		OpHelmholtzRhs(const string re_field_name,const string im_field_name,BlockData &data,CommonData &common_data,bool imim):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useImim(imim),field_name(re_field_name),useTsF(true) {}
	
		Vec F;
		OpHelmholtzRhs(const string re_field_name,const string im_field_name,Vec _F,BlockData &data,CommonData &common_data,bool imim):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useImim(imim),field_name(re_field_name),useTsF(false),F(_F) { }
	
		ublas::vector<double> Nf;
	
		/** \brief calculate helmholtz operator apply on lift.
		  *
		  * F = int diffN^T gard_T0 - N^T k T0 dOmega^2
		  *
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
	
			if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
				PetscFunctionReturn(0);
			}
	
			try {
	
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);
	
				PetscErrorCode ierr;
	
				int nb_row_dofs = data.getIndices().size();
				Nf.resize(nb_row_dofs);
				bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));
	
				double wAvenumber = dAta.aNgularfreq/dAta.sPeed;  //wave number K is the propotional to the frequency of 
				//incident wave and represents number of waves per wave length 2Pi - 2Pi/K 				
				double wAvenUmber = pow(wAvenumber,2.0);
				

				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	
					double val = getVolume()*getGaussPts()(3,gg);
	
					if(getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
					}
					
					//cerr << val << endl;
					//cerr << data.getDiffN() << endl;
					//cerr << data.getIndices() << endl;
					//cerr << commonData.gradAtGaussPts << endl;
					//cblas
					//cblas_dgemv(CblasRowMajor,CblasNoTrans,nb_row_dofs,3,val,
					//&data.getDiffN()(gg,0),3,&commonData.gradAtGaussPts(gg,0),1,
					//1.,&Nf[0],1);
					double const_kT_Re = wAvenUmber*commonData.pressureReAtGaussPts[gg];
					double const_kT_Im = wAvenUmber*commonData.pressureImAtGaussPts[gg];
					std::cout << "\n data.getDiffN(gg,nb_row_dofs) = \n" << data.getDiffN(gg,nb_row_dofs) << std::endl;
					std::cout << "\n commonData.getGradReAtGaussPts(gg) =\n" << commonData.getGradReAtGaussPts(gg) << std::endl;
					std::cout << "\n data.getIndices().size() = \n" << data.getIndices().size() << std::endl;
					//ublas::vector<double> dNdT_Re(data.getN().size2());
					ublas::vector<double> dNdT_Re(data.getIndices().size());
					dNdT_Re.clear();
					ublas::noalias(dNdT_Re) = prod(data.getDiffN(gg,nb_row_dofs),commonData.getGradReAtGaussPts(gg));
					ublas::vector<double> dNdT_Im(data.getIndices().size());
					dNdT_Im.clear();
					ublas::noalias(dNdT_Im) = prod(data.getDiffN(gg,nb_row_dofs),commonData.getGradImAtGaussPts(gg));
					if(field_name.compare(0,6,"rePRES") == 0) {
	                
					
					//std::string wait;
					//std::cout <<"\n is imim true ???? "<< useImim << std::endl;
					//std::cout << "\n commonData.pressureAtGaussPts[gg] = \n" << commonData.pressureAtGaussPts[gg] << std::endl;
					
					//ublas
					//wait to be fixed and improved here.
					
					//if(!useImim) {
					//	std::string wait;
					//	std::cout << "\n commonData.pressureAtGaussPts_re = \n" << commonData.pressureAtGaussPts << std::endl;
					//}
					//else {
					//	std::string wait;
					//	std::cout << "\n commonData.pressureAtGaussPts_im = \n" << commonData.pressureAtGaussPts << std::endl;
					//}
						
					cblas_daxpy(data.getN().size2(), -const_kT_Re, &data.getN(gg)[0], 1, &dNdT_Re[0], 1);
					cblas_daxpy(data.getN().size2(), -const_kT_Im, &data.getN(gg)[0], 1, &dNdT_Im[0], 1);
					ublas::noalias(Nf) += val*(dNdT_Re-dNdT_Im); 
					} else if(field_name.compare(0,6,"imPRES") == 0) {
					cblas_daxpy(data.getN().size2(), -const_kT_Re, &data.getN(gg)[0], 1, &dNdT_Re[0], 1);
					cblas_daxpy(data.getN().size2(), -const_kT_Im, &data.getN(gg)[0], 1, &dNdT_Im[0], 1);
					ublas::noalias(Nf) += val*(dNdT_Re+dNdT_Im); 
					}
					//if(useImim) {
					//ublas::noalias(Nf) += -val*(dNdT); }
					//else{
					//ublas::noalias(Nf) += val*(dNdT); }
					
				}
	
				//cerr << Nf << endl;
				if(useTsF) {
					ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
										&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
				} else {
					ierr = VecSetValues(F,data.getIndices().size(),
										&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
	
				}
	
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
	
			PetscFunctionReturn(0);
		}
	
	};
	
	
	/** \biref operator to calculate right hand side wave source term F
	* \infroup mofem_helmholtz_elem
	*/
	struct OpHelmholtzRhs_F: public TetElementForcesAndSourcesCore::UserDataOperator {
		
		BlockData &dAta;
		CommonData &commonData;
		bool useTsF;
		bool useScalar;
		OpHelmholtzRhs_F(const string re_field_name,const string im_field_name,BlockData &data,CommonData &common_data,bool usescalar):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useTsF(false),useScalar(usescalar) {}  //here useTsF(ture) originally.
		
		Vec F;
		OpHelmholtzRhs_F(const string re_field_name,const string im_field_name,Vec _F,BlockData &data,CommonData &common_data,bool usescalar):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useTsF(false),F(_F),useScalar(usescalar) { }
		
		ublas::vector<double> Nf;
		
		/** \brief calculate Helmholtz RHS source term.
		  *
		  * F = int diffN^T F(x,y,z) dOmega^2
		  *
		  */
				
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			
			if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
				PetscFunctionReturn(0);
			}
			
			
			try {
				
				PetscErrorCode ierr;
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				int nb_row = data.getN().size2();
				Nf.resize(nb_row);
				bzero(&Nf[0],nb_row*sizeof(double));
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
					
					double val = getVolume()*getGaussPts()(3,gg);
					
					if(getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
					}
				 ublas::vector<double> F_S;
				 //Scalar value source term F
				    if(useScalar) {
				    const double f_scalar = -1; //scalar source vector F=1, wait for modification. 
				    ublas::scalar_vector<double> F (data.getN().size1(),f_scalar); //scalar source vector 
					F_S = F;
				    } 
				 //Functional value source term F
				    else {
						//ublas::vector<double> F (data.getN().size1());
						//double amplitude = 1;
						//double wAvenumber = dAta.aNgularfreq/dAta.sPeed;  //wave number K is the propotional to the frequency of 
						//ublas::vector<double> d (3); //incident wave direction. Unit vector.wait to move to .cpp executable file.
						//d (0) = 1;  // x direction
						//d (1) = 0;
						//d (2) = 0;
						//ublas::vector<double> sPatialpsn (3); //spatial position of the nodes.
						//sPatialpsn (0) = getCoordsAtGaussPts()(gg,0);
						//sPatialpsn (1) = getCoordsAtGaussPts()(gg,1);
						//sPatialpsn (2) = getCoordsAtGaussPts()(gg,2);
						//double const1 = inner_prod(sPatialpsn,d);
						//double iNcidentwave;
						//if(d (0) != 0 && d (1) != 0 && d (2) != 0){
						//	ublas::vector<double> const2 (3);
						//	const2 (0) = sPatialpsn (0)*exp(wAvenumber*const1);
						//	const2 (1) = sPatialpsn (1)*exp(wAvenumber*const1);
						//	const2 (2) = sPatialpsn (2)*exp(wAvenumber*const1);
						//double iNcidentwave = inner_prod(const2,getNormal());}
						//else if(d (0) != 0){
						//	ublas::vector<double> const2 (3);
						//	const2 (0) = sPatialpsn (0)*exp(wAvenumber*const1);
						//	const2 (1) = 0;
						//	const2 (2) = 0;
						//iNcidentwave = inner_prod(const2,getNormal());}
						
						////double iNcidentwave = exp(wAvenumber*const1); 
						
						////std::string wait;
						////std::cout << "\n getNormal = \n" << getNormal() << std::endl;
						////std::cout << "\n getNormals_at_GaussPt() = \n" << getNormals_at_GaussPt() << std::endl;
						
						//	F = amplitude*iNcidentwave;  //FluxData.HeatfluxCubitBcData.data.value1 * area
						////cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
						////ublas::noalias(F) += val*flux*data.getN(gg,nb_dofs);
		
				    ublas::vector<double> F (data.getN().size1());
				    for(unsigned int jj = 0;jj<data.getN().size1();jj++) {
				    F (jj) = 1-pow(getCoordsAtGaussPts()(jj,0),2.0);           //wait for modification for F=1-x^2
					
				    }
					F_S =F;
				 //ublas::vector<double> x( y.size() ); 
				 //std::copy( y.begin(), y.end(), x.begin() ); 
					
				}
				val *= F_S[gg];
				ublas::noalias(Nf) += val*data.getN(gg);
				
				//cerr << data.getIndices() << endl;
				//cerr << data.getDiffN() << endl;
				
				//cerr << Nf << endl;
				}
				if(useTsF) {
					ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
										&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
				} else {
					ierr = VecSetValues(F,data.getIndices().size(),
										&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
					
				}
				
				
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
			
			PetscFunctionReturn(0);
		}
		
	};
	
	
	
	/** \biref operator to calculate left hand side of stiffness terms
	  * \infroup mofem_helmholtz_elem
	  */
	struct OpHelmholtzLhs_rereA: public TetElementForcesAndSourcesCore::UserDataOperator {
	
		BlockData &dAta;
		CommonData &commonData;
        bool useTsB;
		OpHelmholtzLhs_rereA(const string re_field_name,const string im_field_name,BlockData &data,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useTsB(true) { }
		
		Mat A;
		OpHelmholtzLhs_rereA(const string re_field_name,const string im_field_name,Mat _A,BlockData &data,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useTsB(false),A(_A) {}
	
		ublas::matrix<double> K,transK;
	
		/** \brief calculate helmholtz stiffness matrix
		  *
		  * K = int diffN^T  diffN^T dOmega^2
		  *
		  */
		PetscErrorCode doWork(
			int row_side,int col_side,
			EntityType row_type,EntityType col_type,
			DataForcesAndSurcesCore::EntData &row_data,
			DataForcesAndSurcesCore::EntData &col_data) {
			PetscFunctionBegin;
	
			if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
				PetscFunctionReturn(0);
			}
	
			try {
	
				if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
				if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
	
				int nb_row = row_data.getN().size2();
				int nb_col = col_data.getN().size2();
				K.resize(nb_row,nb_col);
				bzero(&*K.data().begin(),nb_row*nb_col*sizeof(double));
	
				for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	
					double val = getVolume()*getGaussPts()(3,gg);
					if(getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
					}
	
					//cblas
					//double *diff_N_row,*diff_N_col;
					//diff_N_row = &row_data.getDiffN()(gg,0);
					//diff_N_col = &col_data.getDiffN()(gg,0);
					//cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
					//nb_row,nb_col,3,
					//val,diff_N_row,3,diff_N_col,3,1.,&K(0,0),nb_col);
	
					//ublas
					noalias(K) += val*prod(row_data.getDiffN(gg,nb_row),trans(col_data.getDiffN(gg,nb_col)));
	
				}
	
				PetscErrorCode ierr;
				if(!useTsB) {
					const_cast<FEMethod*>(getFEMethod())->ts_B = A;    //FEMethod does not belong to fieldinterface anymore.
				}
				ierr = MatSetValues(
						   A,  //(getFEMethod()->ts_B), instead in New thermal element. wait
						   nb_row,&row_data.getIndices()[0],
						   nb_col,&col_data.getIndices()[0],
						   &K(0,0),ADD_VALUES); CHKERRQ(ierr);
				if(row_side != col_side || row_type != col_type) {
					transK.resize(nb_col,nb_row);
					noalias(transK) = trans( K );
					ierr = MatSetValues(
							   A,
							   nb_col,&col_data.getIndices()[0],
							   nb_row,&row_data.getIndices()[0],
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

	struct OpHelmholtzLhs_imimA: public TetElementForcesAndSourcesCore::UserDataOperator {
		
		BlockData &dAta;
		CommonData &commonData;
	
        bool useTsB;
		OpHelmholtzLhs_imimA(const string re_field_name,const string im_field_name,BlockData &data,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useTsB(true) { }
		
		Mat A;
		OpHelmholtzLhs_imimA(const string re_field_name,const string im_field_name,Mat _A,BlockData &data,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useTsB(false),A(_A) {}
		
		ublas::matrix<double> K,transK;
		
		/** \brief calculate helmholtz stiffness matrix
		  *
		  * K = int diffN^T  diffN^T dOmega^2
		  *
		  */
		PetscErrorCode doWork(
			int row_side,int col_side,
			EntityType row_type,EntityType col_type,
			DataForcesAndSurcesCore::EntData &row_data,
			DataForcesAndSurcesCore::EntData &col_data) {
			PetscFunctionBegin;
			
			if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
				PetscFunctionReturn(0);
			}
			
			try {
				
				if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
				if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
				
				int nb_row = row_data.getN().size2();
				int nb_col = col_data.getN().size2();
				K.resize(nb_row,nb_col);
				bzero(&*K.data().begin(),nb_row*nb_col*sizeof(double));
				
				for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
					
					double val = getVolume()*getGaussPts()(3,gg);
					if(getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
					}
					
					//cblas
					//double *diff_N_row,*diff_N_col;
					//diff_N_row = &row_data.getDiffN()(gg,0);
					//diff_N_col = &col_data.getDiffN()(gg,0);
					//cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,
					//nb_row,nb_col,3,
					//val,diff_N_row,3,diff_N_col,3,1.,&K(0,0),nb_col);
					
					//ublas
					noalias(K) += val*prod(row_data.getDiffN(gg,nb_row),trans(col_data.getDiffN(gg,nb_col)));
					
				}
				
				PetscErrorCode ierr;
				if(!useTsB) {
					const_cast<FEMethod*>(getFEMethod())->ts_B = A;
				}
				ierr = MatSetValues(
						   A,
						   nb_row,&row_data.getIndices()[0],
						   nb_col,&col_data.getIndices()[0],
						   &K(0,0),ADD_VALUES); CHKERRQ(ierr);
				if(row_side != col_side || row_type != col_type) {
					transK.resize(nb_col,nb_row);
					noalias(transK) = trans( K );
					ierr = MatSetValues(
							   A,
							   nb_col,&col_data.getIndices()[0],
							   nb_row,&row_data.getIndices()[0],
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

	/** \brief operator to calculate left hand side of element mass matrix terms
	  * \infroup mofem_helmholtz_elem * Lumped mass matrix does not applied to Hierarchical shape functions.
	  */
	struct OpMassLsh: public TetElementForcesAndSourcesCore::UserDataOperator {
	
		BlockData &dAta;
		CommonData &commonData;
		bool useTsB;
		Mat A;
		OpMassLsh(const string re_field_name,const string im_field_name,Mat _A,BlockData &data,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),commonData(common_data),useTsB(false),A(_A) {}
	
		ublas::matrix<double> M,transM;
	
		/** \brief calculate Helmholtz Mass matrix
		  *
		  * M = int minus N^T K^2 N dOmega^2
		  *
		  */
		PetscErrorCode doWork(
			int row_side,int col_side,
			EntityType row_type,EntityType col_type,
			DataForcesAndSurcesCore::EntData &row_data,
			DataForcesAndSurcesCore::EntData &col_data) {
			PetscFunctionBegin;
	
			try {
	
				if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
				if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
	
				int nb_row = row_data.getN().size2();
				int nb_col = col_data.getN().size2();
				M.resize(nb_row,nb_col);
				bzero(&*M.data().begin(),nb_row*nb_col*sizeof(double));
	
				for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	                 
					double wAvenumber = dAta.aNgularfreq/dAta.sPeed;  //wave number K is the propotional to the frequency of 
					                                       // wave and represents number of waves per wave length 2Pi - 2Pi/K
					
					
					double wAvenUmber = pow(wAvenumber,2.0);
					
					double val = -wAvenUmber*getVolume()*getGaussPts()(3,gg);
					if(getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
					}
					val *= getFEMethod()->ts_a;
					//cblas
					//double *N_row,*N_col;
					//N_row = &row_data.getN()(gg,0);
					//N_col = &col_data.getN()(gg,0);
					//cblas_dger(CblasRowMajor,
					//  nb_row,nb_col,val,N_row,1,N_col,1,&M(0,0),nb_col);
					//ublas
					noalias(M) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
					
				}
	
				PetscErrorCode ierr;
				if(!useTsB) {
					const_cast<FEMethod*>(getFEMethod())->ts_B = A;
				}
				ierr = MatSetValues(
						   (getFEMethod()->ts_B),
						   nb_row,&row_data.getIndices()[0],
						   nb_col,&col_data.getIndices()[0],
						   &M(0,0),ADD_VALUES); CHKERRQ(ierr);
				if(row_side != col_side || row_type != col_type) {
					transM.resize(nb_col,nb_row);
					noalias(transM) = trans(M);
					ierr = MatSetValues(
							   (getFEMethod()->ts_B),
							   nb_col,&col_data.getIndices()[0],
							   nb_row,&row_data.getIndices()[0],
							   &transM(0,0),ADD_VALUES); CHKERRQ(ierr);
				}
	
	
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
	
	
			PetscFunctionReturn(0);
		}
	
	};
	




	/** \brief operator for calculate Impedance flux and assemble to right hand side
	  * \infroup mofem_helmholtz_elem
	  */
	struct OpHelmholtzFlux:public TriElementForcesAndSurcesCore::UserDataOperator {
	    
		BlockData &datA;
		FluxData &dAta;
		bool ho_geometry;
		bool useTsF;
		OpHelmholtzFlux(const string re_field_name,const string im_field_name,BlockData &DATA,FluxData &data,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			datA(DATA),dAta(data),ho_geometry(_ho_geometry),useTsF(true) { }
	
		Vec F;
		OpHelmholtzFlux(const string re_field_name,const string im_field_name,Vec _F,
				   BlockData &DATA,FluxData &data,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			datA(DATA),dAta(data),ho_geometry(_ho_geometry),useTsF(false),F(_F) { }
	
		ublas::vector<FieldData> Nf;
	
		/** \brief calculate Impedance flux
		  *
		  * F = int_S N^T * g dS  , F=0 represents sound hard wave (mere reflecting no absorption) 
		  *
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
	
			if(data.getIndices().size()==0) PetscFunctionReturn(0);
			if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
	
			PetscErrorCode ierr;
	
			const FENumeredDofMoFEMEntity *dof_ptr;
			ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
			int rank = dof_ptr->get_max_rank();
	
			int nb_dofs = data.getIndices().size()/rank;
	
			Nf.resize(data.getIndices().size());
			//bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));
			fill(Nf.begin(),Nf.end(),0);
	
			//cerr << getNormal() << endl;
			//cerr << getNormals_at_GaussPt() << endl;
	
			for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	
				double val = getGaussPts()(2,gg);
				double flux;
				
				//Torus Example
				//double k = datA.aNgularfreq/datA.sPeed;  //wave number K is the propotional to the frequency of 
				//double nEumanntOrus;
				//ublas::vector<double> const1 (3); //spatial position of the nodes.
				//const1 (0) = cos(k*getCoordsAtGaussPts()(gg,0)/sqrt(3))*sin(k*getCoordsAtGaussPts()(gg,1)/sqrt(3))*sin(k*getCoordsAtGaussPts()(gg,2)/sqrt(3));
				//const1 (1) = sin(k*getCoordsAtGaussPts()(gg,0)/sqrt(3))*cos(k*getCoordsAtGaussPts()(gg,1)/sqrt(3))*sin(k*getCoordsAtGaussPts()(gg,2)/sqrt(3));
				//const1 (2) = sin(k*getCoordsAtGaussPts()(gg,0)/sqrt(3))*sin(k*getCoordsAtGaussPts()(gg,1)/sqrt(3))*cos(k*getCoordsAtGaussPts()(gg,2)/sqrt(3));
                //double const2 = inner_prod(const1,-getNormal());
				//nEumanntOrus = (k/sqrt(3))*const2;
				

				
				if(ho_geometry) {
					double area = norm_2(getNormals_at_GaussPt(gg)); //cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
					flux = dAta.dAta.data.value1*area;  //FluxData.HeatfluxCubitBcData.data.value1 * area
				} else {
					flux = dAta.dAta.data.value1*getArea();
				}
				//cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
				ublas::noalias(Nf) += val*flux*data.getN(gg,nb_dofs);
	
			}
	
			//cerr << "VecSetValues\n";
			//cerr << Nf << endl;
			//cerr << data.getIndices() << endl;
	
			if(useTsF) {
				ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
									&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
			} else {
				ierr = VecSetValues(F,data.getIndices().size(),
									&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
			}
	
			PetscFunctionReturn(0);
		}
	
	};
	
	
	/** \brief operator for calculate Impedance flux and assemble to right hand side
	* \infroup mofem_helmholtz_elem
	*/
	struct OpHelmholtzIncidentWave:public TriElementForcesAndSurcesCore::UserDataOperator {
		
		ublas::matrix<double> &hoCoords;
		BlockData &datA;
		FluxData &dAta;
		bool ho_geometry;
		bool useTsF;
		string re_field;
		string im_field;
		//double (*fUN)(double x,double y,double z);
		OpHelmholtzIncidentWave(const string re_field_name,const string im_field_name,BlockData &DATA,FluxData &data,ublas::matrix<double> &ho_coords,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			datA(DATA),dAta(data),ho_geometry(_ho_geometry),useTsF(true),re_field(re_field_name),im_field(im_field_name),hoCoords(ho_coords) { }
		
		Vec F;
		OpHelmholtzIncidentWave(const string re_field_name,const string im_field_name,Vec _F,
						BlockData &DATA,FluxData &data,ublas::matrix<double> &ho_coords,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			datA(DATA),dAta(data),ho_geometry(_ho_geometry),useTsF(false),F(_F),re_field(re_field_name),im_field(im_field_name),hoCoords(ho_coords) { }
		
		ublas::vector<FieldData> Nf;
		
		/** \brief calculate Impedance flux
		  *
		  * F = int_S N^T * g gradient (exp(i*k*x dot d) ) * n dS  , F=0 represents sound hard wave (mere reflecting no absorption) 
		  *
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			
			if(data.getIndices().size()==0) PetscFunctionReturn(0);
			if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
			
			PetscErrorCode ierr;
			
			const FENumeredDofMoFEMEntity *dof_ptr;
			ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
			int rank = dof_ptr->get_max_rank();
			
			int nb_dofs = data.getIndices().size()/rank;
			
			Nf.resize(data.getIndices().size());
			//bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));
			fill(Nf.begin(),Nf.end(),0);
			
			//cerr << getNormal() << endl;
			//cerr << getNormals_at_GaussPt() << endl;
			
			for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
				
				double val = getGaussPts()(2,gg);
				
				double flux;
				double wAvenumber = datA.aNgularfreq/datA.sPeed;  //wave number K is the propotional to the frequency of 
				double iNcidentwave;
				
				//complex< double > result;
				//ublas::vector< complex< double > > const2 (3);
				//const complex< double > i( 0.0, 1.0 );
				//ublas::vector<double> d (3); //incident wave direction. Unit vector.wait to move to .cpp executable file.
				//d (0) = 1;  // x direction
				//d (1) = 0;
				//d (2) = 0;
				//ublas::vector<double> sPatialpsn (3); //spatial position of the nodes.
				//sPatialpsn (0) = getCoordsAtGaussPts()(gg,0);
				//sPatialpsn (1) = getCoordsAtGaussPts()(gg,1);
				//sPatialpsn (2) = getCoordsAtGaussPts()(gg,2);
				//double const1 = inner_prod(sPatialpsn,d);
				
				//if(d (0) != 0 && d (1) != 0 && d (2) != 0){
				//	//ublas::vector<double> const2 (3);
				//	const2 (0) = i*sPatialpsn (0)*exp(i*wAvenumber*const1);
				//	const2 (1) = i*sPatialpsn (1)*exp(i*wAvenumber*const1);
				//	const2 (2) = i*sPatialpsn (2)*exp(i*wAvenumber*const1);
				//    result = inner_prod(const2,getNormal());}
				//else if(d (0) != 0){
				//	//ublas::vector<double> const2 (3);
				//	const2 (0) = i*sPatialpsn (0)*exp(i*wAvenumber*const1);
				//	const2 (1) = 0;
				//	const2 (2) = 0;
				//    result = inner_prod(const2,getNormal());}
				////////////////////********************************/////////////
				unsigned int nb_gauss_pts = data.getN().size1();
				double x,y,z;
					
					if(hoCoords.size1() == data.getN().size1()) {
					//if(hoCoords.size1() == nb_gauss_pts) {
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
					
					
					

				const double pi = atan( 1.0 ) * 4.0;
				double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
				double theta = atan2(y,x)+2*pi; //the arctan of radians (y/x)
				const double k = wAvenumber;  //Wave number
				const double a = 0.5;         //radius of the sphere,wait to modify by user
				const double const1 = k * a;
				//double const2 = k * R;
				
				const complex< double > i( 0.0, 1.0 );
				
				// magnitude of incident wave
				//const double phi_incident_mag = 1.0;
				
				const double tol = 1.0e-10;
				double max = 0.0;
				double min = 999999.0;
				
				complex< double > result = 0.0;
				complex< double > prev_result;
				
				double error = 100.0;
				unsigned int n = 0; //initialized the infinite series loop
				
				while( error > tol )  //finding the acoustic potential in one single point.
				{
					double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );  //The derivative of bessel function
					//complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
					//complex< double > hn_der = 0.5 * ( sph_hankel_1( n - 1, const1 ) -
					//( sph_hankel_1( n, const1 ) + const1 * sph_hankel_1( n + 1, const1 ) ) / const1 );
					double Pn = legendre_p( n, cos( theta ) );  //Legendre
					//complex< double >hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
					prev_result = result;
					result -= k * pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der * Pn;
					error = abs( abs( result ) - abs( prev_result ) );
					++n;
				}
				
				if(re_field.compare(0,6,"rePRES")==0){
					
				iNcidentwave = std::real(result);
				} else if(re_field.compare(0,6,"imPRES")==0)
					{iNcidentwave = std::imag(result);}
				

				flux = dAta.dAta.data.value1*iNcidentwave;  //FluxData.HeatfluxCubitBcData.data.value1 * area
			
			
				//cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
				ublas::noalias(Nf) += val*flux*data.getN(gg,nb_dofs);
			}
			
			std::cout << "\n Nf = \n" << Nf << std::endl;
			//cerr << "VecSetValues\n";
			//cerr << Nf << endl;
			//cerr << data.getIndices() << endl;
			
			if(useTsF) {
				ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
									&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
			} else {
				ierr = VecSetValues(F,data.getIndices().size(),
									&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
			}
			
			PetscFunctionReturn(0);
		}
		
	};
	

	/// \biref operator to calculate Impedance on body surface and assemble to imaginary lhs of equations 
	struct OpImpedanceLhs_reimC:public TriElementForcesAndSurcesCore::UserDataOperator {
	
		ImpedanceData &dAta;
		BlockData &blockData;
		bool ho_geometry;
		bool useTsB;
	
		OpImpedanceLhs_reimC(const string re_field_name,const string im_field_name,
						ImpedanceData &data,BlockData &block_data,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),ho_geometry(_ho_geometry),useTsB(true),blockData(block_data) {}
	
		Mat A;
		OpImpedanceLhs_reimC(const string re_field_name,const string im_field_name,Mat _A,
						ImpedanceData &data,BlockData &block_data,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),ho_geometry(_ho_geometry),useTsB(false),A(_A),blockData(block_data) {}
	
		ublas::matrix<double> K,transK;
		/** \brief calculate helmholtz Impedance term in the imaginary lhs of equations
		 *
		 * K = intS N^T i sIgma K N dS
		 */
		PetscErrorCode doWork(
			int row_side,int col_side,
			EntityType row_type,EntityType col_type,
			DataForcesAndSurcesCore::EntData &row_data,
			DataForcesAndSurcesCore::EntData &col_data) {
			PetscFunctionBegin;
	
			try {
	
				if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
				if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
	
				int nb_row = row_data.getN().size2();
				int nb_col = col_data.getN().size2();
				K.resize(nb_row,nb_col);
				bzero(&*K.data().begin(),nb_row*nb_col*sizeof(double));
	
				for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {  /* Integrate the shape functions in one element */
					double waveNumber = blockData.aNgularfreq/blockData.sPeed;
                    
				    //std::string wait;
					//std::cout << "\n waveNumber = \n" << waveNumber << std::endl;
        					
					double impedanceConst;
					if(ho_geometry) {
						double area = norm_2(getNormals_at_GaussPt(gg));
						impedanceConst = -dAta.sIgma*waveNumber*area;
					}   else {
						impedanceConst = -dAta.sIgma*waveNumber*getArea();
					}
					double val = getGaussPts()(2,gg)*impedanceConst;
					noalias(K) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
	
				}
	
				PetscErrorCode ierr;
				if(!useTsB) {
					const_cast<FEMethod*>(getFEMethod())->ts_B = A;  //to change the pointer getFEMethod()'s private member ts_B value with A.
				}                                                                    //getFEMethod() belong to class FieldInterface::FEMethod
				ierr = MatSetValues(
						   (getFEMethod()->ts_B),
						   nb_row,&row_data.getIndices()[0],              //MatSetValues(A,i,m,j,n,value,INSERT_VALUES,ierr), i*j the size of input matrix,
						   nb_col,&col_data.getIndices()[0],              // m and n - the position of parent matrix to allowed the input matrix.
						   &K(0,0),ADD_VALUES); CHKERRQ(ierr);
				if(row_side != col_side || row_type != col_type) {
					transK.resize(nb_col,nb_row);
					noalias(transK) = trans( K );
					ierr = MatSetValues(
							   (getFEMethod()->ts_B),
							   nb_col,&col_data.getIndices()[0],
							   nb_row,&row_data.getIndices()[0],
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

	
	/// \biref operator to calculate Impedance on body surface and assemble to imaginary lhs of equations 
	struct OpImpedanceLhs_imreC:public TriElementForcesAndSurcesCore::UserDataOperator {
		
		ImpedanceData &dAta;
		BlockData &blockData;
		bool ho_geometry;
		bool useTsB;
		
		OpImpedanceLhs_imreC(const string re_field_name,const string im_field_name,
					   ImpedanceData &data,BlockData &block_data,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),ho_geometry(_ho_geometry),useTsB(true),blockData(block_data) {}
		
		Mat A;
		OpImpedanceLhs_imreC(const string re_field_name,const string im_field_name,Mat _A,
					   ImpedanceData &data,BlockData &block_data,bool _ho_geometry = false):
			TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
			dAta(data),ho_geometry(_ho_geometry),useTsB(false),A(_A),blockData(block_data) {}
		
		ublas::matrix<double> K,transK;
		/** \brief calculate helmholtz Impedance term in the imaginary lhs of equations
		 *
		 * K = - intS N^T i sIgma K N dS
		 */
		PetscErrorCode doWork(
			int row_side,int col_side,
			EntityType row_type,EntityType col_type,
			DataForcesAndSurcesCore::EntData &row_data,
			DataForcesAndSurcesCore::EntData &col_data) {
			PetscFunctionBegin;
			
			try {
				
				if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
				if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
				
				int nb_row = row_data.getN().size2();
				int nb_col = col_data.getN().size2();
				K.resize(nb_row,nb_col);
				bzero(&*K.data().begin(),nb_row*nb_col*sizeof(double));
				
				for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {  /* Integrate the shape functions in one element */
					
					double waveNumber = blockData.aNgularfreq/blockData.sPeed; //for impedance BC
					//double waveNumber = 1.0;                                   //for neumann BC

					double impedanceConst;
					if(ho_geometry) {
						double area = norm_2(getNormals_at_GaussPt(gg));
						impedanceConst = dAta.sIgma*waveNumber*area;
					}   else {
						impedanceConst = dAta.sIgma*waveNumber*getArea();
					}
					double val = getGaussPts()(2,gg)*impedanceConst;
					noalias(K) += val*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
					
				}
				
				PetscErrorCode ierr;
				if(!useTsB) {
					const_cast<FEMethod*>(getFEMethod())->ts_B = A;  //to change the pointer getFEMethod()'s private member ts_B value with A.
				}                                                                    //getFEMethod() belong to class FieldInterface::FEMethod
				ierr = MatSetValues(
						   (getFEMethod()->ts_B),
						   nb_row,&row_data.getIndices()[0],              //MatSetValues(A,i,m,j,n,value,INSERT_VALUES,ierr), i*j the size of input matrix,
						   nb_col,&col_data.getIndices()[0],              // m and n - the position of parent matrix to allowed the input matrix.
						   &K(0,0),ADD_VALUES); CHKERRQ(ierr);
				if(row_side != col_side || row_type != col_type) {
					transK.resize(nb_col,nb_row);
					noalias(transK) = trans( K );
					ierr = MatSetValues(
							   (getFEMethod()->ts_B),
							   nb_col,&col_data.getIndices()[0],
							   nb_row,&row_data.getIndices()[0],
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
	
	/* The functional load vector
	*/
	//const EntityHandle* conn_face;
	//int num_nodes;
	//getCoordsAtGaussPts()(gg,d);//d=x || y || z 
	
	
	/** \brief add helmholtz element on tets
	  * \infroup mofem_helmholtz_elem
	  *
	  * It get data from block set and define element in moab
	  *w
	  * \param problem name
	  * \param field name
	  * \param name of mesh nodal positions (if not defined nodal coordinates are used)
	  */
	PetscErrorCode addHelmholtzElements(
		const string problem_name,const string re_field_name,const string im_field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
		PetscFunctionBegin;
	
		PetscErrorCode ierr;
		ErrorCode rval;
		
		//to do matrix A and C; add Finite Element "HELMHOLTZ_FE"
		ierr = mField.add_finite_element("HELMHOLTZ_FE",MF_ZERO); CHKERRQ(ierr ); 
		ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_FE",im_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_FE",im_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FE",im_field_name); CHKERRQ(ierr);
		
		if(mField.check_field(mesh_nodals_positions)) {
			ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FE",mesh_nodals_positions); CHKERRQ(ierr);
		}
		ierr = mField.modify_problem_add_finite_element(problem_name,"HELMHOLTZ_FE"); CHKERRQ(ierr);
	
		//takes skin of block of entities
		//Skinner skin(&mField.get_moab());
		// loop over all blocksets and get data which name is FluidPressure
		
		for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
				
			if(it->get_Cubit_name().compare(0,13,"MAT_HELMHOLTZ") == 0) {
				
				//get block attributes
				vector<double> attributes;
				ierr = it->get_Cubit_attributes(attributes); CHKERRQ(ierr);
				if(attributes.size()<2) {
					SETERRQ1(PETSC_COMM_SELF,1,"not enough block attributes to deffine fluid pressure element, attributes.size() = %d ",attributes.size());
				}
				
				blockData.aNgularfreq = attributes[0];
				blockData.sPeed = attributes[1];
				setOfBlocks[it->get_msId()].aNgularfreq = attributes[0];
				setOfBlocks[it->get_msId()].sPeed = attributes[1];
				
				rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
				ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"HELMHOLTZ_FE"); CHKERRQ(ierr);
			}
		}
		
		//for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_HELMHOLTZSET,it)) {
		//	Mat_Helmholtz pres_data;
		//	ierr = it->get_attribute_data_structure(pres_data); CHKERRQ(ierr);
		//	//wait
		//	blockData.aNgularfreq = pres_data.data.AngularFreq;
		//	blockData.sPeed = pres_data.data.Speed;
		//	//wait
		//	setOfBlocks[it->get_msId()].aNgularfreq = pres_data.data.AngularFreq;
		//	setOfBlocks[it->get_msId()].sPeed = pres_data.data.Speed;
		//	rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
		//	ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"HELMHOLTZ_FE"); CHKERRQ(ierr);
	
		//}
	
		PetscFunctionReturn(0);
	}	
	
	/** \brief add helmholtz flux element
	  * \infroup mofem_helmholtz_elem
	  *
	  * It get data from helmholtz flux set and define elemenet in moab. Alternatively
	  * uses block set with name HELMHOLTZ_FLUX.
	  *
	  * \param problem name
	  * \param field name
	  * \param name of mesh nodal positions (if not defined nodal coordinates are used)
	  */
	PetscErrorCode addHelmholtzFluxElement(
		const string problem_name,const string re_field_name,const string im_field_name) {
		PetscFunctionBegin;
	
		PetscErrorCode ierr;
		ErrorCode rval;
	
		ierr = mField.add_finite_element("HELMHOLTZ_FLUX_FE",MF_ZERO); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_FLUX_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_FLUX_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FLUX_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_FLUX_FE",im_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_FLUX_FE",im_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FLUX_FE",im_field_name); CHKERRQ(ierr);
		
		if(mField.check_field("MESH_NODE_POSITIONS")) {
			ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_FLUX_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
		}
		
		ierr = mField.modify_problem_add_finite_element(problem_name,"HELMHOLTZ_FLUX_FE"); CHKERRQ(ierr);
	
		for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|HEATFLUXSET,it)) {
			ierr = it->get_cubit_bc_data_structure(setOfFluxes[it->get_msId()].dAta); CHKERRQ(ierr);
			rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfFluxes[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
			ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluxes[it->get_msId()].tRis,"HELMHOLTZ_FLUX_FE"); CHKERRQ(ierr);
		}
	
		//this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
		//not elegant, but good enough
		for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
			if(it->get_Cubit_name().compare(0,9,"HELMHOLTZ_FLUX") == 0) {
				vector<double> data;
				ierr = it->get_Cubit_attributes(data); CHKERRQ(ierr);
				if(data.size()!=1) {
					SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
				}
				strcpy(setOfFluxes[it->get_msId()].dAta.data.name,"HelmholtzFlux");
				setOfFluxes[it->get_msId()].dAta.data.flag1 = 1;
				setOfFluxes[it->get_msId()].dAta.data.value1 = data[0];
				//cerr << setOfFluxes[it->get_msId()].dAta << endl;
				rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfFluxes[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
				ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluxes[it->get_msId()].tRis,"HELMHOLTZ_FLUX_FE"); CHKERRQ(ierr);
	
			}
		}
	
	
		PetscFunctionReturn(0);
	}


	/** \brief add Impedance element
	* \infroup mofem_helmholtz_elem
	*
	* It get data from convection set and define elemenet in moab. Alternatively
	* uses block set with name IMPEDANCE.
	*
	* \param problem name
	* \param field name
	* \param name of mesh nodal positions (if not defined nodal coordinates are used)
	*/
	PetscErrorCode addHelmholtzImpedanceElement(
		const string problem_name,const string re_field_name,const string im_field_name) {
		PetscFunctionBegin;
	
		PetscErrorCode ierr;
		ErrorCode rval;
	    
		ierr = mField.add_finite_element("HELMHOLTZ_IMPEDANCE_FE",MF_ZERO); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_IMPEDANCE_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_IMPEDANCE_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_IMPEDANCE_FE",re_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_row("HELMHOLTZ_IMPEDANCE_FE",im_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_col("HELMHOLTZ_IMPEDANCE_FE",im_field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_data("HELMHOLTZ_IMPEDANCE_FE",im_field_name); CHKERRQ(ierr);

		ierr = mField.modify_problem_add_finite_element(problem_name,"HELMHOLTZ_IMPEDANCE_FE"); CHKERRQ(ierr);
	
		//this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
		//not elegant, but good enough
		for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
			if(it->get_Cubit_name().compare(0,9,"IMPEDANCE") == 0) {
				
				vector<double> data;
				ierr = it->get_Cubit_attributes(data); CHKERRQ(ierr);
				if(data.size()!=2) {
					SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
				}
				setOfImpedance[it->get_msId()].g = data[0];
				setOfImpedance[it->get_msId()].sIgma = data[1];
				//cerr << setOfFluxes[it->get_msId()].dAta << endl;
				rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfImpedance[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
				ierr = mField.add_ents_to_finite_element_by_TRIs(setOfImpedance[it->get_msId()].tRis,"HELMHOLTZ_IMPEDANCE_FE"); CHKERRQ(ierr);
	
			}
		}
	
	
		PetscFunctionReturn(0);
	}
	
	
	/** \brief this function is used in case of stationary problem to set elements for rhs
	  * \infroup mofem_helmholtz_elem
	  */
	PetscErrorCode setHelmholtzFiniteElementRhs_FOperators(string re_field_name,string im_field_name,Vec &F,bool useScalar) {
		PetscFunctionBegin;
		map<int,BlockData>::iterator sit = setOfBlocks.begin();
		feRhs.get_op_to_do_Rhs().push_back(new OpGetGradReAtGaussPts(re_field_name,commonData));
		for(;sit!=setOfBlocks.end();sit++) {
			//add finite element
			feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs_F(re_field_name,im_field_name,F,sit->second,commonData,useScalar));
		}
		PetscFunctionReturn(0);
	}
	
	PetscErrorCode setHelmholtzFiniteElementRhsOperators_rere(string re_field_name,string im_field_name,Vec &F) {
		PetscFunctionBegin;
		bool imim = false;
		map<int,BlockData>::iterator sit = setOfBlocks.begin();
		
		for(;sit!=setOfBlocks.end();sit++) {
			feRhs.get_op_to_do_Rhs().push_back(new OpGetGradReAtGaussPts(re_field_name,commonData));
			std::cout << "\n commonData.getGradReAtGaussPts(1) = \n" <<  commonData.getGradReAtGaussPts(1) << std::endl;
			feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureReAtGaussPts(re_field_name,commonData));
			std::cout << "\n commonData.pressureReAtGaussPts = \n" <<  commonData.pressureReAtGaussPts << std::endl;
			feRhs.get_op_to_do_Rhs().push_back(new OpGetGradImAtGaussPts(im_field_name,commonData));
			feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureImAtGaussPts(im_field_name,commonData));
			//add finite element
			feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs(re_field_name,re_field_name,F,sit->second,commonData,imim));
			feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs(im_field_name,im_field_name,F,sit->second,commonData,imim));

		}
		
		//imim = true;
		//sit = setOfBlocks.begin();
		//feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(im_field_name,commonData));
		//feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureAtGaussPts(im_field_name,commonData));


		//for(;sit!=setOfBlocks.end();sit++) {
		//	//add finite element
		//	feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs(re_field_name,re_field_name,F,sit->second,commonData,imim));
		//}
		
		PetscFunctionReturn(0);
	}
	
   	//PetscErrorCode setHelmholtzFiniteElementRhsOperators_imim(string re_field_name,string im_field_name,Vec &F) {
	//	bool imim = false;
	//	PetscFunctionBegin;
	//	map<int,BlockData>::iterator sit = setOfBlocks.begin();

	//	for(;sit!=setOfBlocks.end();sit++) {
	//		feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(im_field_name,commonData));
	//		feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureAtGaussPts(im_field_name,commonData)); //get the pressure at gaussian points
	//		//add finite element
	//		feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs(im_field_name,im_field_name,F,sit->second,commonData,imim));
	//	}
	//	
	//	//sit = setOfBlocks.begin();
	//	//feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(re_field_name,commonData));
	//	//feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureAtGaussPts(re_field_name,commonData)); //get the pressure at gaussian points
	//	//for(;sit!=setOfBlocks.end();sit++) {
	//	//	//add finite element
	//	//	feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs(im_field_name,im_field_name,F,sit->second,commonData,imim));
	//	//}
	//	
	//	PetscFunctionReturn(0);
	//}
	
	
	
	
	/** \brief this fucntion is used in case of stationary Helmholtz problem for lhs stiffness term
	  * \infroup mofem_helmholtz_elem
	  */
	PetscErrorCode setHelmholtzFiniteElementLhsOperators(string re_field_name,string im_field_name,Mat A) {
		PetscFunctionBegin;
		map<int,BlockData>::iterator sit = setOfBlocks.begin();
		for(;sit!=setOfBlocks.end();sit++) {
			//add finite elemen
			feLhs.get_op_to_do_Lhs().push_back(new OpHelmholtzLhs_rereA(re_field_name,re_field_name,A,sit->second,commonData));	
			feLhs.get_op_to_do_Lhs().push_back(new OpHelmholtzLhs_imimA(im_field_name,im_field_name,A,sit->second,commonData));

		}
		PetscFunctionReturn(0);
	}
	
	/** \brief this fucntion is used in case of stationary helmholtz problem for lhs Mass term
	* \infroup mofem_helmholtz_elem
	*/
	PetscErrorCode setHelmholtzMassFiniteElementLhsOperators(string re_field_name,string im_field_name,Mat A) {
		PetscFunctionBegin;
		map<int,BlockData>::iterator sit = setOfBlocks.begin();
		for(;sit!=setOfBlocks.end();sit++) {
			//add finite elemen
			feLhs.get_op_to_do_Lhs().push_back(new OpMassLsh(re_field_name,re_field_name,A,sit->second,commonData));
			feLhs.get_op_to_do_Lhs().push_back(new OpMassLsh(im_field_name,im_field_name,A,sit->second,commonData));
		}
		PetscFunctionReturn(0);
	}
	
	/** \brief this function is used in case of statonary problem for helmholtz flux terms
	* \infroup mofem_helmholtz_elem
	*/
	PetscErrorCode setHelmholtzFluxFiniteElementRhsOperators(string re_field_name,string im_field_name,Vec &F,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
		PetscFunctionBegin;
		bool ho_geometry = false;
		//if(mField.check_field(mesh_nodals_positions)) {
		//	ho_geometry = true;
		//}
		map<int,FluxData>::iterator sit = setOfFluxes.begin();
		for(;sit!=setOfFluxes.end();sit++) {
			//add finite element
			feFlux.get_op_to_do_Rhs().push_back(new OpHelmholtzFlux(re_field_name,im_field_name,F,blockData,sit->second,ho_geometry));
		}
		PetscFunctionReturn(0);
	}
	/** \brief this function is used in case of statonary for helmholtz incident wave flux terms
	
	*/
	
	PetscErrorCode setHelmholtzIncidentWaveFiniteElementRhsOperators(string re_field_name,string im_field_name,Vec &F,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
		PetscFunctionBegin;
		bool ho_geometry = false;
		if(mField.check_field(mesh_nodals_positions)) {
		//	ho_geometry = true;
		feIncidentWave.get_op_to_do_Rhs().push_back(new OpHoCoordTri(mesh_nodals_positions,hoCoordsTri));
		}
		
		map<int,FluxData>::iterator sit = setOfFluxes.begin();
		for(;sit!=setOfFluxes.end();sit++) {
			//add finite element
			feIncidentWave.get_op_to_do_Rhs().push_back(new OpHelmholtzIncidentWave(re_field_name,re_field_name,F,blockData,sit->second,hoCoordsTri,ho_geometry));
			feIncidentWave.get_op_to_do_Rhs().push_back(new OpHelmholtzIncidentWave(im_field_name,im_field_name,F,blockData,sit->second,hoCoordsTri,ho_geometry));

		}
		PetscFunctionReturn(0);
	}
	
	/* \brief linear Steady Impedance terms in rhs and lhs
	*/
	//PetscErrorCode setHelmholtzImpedanceFiniteElementRhsOperators(string field_name,Vec &F,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
	//	PetscFunctionBegin;
	//	bool ho_geometry = false;
	//	if(mField.check_field(mesh_nodals_positions)) {
	//		ho_geometry = true;
	//	}
	//	map<int,ImpedanceData>::iterator sit = setOfImpedance.begin();
	//	for(;sit!=setOfImpedance.end();sit++) {
	//		//add finite element
	//		feImpedanceRhs.get_op_to_do_Rhs().push_back(new OpImpedanceRhs(field_name,F,sit->second,ho_geometry));
	//	}
	//	PetscFunctionReturn(0);
	//}
	
	PetscErrorCode setHelmholtzImpedanceFiniteElementLhsOperators(string re_field_name,string im_field_name,Mat A,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
		PetscFunctionBegin;
		bool ho_geometry = false;
		//if(mField.check_field(mesh_nodals_positions)) {
		//	ho_geometry = true;
		//}
		map<int,ImpedanceData>::iterator sit = setOfImpedance.begin();
		for(;sit!=setOfImpedance.end();sit++) {
			//add finite element
			feImpedanceLhs.get_op_to_do_Lhs().push_back(new OpImpedanceLhs_reimC(re_field_name,im_field_name,A,sit->second,blockData,ho_geometry));
			feImpedanceLhs.get_op_to_do_Lhs().push_back(new OpImpedanceLhs_imreC(im_field_name,re_field_name,A,sit->second,blockData,ho_geometry));
		}
		PetscFunctionReturn(0);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	};
	
	}
	
	#endif //__HELMHOLTZ_ELEMENT_HPP
	
	/***************************************************************************//**
	 * \defgroup mofem_helmholtz_elem Helmholtz element
	 * \ingroup mofem_forces_and_sources
	 ******************************************************************************/
