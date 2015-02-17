/** \file L^2NormElement.hpp
 * \brief Operators and data structures for L^2Norm analysis
 *
 * Implementation of L^2Norm element for error analysis
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

#ifndef __L2NORM_ELEMENT_HPP
#define __L2NORM_ELEMENT_HPP

using namespace boost::numeric;
using namespace MoFEM;
#include<moab/Skinner.hpp>
namespace MoFEM {

//ubals::vector<double> analyticalFunction(
//double x,double y,double z) = 0;
	
/** \brief finite element to appeximate analytical solution on surface
  */
struct L2normElement {

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
	
	struct MyTriFE: public TriElementForcesAndSurcesCore {
		MyTriFE(FieldInterface &m_field): TriElementForcesAndSurcesCore(m_field) {}
		int getRule(int order) { return order; };  //why not return ceil(order/2);? wait 
    };
	
	
	FieldInterface &mField;
    L2normElement(
        FieldInterface &m_field):
        feRhs(m_field),feLhs(m_field),mField(m_field) {}
	
	
    //MyTriFE feApprox; 
    //MyTriFE& getLoopFeApprox() { return feApprox; } 
	
    ublas::matrix<double> hoCoords;
	//get the higher order approximation coordinates.
    struct OpHoCoord: public TriElementForcesAndSurcesCore::UserDataOperator {
	
		
		ublas::matrix<double> &hoCoords;
		OpHoCoord(const string field_name,ublas::matrix<double> &ho_coords): 
			TriElementForcesAndSurcesCore::UserDataOperator(field_name),
			hoCoords(ho_coords) {}
	
		/*	
		Cartesian coordinates for gaussian points inside elements
		X^coordinates = DOF dot* N
		*/
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
	
			try {
	
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
	
				hoCoords.resize(data.getN().size1(),3);
				if(type == MBVERTEX) {
					hoCoords.clear();
				}
	
				int nb_dofs = data.getN().size2();
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
					for(int dd = 0;dd<3;dd++) {
						hoCoords(gg,dd) += cblas_ddot(nb_dofs,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
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
	
	/** \brief common data used by volume elements
	* \infroup mofem_helmholtz_elem
	*/
	struct CommonData {
		ublas::vector<double> numValueAtGaussPts;   //Helmholtz Pressure
		ublas::vector<double> numValueRateAtGaussPts;
		ublas::matrix<double> gradAtGaussPts;
		inline ublas::matrix_row<ublas::matrix<double> > getGradAtGaussPts(const int gg) {
			return ublas::matrix_row<ublas::matrix<double> >(gradAtGaussPts,gg);
		}
	};
	CommonData commonData;
	
	/// \brief operator to calculate pressure gradient at Gauss points
	struct OpGetGradAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
		
		CommonData &commonData;
		OpGetGradAtGaussPts(const string field_name,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(field_name),
			commonData(common_data) {}
		
		/** \brief operator calculating pressure gradients
		  *
		  * gradient of numerical solution is calculated multiplying derivatives of shape functions by degrees of freedom.
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			try {
				
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				int nb_dofs = data.getFieldData().size();
				int nb_gauss_pts = data.getN().size1();
				
				//initialize
				commonData.gradAtGaussPts.resize(nb_gauss_pts,3);
				if(type == MBVERTEX) {
					fill(commonData.gradAtGaussPts.data().begin(),commonData.gradAtGaussPts.data().end(),0);
				}
				
				for(int gg = 0;gg<nb_gauss_pts;gg++) {
					ublas::noalias(commonData.getGradAtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
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
	//struct OpGetFieldAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
	
	//	ublas::vector<double> &fieldAtGaussPts;
	//	OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
	//		TetElementForcesAndSourcesCore::UserDataOperator(field_name),
	//		fieldAtGaussPts(field_at_gauss_pts) {}
	
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
	struct OpGetTetNumericalValueAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetNumericalValueAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field_name,common_data.numValueAtGaussPts) {}
	};
	
	/** \brief operator to calculate pressure at Gauss pts
	* \infroup mofem_helmholtz_elem
	*/
	struct OpGetTriNumericalValueAtGaussPts: public OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore> {
		OpGetTriNumericalValueAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore>(field_name,common_data.numValueAtGaussPts) {}
	};
	
	
	/** \brief Rhs operaetar used to build matrix
	*/
	template<typename OP>  //Calculate L^2 norm error in Tets or surface triangles
	struct OpRhsError:public OP::UserDataOperator {
		
		CommonData &commonData;
		Vec F;
		ublas::matrix<double> &hoCoords;
		double (*fUN)(double x,double y,double z);
		bool useTsF;
		OpRhsError(const string re_field_name,const string im_field_name,ublas::matrix<double> &ho_coords,
				   CommonData &common_data,   
				   double (*fun)(double x,double y,double z)): 
			OP::UserDataOperator(re_field_name,im_field_name),
			hoCoords(ho_coords),fUN(fun),commonData(common_data),useTsF(true) {}
		
		OpRhsError(const string re_field_name,const string im_field_name,ublas::matrix<double> &ho_coords,
			  Vec _F,CommonData &common_data,   
			  double (*fun)(double x,double y,double z)): 
			OP::UserDataOperator(re_field_name,im_field_name),
			hoCoords(ho_coords),fUN(fun),commonData(common_data),useTsF(false),F(_F) {}
		
		ublas::vector<FieldData> UmUh;
		//struct OpRhs:public TriElementForcesAndSurcesCore::UserDataOperator {

		//  ublas::matrix<double> &hoCoords;
		//  double (*fUN)(double x,double y,double z);
		//  OpRhs(const string re_field_name,const string im_field_name,ublas::matrix<double> &ho_coords,
		//double (*fun)(double x,double y,double z)): 
		//TriElementForcesAndSurcesCore::UserDataOperator(re_field_name,im_field_name),
		//hoCoords(ho_coords),fUN(fun)  {}

		//  ublas::vector<FieldData> NTf;
		
		/*	
		Rhs force vector merely with Dirichlet values
		F = int_S N^T Fun dS
		*/
		
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			if (TypeIsTet<OP>::useTet) {
			PetscFunctionBegin;
			PetscErrorCode ierr;
			
			try {

				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				//if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);
				
				unsigned int nb_row = data.getN().size2();
				if(nb_row != data.getIndices().size()) {
					SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
							"currently works only for scalar fields, extension to fields with higher rank need to be implemented");
				}
				UmUh.resize(nb_row);

				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

					double x,y,z;

					
					//Integrate over volume
					double val = this->getVolume()*this->getGaussPts()(3,gg); //or OP::UserDataOperator::getVolume()
					
					if(this->getHoGaussPtsDetJac().size()>0) {
						val *= this->getHoGaussPtsDetJac()[gg]; ///< higher order geometry
						x = hoCoords(gg,0);
						y = hoCoords(gg,1);
						z = hoCoords(gg,2);
					} else{
						x = this->getCoordsAtGaussPts()(gg,0);  //Coordinates of Gaussian points for higher order geometry
						y = this->getCoordsAtGaussPts()(gg,1);
						z = this->getCoordsAtGaussPts()(gg,2);
					}
						
					
					double a = fUN(x,y,z);
					double uExact = pow(a,2.0);
					double uAnalyt = pow(commonData.numValueAtGaussPts[gg],2.0);
					
					//for (i=0; i<=arrayset; i++){
					//	a[i]= pi*(pow(static_cast<float>(g[i]),2));
					//}
					
                   //#include <algorithm>
                   //		using namespace std;
                   //		
					//static inline double computeSquare (double x) { return x*x; }
					
					//...
					
					//std::transform(Array1.begin(), Array1.end(), Array1.begin(), (double(*)(double)) sqrt);
					//std::transform(Array2.begin(), Array2.end(), Array2.begin(), computeSquare);
					bool relativeError = false;
					if (relativeError){
						
					ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val)/sqrt(uExact*val);
					} else {
					ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val); }
				}
	
				}


			 catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			 } 
			 PetscFunctionReturn(0);
			 
			} else if (!TypeIsTet<OP>::useTet) {
			
				PetscFunctionBegin;
				
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				//if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
				
				PetscErrorCode ierr;
				
				const FENumeredDofMoFEMEntity *dof_ptr;
				ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
				int rank = dof_ptr->get_max_rank();
				
				int nb_dofs = data.getIndices().size()/rank;
				UmUh.resize(data.getIndices().size();
				
				//bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));
				fill(Nf.begin(),Nf.end(),0);
				
				//cerr << getNormal() << endl;
				//cerr << getNormals_at_GaussPt() << endl;
				
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
					
					
					double x,y,z;
					//Integrate over surface
					double val = getGaussPts()(2,gg);
					if(hoCoords.size1() == data.getN().size1()) {
						double area = norm_2(this->getNormals_at_GaussPt(gg)); 
						val *= area;
						x = hoCoords(gg,0);
						y = hoCoords(gg,1);
						z = hoCoords(gg,2);
					} else {
						val *= this->getArea();
						x = this->getCoordsAtGaussPts()(gg,0);  //Coordinates of Gaussian points for higher order geometry
						y = this->getCoordsAtGaussPts()(gg,1);
						z = this->getCoordsAtGaussPts()(gg,2);
					}
			
					
					double a = fUN(x,y,z);
					double uExact = pow(a,2.0);
					double uAnalyt = pow(commonData.numValueAtGaussPts[gg],2.0);
					
					
					//cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
					bool relativeError = false;
					if (relativeError){
						
						ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val)/sqrt(uExact*val);
					} else {
					ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val) }
					//ublas::noalias(Nf) += val*flux*data.getN(gg,nb_dofs);
					
				}
										
				if(useTsF) {
				ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
									&data.getIndices()[0],&*UmUh.data().begin(),ADD_VALUES); CHKERRQ(ierr);	
				} else {
					ierr = VecSetValues(F,data.getIndices().size(),
										&data.getIndices()[0],&*UmUh.data().begin(),ADD_VALUES); CHKERRQ(ierr);	
				}
				
				//ierr = VecSetValues(getFEMethod()->snes_f,data.getIndices().size(),
				//					&data.getIndices()[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);


			
				catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
					SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
				}
				
				PetscFunctionReturn(0);
			
		
    };
	
}
	};
	
	
	template<typename OP,typename OS = void>
	struct TypeIsTet
	{
		static const bool useTet = false;
	};
	
	template<typename OS>
	struct TypeIsTet<TetElementForcesAndSourcesCore,OS>
	{
		static const bool useTet = true;
	};

/** \brief operator to calculate tempereature at Gauss pts
  * \infroup mofem_thermal_elem
  */
struct OpTetRhs: public OpRhsError<TetElementForcesAndSourcesCore> {
    OpTetRhs(const string re_field_name,const string im_field_name,ublas::matrix<double> &ho_coords,
		Vec _F,CommonData &common_data,   
		double (*fun)(double x,double y,double z)): 
		OpRhsError<TetElementForcesAndSourcesCore>(re_field_name,im_field_name,ho_coords,
				_F,common_data,fun) {}

};

/** \brief operator to calculate tempereature at Gauss pts
  * \infroup mofem_thermal_elem
  */
struct OpTriRhs: public OpRhsError<TriElementForcesAndSurcesCore> {
    OpTriRhs(const string re_field_name,const string im_field_name,ublas::matrix<double> &ho_coords,
		Vec _F,CommonData &common_data,   
		double (*fun)(double x,double y,double z)): 
		OpRhsError<TriElementForcesAndSurcesCore>(re_field_name,im_field_name,ho_coords,
				_F,common_data,fun) {}

};
///** \brief operator to calculate tempereature at Gauss pts
//  * \infroup mofem_thermal_elem
//  */
//struct OpGetTriTemperatureAtGaussPts: public OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore> {
//    OpGetTriTemperatureAtGaussPts(const string field_name,CommonData &common_data):
//		OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore>(field_name,common_data.temperatureAtGaussPts) {}
//};


/*
  Add the L2norm element with same problem and same field as the original problem
 
 */
PetscErrorCode addL2NormElements(
	const string problem_name,string fe,const string re_field_name,const string im_field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
	PetscFunctionBegin;
	PetscErrorCode ierr;
	ErrorCode rval;
	ierr = mfield.add_finite_element(fe,MF_ZERO); CHKERRQ(ierr);
    ierr = mfield.modify_finite_element_add_field_row(fe,re_field_name); CHKERRQ(ierr);
    ierr = mfield.modify_finite_element_add_field_col(fe,re_field_name); CHKERRQ(ierr);
    ierr = mfield.modify_finite_element_add_field_data(fe,re_field_name); CHKERRQ(ierr);
	
	ierr = mfield.modify_finite_element_add_field_row(fe,im_field_name); CHKERRQ(ierr);
    ierr = mfield.modify_finite_element_add_field_col(fe,im_field_name); CHKERRQ(ierr);
    ierr = mfield.modify_finite_element_add_field_data(fe,im_field_name); CHKERRQ(ierr);
	
	if(m_field.check_field(nodals_positions)) {
		ierr = m_field.modify_finite_element_add_field_data(fe,nodals_positions); CHKERRQ(ierr);
    }
	ierr = m_field.modify_problem_add_finite_element(problem,fe); CHKERRQ(ierr);
	
	
	Range tRis;
	for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"L2NROM_ERROR",it)) {
		rval = moab.get_entities_by_type(it->get_meshset(),MBTRI,tRis,true); CHKERR_PETSC(rval);
		ierr = m_field.add_ents_to_finite_element_by_TRIs(tRis,fe); CHKERRQ(ierr);
	}
	
    
    PetscFunctionReturn(0);
	
	}
	
	
	
	//PetscErrorCode setHelmholtzFiniteElementRhsOperators_rere(string re_field_name,string im_field_name,Vec &F) {
	//	PetscFunctionBegin;
	//	map<int,BlockData>::iterator sit = setOfBlocks.begin();
	//	feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(re_field_name,commonData));
	//	feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureAtGaussPts(re_field_name,commonData));
	//	for(;sit!=setOfBlocks.end();sit++) {
	//		//add finite element
	//		feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs(re_field_name,re_field_name,F,sit->second,commonData));
	//	}
	//	PetscFunctionReturn(0);
	//}
	
   	//PetscErrorCode setHelmholtzFiniteElementRhsOperators_imim(string re_field_name,string im_field_name,Vec &F) {
	//	PetscFunctionBegin;
	//	map<int,BlockData>::iterator sit = setOfBlocks.begin();
	//	feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(im_field_name,commonData));
	//	feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureAtGaussPts(im_field_name,commonData)); //get the pressure at gaussian points
	//	for(;sit!=setOfBlocks.end();sit++) {
	//		//add finite element
	//		feRhs.get_op_to_do_Rhs().push_back(new OpHelmholtzRhs(im_field_name,im_field_name,F,sit->second,commonData));
	//	}
	//	PetscFunctionReturn(0);
	//}
	
	


Range tRis;


PetscErrorCode setL2NormRelativeError(FieldInterface &m_field,string re_field_name,string im_field_name,Vec &F,double (*fun)(double x,double y,double z),
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    if(m_field.check_field(nodals_positions)) {
		feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(re_field_name,commonData));
		feRhs.get_op_to_do_Rhs().push_back(new OpGetTetPressureAtGaussPts(re_field_name,commonData));
		feRhs.get_op_to_do_Rhs().push_back(new OpHoCoord(nodals_positions,hoCoords));
    }
    feRhs.get_op_to_do_Rhs().push_back(new OpTetRhs(re_field_name,im_field_name,hoCoords,F,commonData,fun));
    //approxField.getLoopFeApprox().get_op_to_do_Lhs().push_back(new ApproxField::OpLhs(re_field_name,approxField.hoCoords));
	PetscFunctionReturn(0);
}





};




}


#endif //__L2NORM_ELEMENT_HPP


/***************************************************************************//**
 * \defgroup mofem_L2Norm_elem L2Norm element
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/





	
