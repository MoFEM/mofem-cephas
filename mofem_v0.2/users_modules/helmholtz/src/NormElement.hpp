/** \file NormElement.hpp
 * \brief Operators and data structures for L^2Norm analysis
 *
 * Implementation of L^2 and H_1 Norm element for error analysis
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

#ifndef __NORM_ELEMENT_HPP
#define __NORM_ELEMENT_HPP

using namespace boost::numeric;
using namespace MoFEM;
#include<moab/Skinner.hpp>
namespace MoFEM {

//calculate the norm of error for scalar field
	
/** \brief finite element to appeximate analytical solution on surface
  */
struct NormElement {

	struct MyVolumeFE: public TetElementForcesAndSourcesCore {
		MyVolumeFE(FieldInterface &mField): TetElementForcesAndSourcesCore(mField) {}
		
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
		int getRule(int order) { return order; };
	};
	MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
	MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element
	MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
	MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element
	
	struct MyTriFE: public TriElementForcesAndSurcesCore {
		MyTriFE(FieldInterface &m_field): TriElementForcesAndSurcesCore(m_field) {}
		int getRule(int order) { return ceil(order/2); };  
    };
	
	
	FieldInterface &m_field;
    NormElement(
        FieldInterface &mField):
        feRhs(mField),feLhs(mField),m_field(mField) {}
	
	//Field data
	struct CommonData {
		//double fIeld1; // value of scalar in field 1
		//double fIeld2;  // etc.
		
		ublas::vector<double> fieldValue1AtGaussPts;
		ublas::vector<double> fieldValue1RateAtGaussPts;
		ublas::matrix<double> gradField1AtGaussPts;
		
		ublas::vector<double> fieldValue2AtGaussPts;
		ublas::vector<double> fieldValue2RateAtGaussPts;
		ublas::matrix<double> gradField2AtGaussPts;
		inline ublas::matrix_row<ublas::matrix<double> > getGradField1AtGaussPts(const int gg) {
			return ublas::matrix_row<ublas::matrix<double> >(gradField1AtGaussPts,gg);
		};
	
		inline ublas::matrix_row<ublas::matrix<double> > getGradField2AtGaussPts(const int gg) {
			return ublas::matrix_row<ublas::matrix<double> >(gradField2AtGaussPts,gg);		
		};
	
	};
	CommonData commonData;
	
	struct BlockData {
		Range tEts; ///< constatins elements in block set
	};
	map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData
	//BlockData blockData;
	
	//map<int,CommonData> setOfCommons; ///< maps block set id with appropiate BlockData
	
	
	/// \brief operator to calculate temperature gradient at Gauss points
	struct OpGetGradField1AtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
	
		CommonData &commonData;
		OpGetGradField1AtGaussPts(const string field_name,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(field_name),
			commonData(common_data) {}
	
		/** \brief operator calculating temperature gradients
		  *
		  * temerature gradient is calculated multiplying derivatives of shape functions by degrees of freedom.
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			try {
	
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				int nb_dofs = data.getFieldData().size();
				int nb_gauss_pts = data.getN().size1();
	
				//initialize
				commonData.gradField1AtGaussPts.resize(nb_gauss_pts,3);
				if(type == MBVERTEX) {
					fill(commonData.gradField1AtGaussPts.data().begin(),commonData.gradField1AtGaussPts.data().end(),0);
				}
	
				for(int gg = 0;gg<nb_gauss_pts;gg++) {
					ublas::noalias(commonData.getGradField1AtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
				}
	
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
	
			PetscFunctionReturn(0);
		}
	
	};
	
	/// \brief operator to calculate temperature gradient at Gauss points
	struct OpGetGradField2AtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
		
		CommonData &commonData;
		OpGetGradField2AtGaussPts(const string field_name,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(field_name),
			commonData(common_data) {}
		
		/** \brief operator calculating temperature gradients
		  *
		  * temerature gradient is calculated multiplying derivatives of shape functions by degrees of freedom.
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			try {
				
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				int nb_dofs = data.getFieldData().size();
				int nb_gauss_pts = data.getN().size1();
				
				//initialize
				commonData.gradField2AtGaussPts.resize(nb_gauss_pts,3);
				if(type == MBVERTEX) {
					fill(commonData.gradField2AtGaussPts.data().begin(),commonData.gradField2AtGaussPts.data().end(),0);
				}
				
				for(int gg = 0;gg<nb_gauss_pts;gg++) {
					ublas::noalias(commonData.getGradField2AtGaussPts(gg)) += prod( trans(data.getDiffN(gg,nb_dofs)), data.getFieldData() );
				}
				
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
			
			PetscFunctionReturn(0);
		}
		
	};
	
	/** \brief opearator to caulate tempereature  and rate of temperature at Gauss points
    * \infroup mofem_thermal_elem
    */
	template<typename OP>
	struct OpGetFieldAtGaussPts: public OP::UserDataOperator {
	
		ublas::vector<double> &fieldAtGaussPts;
		OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
			OP::UserDataOperator(field_name),
			fieldAtGaussPts(field_at_gauss_pts) {}
	
		/** \brief operator calculating temperature and rate of temperature
		  *
		  * temperature temperature or rate of temperature is calculated multiplyingshape functions by degrees of freedom
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
				//std::cout << "\n data.getFieldData() = \n" << data.getFieldData() << std::endl;
				//std::cout << "\n nb_dofs = \n" << nb_dofs << std::endl;
				//std::cout << "\n nb_gauss_pts = \n" << nb_gauss_pts << std::endl;
				for(int gg = 0;gg<nb_gauss_pts;gg++) {
					fieldAtGaussPts[gg] += inner_prod(data.getN(gg,nb_dofs),data.getFieldData());
					//std::cout << "\n fieldAtGaussPts[gg] \n" << fieldAtGaussPts[gg] << std::endl;
					
				}
	
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
			}
	
			PetscFunctionReturn(0);
		}
	
	};
	
	
	/** \brief operator to calculate tempereature at Gauss pts
    * \infroup mofem_thermal_elem
    */
	struct OpGetTetField1AtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetField1AtGaussPts(const string field1_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field1_name,common_data.fieldValue1AtGaussPts) {}
	};
	
	struct OpGetTetField2AtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetField2AtGaussPts(const string field2_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field2_name,common_data.fieldValue2AtGaussPts) {}
	};
	
	
	/** \brief operator to calculate tempereature at Gauss ptss
	  * \infroup mofem_thermal_elem
	  */
	struct OpGetTriFieldAtGaussPts: public OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore> {
		OpGetTriFieldAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore>(field_name,common_data.fieldValue1AtGaussPts) {}
	};
	
	/** \brief operator to calculate temperature rate at Gauss pts
	  * \infroup mofem_thermal_elem
	  */
	struct OpGetTetFieldRateAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetFieldRateAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field_name,common_data.fieldValue1RateAtGaussPts) {}
	};
	
	
	
	/** \brief Rhs operaetar used loop differences between two fields
	*/
	template<typename OP>  //Calculate L^2 norm error in Tets or surface triangles
	struct OpRhsError:public OP::UserDataOperator {
		
		CommonData &commonData;
		Vec F;
		bool usel2;
        //ublas::vector<double> fieldValue1AtGaussPts;
		//ublas::vector<double> fieldValue2AtGaussPts;
		//ublas::<double> gradFieldValue1AtGaussPts;
		//ublas::matrix<double> gradFieldValue2AtGaussPts;
		ublas::vector<double> Nf;
		
		bool useTsF;
		OpRhsError(const string re_field_name,const string im_field_name,
				   CommonData &common_data,bool useL2
				   ): 
			OP::UserDataOperator(re_field_name,im_field_name),
			commonData(common_data),usel2(useL2),useTsF(true) {}
		
		OpRhsError(const string re_field_name,const string im_field_name,
			  Vec _F,CommonData &common_data,bool useL2
			  ): 
			OP::UserDataOperator(re_field_name,im_field_name),
			commonData(common_data),usel2(useL2),useTsF(false),F(_F) {}
		
		ublas::vector<double> eRror;

		
		/*	
		Rhs force vector merely with Dirichlet values
		F = sqrt[ int_S (F1 - F2)^2 dS ]
		*/
		
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			if (TypeIsTet<OP>::useTet) {
			PetscFunctionBegin;
			PetscErrorCode ierr;
			
			try {
				std::cout << "\n commonData.fieldValue1AtGaussPts =\n " << commonData.fieldValue1AtGaussPts << std::endl;

				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				//if(dAta.tEts.find(getMoFEMFEPtr()->get_ent())==dAta.tEts.end()) PetscFunctionReturn(0);
				
				unsigned int nb_row = data.getN().size2();
				if(nb_row != data.getIndices().size()) {
					SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
							"currently works only for scalar fields, extension to fields with higher rank need to be implemented");
				}
				//resize error in specific elements on each vertex equal to dofs
				eRror.resize(nb_row);
				Nf.resize(nb_row);
				cout << "\n nb_row = \n " << nb_row << std::endl;
				Nf.clear();
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
					
					//Integrate over volume
					double val = this->getVolume()*this->getGaussPts()(3,gg); //or OP::UserDataOperator::getVolume()
					
					if(this->getHoGaussPtsDetJac().size()>0) {
						val *= this->getHoGaussPtsDetJac()[gg]; ///< higher order geometry

					} else{
					}
						
					ublas::vector<double> uAnaly = commonData.fieldValue1AtGaussPts;
					ublas::vector<double> uNumer = commonData.fieldValue2AtGaussPts;
					ublas::matrix_row< ublas::matrix<double> > uAnalyGrad = commonData.getGradField1AtGaussPts(gg);
					ublas::matrix_row< ublas::matrix<double> > uNumerGrad = commonData.getGradField2AtGaussPts(gg);
					
                   //std::cout << "\n commonData.fieldValue1AtGaussPts =\n " << commonData.fieldValue1AtGaussPts << std::endl;
				   //std::cout << "\n commonData.getGradField2AtGaussPts(gg) = \n" << commonData.getGradField2AtGaussPts(gg) << std::endl;
					ublas::vector<double> sqError;
					sqError.resize(1); //resize vector to scalar.
					if(usel2) {
					eRror = uAnaly - uNumer;
					sqError[0] = pow(eRror(gg),2.0);
					//for (i=0; i<=arrayset; i++){
					//	a[i]= pi*(pow(static_cast<float>(g[i]),2));
					//}
					
                   //#include <algorithm>
                   //		using namespace std;
                   //		
					//static inline double computeSquare (double x) { return x*x; }
					
					//std::transform(Array1.begin(), Array1.end(), Array1.begin(), (double(*)(double)) sqrt);
					//std::transform(Array2.begin(), Array2.end(), Array2.begin(), computeSquare);
					} else{
					eRror = uAnaly - uNumer;
					//double sqGradError
					double sqGradError = ublas::inner_prod((commonData.getGradField1AtGaussPts(gg)-commonData.getGradField2AtGaussPts(gg)),(commonData.getGradField1AtGaussPts(gg)-commonData.getGradField2AtGaussPts(gg)));
					sqError[0] = sqGradError + pow(eRror(gg),2.0);
					}
					
					ublas::noalias(Nf) += val*sqError;
					
					
				//	bool relativeError = false;
				//	if (relativeError){
				//		
				//	ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val)/sqrt(uExact*val);
				//	} else {
				//	ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val); }
			    }
				 
				if(useTsF) {
					//ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
					//					&data.getIndices()[0],&*Nf.data().begin(),ADD_VALUES); CHKERRQ(ierr);	
					return 0;
				} else {
					ierr = VecSetValues(F,data.getIndices().size(),
										&data.getIndices()[0],&*Nf.data().begin(),ADD_VALUES); CHKERRQ(ierr);	
				}
	
				}


			 catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			 } 
			 PetscFunctionReturn(0);
			 
			} else if (!TypeIsTet<OP>::useTet) {
			
				//PetscFunctionBegin;
				
				//if(data.getIndices().size()==0) PetscFunctionReturn(0);
				////if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
				
				//PetscErrorCode ierr;
				
				//const FENumeredDofMoFEMEntity *dof_ptr;
				//ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
				//int rank = dof_ptr->get_max_rank();
				
				//int nb_dofs = data.getIndices().size()/rank;
				////UmUh.resize(data.getIndices().size();
				
				////bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));
				//fill(Nf.begin(),Nf.end(),0);
				
				////cerr << getNormal() << endl;
				////cerr << getNormals_at_GaussPt() << endl;
				
				//for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
				//	
				//	//////NEED TO BE FIXED!!!!
				//	double x,y,z;
				//	//Integrate over surface
				//	double val = getGaussPts()(2,gg);
				//	if(hoCoords.size1() == data.getN().size1()) {
				//		double area = norm_2(this->getNormals_at_GaussPt(gg)); 
				//		val *= area;
				//		x = hoCoords(gg,0);
				//		y = hoCoords(gg,1);
				//		z = hoCoords(gg,2);
				//	} else {
				//		val *= this->getArea();
				//		x = this->getCoordsAtGaussPts()(gg,0);  //Coordinates of Gaussian points for higher order geometry
				//		y = this->getCoordsAtGaussPts()(gg,1);
				//		z = this->getCoordsAtGaussPts()(gg,2);
				//	}
			
				//	
				//	double a = fUN(x,y,z);
				//	double uExact = pow(a,2.0);
				//	double uAnalyt = pow(commonData.numValueAtGaussPts[gg],2.0);
				//	
				//	
				//	//cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);
				//	bool relativeError = false;
				//	if (relativeError){
				//		
				//		ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val)/sqrt(uExact*val);
				//	} else {
				//	ublas::noalias(UmUh) += sqrt((uExact-uAnalyt)*val) }
				//	//ublas::noalias(Nf) += val*flux*data.getN(gg,nb_dofs);
				//	
				//}
										
				//if(useTsF) {
				//ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),
				//					&data.getIndices()[0],&*UmUh.data().begin(),ADD_VALUES); CHKERRQ(ierr);	
				//} else {
				//	ierr = VecSetValues(F,data.getIndices().size(),
				//						&data.getIndices()[0],&*UmUh.data().begin(),ADD_VALUES); CHKERRQ(ierr);	
				//}
				
				//ierr = VecSetValues(getFEMethod()->snes_f,data.getIndices().size(),
				//					&data.getIndices()[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);


			
				//catch (const std::exception& ex) {
				//ostringstream ss;
				//ss << "throw in method: " << ex.what() << endl;
				//	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
				//}
				
				//PetscFunctionReturn(0);
			return 0;
		
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

/** \brief operator to calculate error at Gauss pts
  * \infroup mofem_thermal_elem
  */
struct OpTetRhs: public OpRhsError<TetElementForcesAndSourcesCore> {
    OpTetRhs(const string re_field_name,const string im_field_name,
		Vec _F,CommonData &common_data,   
		bool useL2): 
		OpRhsError<TetElementForcesAndSourcesCore>(re_field_name,im_field_name,
				_F,common_data,useL2) {}

};



/*
  Add the error norm element with same problem and same field as the original problem
 
 */
PetscErrorCode addNormElements(
	const string problem,string fe,const string norm_field_name,const string field1_name,
	const string field2_name,
	const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
	PetscFunctionBegin;
	PetscErrorCode ierr;
	ErrorCode rval;
	ierr = m_field.add_finite_element(fe,MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row(fe,norm_field_name); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(fe,norm_field_name); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(fe,norm_field_name); CHKERRQ(ierr);
	
    ierr = m_field.modify_finite_element_add_field_data(fe,field1_name); CHKERRQ(ierr);
	
    ierr = m_field.modify_finite_element_add_field_data(fe,field2_name); CHKERRQ(ierr);
	
	
	if(m_field.check_field(mesh_nodals_positions)) {
		ierr = m_field.modify_finite_element_add_field_data(fe,mesh_nodals_positions); CHKERRQ(ierr);
    }
	ierr = m_field.modify_problem_add_finite_element(problem,fe); CHKERRQ(ierr);
	
	
	//Range tEts;
	for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"MAT_NORM",it)) {
		//rval = m_Field.get_moab().get_entities_by_type(it->get_meshset(),MBTET,tEts,true); CHKERR_PETSC(rval);
		rval = m_field.get_moab().get_entities_by_type(it->get_meshset(),MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);

		ierr = m_field.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,fe); CHKERRQ(ierr);
	}
	
    
    PetscFunctionReturn(0);
	
	}

//Range tEts;


PetscErrorCode setNormFiniteElementRhsOperator(string norm_field_name,string field1_name,
	string field2_name,Vec &F,bool usel2,
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
	//ublas::vector<double> field_Value1AtGaussPts;
	//ublas::vector<double> field_Value2AtGaussPts;
	map<int,BlockData>::iterator sit = setOfBlocks.begin();
	
	for(;sit!=setOfBlocks.end();sit++) {
		
		//Calculate field values at gaussian points for field1 and field2; 
		feRhs.get_op_to_do_Rhs().push_back(new OpGetTetField1AtGaussPts(field1_name,commonData));
		std::cout << "\n commonData.fieldValue1AtGaussPts =\n " << commonData.fieldValue1AtGaussPts << std::endl;
		feRhs.get_op_to_do_Rhs().push_back(new OpGetTetField2AtGaussPts(field2_name,commonData));
		feRhs.get_op_to_do_Rhs().push_back(new OpGetGradField1AtGaussPts(field1_name,commonData));
		feRhs.get_op_to_do_Rhs().push_back(new OpGetGradField2AtGaussPts(field2_name,commonData));
		
    feRhs.get_op_to_do_Rhs().push_back(new OpTetRhs(norm_field_name,norm_field_name,F,commonData,usel2));
}

	PetscFunctionReturn(0);



}




};

}


#endif //__NORM_ELEMENT_HPP


/***************************************************************************//**
 * \defgroup mofem_Norm_elem Norm element
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/





	
