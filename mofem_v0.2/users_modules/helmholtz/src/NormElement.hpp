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
		//two different fields as inputs		
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
	
	// \brief operator to calculate field gradient at Gauss points
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
				if(data.getFieldData().size()==0) PetscFunctionReturn(0);
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
	
	/// \brief operator to calculate field gradient at Gauss points
	struct OpGetGradField2AtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {
		
		CommonData &commonData;
		OpGetGradField2AtGaussPts(const string field_name,CommonData &common_data):
			TetElementForcesAndSourcesCore::UserDataOperator(field_name),
			commonData(common_data) {}
		
		/** \brief operator calculating field gradients
		  *
		  * field gradient is calculated multiplying derivatives of shape functions by degrees of freedom.
		  */
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			try {
				
				if(data.getFieldData().size()==0) PetscFunctionReturn(0);
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
	
	/** \brief opearator to caulate field value and rate of field value at Gauss points
    * \infroup mofem_thermal_elem
    */
	template<typename OP>
	struct OpGetFieldAtGaussPts: public OP::UserDataOperator {
	
		ublas::vector<double> &fieldAtGaussPts;
		OpGetFieldAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
			OP::UserDataOperator(field_name),
			fieldAtGaussPts(field_at_gauss_pts) {}
	
		/** \brief operator calculating field and rate of field
		  *
		  * field value or rate of field value is calculated multiplyingshape functions by degrees of freedom
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
	
	
	/** \brief operator to calculate field at Gauss pts
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
	
	
	/** \brief operator to calculate field at Gauss ptss
	  * \infroup mofem_thermal_elem
	  */
	struct OpGetTriFieldAtGaussPts: public OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore> {
		OpGetTriFieldAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TriElementForcesAndSurcesCore>(field_name,common_data.fieldValue1AtGaussPts) {}
	};
	
	/** \brief operator to calculate field rate at Gauss pts
	  * \infroup mofem_thermal_elem
	  */
	struct OpGetTetFieldRateAtGaussPts: public OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore> {
		OpGetTetFieldRateAtGaussPts(const string field_name,CommonData &common_data):
			OpGetFieldAtGaussPts<TetElementForcesAndSourcesCore>(field_name,common_data.fieldValue1RateAtGaussPts) {}
	};
	
	
	/** \brief Lhs operaetar for tetrahedral used to build matrix
	*/
    struct OpLhs:public TetElementForcesAndSourcesCore::UserDataOperator {
		
		Mat A;
		//bool solveBc;
		OpLhs(const string re_field_name,Mat _A): 
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name),
			A(_A) { }
		
		OpLhs(const string re_field_name): 
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name) { }
		
		ublas::matrix<FieldData> NTN,transNTN;
		
		/*	
		Lhs mass matrix
		A = N^T N
		*/
		
		PetscErrorCode doWork(
			int row_side,int col_side,
			EntityType row_type,EntityType col_type,
			DataForcesAndSurcesCore::EntData &row_data,
			DataForcesAndSurcesCore::EntData &col_data) {
			PetscFunctionBegin;
			
			PetscErrorCode ierr;
			
			try {
				
				if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
				if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
				
				int nb_row = row_data.getIndices().size();
				int nb_col = col_data.getIndices().size();
				
				if(nb_row != row_data.getIndices().size()) {
					SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
							"currently works only for scalar fields, extension to fields with higher rank need to be implemented");
				}
				
				NTN.resize(nb_row,nb_col);
				NTN.clear();
				
				for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
					double val = getVolume()*getGaussPts()(3,gg);
	
					//if(hoCoords.size1() == row_data.getN().size1()) {
					if(this->getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
					} 
					

					cblas_dger(CblasRowMajor,nb_row,nb_col,val,
							   &row_data.getN(gg)[0],1,&col_data.getN(gg)[0],1,
							   &NTN(0,0),nb_col);
					
						
				}
					
					ierr = MatSetValues(
							   A,
							   nb_row,&row_data.getIndices()[0],
							   nb_col,&col_data.getIndices()[0],
							   &NTN(0,0),ADD_VALUES); CHKERRQ(ierr);
					if(row_side != col_side || row_type != col_type) {
						transNTN.resize(nb_col,nb_row);
						noalias(transNTN) = trans(NTN);
						ierr = MatSetValues(
								   A,
								   nb_col,&col_data.getIndices()[0],
								   nb_row,&row_data.getIndices()[0],
								   &transNTN(0,0),ADD_VALUES); CHKERRQ(ierr);
					}
					
					
					
				
			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}
			
			PetscFunctionReturn(0);
		}
    };
	
	
	/** \brief Rhs operaetar used loop differences between two fields
	*/
	struct OpRhs:public TetElementForcesAndSourcesCore::UserDataOperator {
		
		CommonData &commonData;
		Vec F;//norm error
		//Vec D;//relative error
		bool useL2;
		bool useTsF;
		bool useRela;//use relative error
		ublas::vector<double> Nf;
		ublas::vector<double> rElative_error;
		OpRhs(const string re_field_name,const string im_field_name,
				   CommonData &common_data,bool usel2,bool userela
				   ): 
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			commonData(common_data),useL2(usel2),useTsF(true),useRela(userela) {}
		
		OpRhs(const string re_field_name,const string im_field_name,
			  Vec _F,CommonData &common_data,bool usel2,bool userela
			  ): 
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name,im_field_name),
			commonData(common_data),useL2(usel2),useTsF(false),useRela(userela),F(_F) {}
		
		
		/*	
		Rhs force vector merely with field values
		F = sqrt[ int_S (F1 - F2)^2 dS ]
		*/
		
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			
			PetscFunctionBegin;
			PetscErrorCode ierr;
			
			try {
				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				
				
				
				unsigned int nb_row = data.getIndices().size();
				if(nb_row != data.getIndices().size()) {
					SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
							"currently works only for scalar fields, extension to fields with higher rank need to be implemented");
				}
				//resize error in specific elements on each vertex equal to dofs
				Nf.resize(nb_row);
				rElative_error.resize(nb_row);
				Nf.clear();
				rElative_error.clear();
				ublas::vector<double> uAnaly = commonData.fieldValue1AtGaussPts;
				ublas::vector<double> uNumer = commonData.fieldValue2AtGaussPts;
				
				double eRror;
				double sqError;
				
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
					
					//Integrate over volume
					double val = getVolume()*getGaussPts()(3,gg);//this->getGaussPts()(3,gg); 
					
					if(this->getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry

					} else{
					}
					ublas::matrix_row< ublas::matrix<double> > uAnalyGrad = commonData.getGradField1AtGaussPts(gg);
					ublas::matrix_row< ublas::matrix<double> > uNumerGrad = commonData.getGradField2AtGaussPts(gg);
		

					if(useL2) { //case L2 norm
						
						double aa = abs(uAnaly(gg));
						double bb = abs(uNumer(gg));
						eRror = aa - bb;
						//eRror = uAnaly(gg) - uNumer(gg);
						sqError = pow(eRror,2.0);

					} else if(!useL2) { //case H1 norm
					
						double aa = uAnaly(gg);
						double bb = uNumer(gg);
						eRror = aa - bb;
						double sqGradError = ublas::inner_prod((commonData.getGradField1AtGaussPts(gg)-commonData.getGradField2AtGaussPts(gg)),(commonData.getGradField1AtGaussPts(gg)-commonData.getGradField2AtGaussPts(gg)));
					
						sqError = sqGradError + pow(eRror,2.0);
						
					}
					//need to calculate sqrt of norm^2
					if(!useRela) { //case Norm error
						
						ublas::noalias(Nf) += val*sqError*data.getN(gg);
					
					} else if(useRela) { //case relative error
						
						double sqUanaly = pow(norm_inf(uAnaly),2.0);
				
						ublas::noalias(rElative_error) += val*(pow(eRror,2.0)/sqUanaly)*data.getN(gg);

					}

			    }
				
				/*  take sqrt of ||error|| */
				//if(!useRela) {
				//	//std::transform(Nf.begin(), Nf.end(), Nf.begin(), (double(*)(double)) sqrt);
				//} else {
				//	//std::transform(rElative_error.begin(), rElative_error.end(), rElative_error.begin(), (double(*)(double)) sqrt);
				
				//}

				if(!useRela) {
					ierr = VecSetValues(F,data.getIndices().size(),
					&data.getIndices()[0],&*Nf.data().begin(),ADD_VALUES); CHKERRQ(ierr);} else {
					ierr = VecSetValues(F,data.getIndices().size(),
										&data.getIndices()[0],&*rElative_error.data().begin(),ADD_VALUES); CHKERRQ(ierr);	
				}
	            
			}
			
			 catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			 } 
			 PetscFunctionReturn(0);
			 

	
}
	};



/*
  Add the error norm element with same problem and same field as the original problem
 
 */
PetscErrorCode addNormElements(
	const string problem,string fe,const string norm_field_name,
	const string field1_name,const string field2_name,
	//const string field3_name,const string field4_name,
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
	string field2_name,Mat A,Vec &F,bool usel2,bool userela,
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

	map<int,BlockData>::iterator sit = setOfBlocks.begin();
	
	for(;sit!=setOfBlocks.end();sit++) {
		
		//Calculate field values at gaussian points for field1 and field2; 
		feRhs.get_op_to_do_Rhs().push_back(new OpGetTetField1AtGaussPts(field1_name,commonData));
		feRhs.get_op_to_do_Rhs().push_back(new OpGetTetField2AtGaussPts(field2_name,commonData));
		feRhs.get_op_to_do_Rhs().push_back(new OpGetGradField1AtGaussPts(field1_name,commonData));
		feRhs.get_op_to_do_Rhs().push_back(new OpGetGradField2AtGaussPts(field2_name,commonData));
		
		feLhs.get_op_to_do_Lhs().push_back(new OpLhs(norm_field_name,A));

		feRhs.get_op_to_do_Rhs().push_back(new OpRhs(norm_field_name,norm_field_name,F,commonData,usel2,userela));
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





	
