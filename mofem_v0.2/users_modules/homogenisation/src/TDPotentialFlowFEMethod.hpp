/* Copyright (C) 2014, GUOQIANG XUE (XUE AT Glasgow.ac.uk)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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


#ifndef __TDPOTENTIALFLOWFEMETHOD_HPP__
#define __TDPOTENTIALFLOWFEMETHOD_HPP__

#include<moab/Skinner.hpp>

namespace MoFEM {
	
	/** \brief struture grouping operators and data used for 2D potential problems
	 *
	 * In order to assemble matrices and right hand vectors, the loops over
	 * elements, enetities over that elememnts and finally loop over intergration
	 * points are executed.
	 *
	 * Following implementation separte those three cegories of loops and to eeach
	 * loop attach operator.
	 *
	 */
struct TDPotentialElement {	
	
 // brief definition of volume element
 // brefe define surface element, this element is used to intergrate 2d-surface potential flow
	struct MyTriFE: public TriElementForcesAndSurcesCore {
		MyTriFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
		int getRule(int order) { return ceil(order/2); };
	};
	
	MyTriFE fe2DpoetentialRhs; //
	MyTriFE& getLoopFe2DpoetentialRhs() { return fe2DpoetentialRhs; } //
	MyTriFE fe2DpoetentialLhs;
	MyTriFE& getLoopFe2DpoetentialLhs() { return fe2DpoetentialLhs; }
	
	FieldInterface &mField;
	TDPotentialElement(FieldInterface &m_field):
    fe2DpoetentialRhs(m_field),fe2DpoetentialLhs(m_field),
    mField(m_field) {}
	
	/** \brief add 2D Potential element
	 * \infroup mofem_2DPotential_elem
	 *
	 * It get data from 2dPOETNTIAL set and define elemenet in moab. Alternatively
	 * uses block set with name surface101.
	 *
	 * \param problem name
	 * \param field name
	 * \param name of mesh nodal positions (if not defined nodal coordinates are used)
	 */
	
	struct SurfaceData {
		Range tRis;}; 
	map<int,SurfaceData> setOfSurface; //< maps block set id with appropiate data
	

//------------------111----------------------------//	
	
	PetscErrorCode add2DPotential_Element(
											   const string problem_name,const string field_name) {
		PetscFunctionBegin;
		
		PetscErrorCode ierr;
		ErrorCode rval;
		
		ierr = mField.add_finite_element("2DPOTENTIAL_ELEM",MF_ZERO); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_row("2DPOTENTIAL_ELEM",field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_col("2DPOTENTIAL_ELEM",field_name); CHKERRQ(ierr);
		ierr = mField.modify_finite_element_add_field_data("2DPOTENTIAL_ELEM",field_name); CHKERRQ(ierr);
//		if(mField.check_field(mesh_nodals_positions)) {
//			ierr = mField.modify_finite_element_add_field_data("2DPOTENTIAL_ELEM",mesh_nodals_positions); CHKERRQ(ierr);
//		}
		ierr = mField.modify_problem_add_finite_element(problem_name,"2DPOTENTIAL_ELEM"); CHKERRQ(ierr);
		
		/*==================create the meshset for surfaces============================*/

		for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
			if(it->get_Cubit_name().compare(0,11,"Sideset_TRI") == 0) {
				rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfSurface[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
				ierr = mField.add_ents_to_finite_element_by_TRIs(setOfSurface[it->get_msId()].tRis,"2DPOTENTIAL_ELEM"); CHKERRQ(ierr);
				
			}
		}

		
		PetscFunctionReturn(0);
	}
//--------------------111----------------------------//	

	/* \brief finite terms in rhs and lhs
	 */
	
	struct OpTDPotentialLhs:public TriElementForcesAndSurcesCore::UserDataOperator {
		
		SurfaceData &dAta;
		bool useTsB;
		
		OpTDPotentialLhs(const string field_name,
						SurfaceData &data):
		TriElementForcesAndSurcesCore::UserDataOperator(field_name),
		dAta(data),useTsB(true) {}
		
		Mat A;
		OpTDPotentialLhs(const string field_name,Mat _A,
						SurfaceData &data):
		TriElementForcesAndSurcesCore::UserDataOperator(field_name),
		dAta(data),useTsB(false),A(_A) {}
		
		ublas::matrix<double> K,transK;
		/** \brief calculate  stiffness term in the lhs of equations
		 *
		 * K = intS diffN^T alpha diffN dS
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
				
				for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

					double val = getArea()*getGaussPts()(2,gg);
////					if(getHoGaussPtsDetJac().size()>0) {
////						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
////					}
//					cout << "val============"<<val<<endl;
//          cout << "nb_row============"<<nb_row<<endl;
//          cout << "gg============"<<gg<<endl;
//
//          cout << "row_data.getDiffN(gg,nb_row)=="<<row_data.getDiffN()<<endl;

					noalias(K) += val*prod(row_data.getDiffN(),trans(col_data.getDiffN()));  //?

				}

				//cerr << k << endl;
				
				PetscErrorCode ierr;
				if(!useTsB) {
					const_cast<FEMethod*>(getFEMethod())->ts_B = A;
				}
				ierr = MatSetValues(
									(getFEMethod()->ts_B),
									nb_row,&row_data.getIndices()[0],
									nb_col,&col_data.getIndices()[0],
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
	
//--------------------222-----------------------------//
	PetscErrorCode setTDPotentialFiniteElementLhsOperators(string field_name,Mat A,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
		PetscFunctionBegin;
		map<int,SurfaceData>::iterator sit = setOfSurface.begin();
		for(;sit!=setOfSurface.end();sit++) {
			//add finite element
			fe2DpoetentialLhs.get_op_to_do_Lhs().push_back(new OpTDPotentialLhs(field_name,A,sit->second));
		}
		PetscFunctionReturn(0);
	}
//-------------------222----------------------------//
	
	
	
	
	
	
};
}


#endif // __TDPOTENTIALFLOWFEMETHOD_HPP__

