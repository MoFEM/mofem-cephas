/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * DESCRIPTION:Analytical equation for Helmholtz operator and
 * The Norm of approximation error implemented.
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

#ifndef __ANALYTICAL_ELEMENT_HPP__
#define __ANALYTICAL_ELEMENT_HPP__

using namespace boost::numeric;
#include<moab/Skinner.hpp>

#define RND_EPS 1e-4

//Rounding
double roundn(double n)
{
	//break n into fractional part (fract) and integral part (intp)
	double fract, intp;
	fract = modf(n,&intp);
	
//    //round up
//    if (fract>=.5)
//    {
//        n*=10;
//        ceil(n);
//        n/=10;
//    }
//	//round down
//    if (fract<=.5)
//    {
//		n*=10;
//        floor(n);
//        n/=10;
//    }
	// case where n approximates zero, set n to "positive" zero
	if (abs(intp)==0)
	{
		if(abs(fract)<=RND_EPS)
		{
			n=0.0000;
		}
	}
	return n;
}


namespace MoFEM {

struct AnalyticalSolution {

  //ubals::vector<double> analyticalFunction(
    //double x,double y,double z) = 0;

  /** \brief finite element to appeximate analytical solution on surface
    */

	struct MyVolumeFE: public TetElementForcesAndSourcesCore {
		MyVolumeFE(FieldInterface &m_field): TetElementForcesAndSourcesCore(m_field) {}
		
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
		int getRule(int order) { return order; }; //wait for confirmation
	};

	MyVolumeFE FeRhs; ///< cauclate right hand side for tetrahedral elements
	MyVolumeFE& getLoopFeRhs() { return FeRhs; } ///< get rhs volume element
	MyVolumeFE FeLhs; //< calculate left hand side for tetrahedral elements
	MyVolumeFE& getLoopFeLhs() { return FeLhs; } ///< get lhs volume element
	
	FieldInterface &m_field;

	AnalyticalSolution(FieldInterface &mField,Range &tets):
		FeRhs(mField),FeLhs(mField),
		m_field(mField),tEts(tets) {}
	
	//MAT_HELMHOLTZ DATA
	struct BlockData {
		double aNgularfreq; // Angular frequency
		double sPeed;   // Wave Speed
		
		//Range tEts; ///< constatins elements in block set
	};
	map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData
	BlockData blockData;
	
	ublas::matrix<double> hoCoordsTet;
	
	//struct BlockData {
		Range tEts; ///< constatins elements in block set
	//};
	//map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData
	
	  
	struct OpHoCoordTet: public TetElementForcesAndSourcesCore::UserDataOperator {
	
		
		ublas::matrix<double> &hoCoordsTet;
		OpHoCoordTet(const string field_name,ublas::matrix<double> &ho_coords): 
			TetElementForcesAndSourcesCore::UserDataOperator(field_name),
			hoCoordsTet(ho_coords) {}
	
		/*	
		Cartesian coordinates for gaussian points inside elements
		X^coordinates = DOF dot* N
		*/
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
	
			try {
	
				if(data.getFieldData().size()==0) PetscFunctionReturn(0);
	
				hoCoordsTet.resize(data.getN().size1(),3);
				if(type == MBVERTEX) {
					hoCoordsTet.clear();
				}
	
				int nb_dofs = data.getN().size2();
				//calculate x,y,z in each GaussPts
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
					for(int dd = 0;dd<3;dd++) {  
						hoCoordsTet(gg,dd) += cblas_ddot(nb_dofs,&data.getN(gg)[0],1,&data.getFieldData()[dd],3); 
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
	
  
    /** \brief Rhs operaetar for tetrahedral used to build matrix
      */
	
	struct OpRhs:public TetElementForcesAndSourcesCore::UserDataOperator {
		Vec C;
		ublas::matrix<double> &hoCoords;
		double (*fUN)(double x,double y,double z);
		bool solveBc;
		
		OpRhs(const string re_field_name,ublas::matrix<double> &ho_coords,
				 double (*fun)(double x,double y,double z),Vec _C): 
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name),
			hoCoords(ho_coords),fUN(fun),C(_C),solveBc(false)  {}
		
		OpRhs(const string re_field_name,ublas::matrix<double> &ho_coords,
				 double (*fun)(double x,double y,double z)): 
			TetElementForcesAndSourcesCore::UserDataOperator(re_field_name),
			hoCoords(ho_coords),fUN(fun),solveBc(true) {}
		
		ublas::vector<FieldData> NTf;
		
		/*	
		Rhs force vector merely with Dirichlet values
		F = int_S N^T Fun dS
		*/
		
		PetscErrorCode doWork(
			int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
			PetscFunctionBegin;
			PetscErrorCode ierr;
			
			try {

				if(data.getIndices().size()==0) PetscFunctionReturn(0);
				
				int nb_row_dofs = data.getIndices().size();
				unsigned int nb_row = data.getN().size2(); //no. of DOF
				if(nb_row != data.getIndices().size()) {
					SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
							"currently works only for scalar fields, extension to fields with higher rank need to be implemented");
				}
				NTf.resize(nb_row);

				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

					double x,y,z;
					double val = getGaussPts()(3,gg)*getVolume();
					
					if(hoCoords.size1() == data.getN().size1()) {
						
						val *= getHoGaussPtsDetJac()[gg];
						x = hoCoords(gg,0);
						y = hoCoords(gg,1);
						z = hoCoords(gg,2);
					} else {
						
						x = getCoordsAtGaussPts()(gg,0);
						y = getCoordsAtGaussPts()(gg,1);
						z = getCoordsAtGaussPts()(gg,2);
					}
					
					double a = fUN(x,y,z);
					
					
					noalias(NTf) = val*data.getN(gg)*a;
				}
				
				
					if(solveBc){
                    //Set the value for Dirichlet BC
					ierr = VecSetValues(getFEMethod()->snes_f,data.getIndices().size(),
										&data.getIndices()[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);
					} else if(!solveBc){
					//Set the value to calculate the interpolation of Exact Solution 
					ierr = VecSetValues(C,data.getIndices().size(),
										&data.getIndices()[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);
					}
				


			} catch (const std::exception& ex) {
				ostringstream ss;
				ss << "throw in method: " << ex.what() << endl;
				SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
			}

			PetscFunctionReturn(0);
		}

    };
  
	
	PetscErrorCode initializeExactProblem(
		string problem,string fe,string re_field_name,
		string nodals_positions = "MESH_NODE_POSITIONS") {
		PetscFunctionBegin;
		PetscErrorCode ierr;
		ErrorCode rval;
		
		//Add tets elements, wait for correction.
		ierr = m_field.add_finite_element(fe,MF_ZERO); CHKERRQ(ierr);
		ierr = m_field.modify_finite_element_add_field_row(fe,re_field_name); CHKERRQ(ierr);
		ierr = m_field.modify_finite_element_add_field_col(fe,re_field_name); CHKERRQ(ierr);
		ierr = m_field.modify_finite_element_add_field_data(fe,re_field_name); CHKERRQ(ierr);
		
		if(m_field.check_field(nodals_positions)) {
			ierr = m_field.modify_finite_element_add_field_data(fe,nodals_positions); CHKERRQ(ierr);
		}
		//set finite elements for problem
		ierr = m_field.modify_problem_add_finite_element(problem,fe); CHKERRQ(ierr);
		
		//loop over the blockset entity.
		for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
			if(it->get_Cubit_name().compare(0,9,"EXACT_SOL") == 0) {
				
				rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTET,tEts,true); CHKERR_PETSC(rval);
				ierr = m_field.add_ents_to_finite_element_by_TETs(tEts,fe); CHKERRQ(ierr);
		
			}
		}
		
		for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
			
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
				
			}
		}
		
		PetscFunctionReturn(0);
	}
  
  
  PetscErrorCode setExactSolRhsOp(
    string field_name,
    double (*fun)(double x,double y,double z),
	Vec &C,
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
	//loop over tets
    if(m_field.check_field(nodals_positions)) {
		FeRhs.get_op_to_do_Rhs().push_back(new OpHoCoordTet(nodals_positions,hoCoordsTet));
    }
	
	map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
	FeRhs.get_op_to_do_Rhs().push_back(new OpRhs(field_name,hoCoordsTet,fun,C));
		}
	
	PetscFunctionReturn(0);
  }
    
  
  PetscErrorCode solveExactProblem(
	  string problem,string fe,string re_field_name,Vec &C) {
	  PetscFunctionBegin;
	  PetscErrorCode ierr;
  
	  
	  ierr = VecZeroEntries(C); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);


	  
	  ierr = m_field.set_global_VecCreateGhost(problem,ROW,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = m_field.loop_finite_elements(problem,fe,getLoopFeRhs()); CHKERRQ(ierr); //wait; where problem occurs

	  ierr = VecGhostUpdateBegin(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(C,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  
	  ierr = m_field.set_global_VecCreateGhost(problem,ROW,C,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	  //PetscViewer viewer;
	  //ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Vec C.txt",&viewer); CHKERRQ(ierr);
	  //VecView(C,viewer);
	  //ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	  

	  
	  PetscFunctionReturn(0);
  }

};

}

#endif //__ANALYTICAL_ELEMENT_HPP__

