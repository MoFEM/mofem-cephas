/** \file NormElement.hpp
\ingroup mofem_helmholtz_elem
 * \brief Operators and data structures for L^2Norm analysis
 *
 * Implementation of L^2 and H_1 Norm element for error analysis
 *
 * \bug works only for scalar field, in order to implement for vector field, look at field approximation.hpp and helmholtzElement.hpp
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

  double& eRror;
  double& aNaly;
  
  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
    FieldInterface& mField;
    int addToRank; ///< default value 1, i.e. assumes that geometry is approx. by quadratic functions.
    double& eRror;
    double& aNaly;
    MyVolumeFE(FieldInterface &m_field,double &error,double &analy,int add_to_rank): 
        VolumeElementForcesAndSourcesCore(m_field),mField(m_field),addToRank(add_to_rank),
        eRror(error),aNaly(analy) {}
    int getRule(int order) { return order+addToRank; };
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin; 
      PetscErrorCode ierr;
      eRror = 0;
      aNaly = 0;
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode postProcess() {
      PetscFunctionBegin; 
      PetscErrorCode ierr;
      
      int rank;
      MPI_Comm_rank(mField.get_comm(),&rank);
      Vec ghost1;
      Vec ghost2;
      if(!rank) {
        ierr = VecCreateGhostWithArray(mField.get_comm(),1,1,0,NULL,&eRror,&ghost1); CHKERRQ(ierr);
        ierr = VecCreateGhostWithArray(mField.get_comm(),1,1,0,NULL,&aNaly,&ghost2); CHKERRQ(ierr);
      } else {
        int g[] = {0};
        //cout << "\n g = \n" << g << endl;
        //cout << "\n ghost = \n" << ghost << endl;
        ierr = VecCreateGhostWithArray(mField.get_comm(),0,1,1,g,&eRror,&ghost1); CHKERRQ(ierr);
        //cout << "\n g after = \n" << g << endl;
        //cout << "\n ghost after = \n" << ghost << endl;
        ierr = VecCreateGhostWithArray(mField.get_comm(),0,1,1,g,&aNaly,&ghost2); CHKERRQ(ierr);
      }

      ierr = VecGhostUpdateBegin(ghost1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(ghost1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(ghost1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(ghost1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecDestroy(&ghost1); CHKERRQ(ierr);
      
      ierr = VecGhostUpdateBegin(ghost2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(ghost2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(ghost2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(ghost2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecDestroy(&ghost2); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }
    
  };
  
  MyVolumeFE fE; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFe() { return fE; } ///< get lhs volume element
	
  FieldInterface &m_field;
  int addToRank; ///< default value 1, i.e. assumes that geometry is approx. by quadratic functions.
  
  NormElement(
     FieldInterface &mField,double &error,double &analy,int add_to_rank = 1):
     fE(mField,error,analy,add_to_rank),m_field(mField),addToRank(add_to_rank),
     eRror(error),aNaly(analy) {}
	
  //Field data
  struct CommonData {
	
	map<string,ublas::vector<double> > pressureAtGaussPts; 
	map<string,ublas::matrix<double> > gradPressureAtGaussPts;
	ublas::matrix<double> hoCoords;
	 
	map<EntityType, vector< ublas::vector<int> > > imIndices;
	  
	  
  };
  CommonData commonData;
	
  struct VolumeData {
	Range tEts; ///< constatins elements in block set
  };
  map<int,VolumeData> volumeData; ///< maps block set id with appropiate VolumeData	
  
	
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
	
	
	/** \brief Lhs operaetar for tetrahedral used to build matrix
	*/
    struct OpLhs:public VolumeElementForcesAndSourcesCore::UserDataOperator {
		
		Mat A;
		//bool solveBc;
		OpLhs(const string re_field_name,Mat _A): 
			VolumeElementForcesAndSourcesCore::UserDataOperator(re_field_name),
			A(_A) { }
		
		OpLhs(const string re_field_name): 
			VolumeElementForcesAndSourcesCore::UserDataOperator(re_field_name) { }
		
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
	struct OpRhs:public VolumeElementForcesAndSourcesCore::UserDataOperator {
		
		CommonData &commonData;
		Vec F;//norm error
		//Vec D;//relative error
		bool useL2;
		bool useTsF;
		bool useRela;//use relative error
		ublas::vector<double> Nf;
		ublas::vector<double> rElative_error;
		const string normfieldName;
		const string refieldName;
		const string imfieldName;
		
        double& eRror;
        double& aNaly;
        
		OpRhs(const string norm_field_name,const string re_field_name,const string im_field_name,
				   CommonData &common_data,double &error,double &analy,bool usel2,bool userela): 
			VolumeElementForcesAndSourcesCore::UserDataOperator(norm_field_name),
			commonData(common_data),eRror(error),aNaly(analy),
             useL2(usel2),useTsF(true),useRela(userela),normfieldName(norm_field_name)
			,refieldName(re_field_name)
			,imfieldName(im_field_name) {}
		
		OpRhs(const string norm_field_name,const string re_field_name,const string im_field_name,
			  Vec _F,CommonData &common_data,double &error,double &analy,bool usel2,bool userela
			  ): 
			VolumeElementForcesAndSourcesCore::UserDataOperator(norm_field_name),
			commonData(common_data),eRror(error),aNaly(analy),
            useL2(usel2),useTsF(false),useRela(userela),F(_F),normfieldName(norm_field_name)
			,refieldName(re_field_name)
			,imfieldName(im_field_name) {}
		
		
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
				ublas::vector<double> &u_analy = commonData.pressureAtGaussPts[refieldName];
				ublas::vector<double> &u_numer = commonData.pressureAtGaussPts[imfieldName];
				
				ublas::matrix<double> &uAnalyGrad = commonData.gradPressureAtGaussPts[refieldName];
				ublas::matrix<double> &uNumerGrad = commonData.gradPressureAtGaussPts[imfieldName];
				double error;
				
				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
                  
					//Integrate over volume
					double val = getVolume()*getGaussPts()(3,gg);//this->getGaussPts()(3,gg); 
					if(this->getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry

					}
					
					const ublas::matrix_row<ublas::matrix<double> > u_analy_grad(uAnalyGrad,gg);
					const ublas::matrix_row<ublas::matrix<double> > u_numer_grad(uNumerGrad,gg);

                    ublas::vector<double> GradError = u_analy_grad - u_numer_grad;
                    
					if(useL2) { //case L2 norm
                      
						error = u_analy[gg] - u_numer[gg];
						eRror += error*error*val;
                        aNaly += u_analy(gg)*u_analy(gg)*val;
					} else if(!useL2) { //case H1 norm
					
						error = u_analy[gg] - u_numer[gg];
						
                        aNaly += u_analy(gg)*u_analy(gg)*val;
						eRror += (ublas::inner_prod(GradError,GradError) + error*error)*val;
						
					}
					//need to calculate sqrt of norm^2
					if(!useRela) { //case Norm error
						
						ublas::noalias(Nf) += val*error*data.getN(gg,nb_row);
						
					} else if(useRela) { //case relative error

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
					&data.getIndices()[0],&*Nf.data().begin(),ADD_VALUES); CHKERRQ(ierr);} 
				else {
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
      rval = m_field.get_moab().get_entities_by_type(it->get_meshset(),MBTET,volumeData[it->get_msId()].tEts,true); CHKERR_PETSC(rval);

      ierr = m_field.add_ents_to_finite_element_by_TETs(volumeData[it->get_msId()].tEts,fe); CHKERRQ(ierr);
	}
	
    
    PetscFunctionReturn(0);
	
	}


    
    

PetscErrorCode setNormFiniteElementRhsOperator(string norm_field_name,string field1_name,
	string field2_name,Vec &F,bool usel2,bool userela,
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

	fE.getRowOpPtrVector().push_back(new OpGetValueAndGradAtGaussPts(field1_name,commonData));
	fE.getRowOpPtrVector().push_back(new OpGetValueAndGradAtGaussPts(field2_name,commonData));
	
	map<int,VolumeData>::iterator sit = volumeData.begin();
	
	for(;sit!=volumeData.end();sit++) {
		
		//Calculate field values at gaussian points for field1 and field2; 
      fE.getRowOpPtrVector().push_back(new OpRhs(norm_field_name,field1_name,field2_name,F,commonData,eRror,aNaly,usel2,userela));

	}
	
	PetscFunctionReturn(0);
}


PetscErrorCode setNormFiniteElementLhsOperator(string norm_field_name,string field1_name,
	string field2_name,Mat A,bool usel2 = false,bool userela = false,
    string nodals_positions = "MESH_NODE_POSITIONS") {
  PetscFunctionBegin;

  map<int,VolumeData>::iterator sit = volumeData.begin();
  
  for(;sit!=volumeData.end();sit++) {
    
    //Calculate field values at gaussian points for field1 and field2; 

    fE.getRowColOpPtrVector().push_back(new OpLhs(norm_field_name,A));

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





	
