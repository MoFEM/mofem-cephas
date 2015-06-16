/** \file NormElement.hpp
\ingroup mofem_helmholtz_elem
 * \brief Operators and data structures for L^2Norm analysis
 *
 * Implementation of L^2 and H_1 Norm element for error analysis
 *
 * \bug works only for scalar field, in order to implement for vector field,
   look at field approximation.hpp and helmholtzElement.hpp
   In addition, the p error estimator line command option need to be implemented.
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

  /// \brief  Volume element
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
      //PetscErrorCode ierr;
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
        ierr = VecCreateGhostWithArray(mField.get_comm(),0,1,1,g,&eRror,&ghost1); CHKERRQ(ierr);
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
  double& eRror;
  double& aNaly;

  NormElement(
     FieldInterface &mField,double &error,double &analy,int add_to_rank = 1):
     fE(mField,error,analy,add_to_rank),
     m_field(mField),
     addToRank(add_to_rank),
     eRror(error),
     aNaly(analy) {}

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
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
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


  /** \brief Lhs Matrix for Norm Element
    \ingroup mofem_helmholtz_elem

    \f[
    A_{ik} = \int_{\Omega^e} N_i N_k \textrm{d}V
    \f]

  */
    struct OpLhs:public VolumeElementForcesAndSourcesCore::UserDataOperator {

		Mat A;

		OpLhs(const string field_name,Mat _A):
			VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROWCOL),
			A(_A) { }

		OpLhs(const string field_name):
			VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROWCOL) { }

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


    /** \brief Rhs vector for Norm Element
      \ingroup mofem_helmholtz_elem

      \f[
      F_i = \int_{\Omega^e} (\Phi_{ref} - \Phi_{fem}) N_i  \textrm{d}V
      \f]

    */

	struct OpRhs:public VolumeElementForcesAndSourcesCore::UserDataOperator {

		CommonData &commonData;
		Vec F;//norm error
		bool useL2;
		bool useTsF;
		bool useRela;//use relative error
		ublas::vector<double> Nf;
		ublas::vector<double> rElative_error;
		const string normfieldName;
		const string anfieldName1;
    const string anfieldName2;
		const string nufieldName1;
    const string nufieldName2;

    double& eRror;
    double& aNaly;

		OpRhs(const string norm_field_name,const string an_field_name1,const string nu_field_name1,
           const string an_field_name2,const string nu_field_name2,
				   CommonData &common_data,double &error,double &analy,bool usel2):
			VolumeElementForcesAndSourcesCore::UserDataOperator(norm_field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
			commonData(common_data),eRror(error),aNaly(analy),
      useL2(usel2),useTsF(true),normfieldName(norm_field_name)
			,anfieldName1(an_field_name1)
			,nufieldName1(nu_field_name1)
      ,anfieldName2(an_field_name2)
      ,nufieldName2(nu_field_name2) {}

		OpRhs(const string norm_field_name,const string an_field_name1,const string nu_field_name1,
           const string an_field_name2,const string nu_field_name2,
			  Vec _F,CommonData &common_data,double &error,double &analy,bool usel2
			  ):
			VolumeElementForcesAndSourcesCore::UserDataOperator(norm_field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
			commonData(common_data),eRror(error),aNaly(analy),
      useL2(usel2),useTsF(false),F(_F),normfieldName(norm_field_name)
      ,anfieldName1(an_field_name1)
			,nufieldName1(nu_field_name1)
      ,anfieldName2(an_field_name2)
      ,nufieldName2(nu_field_name2) {}


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
				ublas::vector<double> &u_analy1 = commonData.pressureAtGaussPts[anfieldName1];
				ublas::vector<double> &u_numer1 = commonData.pressureAtGaussPts[nufieldName1];

        ublas::vector<double> &u_analy2 = commonData.pressureAtGaussPts[anfieldName2];
				ublas::vector<double> &u_numer2 = commonData.pressureAtGaussPts[nufieldName2];


				ublas::matrix<double> &uAnalyGrad1 = commonData.gradPressureAtGaussPts[anfieldName1];
				ublas::matrix<double> &uNumerGrad1 = commonData.gradPressureAtGaussPts[nufieldName1];

        ublas::matrix<double> &uAnalyGrad2 = commonData.gradPressureAtGaussPts[anfieldName2];
        ublas::matrix<double> &uNumerGrad2 = commonData.gradPressureAtGaussPts[nufieldName2];

				double error;
        double PlainError;

				for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

					//Integrate over volume
					double val = getVolume()*getGaussPts()(3,gg);//this->getGaussPts()(3,gg);
					if(this->getHoGaussPtsDetJac().size()>0) {
						val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry

					}

					const ublas::matrix_row<ublas::matrix<double> > u_analy_grad1(uAnalyGrad1,gg);
					const ublas::matrix_row<ublas::matrix<double> > u_numer_grad1(uNumerGrad1,gg);
          const ublas::matrix_row<ublas::matrix<double> > u_analy_grad2(uAnalyGrad2,gg);
          const ublas::matrix_row<ublas::matrix<double> > u_numer_grad2(uNumerGrad2,gg);

          ublas::vector<double> GradError1 = u_analy_grad1 - u_numer_grad1;
          ublas::vector<double> GradError2 = u_analy_grad2 - u_numer_grad2;

					if(useL2) { //case L2 norm
            /* real and imaginary part of error multiply its complex conjugate */
						error = (u_analy1[gg] - u_numer1[gg])*(u_analy1[gg] - u_numer1[gg])
                    + (u_analy2[gg] - u_numer2[gg])*(u_analy2[gg] - u_numer2[gg]);
						eRror += error*val;
            aNaly += (u_analy1(gg)*u_analy1(gg) + u_analy2(gg)*u_analy2(gg))*val;

					} else if(!useL2) { //case H1 norm

            PlainError = (u_analy1[gg] - u_numer1[gg])*(u_analy1[gg] - u_numer1[gg])
                       + (u_analy2[gg] - u_numer2[gg])*(u_analy2[gg] - u_numer2[gg]);
						error = ublas::inner_prod(GradError1,GradError1) + ublas::inner_prod(GradError2,GradError2)
                    + PlainError;

            aNaly += (u_analy1(gg)*u_analy1(gg) + u_analy2(gg)*u_analy2(gg) +
                     ublas::inner_prod(u_analy_grad1,u_analy_grad1)
                      + ublas::inner_prod(u_analy_grad2,u_analy_grad2)) * val;
						eRror += (ublas::inner_prod(GradError1,GradError1) +
                     ublas::inner_prod(GradError2,GradError2) + PlainError)*val;

					}
					//need to calculate sqrt of norm^2
					ublas::noalias(Nf) += val*error*data.getN(gg,nb_row);


				}

				/*  take sqrt of ||error|| local*/
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



    /** \brief Add Norm Element Problem
      * \ingroup mofem_helmholtz_elem
      *
      * It get data from block set and define element in moab
      *w
      * \param problem name
      * \param field name
      * \param name of mesh nodal positions (if not defined nodal coordinates are used)
      */
  PetscErrorCode addNormElements(
      const string problem,string fe,const string norm_field_name,
      const string an_field_name1,
      const string nu_field_name1,const string an_field_name2,
      const string nu_field_name2,
      const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ErrorCode rval;
      ierr = m_field.add_finite_element(fe,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row(fe,norm_field_name); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col(fe,norm_field_name); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data(fe,norm_field_name); CHKERRQ(ierr);

      ierr = m_field.modify_finite_element_add_field_data(fe,an_field_name1); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data(fe,nu_field_name1); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data(fe,an_field_name2); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data(fe,nu_field_name2); CHKERRQ(ierr);


      if(m_field.check_field(mesh_nodals_positions)) {
        ierr = m_field.modify_finite_element_add_field_data(fe,mesh_nodals_positions); CHKERRQ(ierr);
      }
      ierr = m_field.modify_problem_add_finite_element(problem,fe); CHKERRQ(ierr);


      //Range tEts;
      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"MAT_NORM",it)) {

        rval = m_field.get_moab().get_entities_by_type(it->get_meshset(),MBTET,volumeData[it->get_msId()].tEts,true); CHKERR_PETSC(rval);

        ierr = m_field.add_ents_to_finite_element_by_TETs(volumeData[it->get_msId()].tEts,fe); CHKERRQ(ierr);
      }


      PetscFunctionReturn(0);

     }





  PetscErrorCode setNormFiniteElementRhsOperator(string norm_field_name,string an_field_name1,
      string nu_field_name1,string an_field_name2,
      string nu_field_name2,Vec &F,bool usel2,
      string nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;

      fE.getOpPtrVector().push_back(new OpGetValueAndGradAtGaussPts(an_field_name1,commonData));
      fE.getOpPtrVector().push_back(new OpGetValueAndGradAtGaussPts(nu_field_name1,commonData));
      fE.getOpPtrVector().push_back(new OpGetValueAndGradAtGaussPts(an_field_name2,commonData));
      fE.getOpPtrVector().push_back(new OpGetValueAndGradAtGaussPts(nu_field_name2,commonData));


      map<int,VolumeData>::iterator sit = volumeData.begin();

      for(;sit!=volumeData.end();sit++) {

          //Calculate field values at gaussian points for field1 and field2;
        fE.getOpPtrVector().push_back(new OpRhs(norm_field_name,an_field_name1,nu_field_name1,an_field_name2,nu_field_name2,F,commonData,eRror,aNaly,usel2));

      }

      PetscFunctionReturn(0);
  }


  PetscErrorCode setNormFiniteElementLhsOperator(string norm_field_name,
      Mat A,bool usel2 = false,
      string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    map<int,VolumeData>::iterator sit = volumeData.begin();

    for(;sit!=volumeData.end();sit++) {

      fE.getOpPtrVector().push_back(new OpLhs(norm_field_name,A));

    }

    PetscFunctionReturn(0);
   }


  };

  }


#endif //__NORM_ELEMENT_HPP


/***************************************************************************//**
 *
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
