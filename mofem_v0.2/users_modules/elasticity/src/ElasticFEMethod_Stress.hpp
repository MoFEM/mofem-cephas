/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ELASTICFEMETHOD_STRESS_HPP__
#define __ELASTICFEMETHOD_STRESS_HPP__

namespace ObosleteUsersModules {
  
  struct ElasticFEMethod_Stress: public ElasticFEMethod {
    
    ublas::matrix<double> ptQuery;
    
    ElasticFEMethod_Stress(FieldInterface& _mField,Mat _Aij,Vec _X,Vec _F,double _lambda,double _mu, ublas::matrix<double> _ptquery, string _field_name = "DISPLACEMENT"):
    ElasticFEMethod(_mField,_Aij,_X,_F,_lambda,_mu,_field_name),ptQuery(_ptquery) {}
    
    
    
    
    struct CommonData {
      map<EntityHandle, ublas::vector<ublas::vector<double> > > StressMap;  //map[element, nb_gausspoint[stress]]
    };
    CommonData commonData;

    
    virtual PetscErrorCode Calculate_Stress() {
      PetscFunctionBegin;
      
      try {
        
        double _lambda,_mu;
        ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
        ierr = calculateD(_lambda,_mu); CHKERRQ(ierr);
//        cout<<"D "<< D <<endl;

        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);
        unsigned int g_dim = g_NTET.size()/4;
        
        EntityHandle fe_ent;  fe_ent=fePtr->get_ent();
        commonData.StressMap[fe_ent].resize(g_dim);

//        cout<<"g_dim "<< g_dim <<endl;

        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ublas::matrix< FieldData > GradU = *viit;
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
//            cout<<"VoightStrain "<< VoightStrain <<endl;

            ublas::vector<FieldData> VoightStress = prod(D,VoightStrain);
            cout<<"VoightStress "<< VoightStress <<endl;
            commonData.StressMap[fe_ent](gg).resize(6);
            commonData.StressMap[fe_ent](gg)=VoightStress;
//            cout<<"commonData.StressMap[fe_ent](gg) "<< commonData.StressMap[fe_ent](gg) <<endl;

            

          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
          }
        }
        
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      
      PetscFunctionReturn(0);
    }

    
    
    ublas::matrix<double> gaussPts;
    virtual PetscErrorCode Get_g_NTET() {
      PetscFunctionBegin;
      
      
      
      int nb_gauss_pts=ptQuery.size2();
      gaussPts.resize(4,nb_gauss_pts);
      for(int ii=0; ii<nb_gauss_pts; ii++){
        for(int jj=0; jj<3; jj++){
          gaussPts(jj,ii)=ptQuery(jj,ii);
          gaussPts(3,ii)=0;
        }
      }
      
      g_NTET.resize(4*nb_gauss_pts);
      ierr = ShapeMBTET(&g_NTET[0],&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
      G_TET_W = &gaussPts(3,0);

//      cout<<"gaussPts "<<gaussPts<<endl;
      
      PetscFunctionReturn(0);
    }

    
    
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

//      cout<<"ptQuery.size2() "<<ptQuery.size2()<<endl;
//      cout<<"ptQuery.size1() "<<ptQuery.size1()<<endl;

//      cout<<"Hello from new operator"<<endl;
      ierr = Calculate_Stress(); CHKERRQ(ierr);
      
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    
    
    
    
  };
  
  
}

#endif //__ELASTICFEMETHOD_STRESS_HPP__
