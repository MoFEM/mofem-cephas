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

#ifndef __ELASTICFEMETHOD_STRAIN_CALCULATION_HPP__
#define __ELASTICFEMETHOD_STRAIN_CALCULATION_HPP__

namespace ObosleteUsersModules {
  
  struct ElasticFEMethod_strain_calculation: public ElasticFEMethod {
    
    double wt_Gauss;
    ElasticFEMethod_strain_calculation(FieldInterface& _mField,Mat _Aij,Vec _X,Vec _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
    ElasticFEMethod(_mField,_Aij,_X,_F,_lambda,_mu,_field_name){}
    
    
    
    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;
      
      try {
        
        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        
        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);
        

        unsigned int g_dim = g_NTET.size()/4;
        cout<<"g_dim "<<g_dim<<endl;
        
        ublas::vector<FieldData> DISP;

        DISP.resize(4);
        const EntityHandle* conn;
        int num_nodes;
        rval = mField.get_moab().get_connectivity(fePtr->get_ent(),conn,num_nodes,true); CHKERR_PETSC(rval);
        cout<<"num_nodes  =  "<<num_nodes<<endl;
        for(int nn = 0;nn<num_nodes; nn++) {
          for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,fieldName,conn[nn],iit)) {
            cout<<"iit->get_FieldData()  "<<iit->get_FieldData()<<endl;
          }
        }

//        string wait;
//        cin>>wait;

        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
          try {
            ublas::matrix< FieldData > GradU = *viit;
            if(!invH.empty()) {
              GradU = prod( GradU, invH[gg] );
            }
            
            cout<<"GradU "<<GradU<<endl;

            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            cout<<"VoightStrain "<<VoightStrain<<endl;
            
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

    
    
    
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      
      cout<<"Hi from ElasticFEMethod_strain_calculation "<<endl;
      ierr = Fint(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    
  };
  
  
}

#endif //__ELASTICFEMETHOD_MATRIX_HPP__
