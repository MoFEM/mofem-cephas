/* Copyright (C) 2013, Michel Cortis <mikecortis at gmail.com>
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

#ifndef __ELASTICFEMETHODTRANSISORVE_HPP__
#define __ELASTICFEMETHODTRANSISORVE_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethodTransIso.hpp"

namespace MoFEM {

struct ElasticFEMethodTransIsoRVE: public TranIsotropicFibreDirRotElasticFEMethod{

    ElasticFEMethodTransIsoRVE(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _D,Vec& _F,
                                            double _lambda,double _mu,double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp):
    TranIsotropicFibreDirRotElasticFEMethod(_mField, _dirihlet_ptr,_Aij,_D, _F, _lambda, _mu, _E_p, _E_z, _nu_p, _nu_pz, _G_zp)  {
        
        double def_VAL2[3] = {0,0,0};
        // create a new tag th_fibre_dir
        rval = moab.tag_get_handle( "POT_FLOW_FIBRE_DIR",3,MB_TYPE_DOUBLE,th_fibre_dir,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
};
    
    
    PetscErrorCode Fint() {
        PetscFunctionBegin;
        
        double _lambda,_mu;
        ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
        ierr = calulateD(_lambda,_mu); CHKERRQ(ierr);
        
        //Gradient at Gauss points;
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        ierr = GetGaussDiffDataVector("DISPLACEMENT_FROM_APP_STRAIN",GradU_at_GaussPt); CHKERRQ(ierr);
        
//        cout<<"GradU_at_GaussPt "<< GradU_at_GaussPt[0] << endl;
        
        unsigned int g_dim = g_NTET.size()/4;
        assert(GradU_at_GaussPt.size() == g_dim);
        NOT_USED(g_dim);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
            ublas::matrix< FieldData > GradU = *viit;
            ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
            ublas::vector< FieldData > VoightStrain(6);
            VoightStrain[0] = Strain(0,0);
            VoightStrain[1] = Strain(1,1);
            VoightStrain[2] = Strain(2,2);
            VoightStrain[3] = 2*Strain(0,1);
            VoightStrain[4] = 2*Strain(1,2);
            VoightStrain[5] = 2*Strain(2,0);
            double w = V*G_W_TET[gg];
            ublas::vector<FieldData> VoightStress = prod(w*D_At_GaussPoint[gg],VoightStrain);
            //BT * VoigtStress
            for(int rr = 0;rr<row_mat;rr++) {
                f_int.resize(row_mat);
                ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
                if(gg == 0) {
                    f_int[rr] = prod( trans(B), VoightStress );
                } else {
                    f_int[rr] += prod( trans(B), VoightStress );
                }
            }
        }
        
        PetscFunctionReturn(0);
    }
    
    
    
    
    PetscErrorCode Fint(Vec F_int) {
        PetscFunctionBegin;
        ierr = ElasticFEMethodTransIsoRVE::Fint(); CHKERRQ(ierr);
        for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
            
        }
        PetscFunctionReturn(0);
    }

    
    
};

}
#endif //__ELASTICFEMETHODTRANSISORVE_HPP__

