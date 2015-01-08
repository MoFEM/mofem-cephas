/** \file GroundSurfaceTemerature.hpp 
 * \brief Operators and data structures for thermal analys
 *
 * Implementation of boundary conditions for ground temerature
 *
 */

/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __GROUNDSURFACETEMERATURE_HPP
#define __GROUNDSURFACETEMERATURE_HPP

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif 

/** \brief Implementation of ground surface temperature
  *
  * Ground surface temperature simulation for different land covers
  * William R. Herb *, Ben Janke, Omid Mohseni, Heinz G. Stefan
  * Journal of Hydrology (2008) 356, 327– 343
  */
struct GroundSurfaceTemerature {

  ThermalElement::CommonData &commonData;
  GroundSurfaceTemerature(ThermalElement &termal_elem):
    commonData(termal_elem.commonData) {};

  struct Parameters {	

    double alpha; 	//< Solar albedo
    double Cfc; 	//< Surface heat/moisture transfer coefficient for forced convection
    double Cnc;		//< Coefficient for natural convection
    double CSh;		//< Wind sheltering coefficient
    double eps;		//< Pavement emissivity
    double rhoCp;	//< Density specific heat pavement (J/m3/°C)
    Range tRis;		//< Triangles on which parameters are defined

    Parameters(Range tris): 
      tRis(tris) {}

  };

  struct Asphalt: public Parameters {
    Asphalt(Range tris): Parameters(tris) {
      alpha = 0.12;
      Cfc = 0.0015;
      Cnc = 0.0015;
      CSh = 1.;
      eps = 0.94;
      rhoCp = 2.0e06;
    }
  };

  struct Concrete: public Parameters {
    Concrete(Range tris): Parameters(tris) {
      alpha = 0.20;
      Cfc = 0.0015;
      Cnc = 0.0015;
      CSh = 1.;
      eps = 0.94;
      rhoCp = 2.0e06;
    }
  };

  struct BareSoil: public Parameters {
    BareSoil(Range tris): Parameters(tris) {
      alpha = 0.15;
      Cfc = 0.003;
      Cnc = 0.0015;
      CSh = 1.;
      eps = 0.94;
      rhoCp = 2.0e06;
    }
  };

  struct TmieDependendData {

    double T0; // reference temperature (K)
    double e0; // reference saturation vapor pressure (es at a certain temp, usually 0 deg C) (Pa)
    double Rv; // gas constant for water vapor (J*K/Kg)
    double Lv; // latent heat of vaporization of water (J)

    double u10;		//< wind at high 10m (m/s)	
    double CR; 		//< cloudness factor (0–1, dimensionless)
    double Ta;		//< air temperature (C)
    double Td;		//< dew point temperature (C)
    double P;		//< pressure
    double Rs; 		//< observed solar radiation (W/m2)

    double ea;		//< athmospheric vapour preassure (Pa)
    double phia;	//< athmospheric virtual temperature

    template <typename TYPE> 
    double calulateVapourPressure(TYPE T) { 
      return es0*exp((lv/Rv)*((1./T0)-(1./(T+T0))));
    }
 
    template <typename TYPE>   
    double calulateRH(TYPE T) { 
      return calulateVapourPressure(Td)/calulateVapourPressure(T);
    }

    template <typename TYPE>   
    double calulateMixingRatio(TYPE T,double P) {
      const double c = 0.62197;
      double e = calulateVapourPressure(T);
      return c*e/(P-e);
    }

    template <typename TYPE>   
    double calculateAbsoluteVirtualTempertaure(TYPE T,double P) {
      const double c = 0.379;
      double e = calulateVapourPressure(T);
      return (T+T0)/(1-c*e/P);
    }

    TmieDependendData() {

      T0 = 273.15; // reference temperature (K)
      e0 = 611; // reference saturation vapor pressure (es at a certain temp, usually 0 deg C) (Pa)
      Rv = 461.5; // gas constant for water vapor (J*K/Kg)
      Lv = 2.5e6; // latent heat of vaporization of water (J)

    }

    virtual set() = 0;
  }

  boost::ptr_vector<Parameters> blockData;

  PetscErrorCode addSurfaces() {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;
  
    ierr = mField.add_finite_element("GROUND_SURFACE_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("GROUND_SURFACE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("GROUND_SURFACE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("GROUND_SURFACE_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("GROUND_SURFACE_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->get_Cubit_name().compare(0,9,"ASPHALT") == 0) {
	Range tris;
        rval = mField.get_moab().get_entities_by_type(it->meshset,tris,true); CHKERR_PETSC(rval);
	blockData.push_back(new Asphalt(tris));
      }
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->get_Cubit_name().compare(0,9,"CONCRETE") == 0) {
	Range tris;
        rval = mField.get_moab().get_entities_by_type(it->meshset,tris,true); CHKERR_PETSC(rval);
	blockData.push_back(new Asphalt(tris));
      }
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->get_Cubit_name().compare(0,9,"BARESOIL") == 0) {
	Range tris;
        rval = mField.get_moab().get_entities_by_type(it->meshset,tris,true); CHKERR_PETSC(rval);
	blockData.push_back(new Asphalt(tris));
      }
    }

    PetscFunctionReturn(0);
  }

  struct Op:public TriElementForcesAndSurcesCore::UserDataOperator {
        
    ThermalElement::CommonData &commonData; 
    TmieDependendData &timeData;
    Parameters &pArameters;
    int tAg;
    bool ho_geometry;

    OpRadiationRhs(
      const string field_name,
      TmieDependendData &time_data,
      Parameters &parameters,
      ThermalElement::CommonData &common_data,
      int tag,bool _ho_geometry = false):
	TriElementForcesAndSurcesCore::UserDataOperator(field_name),
	commonData(common_data),
	timeData(time_data),
	pArameters(parameters),
	tAg(tag), ho_geometry(_ho_geometry) { symm = false; }

    PetscErrorCode record(int gg,double &f) {
      PetscFunctionBegin;

      // sigma = 5.67037321×10−8 (J/s) m−2 K−4 
      // sigma = 8165.3×10−8 (J/day) m−2 K−4
      // sigma = 0.081653 (kJ/day) m−2 K−4
      const double sigma = 0.081653
      
      trace_on(tAg);
      {
	adouble T <<= commonData.temperatureAtGaussPts[gg];
	adouble Tk = T+273.15;

	adouble Tk4 = pow(Tk,4);
	adouble sigma_eps = pArameters.eps*pArameters.sigma;

	adouble hlo = sigma_eps*Tk4; // outgoing longwave radiation (W/m2)
	adouble hli = sigma_eps*(CR+0.67*(1-CR)*pow(timeData.ea,0.08))*pow((timeData.Ta+273.15),4); // incoming longwave radiation (W/m2)
	adouble hs = (1-pArameters.alpha)*timeData.Rs; // net solar radiation (W/m2)

	adouble phis = timeData.calculateAbsoluteVirtualTempertaure(T,timeData.P); // surface virtual temperature
	adouble delta_phi = phi_s - timeData.phia;
	adouble us = pArameters.CSh*timeData.u10; // win speed with sheltering coeeficient
	adouble hconv = pArameters.rhoCp*(pArameters.Cfc*pArameters.us+pArameters.Cnc*pow(delta_phi,0.33))*(T-timeData.Ta);

	adouble hrad = hs+hli-hlo;
	adouble hnet = hrad-hconv;

	hnet >>= f;
      }
      trace_off();

      PetscFunctionReturn(0);
    }
  
    ublas::vector<FieldData> Nf;
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
  
      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
  
      PetscErrorCode ierr;
  
      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();
      int nb_row_dofs = data.getIndices().size()/rank;
  
      Nf.resize(data.getIndices().size());
      Nf.clear();
  
      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
  
	double val = getGaussPts()(2,gg);

	if(ho_geometry) {
          val *= norm_2(getNormals_at_GaussPt(gg));
        } else {
          val* = getArea();
        }

	double f;
	if(gg == 0) {
	  ierr = record(gg,f); CHKERRQ(ierr);
	} else {
	  double T = commonData.temperatureAtGaussPts[gg];
	  int r;
	  r = function(tAg,1,1,&T,&f);
	  if(r!=3) { // function is locally analytic
	    SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
	  }
	}
  
        ublas::noalias(Nf) += val*f*data.getN(gg,nb_row_dofs);
  
      }
  
      ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
  
      PetscFunctionReturn(0);
    }

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

        int nb_row = row_data.getN().size2();
        int nb_col = col_data.getN().size2();

        M.resize(nb_row,nb_col);
	bzero(&*M.data().begin(),nb_row*nb_col*sizeof(double));

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  double T = commonData.temperatureAtGaussPts[gg];
	  double Tk = T+273.15;
	  double Tk4 = pow(Tk,4);
	  double Tk3 = pow(Tk,3);

	  double val = getGaussPts()(2,gg);

	  if(ho_geometry) {
	    val *= norm_2(getNormals_at_GaussPt(gg));
	  } else {
	    val* = getArea();
	  }

	  double f;
	  if(gg == 0) {
	    ierr = record(gg,f); CHKERRQ(ierr);
	  }

	  double grad[1][1];
	  //play recorder for jacobians
	  int r;
	  r = jacobian(tAg,1,1,&T,grad);
	  if(r!=3) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	  }

	  noalias(M) += val*grad[0][0]*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col) );
        }

	ierr = MatSetValues(
	  (getFEMethod()->ts_B),
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &M(0,0),ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }
  
  };
  
};

#endif //__GROUNDSURFACETEMERATURE_HPP

