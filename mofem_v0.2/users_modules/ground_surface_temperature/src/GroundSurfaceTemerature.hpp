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

  FieldInterface &mField;

  /** \brief common data used by volume elements
    * \infroup mofem_thermal_elem
    */
  struct CommonData {
    ublas::vector<double> temperatureAtGaussPts;
  };
  CommonData commonData;

  /** \brief define surface element
    *
    * This element is used to integrate heat fluxes; convection and radiation
    */
  struct MyTriFE: public TriElementForcesAndSurcesCore {
    MyTriFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
    int getRule(int order) { return order; };
  };
  MyTriFE feGroundSurfaceRhs; //< radiation element
  MyTriFE feGroundSurfaceLhs;
  MyTriFE& getFeGroundSurfaceRhs() { return feGroundSurfaceRhs; }
  MyTriFE& getFeGroundSurfaceLhs() { return feGroundSurfaceLhs; }

  GroundSurfaceTemerature(FieldInterface &m_field):
    mField(m_field),
    feGroundSurfaceRhs(m_field),
    feGroundSurfaceLhs(m_field) {};

  struct Parameters {	

    double alpha; 	//< Solar albedo
    double d;		//< Constatnt used to clulate albedo for urtain angle
    double Cfc; 	//< Surface heat/moisture transfer coefficient for forced convection
    double Cnc;		//< Coefficient for natural convection
    double CSh;		//< Wind sheltering coefficient
    double eps;		//< Pavement emissivity
    double rhoCp;	//< Density specific heat pavement (J/m3/°C)
    Range tRis;		//< Triangles on which parameters are defined

    Parameters(Range &tris): 
      tRis(tris) {}

  };

  struct Asphalt: public Parameters {
    Asphalt(Range &tris): Parameters(tris) {
      alpha = 0.12;
      d = 0.25; // not estimated, some goods given number
      Cfc = 0.0015;
      Cnc = 0.0015;
      CSh = 1.;
      eps = 0.94;
      rhoCp = 2.0e06;
    }
  };

  struct Concrete: public Parameters {
    Concrete(Range &tris): Parameters(tris) {
      alpha = 0.20;
      d = 0.25; // not estimated, some goods given number
      Cfc = 0.0015;
      Cnc = 0.0015;
      CSh = 1.;
      eps = 0.94;
      rhoCp = 2.0e06;
    }
  };

  struct BareSoil: public Parameters {
    BareSoil(Range &tris): Parameters(tris) {
      alpha = 0.15;
      d = 0.25; // not estimated, some goods given number
      Cfc = 0.003;
      Cnc = 0.0015;
      CSh = 1.;
      eps = 0.94;
      rhoCp = 2.0e06;
    }
  };

  PetscErrorCode addSurfaces(const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
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
      if(it->get_Cubit_name().compare(0,7,"ASPHALT") == 0) {
	Range tris;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
	blockData.push_back(new Asphalt(tris));
	ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"GROUND_SURFACE_FE"); CHKERRQ(ierr);
      }
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->get_Cubit_name().compare(0,8,"CONCRETE") == 0) {
	Range tris;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
	blockData.push_back(new Concrete(tris));
	ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"GROUND_SURFACE_FE"); CHKERRQ(ierr);
      }
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->get_Cubit_name().compare(0,8,"BARESOIL") == 0) {
	Range tris;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
	blockData.push_back(new BareSoil(tris));
	ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"GROUND_SURFACE_FE"); CHKERRQ(ierr);
      }
    }

    PetscFunctionReturn(0);
  }

  static double netSolarRadiation(double alpha,double d,double cos_omega,GenricClimateModel *time_data_ptr) {
    // Parameterizing the Dependence of Surface Albedo on Solar Zenith Angle Using
    // Atmospheric Radiation Measurement Program Observations
    // F. Yang
    // Environmental Modeling Center National Centers for Environmental Prediction Camp Springs, Maryland
    alpha = alpha*(1+d/(1+2*d*cos_omega)); 
    return (1-alpha)*time_data_ptr->Rs; // net solar radiation (W/m2)
  }

  static double incomingLongWaveRadiation(double eps,GenricClimateModel *time_data_ptr) {
    const double sigma = 5.67037321e-8;
    double sigma_eps = eps*sigma;
    double ea = time_data_ptr->calulateVapourPressure(time_data_ptr->calulateVapourPressure(time_data_ptr->Td));
    // incoming longwave radiation (W/m2)
    return sigma_eps*(time_data_ptr->CR+0.67*(1-time_data_ptr->CR)*pow(ea,0.08))*pow((time_data_ptr->Ta+273.15),4); 
  }

  struct PreProcess: public MoFEM::FEMethod {

    GenricClimateModel *timeDataPtr;
    PreProcess(GenricClimateModel *time_data_ptr):
      timeDataPtr(time_data_ptr) {};
      
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ierr = timeDataPtr->set(ts_t); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }


  };

  struct SolarRadiationPreProcessor: public MoFEM::FEMethod {

    FieldInterface &mField;
    GenricClimateModel *timeDataPtr;
    Parameters *pArametersPtr;
    AdaptiveKDTree kdTree;
    double ePs;
    bool iNit;

    SolarRadiationPreProcessor(
      FieldInterface &m_field,
      GenricClimateModel *time_data_ptr,
      Parameters *parameters_ptr,
      double eps = 1e-6):
      mField(m_field),
      timeDataPtr(time_data_ptr),
      pArametersPtr(parameters_ptr),
      kdTree(&m_field.get_moab()),
      ePs(eps),iNit(false) {
      azimuth = zenith = 0;
    };
    ~SolarRadiationPreProcessor() {
      if(kdTree_rootMeshset) {
	mField.get_moab().delete_entities(&kdTree_rootMeshset,1);
      }
    }
 
    Range sKin,skinNodes;
    EntityHandle kdTree_rootMeshset;

    PetscErrorCode getSkin(Range &tets) {
      PetscFunctionBegin;
      ErrorCode rval;
      PetscErrorCode ierr;
      Skinner skin(&mField.get_moab());
      rval = skin.find_skin(0,tets,false,sKin); CHKERR(rval);
      rval = mField.get_moab().create_meshset(MESHSET_SET,kdTree_rootMeshset); CHKERR_PETSC(rval);
      rval = kdTree.build_tree(sKin,&kdTree_rootMeshset); CHKERR_PETSC(rval);
      PetscFunctionReturn(0);
    }
    
    double azimuth,zenith;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
  
      ErrorCode rval;
      PetscErrorCode ierr;
	
      int def_VAL = 0;
      Tag th_solar_exposure;
      rval = mField.get_moab().
	tag_get_handle("SOLAR_EXPOSURE",1,MB_TYPE_INTEGER,th_solar_exposure,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); 

      double zero[3] = {0,0,0};
      Tag th_solar_radiation;
      rval = mField.get_moab().
	tag_get_handle("SOLAR_RADIATION",1,MB_TYPE_DOUBLE,th_solar_radiation,MB_TAG_CREAT|MB_TAG_SPARSE,zero); 

      Tag th_ray_direction;
      rval = mField.get_moab().
	tag_get_handle("SUN_RAY",3,MB_TYPE_DOUBLE,th_ray_direction,MB_TAG_CREAT|MB_TAG_SPARSE,zero); 

      if(!iNit) {
	rval = mField.get_moab().get_connectivity(pArametersPtr->tRis,skinNodes,true); CHKERR_PETSC(rval);
      }

      if(iNit) {
	if(azimuth == timeDataPtr->azimuth && zenith == timeDataPtr->zenith) {
	  PetscFunctionReturn(0);
	}
      }
      iNit = true;

      azimuth = timeDataPtr->azimuth;
      zenith = timeDataPtr->zenith;
      //assume that X pointing to North
      double ray_unit_dir[] = {
	cos(azimuth*M_PI/180)*sin(zenith*M_PI/180), 
	sin(azimuth*M_PI/180)*sin(zenith*M_PI/180), 
	cos(zenith*M_PI/180) };

      vector<EntityHandle> triangles_out;
      vector<double> distance_out;

      double diffN[6];
      ierr = ShapeDiffMBTRI(diffN); CHKERRQ(ierr);
      Range::iterator tit = pArametersPtr->tRis.begin();
      for(;tit!=pArametersPtr->tRis.end();tit++) {

	if(ray_unit_dir[2]<=0) { 
	  rval = mField.get_moab().tag_set_data(th_solar_radiation,&*tit,1,zero); CHKERR_PETSC(rval);
	  continue;
	}

	int num_nodes;
        const EntityHandle* conn;
	rval = mField.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
	double coords[9]; 
	rval = mField.get_moab().get_coords(conn,3,coords); CHKERR_PETSC(rval);

	double normal[3];
	ierr = ShapeFaceNormalMBTRI(diffN,coords,normal); CHKERRQ(ierr);

	for(int nn = 1;nn<3;nn++) {
	  for(int dd = 0;dd<3;dd++) {
	    coords[dd] += coords[3*nn+dd];
	  }
	}
	for(int dd = 0;dd<3;dd++) {
	  coords[dd] /= 3;
	}

	triangles_out.resize(0);
	distance_out.resize(0);
	rval = kdTree.ray_intersect_triangles(kdTree_rootMeshset,
	    1e-12,
	    ray_unit_dir,coords,
	    triangles_out,
	    distance_out); CHKERR(rval);

	double exposed = 0;
	if(triangles_out.size()>0) {
	  for(int nn = 0;nn<triangles_out.size();nn++) {
	    if(exposed<distance_out[nn]) exposed = distance_out[nn];
	  }
	}

	double hsol;
	if(exposed>ePs) {
	  hsol = 0;
	} else {
	  double cos_phi = 0;
	  for(int nn = 0;nn<3;nn++) {
	    cos_phi += normal[nn]*ray_unit_dir[nn];
	  }
	  cos_phi /= sqrt(pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2));
	  cos_phi = fabs(cos_phi);

	  hsol = netSolarRadiation(pArametersPtr->alpha,pArametersPtr->d,cos_phi,timeDataPtr);	
	}
	
	rval = mField.get_moab().tag_set_data(th_solar_radiation,&*tit,1,&hsol); CHKERR_PETSC(rval);

      }

      Range::iterator nit = skinNodes.begin();
      for(;nit!=skinNodes.end();nit++) {

	if(ray_unit_dir[2]<=0) { 
	  int set = 0;
	  rval = mField.get_moab().tag_set_data(th_solar_exposure,&*nit,1,&set); CHKERR_PETSC(rval);
	  rval = mField.get_moab().tag_set_data(th_ray_direction,&*nit,1,zero); CHKERR_PETSC(rval);
	  continue;
	}
    
	double coords[3]; 
	rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);

	triangles_out.resize(0);
	distance_out.resize(0);
	rval = kdTree.ray_intersect_triangles(kdTree_rootMeshset,
	    1e-12,
	    ray_unit_dir,coords,
	    triangles_out,
	    distance_out); CHKERR(rval);

	double exposed = 0;
	if(triangles_out.size()>0) {
	  for(int nn = 0;nn<triangles_out.size();nn++) {
	    if(exposed<distance_out[nn]) exposed = distance_out[nn];
	  }
	}

	int set = 1;
	if(exposed>ePs) {
	  set = 0;
	}
	rval = mField.get_moab().tag_set_data(th_solar_exposure,&*nit,1,&set); CHKERR_PETSC(rval);
	rval = mField.get_moab().tag_set_data(th_ray_direction,&*nit,1,ray_unit_dir); CHKERR_PETSC(rval);

      }

      PetscFunctionReturn(0);
    }

  };

  boost::ptr_vector<Parameters> blockData;
  boost::ptr_vector<SolarRadiationPreProcessor> preProcessShade;

  /** \brief opearator to caulate tempereature at Gauss points
    * \infroup mofem_thermal_elem
    */
  struct OpGetTriTemperatureAtGaussPts: public TriElementForcesAndSurcesCore::UserDataOperator {

    ublas::vector<double> &fieldAtGaussPts;
    OpGetTriTemperatureAtGaussPts(const string field_name,ublas::vector<double> &field_at_gauss_pts):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
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

  struct OpRhs:public TriElementForcesAndSurcesCore::UserDataOperator {
        
    CommonData &commonData; 
    GenricClimateModel* timeDataPtr;
    Parameters *pArametersPtr;
    int tAg;
    bool ho_geometry;

    OpRhs(
      const string field_name,
      GenricClimateModel *time_data_ptr,
      Parameters *parameters_ptr,
      CommonData &common_data,
      int tag,bool _ho_geometry = false):
	TriElementForcesAndSurcesCore::UserDataOperator(field_name),
	commonData(common_data),
	timeDataPtr(time_data_ptr),
	pArametersPtr(parameters_ptr),
	tAg(tag), ho_geometry(_ho_geometry) {

    }

    int nodalExposure[4],eXposure;
    PetscErrorCode getExposure() {
      PetscFunctionBegin;
      ErrorCode rval;
      PetscErrorCode ierr;
      EntityHandle element_ent = getMoFEMFEPtr()->get_ent();
      const EntityHandle *conn;
      int num_nodes;
      rval = getTriElementForcesAndSurcesCore()->mField.get_moab().get_connectivity(element_ent,conn,num_nodes,true); CHKERR_PETSC(rval);
      int def_VAL = 0;
      Tag th_solar_exposure;
      rval = getTriElementForcesAndSurcesCore()->mField.get_moab().tag_get_handle( 
	"SOLAR_EXPOSURE",1,MB_TYPE_INTEGER,th_solar_exposure,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
      ierr = getTriElementForcesAndSurcesCore()->mField.get_moab().tag_get_data(th_solar_exposure,conn,num_nodes,nodalExposure); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    double T;
    adouble aHlo,aHconv,aHevap,aHrad,aHnet,aDeltaPhi,aT,aTk,aTk4;;
    PetscErrorCode record(double &f) {
      PetscFunctionBegin;

      // sigma = 5.67037321×10−8 (J/s) m−2 K−4 
      // sigma = 8165.3×10−8 (J/day) m−2 K−4
      // sigma = 0.081653 (kJ/day) m−2 K−4
      // STephan-Boltzman
      const double sigma = 5.67037321e-8;
      
      trace_on(tAg);
      {
	aT <<= T;

	aTk = aT+273.15;
	aTk4 = pow(aTk,4);
	
	//radiation
	double sigma_eps = pArametersPtr->eps*sigma;
	aHlo = sigma_eps*aTk4; // outgoing longwave radiation (W/m2)

	//convection
	double us = pArametersPtr->CSh*timeDataPtr->u10; // win speed with sheltering coeeficient
	aDeltaPhi = fmax(1e-12,aT-timeDataPtr->Ta); // this is with assumtion that near the near the surface is the same amout of moisture like in the air 
						    // if ground temerature are lower than air there is no convection

	aHconv = pArametersPtr->rhoCp*(pArametersPtr->Cfc*us+pArametersPtr->Cnc*pow(aDeltaPhi,0.33))*(aT-timeDataPtr->Ta);
	//aHevap = 0; // need to be implemented with moisture model

	aHrad = -aHlo;
	aHnet = 0;//aHrad;//-aHconv;//-aHevap;

	aHnet >>= f;
      }
      trace_off();

      PetscFunctionReturn(0);
    }
  
    ublas::vector<FieldData> Nf;
    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
  
      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(pArametersPtr->tRis.find(getMoFEMFEPtr()->get_ent())==pArametersPtr->tRis.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      if(type == MBVERTEX) {
	ierr = getExposure(); CHKERRQ(ierr);
      }
  
      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();
      int nb_row_dofs = data.getIndices().size()/rank;
  
      Nf.resize(data.getIndices().size());
      Nf.clear();
  
      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = getGaussPts()(2,gg);

	ublas::vector<double> normal;
	if(ho_geometry) {
	  normal = getNormals_at_GaussPt(gg);
	} else {
	  normal = getNormal();
	}
	val *= norm_2(normal)*0.5;
	
	T = commonData.temperatureAtGaussPts[gg];

	double azimuth = timeDataPtr->azimuth;
	double zenith = timeDataPtr->zenith;
	//assume that X pointing to North
	double ray_unit_dir[] = {
	  cos(azimuth*M_PI/180)*sin(zenith*M_PI/180), 
	  sin(azimuth*M_PI/180)*sin(zenith*M_PI/180), 
	  cos(zenith*M_PI/180) };
	double cos_phi = 0;
	for(int nn = 0;nn<3;nn++) {
	  cos_phi += normal[nn]*ray_unit_dir[nn];
	}
	cos_phi /= norm_2(normal);
	
	eXposure = 0;
	for(int nn = 0;nn<3;nn++) {
	  eXposure += getTriElementForcesAndSurcesCore()->dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,nn)*nodalExposure[nn];
	}

	double hnet  = 0;
	/*if(gg == 0) {
	  ierr = record(hnet); CHKERRQ(ierr);
	} else {
	  int r;
	  r = function(tAg,1,1,&T,&hnet);
	  if(r<2) { // function is locally analytic
	    SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
	  }
	}*/
	
	if(eXposure>0) {
	  //hnet += netSolarRadiation(pArametersPtr->alpha,pArametersPtr->d,cos_phi,timeDataPtr);
	}	
	hnet += incomingLongWaveRadiation(pArametersPtr->eps,timeDataPtr);
	hnet /= (double)86400; // number of second in the day
  
        //ublas::noalias(Nf) -= val*hnet*data.getN(gg,nb_row_dofs);
  
      }
  
      ierr = VecSetValues(getFEMethod()->ts_F,data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
  
      PetscFunctionReturn(0);
    }


  
  };

  struct OpLhs:public OpRhs {
        
    OpLhs(
      const string field_name,
      GenricClimateModel *time_data_ptr,
      Parameters *parameters_ptr,
      CommonData &common_data,
      int tag,bool _ho_geometry = false):
	OpRhs(field_name,time_data_ptr,parameters_ptr,common_data,tag,_ho_geometry) {}

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    ublas::matrix<double> NN,transNN;

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

        NN.resize(nb_row,nb_col);
	NN.clear();

        for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  double val = getGaussPts()(2,gg);
	  ublas::vector<double> normal;
	  if(ho_geometry) {
	    normal = getNormals_at_GaussPt(gg);
	  } else {
	    normal = getNormal();
	  }
	  val *= norm_2(normal)*0.5;

	  /*double hnet;
	  if(gg == 0) {
	    ierr = record(hnet); CHKERRQ(ierr);
	  }

	  double grad[1];
	  double* grad_ptr[] = { grad };
	  //play recorder for jacobians
	  int r;
	  r = jacobian(tAg,1,1,&T,grad_ptr);
	  if(r<2) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	  }
	  hnet /= (double)86400; // number of second in the day

	  noalias(NN) -= val*(grad[0]*outer_prod( row_data.getN(gg,nb_row),col_data.getN(gg,nb_col)));*/

        }

	ierr = MatSetValues(
	  (getFEMethod()->ts_B),
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &NN(0,0),ADD_VALUES); CHKERRQ(ierr);
        if(row_side != col_side || row_type != col_type) {
          transNN.resize(nb_col,nb_row);
          noalias(transNN) = trans( NN );
          ierr = MatSetValues(
	    (getFEMethod()->ts_B),
            nb_col,&col_data.getIndices()[0],
            nb_row,&row_data.getIndices()[0],
            &transNN(0,0),ADD_VALUES); CHKERRQ(ierr);
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode setOperators(int tag,
    GenricClimateModel *time_data_ptr,string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    bool ho_geometry = false;
    if(mField.check_field(mesh_nodals_positions)) {
      ho_geometry = true;
    }

    {
      boost::ptr_vector<Parameters>::iterator sit = blockData.begin();
      for(;sit!=blockData.end();sit++) {
	// add finite element operator
	feGroundSurfaceRhs.get_op_to_do_Rhs().push_back(new OpGetTriTemperatureAtGaussPts(field_name,commonData.temperatureAtGaussPts));
	feGroundSurfaceRhs.get_op_to_do_Rhs().push_back(new OpRhs(field_name,time_data_ptr,&*sit,commonData,tag,ho_geometry));
	preProcessShade.push_back(new SolarRadiationPreProcessor(mField,time_data_ptr,&*sit));
      }
    }
    {
      boost::ptr_vector<Parameters>::iterator sit = blockData.begin();
      for(;sit!=blockData.end();sit++) {
	// add finite element operator
	feGroundSurfaceRhs.get_op_to_do_Rhs().push_back(new OpGetTriTemperatureAtGaussPts(field_name,commonData.temperatureAtGaussPts));
	feGroundSurfaceLhs.get_op_to_do_Lhs().push_back(new OpLhs(field_name,time_data_ptr,&*sit,commonData,tag,ho_geometry));
      }
    }

    PetscFunctionReturn(0);
  }


};

#endif //__GROUNDSURFACETEMERATURE_HPP

