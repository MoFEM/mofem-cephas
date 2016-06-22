/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <petsctime.h>

#include <SurfacePressure.hpp>
#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>
#include "ElasticFEMethodTransIso.hpp"

#include "ElasticFE_RVELagrange_Disp.hpp"
#include "ElasticFE_RVELagrange_Disp_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp"
#include "RVEVolume.hpp"

using namespace ObosleteUsersModules;
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/extreme_value_distribution.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/gamma_distribution.hpp>

#include <MCS_RNG.hpp>

#include <FE2_ElasticFEMethod.hpp>
#include <Reliability_SurfacePressure.hpp>
#include <FE2_Macro_Solver_MCS.hpp>
#include <FE2_PostProcStressForReliability.hpp>

#include <ImportProbData.hpp>
#include <LimitStateFunction.hpp>

#include <boost/math/distributions/normal.hpp>

using namespace boost::numeric;


//======================================================

/*******************************************************************************
 *                                                                             *
 *             FUNCTIONS FOR STRUCTURAL RELIABILITY ANLYSIS                    *
 *                                                                             *
 /*******************************************************************************/

//------------------------------------------------------------------------------
// Construct data structure <Stochastic_Model> for collecting data representing
//   statistical information of inputs including probability distribution,
//   correlation matrix and etc.

struct Stochastic_Model {
  int    num_vars;                       // Number of variables
  double dist_type;                      // Distribution type index
  double transf_type;                    // Type of joint distribution
  double R0_method;                      // Method for computation of the modified Nataf correlation matrix
  int flag_sens;                         // Flag for computation of sensitivities w.r.t. parameters
  int ExaminedPly;
  vector<string> NameVars;               // Name of random variables
  ublas::matrix<double> correlation;     // Correlation matrix
  ublas::matrix<double> marg;            // Marginal distribution for each random variable
  ublas::matrix<double> mod_correlation; // modified correlation matrix
  ublas::matrix<double> Lo;              // Chelosky decomposition
  ublas::matrix<double> inv_Lo;          // inverse of matrix Lo
  ublas::vector<double> MatStrength;     // Material strength
  ublas::vector<double> PlyAngle;        // Angle of orientation
};


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

//const double young_modulus = 1;
//const double poisson_ratio = 0.0;

//------------------------------------------------------------------------------
// To transform stress
//

void REL_Stress_Transformation(double theta, ublas::matrix<double> Stress_xyz, ublas::matrix<double> &Stress_123) {
  
  double l1, l2, l3, m1, m2, m3, n1, n2, n3;
  l1=cos(theta);
  m1=sin(theta);
  n1=0.0;
  l2=-sin(theta);
  m2=cos(theta);
  n2=0.0;
  l3=0.0;
  m3=0.0;
  n3=1.0;
  
  ublas::vector<double> Stress_VectorNotation_xyz(6);
  Stress_VectorNotation_xyz(0) = Stress_xyz(0,0);
  Stress_VectorNotation_xyz(1) = Stress_xyz(1,1);
  Stress_VectorNotation_xyz(2) = Stress_xyz(2,2);
  Stress_VectorNotation_xyz(3) = Stress_xyz(0,1);
  Stress_VectorNotation_xyz(4) = Stress_xyz(1,2);
  Stress_VectorNotation_xyz(5) = Stress_xyz(2,0);
  
  ublas::matrix<FieldData> T_strain_inv_T;    T_strain_inv_T.resize(6,6);
  T_strain_inv_T(0,0) =   l1*l1;
  T_strain_inv_T(0,1) =   m1*m1;
  T_strain_inv_T(0,2) =   n1*n1;
  T_strain_inv_T(0,3) = 2*l1*m1;
  T_strain_inv_T(0,4) = 2*m1*n1;
  T_strain_inv_T(0,5) = 2*l1*n1;
  
  T_strain_inv_T(1,0) =   l2*l2;
  T_strain_inv_T(1,1) =   m2*m2;
  T_strain_inv_T(1,2) =   n2*n2;
  T_strain_inv_T(1,3) = 2*l2*m2;
  T_strain_inv_T(1,4) = 2*m2*n2;
  T_strain_inv_T(1,5) = 2*l2*n2;
  
  T_strain_inv_T(2,0) =   l3*l3;
  T_strain_inv_T(2,1) =   m3*m3;
  T_strain_inv_T(2,2) =   n3*n3;
  T_strain_inv_T(2,3) = 2*l3*m3;
  T_strain_inv_T(2,4) = 2*m3*n3;
  T_strain_inv_T(2,5) = 2*l3*n3;
  
  T_strain_inv_T(3,0) = l1*l2;
  T_strain_inv_T(3,1) = m1*m2;
  T_strain_inv_T(3,2) = n1*n2;
  T_strain_inv_T(3,3) = l1*m2+m1*l2;
  T_strain_inv_T(3,4) = m1*n2+n1*m2;
  T_strain_inv_T(3,5) = l1*n2+n1*l2;
  
  T_strain_inv_T(4,0) = l2*l3;
  T_strain_inv_T(4,1) = m2*m3;
  T_strain_inv_T(4,2) = n2*n3;
  T_strain_inv_T(4,3) = l2*m3+m2*l3;
  T_strain_inv_T(4,4) = m2*n3+n2*m3;
  T_strain_inv_T(4,5) = l2*n3+n2*l3;
  
  T_strain_inv_T(5,0) = l1*l3;
  T_strain_inv_T(5,1) = m1*m3;
  T_strain_inv_T(5,2) = n1*n3;
  T_strain_inv_T(5,3) = l1*m3+m1*l3;
  T_strain_inv_T(5,4) = m1*n3+n1*m3;
  T_strain_inv_T(5,5) = l1*n3+n1*l3;
  
  ublas::vector<double> Stress_VectorNotation_123 = prod(T_strain_inv_T, Stress_VectorNotation_xyz);
  
  Stress_123.resize(3,3); Stress_123.clear();
  Stress_123(0,0) = Stress_VectorNotation_123(0);
  Stress_123(1,1) = Stress_VectorNotation_123(1);
  Stress_123(2,2) = Stress_VectorNotation_123(2);
  Stress_123(0,1) = Stress_123(1,0) = Stress_VectorNotation_123(3);
  Stress_123(1,2) = Stress_123(2,1) = Stress_VectorNotation_123(4);
  Stress_123(2,0) = Stress_123(0,2) = Stress_VectorNotation_123(5);
  
}

void REL_Stress_Transformation_Theta(double theta, ublas::matrix<double> Stress_xyz,
                                     ublas::matrix<double> Stress_xyz_r,
                                     ublas::matrix<double> &Stress_123_r) {
  
  double l1, l2, l3, m1, m2, m3, n1, n2, n3;
  double l1_r, l2_r, l3_r, m1_r, m2_r, m3_r, n1_r, n2_r, n3_r;
  l1   =  cos(theta); l1_r = -sin(theta);
  m1   =  sin(theta); m1_r =  cos(theta);
  n1   =  0.0;        n1_r =  0.0;
  l2   = -sin(theta); l2_r = -cos(theta);
  m2   =  cos(theta); m2_r = -sin(theta);
  n2   =  0.0;        n2_r =  0.0;
  l3   =  0.0;        l3_r =  0.0;
  m3   =  0.0;        m3_r =  0.0;
  n3   =  1.0;        n3_r =  0.0;
  
  ublas::vector<double> Stress_VectorNotation_xyz(6);
  Stress_VectorNotation_xyz(0) = Stress_xyz(0,0);
  Stress_VectorNotation_xyz(1) = Stress_xyz(1,1);
  Stress_VectorNotation_xyz(2) = Stress_xyz(2,2);
  Stress_VectorNotation_xyz(3) = Stress_xyz(0,1);
  Stress_VectorNotation_xyz(4) = Stress_xyz(1,2);
  Stress_VectorNotation_xyz(5) = Stress_xyz(2,0);
  
  ublas::vector<double> Stress_VectorNotation_xyz_r(6);
  Stress_VectorNotation_xyz_r(0) = Stress_xyz_r(0,0);
  Stress_VectorNotation_xyz_r(1) = Stress_xyz_r(1,1);
  Stress_VectorNotation_xyz_r(2) = Stress_xyz_r(2,2);
  Stress_VectorNotation_xyz_r(3) = Stress_xyz_r(0,1);
  Stress_VectorNotation_xyz_r(4) = Stress_xyz_r(1,2);
  Stress_VectorNotation_xyz_r(5) = Stress_xyz_r(2,0);
  
  ublas::matrix<FieldData> T_strain_inv_T;    T_strain_inv_T.resize(6,6);
  T_strain_inv_T(0,0) =   l1*l1;
  T_strain_inv_T(0,1) =   m1*m1;
  T_strain_inv_T(0,2) =   n1*n1;
  T_strain_inv_T(0,3) = 2*l1*m1;
  T_strain_inv_T(0,4) = 2*m1*n1;
  T_strain_inv_T(0,5) = 2*l1*n1;
  
  T_strain_inv_T(1,0) =   l2*l2;
  T_strain_inv_T(1,1) =   m2*m2;
  T_strain_inv_T(1,2) =   n2*n2;
  T_strain_inv_T(1,3) = 2*l2*m2;
  T_strain_inv_T(1,4) = 2*m2*n2;
  T_strain_inv_T(1,5) = 2*l2*n2;
  
  T_strain_inv_T(2,0) =   l3*l3;
  T_strain_inv_T(2,1) =   m3*m3;
  T_strain_inv_T(2,2) =   n3*n3;
  T_strain_inv_T(2,3) = 2*l3*m3;
  T_strain_inv_T(2,4) = 2*m3*n3;
  T_strain_inv_T(2,5) = 2*l3*n3;
  
  T_strain_inv_T(3,0) = l1*l2;
  T_strain_inv_T(3,1) = m1*m2;
  T_strain_inv_T(3,2) = n1*n2;
  T_strain_inv_T(3,3) = l1*m2+m1*l2;
  T_strain_inv_T(3,4) = m1*n2+n1*m2;
  T_strain_inv_T(3,5) = l1*n2+n1*l2;
  
  T_strain_inv_T(4,0) = l2*l3;
  T_strain_inv_T(4,1) = m2*m3;
  T_strain_inv_T(4,2) = n2*n3;
  T_strain_inv_T(4,3) = l2*m3+m2*l3;
  T_strain_inv_T(4,4) = m2*n3+n2*m3;
  T_strain_inv_T(4,5) = l2*n3+n2*l3;
  
  T_strain_inv_T(5,0) = l1*l3;
  T_strain_inv_T(5,1) = m1*m3;
  T_strain_inv_T(5,2) = n1*n3;
  T_strain_inv_T(5,3) = l1*m3+m1*l3;
  T_strain_inv_T(5,4) = m1*n3+n1*m3;
  T_strain_inv_T(5,5) = l1*n3+n1*l3;
  
  ublas::matrix<FieldData> T_strain_inv_T_r_Theta;    T_strain_inv_T_r_Theta.resize(6,6);
  T_strain_inv_T_r_Theta(0,0) =   2*l1*l1_r;
  T_strain_inv_T_r_Theta(0,1) =   2*m1*m1_r;
  T_strain_inv_T_r_Theta(0,2) =   2*n1*n1_r;
  T_strain_inv_T_r_Theta(0,3) = 2*(l1_r*m1 + l1*m1_r);
  T_strain_inv_T_r_Theta(0,4) = 2*(m1_r*n1 + m1*n1_r);
  T_strain_inv_T_r_Theta(0,5) = 2*(l1_r*n1 + l1*n1_r);
  
  T_strain_inv_T_r_Theta(1,0) =   2*l2*l2_r;
  T_strain_inv_T_r_Theta(1,1) =   2*m2*m2_r;
  T_strain_inv_T_r_Theta(1,2) =   2*n2*n2_r;
  T_strain_inv_T_r_Theta(1,3) = 2*(l2_r*m2 + l2*m2_r);
  T_strain_inv_T_r_Theta(1,4) = 2*(m2_r*n2 + m2*n2_r);
  T_strain_inv_T_r_Theta(1,5) = 2*(l2_r*n2 + l2*n2_r);
  
  T_strain_inv_T_r_Theta(2,0) =   2*l3*l3_r;
  T_strain_inv_T_r_Theta(2,1) =   2*m3*m3_r;
  T_strain_inv_T_r_Theta(2,2) =   2*n3*n3_r;
  T_strain_inv_T_r_Theta(2,3) = 2*(l3_r*m3 + l3*m3_r);
  T_strain_inv_T_r_Theta(2,4) = 2*(m3_r*n3 + m3*n3_r);
  T_strain_inv_T_r_Theta(2,5) = 2*(l3_r*n3 + l3*n3_r);
  
  T_strain_inv_T_r_Theta(3,0) = l1_r*l2 + l1*l2_r;
  T_strain_inv_T_r_Theta(3,1) = m1_r*m2 + m1*m2_r;
  T_strain_inv_T_r_Theta(3,2) = n1_r*n2 + n1*n2_r;
  T_strain_inv_T_r_Theta(3,3) = l1_r*m2 + l1*m2_r + m1_r*l2 + m1*l2_r;
  T_strain_inv_T_r_Theta(3,4) = m1_r*n2 + m1*n2_r + n1_r*m2 + n1*m2_r;
  T_strain_inv_T_r_Theta(3,5) = l1_r*n2 + l1*n2_r + n1_r*l2 + n1*l2_r;
  
  T_strain_inv_T_r_Theta(4,0) = l2_r*l3 + l2*l3_r;
  T_strain_inv_T_r_Theta(4,1) = m2_r*m3 + m2*m3_r;
  T_strain_inv_T_r_Theta(4,2) = n2_r*n3 + n2*n3_r;
  T_strain_inv_T_r_Theta(4,3) = l2_r*m3 + l2*m3_r + m2_r*l3 + m2*l3_r;
  T_strain_inv_T_r_Theta(4,4) = m2_r*n3 + m2*n3_r + n2_r*m3 + n2*m3_r;
  T_strain_inv_T_r_Theta(4,5) = l2_r*n3 + l2*n3_r + n2_r*l3 + n2*l3_r;
  
  T_strain_inv_T_r_Theta(5,0) = l1_r*l3 + l1*l3_r;
  T_strain_inv_T_r_Theta(5,1) = m1_r*m3 + m1*m3_r;
  T_strain_inv_T_r_Theta(5,2) = n1_r*n3 + n1*n3_r;
  T_strain_inv_T_r_Theta(5,3) = l1_r*m3 + l1*m3_r + m1_r*l3 + m1*l3_r;
  T_strain_inv_T_r_Theta(5,4) = m1_r*n3 + m1*n3_r + n1_r*m3 + n1*m3_r;
  T_strain_inv_T_r_Theta(5,5) = l1_r*n3 + l1*n3_r + n1_r*l3 + n1*l3_r;
  
  T_strain_inv_T_r_Theta = T_strain_inv_T_r_Theta*(M_PI/180);
  
  ublas::vector<double> StressVector_r_1 = prod(T_strain_inv_T, Stress_VectorNotation_xyz_r);
  ublas::vector<double> StressVector_r_2 = prod(T_strain_inv_T_r_Theta, Stress_VectorNotation_xyz);
  ublas::vector<double> Stress_VectorNotation_123_r = StressVector_r_1 + StressVector_r_2;
  
  Stress_123_r.resize(3,3); Stress_123_r.clear();
  Stress_123_r(0,0) = Stress_VectorNotation_123_r(0);
  Stress_123_r(1,1) = Stress_VectorNotation_123_r(1);
  Stress_123_r(2,2) = Stress_VectorNotation_123_r(2);
  Stress_123_r(0,1) = Stress_123_r(1,0) = Stress_VectorNotation_123_r(3);
  Stress_123_r(1,2) = Stress_123_r(2,1) = Stress_VectorNotation_123_r(4);
  Stress_123_r(2,0) = Stress_123_r(0,2) = Stress_VectorNotation_123_r(5);
  
}

int main(int argc, char *argv[]) {
  
  clock_t start_time, finish_time;
  double total_time;
  start_time = clock();
  
  //============================================================================
  //
  //  A. Micro (RVE) Problem
  //
  //============================================================================
  
  ierr = PetscInitialize(&argc,&argv,(char *)0,help); CHKERRQ(ierr);
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  moab::Core mb_instance_RVE;
  Interface& moab_RVE = mb_instance_RVE;
  
  // ===========================================================================
  //
  //  A.I. READ MESH DATA AND FEA COMPUTATION PARAMETERS FROM FILE
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_RVE",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_RVE (MESH FILE NEEDED)");
  }
  
  PetscInt order_RVE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_RVE",&order_RVE,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_RVE = 1;
  }

  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab_RVE.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab_RVE,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab_RVE,PETSC_COMM_WORLD);
  
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_RVE(moab_RVE);
  FieldInterface& m_field_RVE = core_RVE;
  
  
  /*****************************************************************************
   *
   * Get fibre direction information from the potential-flow calculation
   *
   ****************************************************************************/
  
  Tag th_phi;
  //    double def_val  = 0;
  rval = moab_RVE.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi); CHKERR_PETSC(rval);
	
  Tag th_meshset_info;
  int def_meshset_info[2] = {0,0};
  rval = moab_RVE.tag_get_handle("MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info);
  
  int meshset_data[2];
  EntityHandle root = moab_RVE.get_root_set();
  rval = moab_RVE.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);
  
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(meshset_data[0]-1));
  
  ierr = m_field_RVE.build_fields(); CHKERRQ(ierr);
  ierr = m_field_RVE.build_finite_elements(); CHKERRQ(ierr);
  ierr = m_field_RVE.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
  ierr = m_field_RVE.build_problems(); CHKERRQ(ierr);

  /*****************************************************************************
   *
   * Group element into various mesh-set
   * meshset_level0: all element
   * meshset_Matrix: element set for matrix
   * meshset_Inclusion: element set for inclusion/fibre
   *
   ****************************************************************************/
  
  EntityHandle out_meshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  //    ierr = m_field_RVE.get_problem_finite_elements_entities("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
  ierr = m_field_RVE.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  Range LatestRefinedTets;
  rval = moab_RVE.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);
  
  Range LatestRefinedPrisms;
  rval = moab_RVE.get_entities_by_type(out_meshset, MBPRISM,LatestRefinedPrisms,true); CHKERR_PETSC(rval);
	
  cout<<"No of Prisms/Interfaces = "<<LatestRefinedPrisms.size()<<endl;
	
  BitRefLevel problem_bit_level_RVE = bit_levels.back();
  
  EntityHandle meshset_Elastic, meshset_Trans_ISO;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_Trans_ISO); CHKERR_PETSC(rval);
	
	///Getting No. of Fibres to be used for Potential Flow Problem
	int noOfFibres=0;
	for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET|UNKNOWNCUBITNAME,it)) {
		
		std::size_t found=it->get_name().find("PotentialFlow");
		if (found==std::string::npos) continue;
		noOfFibres += 1;
	}
	cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;
	
//	vector<int> fibreList(noOfFibres,0);
//	for (int aa=0; aa<noOfFibres; aa++) {
//		fibreList[aa] = aa + 1;
//	}
  
	vector<Range> RangeFibre(noOfFibres);
	vector<EntityHandle> fibre_meshset(noOfFibres);
	
	for (int ii=0; ii<noOfFibres; ii++) {
		ostringstream sss;
		sss << "POTENTIAL_ELEM" << ii+1;
		for(_IT_GET_FES_BY_NAME_FOR_LOOP_(m_field_RVE, sss.str().c_str() ,it)){
			RangeFibre[ii].insert(it->get_ent());
			rval = moab_RVE.create_meshset(MESHSET_SET,fibre_meshset[ii]); CHKERR_PETSC(rval);
			rval = moab_RVE.add_entities(fibre_meshset[ii],RangeFibre[ii]); CHKERR_PETSC(rval);
			rval = moab_RVE.unite_meshset(meshset_Trans_ISO,fibre_meshset[ii]); CHKERR_PETSC(rval);
		}
	}
  
  rval = moab_RVE.write_file("meshset_Trans_ISO.vtk","VTK","",&meshset_Trans_ISO,1); CHKERR_PETSC(rval);
  
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)){
    
		if(it->get_name() == "MAT_ELASTIC_1") {
			Range TetsInBlock;
			rval = moab_RVE.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
			Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);
      
      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;
      
			rval = moab_RVE.add_entities(meshset_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);
      
		}
	}
  ierr = m_field_RVE.seed_finite_elements(meshset_Elastic); CHKERRQ(ierr);
  
	Range prims_on_problem_bit_level_RVE;
	ierr = m_field_RVE.get_entities_by_type_and_ref_level(problem_bit_level_RVE,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level_RVE); CHKERRQ(ierr);
  //to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level_RVE;
  rval = moab_RVE.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
	rval = moab_RVE.add_entities(meshset_prims_on_problem_bit_level_RVE,prims_on_problem_bit_level_RVE); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level_RVE,BitRefLevel().set()); CHKERRQ(ierr);

  
  
  // ===========================================================================
  //
  // A.II. DEFINE PROBLEM
  //
  // ===========================================================================
  
  /*****************************************************************************
   *
   * Add field
   *  (1) Deterministic fields
   *  (2) Stochastic fields
   *       (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  int field_rank=3;
  ierr = m_field_RVE.add_field("DISP_RVE",H1,field_rank); CHKERRQ(ierr);
  ierr = m_field_RVE.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_finite_element("ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.add_finite_element("TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.add_finite_element("Lagrange_FE"); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   * set field data which finite element use
   * set field row which finite element use
   *
   ****************************************************************************/
  //Define rows/cols and element data
  ierr = m_field_RVE.modify_finite_element_add_field_row("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("ELASTIC_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  //FE Transverse Isotropic
  ierr = m_field_RVE.modify_finite_element_add_field_row("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","DISP_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("TRAN_ISO_FE_RVE","POTENTIAL_FIELD"); CHKERRQ(ierr);
  
  //C and CT
  //======================================================================================================
  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = m_field_RVE.modify_finite_element_add_field_col("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_row("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  
  //data
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_finite_element_add_field_data("Lagrange_FE","DISP_RVE"); CHKERRQ(ierr);
  //======================================================================================================
  
  
  //define problems
  ierr = m_field_RVE.add_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  
  //set finite elements for problem
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.modify_problem_add_finite_element("ELASTIC_PROBLEM_RVE","Lagrange_FE"); CHKERRQ(ierr);
  
  
  //set refinment level for problem
  ierr = m_field_RVE.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_RVE",problem_bit_level_RVE); CHKERRQ(ierr); // problem_bit_level
  
  // ===========================================================================
  //
  // A.III. DECLARE PROBLEM
  //
  // ===========================================================================
  
  
  /*****************************************************************************
   *
   * Add entitities (by tets) to the field
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_ents_to_field_by_TETs(0,"DISP_RVE"); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC_FE_RVE",true); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TETs(meshset_Trans_ISO,"TRAN_ISO_FE_RVE",true); CHKERRQ(ierr);

  Range SurfacesFaces;
  ierr = m_field_RVE.get_cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_FE"); CHKERRQ(ierr);
  
  
  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab_RVE.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
  rval = moab_RVE.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = m_field_RVE.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = m_field_RVE.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_RVE.set_field_order(0,MBTET,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBTRI,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"DISP_RVE",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"DISP_RVE",1); CHKERRQ(ierr);
  
  ierr = m_field_RVE.set_field_order(0,MBTRI,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order_RVE); CHKERRQ(ierr);
  ierr = m_field_RVE.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);
  
  
  
  // ===========================================================================
  //
  //  A.IV. BUILD DATABASE
  //
  // ===========================================================================
  
  //build field
  ierr = m_field_RVE.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field_RVE.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field_RVE.build_adjacencies(problem_bit_level_RVE); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_RVE.build_problems(); CHKERRQ(ierr);
  
  
  
  // ===========================================================================
  //
  //  A.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_RVE.partition_problem("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  ierr = m_field_RVE.partition_finite_elements("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field_RVE.partition_ghost_dofs("ELASTIC_PROBLEM_RVE"); CHKERRQ(ierr);
  
  //print bcs
  ierr = m_field_RVE.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_RVE.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_RVE.print_cubit_materials_set(); CHKERRQ(ierr);
  
  
  
  
  //============================================================================
  //
  //  B. Macro Problem
  //
  //============================================================================
  
  // ===========================================================================
  //
  // B. I. READ MESH DATA AND FINITE ELEMENT ANALYSIS CONTROL PARAMETERS FROM FILE
  //
  // ===========================================================================
  moab::Core mb_instance_Macro;
  Interface& moab_Macro = mb_instance_Macro;
  
  ParallelComm* pcomm_Macro = ParallelComm::get_pcomm(&moab_Macro,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm_Macro =  new ParallelComm(&moab_Macro,PETSC_COMM_WORLD);
  
  /*****************************************************************************
   *
   * Read parameters from line command
   *
   ****************************************************************************/
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file_macro",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file_macro (MESH FILE NEEDED)");
  }
  
  PetscInt order_Macro;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order_Macro",&order_Macro,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order_Macro = 1;
  }
  
  PetscInt NO_Layers;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_NO_Layers",&NO_Layers,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    NO_Layers = 1;
  }
  
  /*****************************************************************************
   *
   * Transfer mesh data to MOAB database
   *
   ****************************************************************************/
  /*
   * Meshing is divided into two groups, one for post-processing, and the rest 
   * in another group.
   */
  rval = moab_Macro.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core_Macro(moab_Macro);
  FieldInterface& m_field_Macro = core_Macro;
  
  //ref meshset ref level 0
  ierr = m_field_Macro.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0_Macro;
  bit_level0_Macro.set(0);
  
  EntityHandle meshset_level0;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field_Macro.seed_ref_level_3D(0,bit_level0_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.get_entities_by_ref_level(bit_level0_Macro,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  EntityHandle meshset_macro_Elastic,meshset_Reliability;
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_macro_Elastic); CHKERR_PETSC(rval);
  rval = moab_Macro.create_meshset(MESHSET_SET,meshset_Reliability); CHKERR_PETSC(rval);
  
  Range TetsInBlock_1st_Ply, TetsInBlock_2nd_Ply, TetsInBlock_3rd_Ply, TetsInBlock_4th_Ply;
  Range TetsInBlock_Reliability_1st,TetsInBlock_Reliability_2nd, TetsInBlock_Reliability_3rd,TetsInBlock_Reliability_4th;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_Macro,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_First") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_1st_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Second") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_2nd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Third") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_3rd_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Fourth") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_4th_Ply,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "RELIABILITY_First") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_Reliability_1st,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "RELIABILITY_Second") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_Reliability_2nd,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "RELIABILITY_Third") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_Reliability_3rd,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "RELIABILITY_Fourth") {
      rval = moab_Macro.get_entities_by_type(it->meshset, MBTET,TetsInBlock_Reliability_4th,true); CHKERR_PETSC(rval);
    }
  }
  
  // ===========================================================================
  //
  // B.II. DEFINE PROBLEM
  //
  // ===========================================================================
  

  /*****************************************************************************
   *
   * Add stochastic field
   * (total 14 field for 1st and 2nd order stochastic PSFEM)
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_field("DISP_MACRO",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field_Macro.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  
  /*****************************************************************************
   *
   * Create finite element for the defined fields
   *
   ****************************************************************************/
  // First layer
  ierr = m_field_Macro.add_finite_element("ELASTIC_1st_Ply"); CHKERRQ(ierr);
  // Second layer
  if (NO_Layers > 1) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    ierr = m_field_Macro.add_finite_element("ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO_REL_1st"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO_REL_2nd"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO_REL_3rd"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_finite_element("ELASTIC_FE_MACRO_REL_4th"); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * set field data which finite element use
   *
   ****************************************************************************/
  // First layer
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_1st_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  // Second layer
  if (NO_Layers > 1) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_2nd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_3rd_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    //Define rows/cols and element data
    ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","DISP_MACRO"); CHKERRQ(ierr);
    ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_4th_Ply","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  }
  
  //Define rows/cols and element data
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL_1st","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL_1st","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_1st","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_1st","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL_2nd","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL_2nd","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_2nd","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_2nd","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL_3rd","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL_3rd","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_3rd","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_3rd","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  ierr = m_field_Macro.modify_finite_element_add_field_row("ELASTIC_FE_MACRO_REL_4th","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_col("ELASTIC_FE_MACRO_REL_4th","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_4th","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO_REL_4th","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //define problems
  ierr = m_field_Macro.add_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //set finite elements for problem
  // First layer
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply"); CHKERRQ(ierr);
  // Second layer
  if (NO_Layers > 1) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  // Third layer
  if (NO_Layers > 2) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  // Fourth layer
  if (NO_Layers > 3) {
    ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_1st"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_2nd"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_3rd"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_4th"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field_Macro.modify_problem_ref_level_add_bit("ELASTIC_PROBLEM_MACRO",bit_level0_Macro); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  // B.III. DECLARE PROBLEM
  //
  // ===========================================================================

  /*****************************************************************************
   *
   * Add entitities (by tets) to the field
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  

  /*****************************************************************************
   *
   * Add finite elements entities
   *
   ****************************************************************************/
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_1st_Ply, "ELASTIC_1st_Ply"); CHKERRQ(ierr);
  if (NO_Layers > 1) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_2nd_Ply,"ELASTIC_2nd_Ply"); CHKERRQ(ierr);
  }
  if (NO_Layers > 2) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_3rd_Ply,"ELASTIC_3rd_Ply"); CHKERRQ(ierr);
  }
  if (NO_Layers > 3) {
    ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_4th_Ply,"ELASTIC_4th_Ply"); CHKERRQ(ierr);
  }
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_Reliability_1st,"ELASTIC_FE_MACRO_REL_1st"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_Reliability_2nd,"ELASTIC_FE_MACRO_REL_2nd"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_Reliability_3rd,"ELASTIC_FE_MACRO_REL_3rd"); CHKERRQ(ierr);
  ierr = m_field_Macro.add_ents_to_finite_element_by_TETs(TetsInBlock_Reliability_4th,"ELASTIC_FE_MACRO_REL_4th"); CHKERRQ(ierr);
  
  /*****************************************************************************
   *
   * Set applied order
   * See reference for detals:
   *   Ainsworth M. and Coyle J. (2003) Hierarchic finite element bases on
   *      unstructured tetrahedral meshes. IJNME, 58(14). pp.2103-2130.
   ****************************************************************************/
  ierr = m_field_Macro.set_field_order(0,MBTET,   "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "DISP_MACRO",order_Macro); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"DISP_MACRO",1); CHKERRQ(ierr);
  //
  ierr = m_field_Macro.set_field_order(0,MBTET,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBTRI,   "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBEDGE,  "MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field_Macro.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
 
 //***************************************************************************
   //Calculate Dmat for each Guass point here
  //Calculate_RVE_Dmat_TransIso_Disp calculate_rve_dmat_TransIso(m_field_Macro);
  
  //ierr = calculate_rve_dmat_TransIso.addElasticElements("DISP_MACRO"); CHKERRQ(ierr);
  //ierr = m_field_Macro.modify_finite_element_add_field_data("ELASTIC_FE_MACRO","Wt"); CHKERRQ(ierr);
  //ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO"); CHKERRQ(ierr);
 //***************************************************************************

  //ierr = MetaNeummanForces::addNeumannBCElements(m_field_Macro,"ELASTIC_PROB","DISP_MACRO"); CHKERRQ(ierr);
  ierr = MetaNeummanForces::addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field_Macro.modify_problem_add_finite_element("ELASTIC_PROBLEM_MACRO","FORCE_FE"); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  //  B.IV. BUILD DATABASE
  //
  // ===========================================================================
  
  //build field
  ierr = m_field_Macro.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = m_field_Macro.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = m_field_Macro.build_adjacencies(bit_level0_Macro); CHKERRQ(ierr);
  
  //build problem
  ierr = m_field_Macro.build_problems(); CHKERRQ(ierr);
  

  // ===========================================================================
  //
  //  B.V. MESH PARTITION
  //
  // ===========================================================================
  
  //partition
  ierr = m_field_Macro.partition_problem("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  //PetscBarrier(PETSC_NULL);
  ierr = m_field_Macro.partition_finite_elements("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field_Macro.partition_ghost_dofs("ELASTIC_PROBLEM_MACRO"); CHKERRQ(ierr);
  
  //m_field_Macro.list_dofs_by_field_name("DISP_MACRO",true);
  
  //print bcs
  ierr = m_field_Macro.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field_Macro.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field_Macro.print_cubit_materials_set(); CHKERRQ(ierr);
  
  // ===========================================================================
  //
  //  C. RELIABILITY ANLYSIS
  //
  // ===========================================================================
  
  cout<<"\n\n";
  cout<<"///////////////////////////////////////////////////////////////////\n";
  cout<<"//                                                               //\n";
  cout<<"//           Reliability calculation starts from here!           //\n";
  cout<<"//           Monte Carlo simultion method!                       //\n";
  cout<<"//                                                               //\n";
  cout<<"/////////////////////////////////////////////////////////////////\n\n";
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 1: read inputs from data file                                //
  //                a) probability data                                       //
  //                 .1) marginal distribution type                           //
  //                 .2) parameters                                           //
  //                 .3) correlation matrix                                   //
  //                b) limit state function                                   //
  //                 .1) indicator of how limit-state function is given       //
  //                     1:                                                   //
  //                c) analysis option                                        //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  
  /*
   *  Read inputs' statistical properties from file to insert into <probdata>
   */
  Stochastic_Model probdata;
  
  // Import data from textile file
  ImportProbData readprobdata;
  ierr = readprobdata.ProbdataFileIn(); CHKERRQ(ierr);
  probdata.marg        = readprobdata.MargProb;
  probdata.correlation = readprobdata.CorrMat;
  probdata.num_vars    = readprobdata.NumVars;
  probdata.MatStrength = readprobdata.MatStrength;
  probdata.NameVars    = readprobdata.NameVars;
  probdata.PlyAngle    = readprobdata.PlyAngle;
  probdata.ExaminedPly = readprobdata.ExaminedLayer;
  
  
  int FailureCriterion; // 1: Tsai-Wu, 2: Tsai-Hill
  string NameOfFailureCriterion;
  FailureCriterion = readprobdata.FailureCriterion;
  
  FE2_Macro_Solver_MCS Solve_FE2_Problem;
  
  LimitStateFunction_MCS TheLSF;
  
  ofstream MCSFile_MS;
  MCSFile_MS.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_MCS_MS.txt",ofstream::out);
  
  ofstream MCSFile_Hashin;
  MCSFile_Hashin.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_MCS_Hashin.txt",ofstream::out);
  
  ofstream MCSFile_TW;
  MCSFile_TW.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_MCS_TW.txt",ofstream::out);
  
  ofstream MCSFile_TH;
  MCSFile_TH.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_MCS_TH.txt",ofstream::out);
  
  ofstream MCSFile_RC;
  MCSFile_RC.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_MCS_RC.txt",ofstream::out);
  
  ofstream MCSFile_Hoffman;
  MCSFile_Hoffman.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Result_MCS_Hoffman.txt",ofstream::out);
  
  double val_G;
  double val_G_MS_LD, val_G_MS_TD, val_G_MS_Shear;
  double val_G_HF, val_G_HM;
  double val_G_TW_2D, val_G_TW;
  double val_G_TH_2D, val_G_TH;
  double val_G_RCF, val_G_RCM;
  double val_G_Hoffman_2D, val_G_Hoffman;
  
  int no_fail = 0;
  int no_mcs = 120000;
  ublas::vector<double> x;
  ublas::matrix<double> StressGP;
  ublas::matrix<double> StressGP_1st, StressGP_2nd, StressGP_3rd, StressGP_4th;
  
  double theta_angle;
  ublas::vector<double> PlyAngle_new;
  PlyAngle_new = probdata.PlyAngle;
  
  /*
   *  Start iteration
   */
  
  clock_t rel_t1, rel_t2;
  double rel_calc_time;
  rel_t1 = clock();
  
  ublas::matrix<double> Dmat;
  ublas::matrix<double> Dmat_1st, Dmat_2nd, Dmat_3rd, Dmat_4th;
  
  MCS_RNG sampling_values;
  
  /*for (int iims = 1; iims<=10000;iims++) {
    ierr = sampling_values.myRNG(probdata.num_vars,probdata.marg,x); CHKERRQ(ierr);
  }*/
  
  for (int imcs=1; imcs<=no_mcs; imcs++) {
    cout<<"\n\n*************************************************\n*\n";
    cout<<"*    This is "<<imcs<<" step!"<<"\n*\n";
    cout<<"*************************************************\n";
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 2: generate samples                                        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////

    //MCS_RNG sampling_values;
    ierr = sampling_values.myRNG(probdata.num_vars,probdata.marg,x); CHKERRQ(ierr);
    cout<<"\nThe random variables: "<<x<<endl;
    /*x(0) = 5800; x(1) = 0.38; x(2) = 0.0841; x(3) = 0.0144; x(4) = 14.96e3;
    x(5) = 208.34e3; x(6) = 14.96e3; x(7) = 1969; x(8) = 1480; x(9) = 38; x(10) = 200; x(11) = 79;
    cout<<"\nThe updated random variables: "<<x<<endl;*/
    
    // ----------------
    // Update ply angle
    // -----------------
    for (int i=1; i<=x.size();i++) {
      if (probdata.NameVars[i].compare(0,11,"orientation") == 0) {
        cout<<"\nAngle "<<x(i-1)<<endl;
        for (int j=0; j<probdata.PlyAngle.size(); j++) {
          cout << "The original ply angle "<<probdata.PlyAngle(j)<<"\t delta x: "<<x(i-1)<<endl;
          PlyAngle_new(j) = probdata.PlyAngle(j) + x(i-1);
          cout << "The modified ply angle "<<PlyAngle_new(j)<<endl;
        }
      }
      else if (probdata.NameVars[i].compare(0,6,"theta1") == 0) {
        PlyAngle_new(0) = x(i-1);
      }
      else if (probdata.NameVars[i].compare(0,6,"theta2") == 0) {
        PlyAngle_new(1) = x(i-1);
      }
      else if (probdata.NameVars[i].compare(0,6,"theta3") == 0) {
        PlyAngle_new(2) = x(i-1);
      }
      else if (probdata.NameVars[i].compare(0,6,"theta4") == 0) {
        PlyAngle_new(3) = x(i-1);
      }
    }
  
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 3: call FE program to calculate structural response        //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    //
    // Evaluate limit-state function and its gradient
    //
    ierr = Solve_FE2_Problem.Calculate_RVEDmat(m_field_RVE,x,probdata.num_vars,probdata.NameVars); CHKERRQ(ierr);
    ierr = Solve_FE2_Problem.Macro_FE_Laminate(m_field_Macro,x,probdata.num_vars,
                                               probdata.NameVars,
                                               PlyAngle_new,
                                               NO_Layers); CHKERRQ(ierr);
    switch (probdata.ExaminedPly) {
      case 0:{cout<<"\n\nThe 1st layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 1st ply
        //--------------------------
        Dmat_1st = Solve_FE2_Problem.Dmat_1st_Ply;
        Dmat_2nd = Solve_FE2_Problem.Dmat_2nd_Ply;
        Dmat_3rd = Solve_FE2_Problem.Dmat_3rd_Ply;
        Dmat_4th = Solve_FE2_Problem.Dmat_4th_Ply;
        break;
      }
      case 1:{cout<<"\n\nThe 1st layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 1st ply
        //--------------------------
        Dmat = Solve_FE2_Problem.Dmat_1st_Ply;
        break;
      }
      case 2: {cout<<"\n\nThe 2nd layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 2nd ply
        //--------------------------
        Dmat = Solve_FE2_Problem.Dmat_2nd_Ply;
        break;
      }
      case 3: {cout<<"\n\nThe 3rd layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 3rd ply
        //--------------------------
        Dmat = Solve_FE2_Problem.Dmat_3rd_Ply;
        break;
      }
      case 4: {cout<<"\n\nThe 4th layer is under examination.\n\n";
        //--------------------------
        // D matrix for the 4th ply
        //--------------------------
        Dmat = Solve_FE2_Problem.Dmat_4th_Ply;
        break;
      }
    }
    
    //theta_angle = PlyAngle_new(probdata.ExaminedPly - 1)*(M_PI/180.0);
  
    //
    // Calculate the stress at Gauss points for specific element(s)
    //
    FE2_PostProcStressForReliability_Zeroth Calc_Stress_1st(m_field_Macro,"DISP_MACRO",Dmat_1st);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_1st",Calc_Stress_1st);  CHKERRQ(ierr);
    
    FE2_PostProcStressForReliability_Zeroth Calc_Stress_2nd(m_field_Macro,"DISP_MACRO",Dmat_2nd);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_2nd",Calc_Stress_2nd);  CHKERRQ(ierr);
    
    FE2_PostProcStressForReliability_Zeroth Calc_Stress_3rd(m_field_Macro,"DISP_MACRO",Dmat_3rd);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_3rd",Calc_Stress_3rd);  CHKERRQ(ierr);
    
    FE2_PostProcStressForReliability_Zeroth Calc_Stress_4th(m_field_Macro,"DISP_MACRO",Dmat_4th);
    ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL_4th",Calc_Stress_4th);  CHKERRQ(ierr);
    
    //cout<<"Stress at GP: "<<Calc_Stress.StressGP<<endl;
    
    theta_angle = PlyAngle_new(0)*(M_PI/180.0);
    StressGP_1st.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_1st.StressGP,StressGP_1st);
    cout<<"1st ply stress at GP: "<<StressGP_1st<<endl;
    theta_angle = PlyAngle_new(1)*(M_PI/180.0);
    StressGP_2nd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_2nd.StressGP,StressGP_2nd);
    cout<<"2nd ply stress at GP: "<<StressGP_2nd<<endl;
    theta_angle = PlyAngle_new(2)*(M_PI/180.0);
    StressGP_3rd.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_3rd.StressGP,StressGP_3rd);
    cout<<"3rd ply stress at GP: "<<StressGP_3rd<<endl;
    theta_angle = PlyAngle_new(3)*(M_PI/180.0);
    StressGP_4th.clear(); REL_Stress_Transformation(theta_angle,Calc_Stress_4th.StressGP,StressGP_4th);
    cout<<"4th ply stress at GP: "<<StressGP_4th<<endl;
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 4: evaluate limit state function                           //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    switch (FailureCriterion) {
      case 0: {
        //
        // Tsai-Wu failure criteria
        //
        NameOfFailureCriterion = "All failure criteria";
        cout<<"\nAll failure criteria are considered!\n";
        /*
        // Maximum stress
        ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_LD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_TD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_MS_Shear); CHKERRQ(ierr);
        // Hashin
        ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_HF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_HM); CHKERRQ(ierr);
        // Tsai Wu
        ierr = TheLSF.gfun_ply_Tsai_Wu_0(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TW_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TW); CHKERRQ(ierr);
        // Tsai Hill
        ierr = TheLSF.gfun_ply_Tsai_Hill_2D(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TH_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TH); CHKERRQ(ierr);
        // Richard Christensen
        ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_RCF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_RCM); CHKERRQ(ierr);
        // Hoffman
        ierr = TheLSF.gfun_ply_Hoffman_2D(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_Hoffman_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_Hoffman); CHKERRQ(ierr);
         */
        //
        // First ply
        //
        // Maximum stress
        ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_MS_LD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_MS_TD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_MS_Shear); CHKERRQ(ierr);
        // Hashin
        ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_HF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_HM); CHKERRQ(ierr);
        // Tsai Wu
        ierr = TheLSF.gfun_ply_Tsai_Wu_0(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_TW_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_TW); CHKERRQ(ierr);
        // Tsai Hill
        ierr = TheLSF.gfun_ply_Tsai_Hill_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_TH_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_TH); CHKERRQ(ierr);
        // Richard Christensen
        ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_RCF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_RCM); CHKERRQ(ierr);
        // Hoffman
        ierr = TheLSF.gfun_ply_Hoffman_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_Hoffman_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,StressGP_1st,val_G_Hoffman); CHKERRQ(ierr);
        
        MCSFile_MS<<setprecision(15)<<val_G_MS_LD<<"\t"<<val_G_MS_TD<<"\t"<<val_G_MS_Shear<<"\t";
        MCSFile_Hashin<<setprecision(15)<<val_G_HF<<"\t"<<val_G_HM<<"\t";
        MCSFile_TW<<setprecision(15)<<val_G_TW_2D<<"\t"<<val_G_TW<<"\t";
        MCSFile_TH<<setprecision(15)<<val_G_TH_2D<<"\t"<<val_G_TH<<"\t";
        MCSFile_RC<<setprecision(15)<<val_G_RCF<<"\t"<<val_G_RCM<<"\t";
        MCSFile_Hoffman<<setprecision(15)<<val_G_Hoffman_2D<<"\t"<<val_G_Hoffman<<"\t";
        
        //
        // Second ply
        //
        // Maximum stress
        ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_MS_LD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_MS_TD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_MS_Shear); CHKERRQ(ierr);
        // Hashin
        ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_HF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_HM); CHKERRQ(ierr);
        // Tsai Wu
        ierr = TheLSF.gfun_ply_Tsai_Wu_0(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_TW_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_TW); CHKERRQ(ierr);
        // Tsai Hill
        ierr = TheLSF.gfun_ply_Tsai_Hill_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_TH_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_TH); CHKERRQ(ierr);
        // Richard Christensen
        ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_RCF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_RCM); CHKERRQ(ierr);
        // Hoffman
        ierr = TheLSF.gfun_ply_Hoffman_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_Hoffman_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,StressGP_2nd,val_G_Hoffman); CHKERRQ(ierr);
        
        MCSFile_MS<<setprecision(15)<<val_G_MS_LD<<"\t"<<val_G_MS_TD<<"\t"<<val_G_MS_Shear<<"\t";
        MCSFile_Hashin<<setprecision(15)<<val_G_HF<<"\t"<<val_G_HM<<"\t";
        MCSFile_TW<<setprecision(15)<<val_G_TW_2D<<"\t"<<val_G_TW<<"\t";
        MCSFile_TH<<setprecision(15)<<val_G_TH_2D<<"\t"<<val_G_TH<<"\t";
        MCSFile_RC<<setprecision(15)<<val_G_RCF<<"\t"<<val_G_RCM<<"\t";
        MCSFile_Hoffman<<setprecision(15)<<val_G_Hoffman_2D<<"\t"<<val_G_Hoffman<<"\t";
        
        
        //
        // Third ply
        //
        // Maximum stress
        ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_MS_LD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_MS_TD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_MS_Shear); CHKERRQ(ierr);
        // Hashin
        ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_HF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_HM); CHKERRQ(ierr);
        // Tsai Wu
        ierr = TheLSF.gfun_ply_Tsai_Wu_0(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_TW_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_TW); CHKERRQ(ierr);
        // Tsai Hill
        ierr = TheLSF.gfun_ply_Tsai_Hill_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_TH_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_TH); CHKERRQ(ierr);
        // Richard Christensen
        ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_RCF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_RCM); CHKERRQ(ierr);
        // Hoffman
        ierr = TheLSF.gfun_ply_Hoffman_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_Hoffman_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,StressGP_3rd,val_G_Hoffman); CHKERRQ(ierr);
        
        MCSFile_MS<<setprecision(15)<<val_G_MS_LD<<"\t"<<val_G_MS_TD<<"\t"<<val_G_MS_Shear<<"\t";
        MCSFile_Hashin<<setprecision(15)<<val_G_HF<<"\t"<<val_G_HM<<"\t";
        MCSFile_TW<<setprecision(15)<<val_G_TW_2D<<"\t"<<val_G_TW<<"\t";
        MCSFile_TH<<setprecision(15)<<val_G_TH_2D<<"\t"<<val_G_TH<<"\t";
        MCSFile_RC<<setprecision(15)<<val_G_RCF<<"\t"<<val_G_RCM<<"\t";
        MCSFile_Hoffman<<setprecision(15)<<val_G_Hoffman_2D<<"\t"<<val_G_Hoffman<<"\t";
        
        
        //
        // Fourth ply
        //
        // Maximum stress
        ierr = TheLSF.gfun_ply_MS_LD(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_MS_LD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_TD(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_MS_TD); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_MS_Shear(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_MS_Shear); CHKERRQ(ierr);
        // Hashin
        ierr = TheLSF.gfun_ply_HF(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_HF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_HM(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_HM); CHKERRQ(ierr);
        // Tsai Wu
        ierr = TheLSF.gfun_ply_Tsai_Wu_0(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_TW_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_TW); CHKERRQ(ierr);
        // Tsai Hill
        ierr = TheLSF.gfun_ply_Tsai_Hill_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_TH_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_TH); CHKERRQ(ierr);
        // Richard Christensen
        ierr = TheLSF.gfun_ply_RCF(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_RCF); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_RCM(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_RCM); CHKERRQ(ierr);
        // Hoffman
        ierr = TheLSF.gfun_ply_Hoffman_2D(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_Hoffman_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Hoffman(x,probdata.NameVars,probdata.MatStrength,StressGP_4th,val_G_Hoffman); CHKERRQ(ierr);
        
        MCSFile_MS<<setprecision(15)<<val_G_MS_LD<<"\t"<<val_G_MS_TD<<"\t"<<val_G_MS_Shear<<"\n";
        MCSFile_Hashin<<setprecision(15)<<val_G_HF<<"\t"<<val_G_HM<<"\n";
        MCSFile_TW<<setprecision(15)<<val_G_TW_2D<<"\t"<<val_G_TW<<"\n";
        MCSFile_TH<<setprecision(15)<<val_G_TH_2D<<"\t"<<val_G_TH<<"\n";
        MCSFile_RC<<setprecision(15)<<val_G_RCF<<"\t"<<val_G_RCM<<"\n";
        MCSFile_Hoffman<<setprecision(15)<<val_G_Hoffman_2D<<"\t"<<val_G_Hoffman<<"\n";
        
        
        if (val_G<0) { no_fail++;}
        break;
      }
      case 1:
        //
        // Tsai-Wu failure criteria
        //
        NameOfFailureCriterion = "Tsai-Wu";
        // Tsai Wu
        ierr = TheLSF.gfun_ply_Tsai_Wu_0(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu_2D(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TW_2D); CHKERRQ(ierr);
        ierr = TheLSF.gfun_ply_Tsai_Wu(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G_TW); CHKERRQ(ierr);

        if (val_G<0) { no_fail++;}
        break;
      case 2:
        //
        // Tsai-Hill failure criteria
        //
        NameOfFailureCriterion = "Tsai-Hill";
        ierr = TheLSF.gfun_ply_Tsai_Hill(x,probdata.NameVars,probdata.MatStrength,StressGP,val_G); CHKERRQ(ierr);
        break;
    }
    
    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    //        STEP 4: output limit state function                             //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////
    
    /*MCSFile_MS<<setprecision(15)<<val_G_MS_LD<<"\t"<<val_G_MS_LD<<"\t"<<val_G_MS_LD<<"\n";
    MCSFile_Hashin<<setprecision(15)<<val_G_HF<<"\t"<<val_G_HF<<"\n";
    MCSFile_TW<<setprecision(15)<<val_G_TW_2D<<"\t"<<val_G_TW<<"\n";
    MCSFile_TH<<setprecision(15)<<val_G_TH_2D<<"\t"<<val_G_TH<<"\n";
    MCSFile_RC<<setprecision(15)<<val_G_RCF<<"\t"<<val_G_RCM<<"\n";
    MCSFile_Hoffman<<setprecision(15)<<val_G_Hoffman_2D<<"\t"<<val_G_Hoffman<<"\n";*/
    /*
    MCSFile<<setprecision(15)<<val_G<<"\t"<<val_G_TH;
    for (int i=0;i<probdata.num_vars;i++) {
      MCSFile<<"\t"<<x(i);
    }
    MCSFile<<"\n";*/
    
    rel_t2 = clock();
    rel_calc_time  = (double)(rel_t2 - rel_t1)/CLOCKS_PER_SEC;
    rel_t1 = rel_t2;
    cout<<"Elapsed time at this step is: "<<rel_calc_time<<" seconds.\n";
  }
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //        STEP 5: Calculate the probability of failure and corresponding    //
  //                estimate of reliability index                             //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  double prob_failure;
  prob_failure = no_fail/no_mcs; cout<<"\n The number of failure is "<<no_fail<<endl;
  double beta;
  using boost::math::normal_distribution;
  normal_distribution<> snorm(0,1);
  if (no_fail == 0) {
    beta = 0.0;
  } else {
    beta = quantile(snorm,prob_failure);
  }
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  //                               FINISH                                     //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  
  finish_time = clock();
  total_time  = (double)(finish_time - start_time)/CLOCKS_PER_SEC;
  
   
  cout<<"\n\n************************************************************\n*\n";
  cout<<"*  The estimate of reliability index (beta) is: "<<beta<<endl;
  cout<<"*  The probability of failure is:               "<<prob_failure<<endl;
  cout<<"*  The failure criterion used is:               "<<NameOfFailureCriterion<<endl;
  cout<<"*  Elapsed time is "<<total_time<<" seconds.\n";
  cout<<"*  The program finishes !!! \n*\n";
  cout<<"****************************************************************"<<endl;
  
  
  ierr = PetscFinalize(); CHKERRQ(ierr);
  
}
