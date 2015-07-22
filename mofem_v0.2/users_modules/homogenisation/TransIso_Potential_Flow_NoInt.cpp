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

//#include "FunctionsForFieldData.hpp"
//#include "cholesky.hpp"

extern "C" {
#include <gm_rule.h>
}

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>
#include <petsctime.h>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

using namespace ObosleteUsersModules;

#include <MethodForForceScaling.hpp>
#include "PotentialFlowFEMethod.hpp"
#include "SurfacePressure.hpp"


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";


int main(int argc, char *argv[]) {

  try {

    PetscInitialize(&argc,&argv,(char *)0,help);

    moab::Core mb_instance;
    Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
    }
    PetscInt order;
    ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      order = 1;
    }

    PetscInt mesh_refinement_level;
    ierr = PetscOptionsGetInt(PETSC_NULL,"-my_mesh_ref_level",&mesh_refinement_level,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      mesh_refinement_level = 0;
    }


    const char *option;
    option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);


    MoFEM::Core core(moab);
    FieldInterface& mField = core;


    //=======================================================================================================
    //Seting nodal coordinates on the surface to make sure they are periodic
    //=======================================================================================================

    Range SurTrisNeg, SurTrisPos;
    ierr = mField.get_cubit_msId_entities_by_dimension(101,SIDESET,2,SurTrisNeg,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 101 = %d\n",SurTrisNeg.size()); CHKERRQ(ierr);
    ierr = mField.get_cubit_msId_entities_by_dimension(102,SIDESET,2,SurTrisPos,true); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SideSet 102 = %d\n",SurTrisPos.size()); CHKERRQ(ierr);

    Range SurNodesNeg,SurNodesPos;
    rval = moab.get_connectivity(SurTrisNeg,SurNodesNeg,true); CHKERR_PETSC(rval);
    cout<<" All nodes on negative surfaces " << SurNodesNeg.size()<<endl;
    rval = moab.get_connectivity(SurTrisPos,SurNodesPos,true); CHKERR_PETSC(rval);
    cout<<" All nodes on positive surfaces " << SurNodesPos.size()<<endl;


    double roundfact=1000.0;   double coords_nodes[3];
    //Populating the Multi-index container with nodes on -ve faces
    for(Range::iterator nit = SurNodesNeg.begin(); nit!=SurNodesNeg.end();  nit++) {
      rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
      //round values to 3 disimal places
      if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
      if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
      if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
      rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
      //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
    }

    ///Populating the Multi-index container with nodes on +ve faces
    for(Range::iterator nit = SurNodesPos.begin(); nit!=SurNodesPos.end();  nit++) {
      rval = moab.get_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
      //round values to 3 disimal places
      if(coords_nodes[0]>=0) coords_nodes[0]=double(int(coords_nodes[0]*roundfact+0.5))/roundfact;  else coords_nodes[0]=double(int(coords_nodes[0]*roundfact-0.5))/roundfact;
      if(coords_nodes[1]>=0) coords_nodes[1]=double(int(coords_nodes[1]*roundfact+0.5))/roundfact;  else coords_nodes[1]=double(int(coords_nodes[1]*roundfact-0.5))/roundfact;
      if(coords_nodes[2]>=0) coords_nodes[2]=double(int(coords_nodes[2]*roundfact+0.5))/roundfact;  else coords_nodes[2]=double(int(coords_nodes[2]*roundfact-0.5))/roundfact;
      rval = moab.set_coords(&*nit,1,coords_nodes);  CHKERR_PETSC(rval);
      //      cout<<"   coords_nodes[0]= "<<coords_nodes[0] << "   coords_nodes[1]= "<< coords_nodes[1] << "   coords_nodes[2]= "<< coords_nodes[2] <<endl;
    }
    //=======================================================================================================


    //add fields
    ierr = mField.add_field("POTENTIAL_FIELD",H1,1); CHKERRQ(ierr);
    ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

    ///Getting No. of Fibres and their index to be used for Potential Flow Problem
    int noOfFibres=0;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|UNKNOWNCUBITNAME,it)) {

      std::size_t found=it->get_name().find("PotentialFlow");
      if (found==std::string::npos) continue;
      noOfFibres += 1;
    }
		cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;

    vector<int> fibreList(noOfFibres,0);
    int aaa=0;

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|UNKNOWNCUBITNAME,it)) {

      std::size_t interfaceFound=it->get_name().find("PotentialFlow_");
      if (interfaceFound==std::string::npos) continue;

      std::string str2 = it->get_name().substr (14,50);

      fibreList[aaa] = atoi(str2.c_str());
      aaa += 1;
    }

    //************************************************************//

    for (int cc = 0; cc < noOfFibres; cc++) {

      ostringstream sss, rrr;

      //add finite elements
      sss << "POTENTIAL_ELEM" << fibreList[cc];
      cout<<sss.str().c_str()<<endl;
      ierr = mField.add_finite_element( sss.str().c_str() ); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row( sss.str().c_str() ,"POTENTIAL_FIELD"); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col( sss.str().c_str() ,"POTENTIAL_FIELD"); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data( sss.str().c_str() ,"POTENTIAL_FIELD"); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data( sss.str().c_str() ,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

      //add problems
      rrr << "POTENTIAL_PROBLEM" << fibreList[cc];
      ierr = mField.add_problem( rrr.str().c_str() ); CHKERRQ(ierr);
      //define problems and finite elements
      ierr = mField.modify_problem_add_finite_element( rrr.str().c_str() , sss.str().c_str() ); CHKERRQ(ierr);

    }

		Tag th_meshset_info;
		int def_meshset_info[2] = {0,0};
		rval = moab.tag_get_handle("MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info);

		int meshset_data[2];
  	EntityHandle root = moab.get_root_set();
		rval = moab.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);

    ierr = mField.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
    vector<BitRefLevel> bit_levels;
    bit_levels.push_back(BitRefLevel().set(0));

		//MESH REFINEMENT
		int ll = 1;


		//End of refinement, save level of refinement
		int meshset_data_root[2]={ll,0};
		rval = moab.tag_set_data(th_meshset_info,&root,1,meshset_data_root); CHKERR_PETSC(rval);

		/******************TETS TO MESHSET AND SAVING TETS ENTITIES******************/
		EntityHandle out_tet_meshset;
		rval = moab.create_meshset(MESHSET_SET,out_tet_meshset); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTET,out_tet_meshset); CHKERRQ(ierr);
		rval = moab.write_file("out_tets.vtk","VTK","",&out_tet_meshset,1); CHKERR_PETSC(rval);
		/*******************************************************/

		/******************PRISMS TO MESHSET AND SAVING PRISMS ENTITIES******************/
		EntityHandle out_meshset_prism;
		rval = moab.create_meshset(MESHSET_SET,out_meshset_prism); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBPRISM,out_meshset_prism); CHKERRQ(ierr);
		rval = moab.write_file("out_prism.vtk","VTK","",&out_meshset_prism,1); CHKERR_PETSC(rval);
		/*******************************************************/

		EntityHandle out_meshset;
		rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
		rval = moab.write_file("out_all_mesh.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
		Range LatestRefinedTets;
		rval = moab.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);


    BitRefLevel problem_bit_level = bit_levels.back();


    ///Adding entities to Field and FE for Potential Flow Problem
    for (int cc = 0; cc < noOfFibres; cc++) {

      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|UNKNOWNCUBITNAME,it)) {

        //      std::size_t found=it->get_name().find("PotentialFlow");
        //      if (found==std::string::npos) continue;
        //          cout<<it->get_name()<<endl;

        ostringstream sss,rrr;
        //set problem level
        sss << "POTENTIAL_ELEM" << fibreList[cc];
        rrr << "PotentialFlow_" << fibreList[cc];

        if(it->get_name() ==  rrr.str().c_str() ) {
          Range TetsInBlock;
          rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
          Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);

          ierr = mField.add_ents_to_field_by_TETs(0,"POTENTIAL_FIELD"); CHKERRQ(ierr);
          ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

          //add finite elements entities
          ierr = mField.add_ents_to_finite_element_by_TETs(block_rope_bit_level, sss.str().c_str()); CHKERRQ(ierr);

        }
      }
    }

    ierr = mField.set_field_order(0,MBVERTEX,"POTENTIAL_FIELD",1); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTET,"POTENTIAL_FIELD",order); CHKERRQ(ierr);

    ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);


    struct myMetaNeummanForces{

      static PetscErrorCode addNeumannFluxBCElements(
                                                     FieldInterface &mField,
                                                     const string problem_name,
                                                     const string field_name,
                                                     const int fibre_id,
                                                     const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {

        PetscFunctionBegin;
        PetscErrorCode ierr;
        ErrorCode rval;

        ostringstream sss,rrr,ppp;
        ppp << "PressureIO_" << fibre_id <<"_1";
        sss << "PressureIO_" << fibre_id <<"_2";
        rrr << "FLUX_FE" <<fibre_id;
        ierr = mField.add_finite_element( rrr.str().c_str() ,MF_ZERO); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_row(rrr.str().c_str(),field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_col(rrr.str().c_str(),field_name); CHKERRQ(ierr);
        ierr = mField.modify_finite_element_add_field_data(rrr.str().c_str(),field_name); CHKERRQ(ierr);
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {

          //          std::size_t PressureFound=it->get_name().find(sss.str().c_str());
          //          if (PressureFound==std::string::npos) continue;
          if (ppp.str().c_str()==it->get_name() || sss.str().c_str()==it->get_name()){

            if(mField.check_field(mesh_nodals_positions)) {
              ierr = mField.modify_finite_element_add_field_data( rrr.str().c_str() ,mesh_nodals_positions); CHKERRQ(ierr);
            }
            ierr = mField.modify_problem_add_finite_element(problem_name, rrr.str().c_str() ); CHKERRQ(ierr);
            Range tris;
            rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
            ierr = mField.add_ents_to_finite_element_by_TRIs(tris, rrr.str().c_str() ); CHKERRQ(ierr);
          }
        }

        PetscFunctionReturn(0);
      }

      static PetscErrorCode setNeumannFluxFiniteElementOperators(
                                                                 FieldInterface &mField,
                                                                 boost::ptr_map<string,NeummanForcesSurface> &neumann_forces,
                                                                 Vec &F,const string field_name,const int fibre_id, const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {

        PetscFunctionBegin;
        PetscErrorCode ierr;
        string fe_name;
        //        fe_name = "FLUX_FE";
        ostringstream sss,rrr,ppp;
        ppp << "PressureIO_" << fibre_id<<"_1";
        sss << "PressureIO_" << fibre_id<<"_2";
        rrr << "FLUX_FE" <<fibre_id;
        fe_name = rrr.str().c_str();
        neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {

          //          std::size_t PressureFound=it->get_name().find(sss.str().c_str());
          //          if (PressureFound==std::string::npos) continue;
          if (ppp.str().c_str()==it->get_name() || sss.str().c_str()==it->get_name()){
            bool ho_geometry = mField.check_field(mesh_nodals_positions);
            ierr = neumann_forces.at(fe_name).addFlux(field_name,F,it->get_msId(),ho_geometry); CHKERRQ(ierr);
            /*pressure_cubit_bc_data data;
             ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
             my_split << *it << endl;
             my_split << data << endl;*/
          }
        }
        PetscFunctionReturn(0);
      }

    };



    for (int cc = 0; cc < noOfFibres; cc++) {
      ostringstream sss;
      sss << "POTENTIAL_PROBLEM" << fibreList[cc];

      //flux boundary conditions
      ierr = myMetaNeummanForces::addNeumannFluxBCElements(mField,sss.str().c_str(),"POTENTIAL_FIELD",fibreList[cc]); CHKERRQ(ierr);
      //set problem level
      ierr = mField.modify_problem_ref_level_add_bit( sss.str().c_str() ,problem_bit_level); CHKERRQ(ierr);
    }

    //build fields
    ierr = mField.build_fields(); CHKERRQ(ierr);
    //build finite elements
    ierr = mField.build_finite_elements(); CHKERRQ(ierr);
    //build adjacencies
    ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);
    //build problem
    ierr = mField.build_problems(); CHKERRQ(ierr);


    //partition problems
    for (int cc = 0; cc < noOfFibres; cc++) {
      ostringstream sss;
      sss << "POTENTIAL_PROBLEM" << fibreList[cc];
      ierr = mField.partition_problem( sss.str().c_str() ); CHKERRQ(ierr);
      ierr = mField.partition_finite_elements( sss.str().c_str() ); CHKERRQ(ierr);
      ierr = mField.partition_ghost_dofs( sss.str().c_str() ); CHKERRQ(ierr);
    }

    //print bcs
    ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
    ierr = mField.print_cubit_pressure_set(); CHKERRQ(ierr);

    //create matrices and vectors
		vector<Vec> F(noOfFibres);
		vector<Vec> D(noOfFibres);
		vector<Mat> A(noOfFibres);

    for (int cc = 0; cc < noOfFibres; cc++) {
      ostringstream sss,rrr,ttt;
      sss << "POTENTIAL_PROBLEM" << fibreList[cc];
      rrr << "POTENTIAL_ELEM" << fibreList[cc];
      ttt << "ZeroPressure_" << fibreList[cc];
      ierr = mField.VecCreateGhost( sss.str().c_str() ,ROW,&F[cc]); CHKERRQ(ierr);
      ierr = mField.VecCreateGhost( sss.str().c_str() ,COL,&D[cc]); CHKERRQ(ierr);
      ierr = mField.MatCreateMPIAIJWithArrays( sss.str().c_str() ,&A[cc]); CHKERRQ(ierr);

      Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
      ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

      //get nodes and other entities to fix
      Range fix_nodes;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|UNKNOWNCUBITNAME,it)) {
        //        std::size_t zeroPressureFound=it->get_name().find(ttt.str().c_str());
        //        if (zeroPressureFound==std::string::npos) continue;
        if (ttt.str().c_str()==it->get_name()){
          rval = moab.get_entities_by_type(it->meshset,MBVERTEX,fix_nodes,true); CHKERR_PETSC(rval);
          Range edges;
          rval = moab.get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
          Range tris;
          rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
          Range adj;
          rval = moab.get_connectivity(tris,adj,true); CHKERR_PETSC(rval);
          fix_nodes.insert(adj.begin(),adj.end());
          rval = moab.get_connectivity(edges,adj,true); CHKERR_PETSC(rval);
          fix_nodes.insert(adj.begin(),adj.end());
          rval = moab.get_adjacencies(tris,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
        }
      }
      FixBcAtEntities fix_dofs(mField,"POTENTIAL_FIELD",A[cc],D[cc],F[cc],fix_nodes);
      //initialize data structure
      ierr = mField.problem_basic_method_preProcess( sss.str().c_str() ,fix_dofs); CHKERRQ(ierr);

      //neuman flux bc elements
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      ierr = myMetaNeummanForces::setNeumannFluxFiniteElementOperators(mField,neumann_forces,F[cc],"POTENTIAL_FIELD",fibreList[cc]); CHKERRQ(ierr);
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = mField.loop_finite_elements( sss.str().c_str() ,mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
      }

      LaplacianElem elem(mField,A[cc],F[cc]);

      ierr = MatZeroEntries(A[cc]); CHKERRQ(ierr);
      ierr = mField.loop_finite_elements( sss.str().c_str() , rrr.str().c_str() ,elem);  CHKERRQ(ierr);

      //post proces fix boundary conditiond
      ierr = mField.problem_basic_method_postProcess( sss.str().c_str() ,fix_dofs); CHKERRQ(ierr);

      ierr = MatAssemblyBegin(A[cc],MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A[cc],MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

      ierr = VecAssemblyBegin(F[cc]); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F[cc]); CHKERRQ(ierr);

      //		VecView(F[0],PETSC_VIEWER_STDOUT_WORLD);
      //set matrix possitives define and symetric for cholesky and icc preceonditionser
      ierr = MatSetOption(A[cc],MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);

      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,A[cc],A[cc]); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);

      ierr = KSPSolve(solver,F[cc],D[cc]); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D[cc],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D[cc],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = mField.set_global_ghost_vector( sss.str().c_str() ,ROW,D[cc],INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      ierr = VecDestroy(&F[cc]); CHKERRQ(ierr);
      ierr = VecDestroy(&D[cc]); CHKERRQ(ierr);
      ierr = MatDestroy(&A[cc]); CHKERRQ(ierr);
    }

    Tag th_phi;
    double def_val = 0;
    rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_FIELD",dof)) {
      EntityHandle ent = dof->get_ent();
      double val = dof->get_FieldData();
      rval = moab.tag_set_data(th_phi,&ent,1,&val); CHKERR_PETSC(rval);
    }

    ProjectionFieldOn10NodeTet ent_method_phi_on_10nodeTet(mField,"POTENTIAL_FIELD",true,false,"PHI");
    ierr = mField.loop_dofs("POTENTIAL_FIELD",ent_method_phi_on_10nodeTet); CHKERRQ(ierr);
    ent_method_phi_on_10nodeTet.set_nodes = false;
    ierr = mField.loop_dofs("POTENTIAL_FIELD",ent_method_phi_on_10nodeTet); CHKERRQ(ierr);

    if(pcomm->rank()==0) {
      rval = moab.write_file("solution_RVE.h5m"); CHKERR_PETSC(rval);
    }

		EntityHandle out_meshset1;
		rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
		ierr = mField.get_entities_by_type_and_ref_level(bit_levels[0],BitRefLevel().set(),MBTET,out_meshset1); CHKERRQ(ierr);
		rval = moab.write_file("solution2.vtk","VTK","",&out_meshset1,1); CHKERR_PETSC(rval);

    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);

      for (int cc = 0; cc < noOfFibres; cc++) {
        EntityHandle out_meshset1;
        rval = moab.create_meshset(MESHSET_SET,out_meshset1); CHKERR_PETSC(rval);
        ostringstream sss,rrr,ttt;
        sss << "POTENTIAL_PROBLEM" << fibreList[cc];
        rrr << "POTENTIAL_ELEM" << fibreList[cc];
        ttt << "out_potential_flow" << fibreList[cc] <<".vtk";
        ierr = mField.get_problem_finite_elements_entities( sss.str().c_str() , rrr.str().c_str() ,out_meshset); CHKERRQ(ierr);
        ierr = mField.get_problem_finite_elements_entities( sss.str().c_str() , rrr.str().c_str() ,out_meshset1); CHKERRQ(ierr);

        rval = moab.write_file( ttt.str().c_str() ,"VTK","",&out_meshset1,1); CHKERR_PETSC(rval);
        rval = moab.delete_entities(&out_meshset1,1); CHKERR_PETSC(rval);
      }

      rval = moab.write_file("out_potential_flow.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }

    PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

}
