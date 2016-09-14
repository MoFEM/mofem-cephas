#include <MoFEM.hpp>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    //Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    MeshsetsManager *meshsets_manager_ptr;
    ierr = m_field.query_interface(meshsets_manager_ptr); CHKERRQ(ierr);

    std::cout << "<<<< SIDESETs >>>>>" << std::endl;

    bool add_block_is_there = false;
    ierr = meshsets_manager_ptr->addMeshset(SIDESET,1002); CHKERRQ(ierr);
    {
      PressureCubitBcData mybc;
      strncpy(mybc.data.name,"Pressure",8);
      mybc.data.flag1 = 0;
      mybc.data.flag2 = 0;
      mybc.data.value1 = 1;
      ierr = meshsets_manager_ptr->setBcData(SIDESET,1002,mybc); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      if(it->getMeshsetId()!=1002) continue;
      add_block_is_there = true;
      PressureCubitBcData mydata;
      ierr = it->getBcDataStructure(mydata); CHKERRQ(ierr);
      //Print data
      std::cout << mydata;
    }
    if(!add_block_is_there) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_OPERATION_UNSUCCESSFUL,"no added block set");
    }

    std::cout << "<<<< BLOCKSETs >>>>>" << std::endl;

    add_block_is_there = false;
    ierr = meshsets_manager_ptr->addMeshset(BLOCKSET,1000,"ADD_BLOCK_SET"); CHKERRQ(ierr);
    std::vector<double> attr(3);
    attr[0] = 0;
    attr[1] = 1;
    attr[2] = 2;
    ierr = meshsets_manager_ptr->setAttribites(BLOCKSET,1000,attr); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      //Get block name
      std::string name = it->getName();
      if (name.compare(0,13,"ADD_BLOCK_SET") == 0) {
        add_block_is_there = true;
        std::vector<double> attributes;
        it->getAttributes(attributes);
        if(attributes.size()!=3) {
          SETERRQ1(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"should be 3 attributes but is %d",attributes.size());
        }
        if(attributes[0]!=0 || attributes[1]!=1 || attributes[2]!=2) {
          SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"wrong values of attributes");
        }
      }
    }
    if(!add_block_is_there) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_OPERATION_UNSUCCESSFUL,"no added block set");
    }
    add_block_is_there = false;
    ierr = meshsets_manager_ptr->addMeshset(BLOCKSET,1001,"MAT_ELASTIC"); CHKERRQ(ierr);
    {
      Mat_Elastic mydata;
      mydata.data.Young = 1;
      mydata.data.Poisson = 0.25;
      ierr = meshsets_manager_ptr->setAttribitesByDataStructure(BLOCKSET,1001,mydata); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      if(it->getMeshsetId()!=1001) continue;
      //Get block name
      std::string name = it->getName();
      if(name.compare(0,13,"MAT_ELASTIC") == 0 && (it->getBcType()&CubitBCType(MAT_ELASTICSET)).any()) {
        add_block_is_there = true;
        Mat_Elastic mydata;
        ierr = it->getAttributeDataStructure(mydata); CHKERRQ(ierr);
        //Print data
        std::cout << mydata;
        if(mydata.data.Young != 1 || mydata.data.Poisson != 0.25) {
          SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"wrong values of attributes");
        }
      }
    }
    if(!add_block_is_there) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_OPERATION_UNSUCCESSFUL,"no added block set");
    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

}
