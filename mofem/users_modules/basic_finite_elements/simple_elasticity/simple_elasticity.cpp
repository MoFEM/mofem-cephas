/** \file elasticity.cpp
 * \ingroup nonlinear_elastic_elem
 * \example elasticity.cpp

 The example shows how to solve the linear elastic problem. An example can read
 file with temperature field, then thermal stresses are included.

 What example can do:
 - take into account temperature field, i.e. calculate thermal stresses and
deformation
 - stationary and time depend field is considered
 - take into account gravitational body forces
 - take in account fluid pressure
 - can work with higher order geometry definition
 - works on distributed meshes
 - multi-grid solver where each grid level is approximation order level
 - each mesh block can have different material parameters and approximation
order

See example how code can be used \cite jordi:2017,
 \image html SquelaDamExampleByJordi.png "Example what you can do with this
code. Analysis of the arch dam of Susqueda, located in Catalonia (Spain)"
width=800px

 This is an example of application code; it does not show how elements are
implemented. Example presents how to:
 - read mesh
 - set-up problem
 - run finite elements on the problem
 - assemble matrices and vectors
 - solve the linear system of equations
 - save results


 If you like to see how to implement finite elements, material, are other parts
of the code, look here;
 - Hooke material, see \ref Hooke
 - Thermal-stress assembly, see \ref  ThermalElement
 - Body forces element, see \ref BodyFroceConstantField
 - Fluid pressure element, see \ref FluidPressure
 - The general implementation of an element for arbitrary Lagrangian-Eulerian
 aformulation for a nonlinear elastic problem is here \ref
 NonlinearElasticElement. Here we limit ourselves to Hooke equation and fix
 mesh, so the problem becomes linear. Not that elastic element is implemented
 with automatic differentiation.

*/

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#include <BasicFiniteElements.hpp>
#include <boost/program_options.hpp>

#include "../poisson/src/AuxPoissonFunctions.hpp"
#include "../poisson/src/PoissonOperators.hpp"
#include <ElasticityNewOperator.hpp>
#include <Hooke.hpp>

using namespace boost::numeric;
using namespace MoFEM;
using namespace std;
namespace po = boost::program_options;

static char help[] = "-my_block_config set block data\n"
                     "\n";

// const double young_modulus = 1;
// const double poisson_ratio = 0.0;

struct BlockOptionData {
  int oRder;
  double yOung;
  double pOisson;
  double initTemp;
  BlockOptionData() : oRder(-1), yOung(-1), pOisson(-2), initTemp(0) {}
};

struct OpK : public VolumeElementForcesAndSourcesCore::UserDataOperator {

  MatrixDouble rowB;
  MatrixDouble colB;
  MatrixDouble CB;
  MatrixDouble K, transK;

  MatrixDouble D;
  double yOung;
  double pOisson;
  double coefficient;

  OpK(bool symm = true)
      : VolumeElementForcesAndSourcesCore::UserDataOperator("U", "U", OPROWCOL,
                                                            symm) {

    pOisson     = 0.1;
    yOung       = 10;
    coefficient = yOung / ((1 + pOisson) * (1 - 2 * pOisson));
    D.resize(6, 6, false);
    D.clear();

    D(0, 0) = 1 - pOisson;
    D(1, 1) = 1 - pOisson;
    D(2, 2) = 1 - pOisson;
    D(3, 3) = 0.5 * (1 - 2 * pOisson);
    D(4, 4) = 0.5 * (1 - 2 * pOisson);
    D(5, 5) = 0.5 * (1 - 2 * pOisson);

    D(0, 1) = pOisson;
    D(0, 2) = pOisson;
    D(1, 0) = pOisson;

    D(1, 2) = pOisson;
    D(2, 0) = pOisson;
    D(2, 1) = pOisson;

    D *= coefficient;
  }

  MoFEMErrorCode makeB(const MatrixAdaptor &diffN, MatrixDouble &B) {

    MoFEMFunctionBeginHot;
    unsigned int nb_dofs = diffN.size1();
    B.resize(6, 3 * nb_dofs, false);
    B.clear();
    for (unsigned int dd = 0; dd < nb_dofs; dd++) {
      const double diff[] = {diffN(dd, 0), diffN(dd, 1), diffN(dd, 2)};
      const int dd3       = 3 * dd;
      for (int rr = 0; rr < 3; rr++) {
        B(rr, dd3 + rr) = diff[rr];
      }
      // gamma_xy
      B(3, dd3 + 0) = diff[1];
      B(3, dd3 + 1) = diff[0];
      // gamma_yz
      B(4, dd3 + 1) = diff[2];
      B(4, dd3 + 2) = diff[1];
      // gamma_xz
      B(5, dd3 + 0) = diff[2];
      B(5, dd3 + 2) = diff[0];
    }

    MoFEMFunctionReturnHot(0);
  }

  /**
     * \brief Do calculations for give operator
     * @param  row_side row side number (local number) of entity on element
     * @param  col_side column side number (local number) of entity on element
     * @param  row_type type of row entity MBVERTEX, MBEDGE, MBTRI or MBTET
     * @param  col_type type of column entity MBVERTEX, MBEDGE, MBTRI or MBTET
     * @param  row_data data for row
     * @param  col_data data for column
     * @return          error code
     */
    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        DataForcesAndSourcesCore::EntData &row_data,
                        DataForcesAndSourcesCore::EntData &col_data) {
    
    MoFEMFunctionBeginHot;
    // get number of dofs on row
    nbRows = row_data.getIndices().size();
    // if no dofs on row, exit that work, nothing to do here
    if (!nbRows)
      MoFEMFunctionReturnHot(0);
    // get number of dofs on column
    nbCols = col_data.getIndices().size();
    // if no dofs on Columbia, exit nothing to do here
    if (!nbCols)
      MoFEMFunctionReturnHot(0);
    // get number of integration points
    nbIntegrationPts = getGaussPts().size2();
    // chekk if entity block is on matrix diagonal
    if (row_side == col_side && row_type == col_type) { // ASK
      isDiag = true;                                    // yes, it is
    } else {
      isDiag = false;
    }
    // integrate local matrix for entity block
    CHKERR iNtegrate(row_data, col_data);

    // asseble local matrix
    CHKERR aSsemble(row_data, col_data);

    MoFEMFunctionReturnHot(0);
  }

protected:
  ///< error code

  int nbRows;           ///< number of dofs on rows
  int nbCols;           ///< number if dof on column
  int nbIntegrationPts; ///< number of integration points
  bool isDiag;          ///< true if this block is on diagonal

  /**
     * \brief Integrate grad-grad operator
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
  virtual MoFEMErrorCode
  iNtegrate(DataForcesAndSourcesCore::EntData &row_data,
            DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBeginHot;

    int nb_dofs_row = row_data.getFieldData().size();
    if (nb_dofs_row == 0)
      MoFEMFunctionReturnHot(0);
    int nb_dofs_col = col_data.getFieldData().size();
    if (nb_dofs_col == 0)
      MoFEMFunctionReturnHot(0);

    K.resize(nb_dofs_row, nb_dofs_col, false);
    K.clear();
    CB.resize(6, nb_dofs_col, false);

    for (int gg = 0; gg != nbIntegrationPts; gg++) {
      // get element volume
      // get integration weight
      double val = getVolume() * getGaussPts()(3, gg);
      // ASK
      const MatrixAdaptor &diffN_row = row_data.getDiffN(gg, nb_dofs_row / 3);
      const MatrixAdaptor &diffN_col = col_data.getDiffN(gg, nb_dofs_col / 3);
      //     printf("\nGauss Point %d at row:\n", gg);
      CHKERR makeB(diffN_row, rowB);

      // printf("\nGauss Point %d at col:\n", gg);
      CHKERR makeB(diffN_col, colB);

      noalias(CB) = prod(D, colB);
      noalias(K) += val * prod(trans(rowB), CB);
    }
    MoFEMFunctionReturnHot(0);
  }

  /**
     * \brief Assemble local entity block matrix
     * @param  row_data row data (consist base functions on row entity)
     * @param  col_data column data (consist base functions on column entity)
     * @return          error code
     */
  virtual MoFEMErrorCode aSsemble(DataForcesAndSourcesCore::EntData &row_data,
                                  DataForcesAndSourcesCore::EntData &col_data) {
    MoFEMFunctionBeginHot;
    // get pointer to first global index on row
    const int *row_indices = &*row_data.getIndices().data().begin();
    // get pointer to first global index on column
    const int *col_indices = &*col_data.getIndices().data().begin();
    Mat B = getFEMethod()->ksp_B != PETSC_NULL ? getFEMethod()->ksp_B
                                               : getFEMethod()->snes_B;
    // assemble local matrix
    CHKERR MatSetValues(B, nbRows, row_indices, nbCols, col_indices,
                        &*K.data().begin(), ADD_VALUES);

    if (!isDiag && sYmm) {
      // if not diagonal term and since global matrix is symmetric assemble
      // transpose term.
      K = trans(K);
      CHKERR MatSetValues(B, nbCols, col_indices, nbRows, row_indices,
                          &*K.data().begin(), ADD_VALUES);
    }
    MoFEMFunctionReturnHot(0);
  }
};

struct ApplyDirichletBc : public MoFEM::FEMethod {

  Range fixFaces, fixNodes, fixSecondNode;

  ApplyDirichletBc(const Range &fix_faces, const Range &fix_nodes,
                   const Range &fix_second_node)
      : MoFEM::FEMethod(), fixFaces(fix_faces), fixNodes(fix_nodes),
        fixSecondNode(fix_second_node) {
    // constructor
  }

  MoFEMErrorCode postProcess() {

    MoFEMFunctionBeginHot;
    std::set<int> set_fix_dofs;
    cerr << fixFaces << endl;
    cerr << fixNodes << endl;
    cerr << fixSecondNode << endl;

    for (_IT_NUMEREDDOF_ROW_FOR_LOOP_(problemPtr, dit)) {
      if (dit->get()->getDofCoeffIdx() == 2) {
        if (fixFaces.find(dit->get()->getEnt()) != fixFaces.end()) {
          set_fix_dofs.insert(dit->get()->getPetscGlobalDofIdx());
        }
      }

      if (fixSecondNode.find(dit->get()->getEnt()) != fixSecondNode.end()) {
        if (dit->get()->getDofCoeffIdx() == 1) {
          printf("The extra node \n");
          set_fix_dofs.insert(dit->get()->getPetscGlobalDofIdx());
        }
      }

      if (fixNodes.find(dit->get()->getEnt()) != fixNodes.end()) {
        set_fix_dofs.insert(dit->get()->getPetscGlobalDofIdx());
      }
    }

    std::vector<int> fix_dofs(set_fix_dofs.size());

    std::copy(set_fix_dofs.begin(), set_fix_dofs.end(), fix_dofs.begin());

    CHKERR MatAssemblyBegin(ksp_B, MAT_FINAL_ASSEMBLY);

    CHKERR MatAssemblyEnd(ksp_B, MAT_FINAL_ASSEMBLY);

    CHKERR VecAssemblyBegin(ksp_f);

    CHKERR VecAssemblyEnd(ksp_f);

    Vec x;
    // ASK
    CHKERR VecDuplicate(ksp_f, &x);

    CHKERR VecZeroEntries(x);

    CHKERR MatZeroRowsColumns(ksp_B, fix_dofs.size(), &fix_dofs[0], 1, x,
                              ksp_f);

    CHKERR VecDestroy(&x);

   MoFEMFunctionReturnHot(0);
  }
};

struct OpPressure : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

  double pressureVal;

  OpPressure(const double pressure_val = 1)
      :

        MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator("U", OPROW),
        pressureVal(pressure_val) {}
  //  virtual ~OpPressure();

  VectorDouble nF;

  FTensor::Index<'i', 3> i;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBeginHot;

    const int nb_dofs = data.getIndices().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);

    nF.resize(nb_dofs, false);
    nF.clear();

    const int nb_gauss_pts = data.getN().size1();
    FTensor::Tensor1<double *, 3> t_normal = getTensor1Normal();
    FTensor::Tensor0<double *> t_base = data.getFTensor0N();

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      double w = 0.5 * getGaussPts()(2, gg); // ASK
      FTensor::Tensor1<double *, 3> t_nf(&nF[0], &nF[1], &nF[2], 3);
      for (int bb = 0; bb != nb_dofs / 3; bb++) {
        t_nf(i) += (w * pressureVal * t_base) * t_normal(i);
        ++t_nf;
        ++t_base;
      }
    }

    CHKERR VecSetValues(getFEMethod()->ksp_f, nb_dofs, &data.getIndices()[0],
                        &nF[0], ADD_VALUES);

    MoFEMFunctionReturnHot(0);
  }
};

int main(int argc, char *argv[]) {

  // Initialize PETSCc
  PetscInitialize(&argc, &argv, (char *)0, help);

  // Create mesh databse
  moab::Core mb_instance;              // create database
  moab::Interface &moab = mb_instance; // create interface to database

  try {
    // Create MoFEM database and link it to MoAB
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    CHKERR DMRegister_MoFEM("DMMOFEM");

    // Get command line options
    int order          = 3;           // default approximation order
    PetscBool flg_test = PETSC_FALSE; // true check if error is numerical error
    CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "", "SimpleElasticProblem",
                             "none");

    // Set approximation order
    CHKERR PetscOptionsInt("-order", "approximation order", "", order, &order,
                           PETSC_NULL);

    // Set testing (used by CTest)
    CHKERR PetscOptionsBool("-test", "if true is ctest", "", flg_test,
                            &flg_test, PETSC_NULL);

    ierr = PetscOptionsEnd();
CHKERRQ(ierr);

    Vec global_error;
    CHKERR PoissonExample::AuxFunctions(m_field).createGhostVec(&global_error);

    boost::shared_ptr<ForcesAndSourcesCore>
        domain_lhs_fe; ///< Volume element for the matrix
    boost::shared_ptr<ForcesAndSourcesCore>
        domain_rhs_fe; ///< Volume element to assemble vector
    boost::shared_ptr<ForcesAndSourcesCore>
        domain_error; ///< Volume element evaluate error
    boost::shared_ptr<ForcesAndSourcesCore>
        post_proc_volume; ///< Volume element to Post-process results
    boost::shared_ptr<ForcesAndSourcesCore> null;

    Simple *simple_interface = m_field.getInterface<MoFEM::Simple>();
//    CHKERR m_field.query_interface(simple_interface);

    CHKERR simple_interface->getOptions();

    CHKERR simple_interface->loadFile();

    CHKERR simple_interface->addDomainField("U", H1, AINSWORTH_LEGENDRE_BASE,
                                            3);

    CHKERR simple_interface->setFieldOrder("U", order);

    Range fix_faces, pressure_faces, fix_nodes, fix_second_node;

    for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field, BLOCKSET, bit)) {
      EntityHandle meshset = bit->getMeshset();
      int id               = bit->getMeshsetId();

      if (id == 1) { // brick-faces

        CHKERR m_field.get_moab().get_entities_by_dimension(meshset, 2,
                                                            fix_faces, true);

        Range adj_ents;
        CHKERR m_field.get_moab().get_adjacencies(fix_faces, 0, false, adj_ents,
                                                  moab::Interface::UNION);

        CHKERR m_field.get_moab().get_adjacencies(fix_faces, 1, false, adj_ents,
                                                  moab::Interface::UNION);
        fix_faces.merge(adj_ents);
      } else if (id == 2) { // node(s)

        CHKERR m_field.get_moab().get_entities_by_dimension(meshset, 0,
                                                            fix_nodes, true);
        
      } else if (id == 3) { // brick pressure faces
        CHKERR m_field.get_moab().get_entities_by_dimension(
            meshset, 2, pressure_faces, true);
        
      } else if (id == 4) { // restrained second node in y direction
        CHKERR m_field.get_moab().get_entities_by_dimension(
            meshset, 0, fix_second_node, true);

      } else {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "Unknown blockset");
      }
    }

    CHKERR simple_interface->addDomainField("U", H1, AINSWORTH_LEGENDRE_BASE,
                                            3);

    CHKERR simple_interface->setFieldOrder("U", order);

    CHKERR simple_interface->defineFiniteElements();

    // Add pressure element
    CHKERR m_field.add_finite_element("PRESSURE");

    CHKERR m_field.modify_finite_element_add_field_row("PRESSURE", "U");

    CHKERR m_field.modify_finite_element_add_field_col("PRESSURE", "U");

    CHKERR m_field.modify_finite_element_add_field_data("PRESSURE", "U");

    CHKERR simple_interface->defineProblem();

    DM dm;
    CHKERR simple_interface->getDM(&dm);

    CHKERR DMMoFEMAddElement(dm, "PRESSURE");

    CHKERR DMMoFEMSetIsPartitioned(dm, PETSC_TRUE);

    CHKERR simple_interface->buildFields();

    CHKERR simple_interface->buildFiniteElements();

    CHKERR m_field.add_ents_to_finite_element_by_dim(
        0, simple_interface->getDim(), simple_interface->getDomainFEName(),
        true);

    CHKERR m_field.build_finite_elements(simple_interface->getDomainFEName());

    CHKERR m_field.add_ents_to_finite_element_by_dim(pressure_faces, 2,
                                                     "PRESSURE");

    CHKERR m_field.build_finite_elements("PRESSURE", &pressure_faces);

    CHKERR simple_interface->buildProblem();

    boost::shared_ptr<VolumeElementForcesAndSourcesCore> elastic_fe(
        new VolumeElementForcesAndSourcesCore(m_field));

    elastic_fe->getOpPtrVector().push_back(new OpK());

    // push operators to elastic_fe
    boost::shared_ptr<FaceElementForcesAndSourcesCore> pressure_fe(
        new FaceElementForcesAndSourcesCore(m_field));
    pressure_fe->getOpPtrVector().push_back(new OpPressure());

    boost::shared_ptr<FEMethod> fix_dofs_fe(
        new ApplyDirichletBc(fix_faces, fix_nodes, fix_second_node));

    boost::shared_ptr<FEMethod> nullFE;

    // Set operators for KSP solver
    CHKERR DMMoFEMKSPSetComputeOperators(
        dm, simple_interface->getDomainFEName(), elastic_fe, nullFE, nullFE);

    CHKERR DMMoFEMKSPSetComputeRHS(dm, "PRESSURE", pressure_fe, nullFE, nullFE);

    Mat A;
    Vec x, f; // ASK

    CHKERR DMCreateMatrix(dm, &A);

    CHKERR DMCreateGlobalVector(dm, &x);

    CHKERR VecDuplicate(x, &f);

    fix_dofs_fe->ksp_B = A;
    fix_dofs_fe->ksp_f = f;
    elastic_fe->ksp_B  = A;
    pressure_fe->ksp_f = f;

    CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getDomainFEName(),
                                    elastic_fe);

    CHKERR DMoFEMLoopFiniteElements(dm, "PRESSURE", pressure_fe);

    // This is done because only post processor is implemented in the
    // ApplyDirichletBc struct
    CHKERR DMoFEMPostProcessFiniteElements(dm, fix_dofs_fe.get());

    KSP solver;

    CHKERR KSPCreate(PETSC_COMM_WORLD, &solver);

    CHKERR KSPSetFromOptions(solver);

    CHKERR KSPSetOperators(solver, A, A);

    CHKERR KSPSetUp(solver);

    CHKERR KSPSolve(solver, f, x);

    VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    // save solution in vector x on mesh
    CHKERR DMoFEMMeshToGlobalVector(dm, x, INSERT_VALUES, SCATTER_REVERSE);

    // Set up post-procesor. It is some generic implementation of finite
    // element.
    PostProcVolumeOnRefinedMesh post_proc(m_field);
    // Add operators to the elements, starting with some generic
    CHKERR post_proc.generateReferenceElementMesh();

    CHKERR post_proc.addFieldValuesPostProc("U");

    CHKERR post_proc.addFieldValuesGradientPostProc("U");

    CHKERR DMoFEMLoopFiniteElements(
        dm, simple_interface->getDomainFEName().c_str(), &post_proc);

    CHKERR post_proc.writeFile("out.h5m");

    CHKERR MatDestroy(&A);

    CHKERR VecDestroy(&x);

    CHKERR VecDestroy(&f);

    CHKERR DMDestroy(&dm);

    // This is a good reference for the future
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF, e.errorCode, e.errorMessage);
  }

  // finish work cleaning memory, getting statistics, etc
   PetscFinalize();

  return 0;
}
