/** \file field_evaluator.cpp
 * \example field_evaluator.cpp
 * \brief test segments distance
 *
 * \ingroup field_evaluator
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

static char help[] = "testing field evaluator\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // Read mesh and create moab and mofem data structures
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Simple interface
    Simple *simple_interface;
    CHKERR m_field.getInterface(simple_interface);
    {

      // get options from command line
      CHKERR simple_interface->getOptions();
      // load mesh file
      CHKERR simple_interface->loadFile();
      // add fields
      CHKERR simple_interface->addDomainField("FIELD1", H1,
                                              AINSWORTH_LEGENDRE_BASE, 3);
      // set fields order
      CHKERR simple_interface->setFieldOrder("FIELD1", 1);
      // setup problem
      CHKERR simple_interface->setUp();

      FieldEvaluatorInterface *field_eval_ptr;
      CHKERR m_field.getInterface(field_eval_ptr);

      std::array<double, 3> point = {0, 0, 0};
      const double dist = 0.1;

      std::array<double, 12> eval_points = {0.0,  0.0, 0.0, 0.0, 0.01, 0.0,
                                            -0.1, 0.0, 0.0, 0.1, 0.0,  0.0};

      auto dm = simple_interface->getDM();
      const MoFEM::Problem *prb_ptr;
      CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

      using VolEle = VolumeElementForcesAndSourcesCore;
      using VolOp = VolumeElementForcesAndSourcesCore::UserDataOperator;

      /**
       * @brief Operator used to check consistency between local coordinates and
       * global cooridnates for integrated points and evaluated points
       *
       */
      struct MyOp : public VolOp {

        MatrixShallowArrayAdaptor<double> evalPoints;

        MyOp(std::array<double, 12> &eval_points)
            : VolOp("FIELD1", OPROW),
              evalPoints(getMatrixAdaptor(eval_points.data(),
                                          eval_points.size() / 3, 3)) {}

        MoFEMErrorCode doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {
          MoFEMFunctionBegin;
          if (type == MBVERTEX) {

            std::cout << "Integration pts" << std::endl;
            std::cout << getGaussPts() << endl;

            std::cout << "Global coordinates " << endl;
            std::cout << getCoordsAtGaussPts() << std::endl;

            for (int gg = 0; gg != getCoordsAtGaussPts().size1(); ++gg) {
              int pt_number = getGaussPts()(3, gg);

              std::cout << "gg " << gg << std::endl;
              std::cout << "pt " << pt_number << std::endl;

              ublas::matrix_row<MatrixDouble> coord_at_gauss_pt(
                  getCoordsAtGaussPts(), gg);
              ublas::matrix_row<MatrixShallowArrayAdaptor<double>> eval_coord(
                  evalPoints, pt_number);

              std::cout << "coord_at_gauss_pt ";
              std::cout << coord_at_gauss_pt << std::endl;

              std::cout << "eval_coord ";
              std::cout << eval_coord << std::endl;

              double error = norm_2(coord_at_gauss_pt - eval_coord);

              if (error > 1e-12)
                SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                         "Difference at %d error = %3.4e", pt_number, error);
            }
          }
          MoFEMFunctionReturn(0);
        }
      };

      // Lambda function is used to set integration rule to -1, that indicates
      // that finite element instace will use non-standard integration points.
      auto get_rule = [&](int order_row, int order_col, int order_data) {
        return -1;
      };

      // Get pointer of FieldEvaluator data. Note finite element and method
      // set integrating points is destroyed when this pointer is releases
      auto data = field_eval_ptr->getData<VolEle>();
      
      // Set operators and integration rule
      if (auto fe_method = data->feMethodPtr.lock()) {
        fe_method->getRuleHook = get_rule;
        fe_method->getOpPtrVector().push_back(new MyOp(eval_points));
      } else
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Pointer to element does not exists");

      // Build tree for particular element
      CHKERR field_eval_ptr->buildTree3D(data,
                                         simple_interface->getDomainFEName());
      // Set points to set on finite elements
      data->setEvalPoints(eval_points.data(), eval_points.size() / 3);

      // Evaluate points on finite elements
      CHKERR field_eval_ptr->evalFEAtThePoint3D(
          &point[0], dist, prb_ptr->getName(),
          simple_interface->getDomainFEName(), data, m_field.get_comm_rank(),
          m_field.get_comm_rank(), nullptr, MF_EXIST, VERY_NOISY);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
