/** \file segments_distance.cpp
 * \brief test segments distance
 *
 * \ingroup mesh_cut
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "testing mesh cut test\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    auto run_test = [](auto w, auto v, auto k, auto l, auto expected_result,
                       double expected_t0, double expected_t1) {
      MoFEMFunctionBegin;
      double t0 = 0, t1 = 0;
      auto result = Tools::minDistanceFromSegments(w, v, k, l, &t0, &t1);
      FTensor::Tensor1<const double *, 3> t_w(w, &w[1], &w[2]);
      FTensor::Tensor1<const double *, 3> t_v(v, &v[1], &v[2]);
      FTensor::Tensor1<const double *, 3> t_k(k, &k[1], &k[2]);
      FTensor::Tensor1<const double *, 3> t_l(l, &l[1], &l[2]);
      FTensor::Tensor1<double, 3> t_delta;
      FTensor::Index<'i',3> i;
      t_delta(i) =
          (t_w(i) + t0 * (t_v(i) - t_w(i))) - (t_k(i) + t1 * (t_l(i) - t_k(i)));
      std::cout << "Result " << result << " : " << t0 << " " << t1 << " dist "
                << sqrt(t_delta(i) * t_delta(i)) << std::endl;

      if (result != expected_result)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Expected solution");
      if (fabs(t0 - expected_t0) > 1e-12)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong value of t0 %3.4e != %3.4e", t0, expected_t0);
      if (fabs(t1 - expected_t1) > 1e-12)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong value of t1 %3.4e != %3.4e", t1, expected_t1);
      MoFEMFunctionReturn(0);
    };

    {
      const double w[] = {-1, 0, 0};
      const double v[] = {1, 0, 0};
      const double k[] = {0, -1, 1};
      const double l[] = {0, 1, 1};
      CHKERR run_test(w, v, k, l, Tools::SOLUTION_EXIST, 0.5, 0.5);
    }

    {
      const double w[] = {-1, 0, 0};
      const double v[] = {1, 0, 0};
      const double k[] = {0, 1, 0};
      const double l[] = {0, 2, 0};
      CHKERR run_test(w, v, k, l, Tools::SOLUTION_EXIST, 0.5, -1);
    }

    {
      const double w[] = {-1, 0, 0};
      const double v[] = {1, 0, 0};
      const double k[] = {-1, -1, 0};
      const double l[] = {1, 1, 0};
      CHKERR run_test(w, v, k, l, Tools::SOLUTION_EXIST, 0.5, 0.5);
    }

    {
      const double w[] = {-1, 0, 1};
      const double v[] = {1, 0, 1};
      const double k[] = {-1, -1, 0};
      const double l[] = {1, 1, 0};
      CHKERR run_test(w, v, k, l, Tools::SOLUTION_EXIST, 0.5, 0.5);
    }

    {
      const double w[] = {0, 0, 0};
      const double v[] = {0, 0, 0};
      const double k[] = {0, 0, 0};
      const double l[] = {0, 1, 0};
      CHKERR run_test(w, v, k, l, Tools::SEGMENT_ONE_IS_POINT, 0, 0);
    }

    {
      const double w[] = {1, 0, 0};
      const double v[] = {1, 0, 0};
      const double k[] = {0, 0, 0};
      const double l[] = {1, 0, 0};
      CHKERR run_test(w, v, k, l, Tools::SEGMENT_ONE_IS_POINT, 0, 1);
    }

    {
      const double w[] = {0, 0, 0};
      const double v[] = {1, 0, 0};
      const double k[] = {1, 0, 0};
      const double l[] = {1, 0, 0};
      CHKERR run_test(w, v, k, l, Tools::SEGMENT_TWO_IS_POINT, 1, 0);
    }

    {
      const double w[] = {-1, 0, 0};
      const double v[] = {1, 0, 0};
      const double k[] = {0, 1, 0};
      const double l[] = {0, 1, 0};
      CHKERR run_test(w, v, k, l, Tools::SEGMENT_TWO_IS_POINT, 0.5, 0);
    }

    {
      const double w[] = {0, 0, 0};
      const double v[] = {0, 0, 0};
      const double k[] = {1, 1, -1};
      const double l[] = {1, 1, 1};
      CHKERR run_test(w, v, k, l, Tools::SEGMENT_ONE_IS_POINT, 0, 0.5);
    }

    {
      const double w[] = {0, 0, 0};
      const double v[] = {0, 0, 1};
      const double k[] = {0, 0, 0};
      const double l[] = {0, 0, 1};
      CHKERR run_test(w, v, k, l, Tools::NO_SOLUTION, 0, 0);
    }

    {
      const double w[] = {1, 0, 0};
      const double v[] = {1, 0, 1};
      const double k[] = {0, 0, 0};
      const double l[] = {0, 0, 1};
      CHKERR run_test(w, v, k, l, Tools::NO_SOLUTION, 0, 0);
    }

    {
      const double w[] = {1, -1, 0};
      const double v[] = {1, 1, 1};
      const double k[] = {0, 1, 0};
      const double l[] = {0, -1, 1};
      CHKERR run_test(w, v, k, l, Tools::SOLUTION_EXIST, 0.5, 0.5);
    }

    {
      const double w[] = {0, 1, 0};
      const double v[] = {0, 2, 0};
      const double k[] = {0, 0, 0};
      const double l[] = {1, 0, 0};
      CHKERR run_test(w, v, k, l, Tools::SOLUTION_EXIST, -1, 0);
    }

    {
      const double w[] = {0, 0, 0};
      const double v[] = {1, 0, 0};
      const double k[] = {0, 1, 0};
      const double l[] = {0, 2, 0};
      CHKERR run_test(w, v, k, l, Tools::SOLUTION_EXIST, 0, -1);
    }


  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
