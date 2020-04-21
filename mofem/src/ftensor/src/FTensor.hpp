/** \file FTensor.hpp
* \brief Tensors class implemented by Walter Landry.
* \ingroup ftensor

*
* The main repository is available under link
<https://bitbucket.org/wlandry/ftensor>.
* For more details look at reference \cite landry2003implementing.
* See link <http://www.wlandry.net/Presentations/FTensor.pdf>.
*
* Note this implementation is modified for proposes of MoFEM.
*
*/

/* Include file for the Fast Tensor classes (FTensor).  Everything is
   in namespace FTensor. */

#pragma once

#include <cmath>
#include <complex>
#include <type_traits>

#ifdef FTENSOR_DEBUG
#include <iostream>
#include <sstream>
#include <stdexcept>
#endif
#include "FTensor/Layout.hpp"

#ifdef WITH_ADOL_C
#include <adolc/adolc.h>
#include <adolc/adtl.h>
#endif

/**
 * \brief Tensors class implemented by Walter Landry.
 * \ingroup ftensor


 * The main repository is available under link
 <https://bitbucket.org/wlandry/ftensor>.
 * For more details look at reference \cite landry2003implementing.
 * See link <http://www.wlandry.net/Presentations/FTensor.pdf>.
 *
 * Note this implementation is modified for proposes of MoFEM.

 */
namespace FTensor
{
  template <class T> class Tensor0;

  template <class T, const int I> class PackPtr
  {};

  template <class T, int Dim> class Tensor1;
  template <class A, class T, int Dim, char i> class Tensor1_Expr;

  template <class T, int Dim1, int Dim2>
  class Tensor2;
  template <class A, class T, int Dim1, int Dim2, char i, char j>
  class Tensor2_Expr;
  template <class A, class T, int N> class Tensor2_number_rhs_0;
  template <class A, class T, int N> class Tensor2_number_rhs_1;

  template <class T, int Dim> class Tensor2_symmetric;
  template <class A, class T, int Dim, char i, char j>
  class Tensor2_symmetric_Expr;

  template <class T, int Dim> class Tensor2_antisymmetric;
  template <class A, class T, int Dim, char i, char j>
  class Tensor2_antisymmetric_Expr;

  template <class T, int Dim0, int Dim1, int Dim2> class Tensor3;
  template <class A, class T, int Dim0, int Dim1, int Dim2, char i, char j,
            char k>
  class Tensor3_Expr;
  template <class A, class T, int N> class Tensor3_number_rhs_0;
  template <class A, class T, int N> class Tensor3_number_rhs_2;
  template <class A, class T, int N> class Tensor3_number_rhs_2;
  template <class A, class T, int N1, int N2> class Tensor3_number_rhs_01;
  template <class A, class T, int N1, int N2> class Tensor3_number_rhs_02;
  template <class A, class T, int N1, int N2> class Tensor3_number_rhs_12;

  template <class T, int Dim01, int Dim2> class Dg;
  template <class A, class T, int Dim01, int Dim2, char i, char j, char k>
  class Dg_Expr;
  template <class A, class T, int N> class Dg_number_rhs_0;
  template <class A, class T, int N> class Dg_number_rhs_2;
  template <class A, class T, int N1, int N2> class Dg_number_rhs_01;
  template <class A, class T, int N1, int N2> class Dg_number_rhs_12;

  template <class T, int Dim0, int Dim12> class Christof;
  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  class Christof_Expr;

  template <class T, int Dim0, int Dim12> class Tensor3_antisymmetric;
  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  class Tensor3_antisymmetric_Expr;

  template <class T, int Dim0, int Dim1, int Dim2, int Dim3> class Tensor4;
  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3, char i,
            char j, char k, char l>
  class Tensor4_Expr;

  template <class T, int Dim> class Riemann;
  template <class A, class T, int Dim, char i, char j, char k, char l>
  class Riemann_Expr;

  template <class T, int Dim01, int Dim23> class Ddg;
  template <class A, class T, int Dim01, int Dim23, char i, char j, char k,
            char l>
  class Ddg_Expr;
  template <class A, class T, int N0, int N1> class Ddg_number_rhs_01;
  template <class A, class T, int N0> class Ddg_number_rhs_0;

}

#include "FTensor/Index.hpp"
#include "FTensor/Number.hpp"
#include "FTensor/promote.hpp"

#include "FTensor/Levi_Civita.hpp"
#include "FTensor/cross.hpp"
#include "FTensor/Kronecker_Delta.hpp"

#include "FTensor/Christof.hpp"
#include "FTensor/Ddg.hpp"
#include "FTensor/Dg.hpp"
#include "FTensor/Riemann.hpp"
#include "FTensor/Tensor0.hpp"
#include "FTensor/Tensor1.hpp"
#include "FTensor/Tensor2.hpp"
#include "FTensor/Tensor2_antisymmetric.hpp"
#include "FTensor/Tensor2_symmetric.hpp"
#include "FTensor/Tensor3.hpp"
#include "FTensor/Tensor3/Tensor3_contracted.hpp"
#include "FTensor/Tensor3_antisymmetric.hpp"
#include "FTensor/Tensor4.hpp"
#include "FTensor/Tensor4/Tensor4_contracted.hpp"

/***************************************************************************/ /**
* \defgroup ftensor Tensor template library
* \brief Efficient Template Tensor library
******************************************************************************/
