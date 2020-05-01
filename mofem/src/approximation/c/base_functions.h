/** \file base_functions.c

*/

/**
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __BASE_FUNCTIONS_H__
#define __BASE_FUNCTIONS_H__

#ifdef __cplusplus
extern "C" {
#endif

/**
\brief Calculate Legendre approximation basis

\ingroup mofem_base_functions

Lagrange polynomial is given by
\f[
L_0(s)=1;\quad L_1(s) = s
\f]
and following terms are generated inductively
\f[
L_{l+1}=\frac{2l+1}{l+1}sL_l(s)-\frac{l}{l+1}L_{l-1}(s)
\f]

Note that:
\f[
s\in[-1,1] \quad \textrm{and}\; s=s(\xi_0,\xi_1,\xi_2)
\f]
where \f$\xi_i\f$ are barycentric coordinates of element.

\param p is approximation order
\param s is position \f$s\in[-1,1]\f$
\param diff_s derivatives of shape functions, i.e. \f$\frac{\partial s}{\partial
\xi_i}\f$ \retval L approximation functions \retval diffL derivatives, i.e.
\f$\frac{\partial L}{\partial \xi_i}\f$ \param dim dimension \return error code

*/
PetscErrorCode Legendre_polynomials(int p, double s, double *diff_s, double *L,
                                    double *diffL, const int dim);

/**
\brief Calculate Jacobi approximation basis

For more details see \cite fuentes2015orientation

\param p is approximation order
\param alpha polynomial parameter
\param x is position \f$s\in[0,t]\f$
\param t range of polynomial
\param diff_x derivatives of shape functions, i.e. \f$\frac{\partial x}{\partial
\xi_i}\f$ \param diff_t derivatives of shape functions, i.e. \f$\frac{\partial
t}{\partial \xi_i}\f$ \retval L approximation functions \retval diffL
derivatives, i.e. \f$\frac{\partial L}{\partial \xi_i}\f$ \param dim dimension
\return error code

*/
PetscErrorCode Jacobi_polynomials(int p, double alpha, double x, double t,
                                  double *diff_x, double *diff_t, double *L,
                                  double *diffL, const int dim);

/**
\brief Calculate integrated Jacobi approximation basis

For more details see \cite fuentes2015orientation

\param p is approximation order
\param alpha polynomial parameter
\param x is position \f$s\in[0,t]\f$
\param t range of polynomial
\param diff_x derivatives of shape functions, i.e. \f$\frac{\partial x}{\partial
\xi_i}\f$ \param diff_t derivatives of shape functions, i.e. \f$\frac{\partial
t}{\partial \xi_i}\f$ \retval L approximation functions \retval diffL
derivatives, i.e. \f$\frac{\partial L}{\partial \xi_i}\f$ \param dim dimension
\return error code

*/
PetscErrorCode IntegratedJacobi_polynomials(int p, double alpha, double x,
                                            double t, double *diff_x,
                                            double *diff_t, double *L,
                                            double *diffL, const int dim);

/**
 \brief Calculate Lobatto base functions \cite FUENTES2015353.

 \ingroup mofem_base_functions

 Order of first function is 2.

 \param p is approximation order
 \param s is a parametrisation of position on segment
 \f$s(\xi_1,\cdots, \xi_{dim})\in[-1,1]\f$ \param diff_s derivatives of shape
 functions, i.e. \f$\frac{\partial s}{\partial \xi_i}\f$ \retval L approximation
 functions \retval diffL derivatives, i.e. \f$\frac{\partial L}{\partial
 \xi_i}\f$
 \param dim dimension of ambient space where the segment is embedded
 \return error code
 - output
 \param L values of basis functions
 \param diffL derivatives of basis functions

*/
PetscErrorCode Lobatto_polynomials(int p, double s, double *diff_s, double *L,
                                   double *diffL, const int dim);

/**
 \brief Calculate Kernel Lobatto base functions.

 \ingroup mofem_base_functions

 This is implemented using definitions from Hermes2d
 <https://github.com/hpfem/hermes> following book by Pavel Solin et al \cite
 solin2003higher.

 \param p is approximation order
 \param s is position \f$s\in[-1,1]\f$
 \param diff_s derivatives of shape functions, i.e. \f$\frac{\partial
 s}{\partial \xi_i}\f$ \retval L approximation functions \retval diffL
 derivatives, i.e. \f$\frac{\partial L}{\partial \xi_i}\f$ \param dim dimension
 \return error code

*/
PetscErrorCode LobattoKernel_polynomials(int p, double s, double *diff_s,
                                         double *L, double *diffL,
                                         const int dim);

/// Definitions taken from Hermes2d code

/// kernel functions for Lobatto base
#define LOBATTO_PHI0(x) (-2.0 * 1.22474487139158904909864203735)
#define LOBATTO_PHI1(x) (-2.0 * 1.58113883008418966599944677222 * (x))
#define LOBATTO_PHI2(x)                                                        \
  (-1.0 / 2.0 * 1.87082869338697069279187436616 * (5 * (x) * (x)-1))
#define LOBATTO_PHI3(x)                                                        \
  (-1.0 / 2.0 * 2.12132034355964257320253308631 * (7 * (x) * (x)-3) * (x))
#define LOBATTO_PHI4(x)                                                        \
  (-1.0 / 4.0 * 2.34520787991171477728281505677 *                              \
   (21 * (x) * (x) * (x) * (x)-14 * (x) * (x) + 1))
#define LOBATTO_PHI5(x)                                                        \
  (-1.0 / 4.0 * 2.54950975679639241501411205451 *                              \
   ((33 * (x) * (x)-30) * (x) * (x) + 5) * (x))
#define LOBATTO_PHI6(x)                                                        \
  (-1.0 / 32.0 * 2.73861278752583056728484891400 *                             \
   (((429 * (x) * (x)-495) * (x) * (x) + 135) * (x) * (x)-5))
#define LOBATTO_PHI7(x)                                                        \
  (-1.0 / 32.0 * 2.91547594742265023543707643877 *                             \
   (((715 * (x) * (x)-1001) * (x) * (x) + 385) * (x) * (x)-35) * (x))
#define LOBATTO_PHI8(x)                                                        \
  (-1.0 / 64.0 * 3.08220700148448822512509619073 *                             \
   ((((2431 * (x) * (x)-4004) * (x) * (x) + 2002) * (x) * (x)-308) * (x) *     \
        (x) +                                                                  \
    7))
#define LOBATTO_PHI9(x)                                                        \
  (-1.0 / 128.0 * 6.4807406984078603784382721642 *                             \
   ((((4199 * (x) * (x)-7956) * (x) * (x) + 4914) * (x) * (x)-1092) * (x) *    \
        (x) +                                                                  \
    63) *                                                                      \
   (x))

/// Derivatives of kernel functions for Lobbatto base
#define LOBATTO_PHI0X(x) (0)
#define LOBATTO_PHI1X(x) (-2.0 * 1.58113883008418966599944677222)
#define LOBATTO_PHI2X(x)                                                       \
  (-1.0 / 2.0 * 1.87082869338697069279187436616 * (10 * (x)))
#define LOBATTO_PHI3X(x)                                                       \
  (-1.0 / 2.0 * 2.12132034355964257320253308631 * (21.0 * (x) * (x)-3.0))
#define LOBATTO_PHI4X(x)                                                       \
  (-1.0 / 4.0 * 2.34520787991171477728281505677 *                              \
   ((84.0 * (x) * (x)-28.0) * (x)))
#define LOBATTO_PHI5X(x)                                                       \
  (-1.0 / 4.0 * 2.54950975679639241501411205451 *                              \
   ((165.0 * (x) * (x)-90.0) * (x) * (x) + 5.0))
#define LOBATTO_PHI6X(x)                                                       \
  (-1.0 / 32.0 * 2.73861278752583056728484891400 *                             \
   (((2574.0 * (x) * (x)-1980.0) * (x) * (x) + 270.0) * (x)))
#define LOBATTO_PHI7X(x)                                                       \
  (-1.0 / 32.0 * 2.91547594742265023543707643877 *                             \
   (((5005.0 * (x) * (x)-5005.0) * (x) * (x) + 1155.0) * (x) * (x)-35.0))
#define LOBATTO_PHI8X(x)                                                       \
  (-1.0 / 64.0 * 3.08220700148448822512509619073 *                             \
   ((((19448.0 * (x) * (x)-24024.0) * (x) * (x) + 8008.0) * (x) * (x)-616.0) * \
    (x)))
#define LOBATTO_PHI9X(x)                                                       \
  (-1.0 / 128.0 * 6.4807406984078603784382721642 *                             \
   ((((37791.0 * (x) * (x)-55692.0) * (x) * (x) + 24570.0) * (x) *             \
     (x)-3276.0) *                                                             \
    (x) * (x)-63.0))

#ifdef __cplusplus
}
#endif

#endif //__BASE_FUNCTIONS_H__
