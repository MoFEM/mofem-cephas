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
\param diff_s derivatives of shape functions, i.e. \f$\frac{\partial s}{\partial \xi_i}\f$
\retval L approximation functions
\retval diffL derivatives, i.e. \f$\frac{\partial L}{\partial \xi_i}\f$
\param dim dimension
\return error code

*/
PetscErrorCode Legendre_polynomials(
  int p,double s,double *diff_s,double *L,double *diffL,const int dim
);


/**
 * \brief Calculate Lobatto base functions

 \param p is approximation order
 \param s is position \f$s\in[-1,1]\f$
 \param diff_s derivatives of shape functions, i.e. \f$\frac{\partial s}{\partial \xi_i}\f$
 \retval L approximation functions
 \retval diffL derivatives, i.e. \f$\frac{\partial L}{\partial \xi_i}\f$
 \param dim dimension
 \return error code

*/
PetscErrorCode Lobatto_polynomials(
  int p,double s,double *diff_s,double *L,double *diffL,const int dim
);

#ifdef __cplusplus
}
#endif

#endif //__BASE_FUNCTIONS_H__
