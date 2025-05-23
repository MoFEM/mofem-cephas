# include "tetrahedron_ncc_rule.h"

/******************************************************************************/

double r8mat_det_4d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[4*4], the matrix whose determinant is desired.

    Output, double R8MAT_DET_4D, the determinant of the matrix.
*/
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
/******************************************************************************/

void reference_to_physical_t4 ( double t[], int n, double ref[], double phy[] )

/******************************************************************************/
/*
  Purpose:

    REFERENCE_TO_PHYSICAL_T4 maps T4 reference points to physical points.

  Discussion:

    Given the vertices of an order 4 physical tetrahedron and a point
    (R,S,T) in the reference tetrahedron, the routine computes the value
    of the corresponding image point (X,Y,Z) in physical space.

    This routine will also be correct for an order 10 tetrahedron,
    if the mapping between reference and physical space
    is linear.  This implies, in particular, that the sides of the
    image tetrahedron are straight, the faces are flat, and
    the "midside" nodes in the physical tetrahedron are
    halfway along the edges of the physical tetrahedron.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, double T[3*4], the coordinates of the vertices.
    The vertices are assumed to be the images of (0,0,0), (1,0,0),
    (0,1,0) and (0,0,1) respectively.

    Input, int N, the number of objects to transform.

    Input, double REF[3*N], points in the reference tetrahedron.

    Output, double PHY[3*N], corresponding points in the
    physical tetrahedron.
*/
{
  int i;
  int j;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*3] = t[i+0*3] * ( 1.0 - ref[0+j*3] - ref[1+j*3] - ref[2+j*3] )
                 + t[i+1*3] *       + ref[0+j*3]
                 + t[i+2*3] *                    + ref[1+j*3]
                 + t[i+3*3] *                                 + ref[2+j*3];
    }
  }

  return;
}
/******************************************************************************/

int tetrahedron_ncc_degree ( int rule )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_DEGREE: degree of an NCC rule for the tetrahedron.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int RULE, the index of the rule.

    Output, int TETRAHEDRON_NCC_DEGREE, the polynomial degree of exactness of
    the rule.
*/
{
  int degree;

  if ( 1 <= rule && rule <= 7 )
  {
    degree = rule - 1;
  }
  else
  {
    degree = -1;
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TETRAHEDRON_NCC_DEGREE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal RULE = %d\n", rule );
    exit ( 1 );
  }

  return degree;
}
/******************************************************************************/

int tetrahedron_ncc_order_num ( int rule )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_ORDER_NUM: order of an NCC rule for the tetrahedron.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int RULE, the index of the rule.

    Output, int TETRAHEDRON_NCC_ORDER_NUM, the order (number of points)
    of the rule.
*/
{
  int order;
  int order_num;
  int *suborder;
  int suborder_num;

  suborder_num = tetrahedron_ncc_suborder_num ( rule );

  suborder = tetrahedron_ncc_suborder ( rule, suborder_num );

  order_num = 0;
  for ( order = 0; order < suborder_num; order++ )
  {
    order_num = order_num + suborder[order];
  }

  free ( suborder );

  return order_num;
}
/******************************************************************************/

void tetrahedron_ncc_rule ( int rule, int order_num, double xyz[], double w[] )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_RULE returns the points and weights of an NCC rule.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int RULE, the index of the rule.

    Input, int ORDER_NUM, the order (number of points) of the rule.

    Output, double XYZ[3*ORDER_NUM], the points of the rule.

    Output, double W[ORDER_NUM], the weights of the rule.
*/
{
  int o;
  int s;
  int *suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyz;
/*
  Get the suborder information.
*/
  suborder_num = tetrahedron_ncc_suborder_num ( rule );

  suborder_xyz = ( double * ) malloc ( 4 * suborder_num * sizeof ( double ) );
  suborder_w = ( double * ) malloc ( suborder_num * sizeof ( double ) );

  suborder = tetrahedron_ncc_suborder ( rule, suborder_num );

  tetrahedron_ncc_subrule ( rule, suborder_num, suborder_xyz, suborder_w );
/*
  Expand the suborder information to a full order rule.
*/
  o = 0;

  for ( s = 0; s < suborder_num; s++ )
  {
    if ( suborder[s] == 1 )
    {
      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;
    }
/*
  Fourfold symmetry on (A,A,A,B)

    123 AAA
    124 AAB
    142 ABA
    412 BAA
*/
    else if ( suborder[s] == 4 )
    {
      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;
    }
/*
  Sixfold symmetry on (A,A,B,B):

    123 (A,A,B)
    132 (A,B,A),
    134 (A,B,B)
    312 (B,A,A)
    314 (B,A,B)
    341 (B,B,A)
*/
    else if ( suborder[s] == 6 )
    {
      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;
    }
/*
  Twelvefold symmetry on (A,A,B,C):

    123 (A,A,B)
    124 (A,A,C)
    132 (A,B,A)
    134 (A,B,C)
    142 (A,C,A)
    143 (A,C,B)
    312 (B,A,A)
    314 (B,A,C)
    341 (B,C,A)
    412 (C,A,A)
    413 (C,A,B)
    431 (C,B,A)
*/
    else if ( suborder[s] == 12 )
    {
      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;
    }
/*
  24 fold symmetry on (A,B,C,D):

    123 (A,B,C)
    124 (A,B,D)
    132 (A,C,B)
    134 (A,C,D)
    142 (A,D,B)
    143 (A,D,C)
    213 (B,A,C)
    214 (B,A,D)
    231 (B,C,A)
    234 (B,C,D)
    241 (B,D,A)
    243 (B,D,C)
    312 (C,A,B)
    314 (C,A,D)
    321 (C,B,A)
    324 (C,B,D)
    341 (C,D,A)
    342 (C,D,B)
    412 (D,A,B)
    413 (D,A,C)
    421 (D,B,A)
    423 (D,B,C)
    431 (D,C,A)
    432 (D,C,B)
*/
    else if ( suborder[s] == 24 )
    {
      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[0+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[1+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[1+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[4+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[1+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[1+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[1+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[1+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[3+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[2+s*4];
      xyz[1+o*3] = suborder_xyz[3+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[0+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[1+s*4];
      xyz[2+o*3] = suborder_xyz[2+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[0+s*4];
      w[o] = suborder_w[s];
      o = o + 1;

      xyz[0+o*3] = suborder_xyz[3+s*4];
      xyz[1+o*3] = suborder_xyz[2+s*4];
      xyz[2+o*3] = suborder_xyz[1+s*4];
      w[o] = suborder_w[s];
      o = o + 1;
    }
    else
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "TETRAHEDRON_NCC_RULE - Fatal error!\n" );
      fprintf ( stderr, "  Illegal SUBORDER(%d) = %d\n", s, suborder[s] );
      exit ( 1 );
    }
  }

  free ( suborder );
  free ( suborder_xyz );
  free ( suborder_w );

  return;
}
/******************************************************************************/

int tetrahedron_ncc_rule_num ( )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_RULE_NUM returns the number of NCC rules available.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Output, int TETRAHEDRON_NCC_RULE_NUM, the number of rules available.
*/
{
  int rule_num;

  rule_num = 7;

  return rule_num;
}
/******************************************************************************/

int *tetrahedron_ncc_suborder ( int rule, int suborder_num )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBORDER returns the suborders for an NCC rule.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int RULE, the index of the rule.

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int TETRAHEDRON_NCC_SUBORDER[SUBORDER_NUM],
    the suborders of the rule.
*/
{
  int *suborder;

  suborder = ( int * ) malloc ( suborder_num * sizeof ( int ) );

  if ( rule == 1 )
  {
    suborder[0] = 1;
  }
  else if ( rule == 2 )
  {
    suborder[0] = 4;
  }
  else if ( rule == 3 )
  {
    suborder[0] = 4;
    suborder[1] = 6;
  }
  else if ( rule == 4 )
  {
    suborder[0] =  4;
    suborder[1] = 12;
    suborder[2] =  4;
  }
  else if ( rule == 5 )
  {
    suborder[0] =  4;
    suborder[1] = 12;
    suborder[2] =  6;
    suborder[3] = 12;
    suborder[4] =  1;
  }
  else if ( rule == 6 )
  {
    suborder[0] =  4;
    suborder[1] = 12;
    suborder[2] = 12;
    suborder[3] = 12;
    suborder[4] = 12;
    suborder[5] =  4;
  }
  else if ( rule == 7 )
  {
    suborder[0] =  4;
    suborder[1] = 12;
    suborder[2] = 12;
    suborder[3] = 12;
    suborder[4] =  6;
    suborder[5] = 24;
    suborder[6] =  4;
    suborder[7] =  4;
    suborder[8] =  6;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TETRAHEDRON_NCC_SUBORDER - Fatal error!\n" );
    fprintf ( stderr, "  Illegal RULE = %d\n", rule );
    exit ( 1 );
  }

  return suborder;
}
/******************************************************************************/

int tetrahedron_ncc_suborder_num ( int rule )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBORDER_NUM returns the number of suborders for an NCC rule.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int RULE, the index of the rule.

    Output, int TETRAHEDRON_NCC_SUBORDER_NUM, the number of suborders
    of the rule.
*/
{
  int suborder_num;

  if ( rule == 1 )
  {
    suborder_num = 1;
  }
  else if ( rule == 2 )
  {
    suborder_num = 1;
  }
  else if ( rule == 3 )
  {
    suborder_num = 2;
  }
  else if ( rule == 4 )
  {
    suborder_num = 3;
  }
  else if ( rule == 5 )
  {
    suborder_num = 5;
  }
  else if ( rule == 6 )
  {
    suborder_num = 6;
  }
  else if ( rule == 7 )
  {
    suborder_num = 9;
  }
  else
  {
    suborder_num = -1;
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TETRAHEDRON_NCC_SUBORDER_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal RULE = %d\n", rule );
    exit ( 1 );
  }

  return suborder_num;
}
/******************************************************************************/

void tetrahedron_ncc_subrule ( int rule, int suborder_num,
  double suborder_xyz[], double suborder_w[] )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE returns a compressed NCC rule.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int RULE, the index of the rule.

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, double SUBORDER_XYZ[4*SUBORDER_NUM],
    the barycentric coordinates of the abscissas.

    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
*/
{
  int i;
  int s;
  int suborder_w_d;
  int *suborder_w_n;
  int suborder_xyz_d;
  int *suborder_xyz_n;

  suborder_xyz_n = ( int * ) malloc ( 4 * suborder_num * sizeof ( int ) );
  suborder_w_n = ( int * ) malloc ( suborder_num * sizeof ( int ) );

  if ( rule == 1 )
  {
    tetrahedron_ncc_subrule_01 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 2 )
  {
    tetrahedron_ncc_subrule_02 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 3 )
  {
    tetrahedron_ncc_subrule_03 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 4 )
  {
    tetrahedron_ncc_subrule_04 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 5 )
  {
    tetrahedron_ncc_subrule_05 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 6 )
  {
    tetrahedron_ncc_subrule_06 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else if ( rule == 7 )
  {
    tetrahedron_ncc_subrule_07 ( suborder_num, suborder_xyz_n, &suborder_xyz_d,
    suborder_w_n, &suborder_w_d );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "TETRAHEDRON_NCC_SUBRULE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal RULE = %d\n", rule );
    exit ( 1 );
  }

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz[i+s*4] =
          ( double ) ( suborder_xyz_n[i+s*4] )
        / ( double ) ( suborder_xyz_d );
    }
  }
  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w[s] = ( double ) suborder_w_n[s] / ( double ) suborder_w_d;
  }

  free ( suborder_w_n );
  free ( suborder_xyz_n );

  return;
}
/******************************************************************************/

void tetrahedron_ncc_subrule_01 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE_01 returns a compressed NCC rule 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.

    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.

    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.

    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
  int i;
  int s;
  int suborder_xyz_n_01[4*1] = {
      1, 1, 1, 1
  };
  int suborder_xyz_d_01 = 4;
  int suborder_w_n_01[1] = { 1 };
  int suborder_w_d_01 = 1;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz_n[i+s*4] = suborder_xyz_n_01[i+s*4];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_01;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_01[s];
  }
  *suborder_w_d = suborder_w_d_01;

  return;
}
/******************************************************************************/

void tetrahedron_ncc_subrule_02 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE_02 returns a compressed NCC rule 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.

    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.

    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.

    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
  int i;
  int s;
  int suborder_xyz_n_02[4*1] = {
      0, 0, 0, 1
  };
  int suborder_xyz_d_02 = 1;
  int suborder_w_n_02[1] = { 1 };
  int suborder_w_d_02 = 4;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz_n[i+s*4] = suborder_xyz_n_02[i+s*4];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_02;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_02[s];
  }
  *suborder_w_d = suborder_w_d_02;

  return;
}
/******************************************************************************/

void tetrahedron_ncc_subrule_03 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE_03 returns a compressed NCC rule 3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.

    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.

    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.

    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
  int i;
  int s;
  int suborder_xyz_n_03[4*2] = {
      0, 0, 0, 2,
      1, 1, 0, 0
  };
  int suborder_xyz_d_03 = 2;
  int suborder_w_n_03[2] = { -1, 4 };
  int suborder_w_d_03 = 20;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz_n[i+s*4] = suborder_xyz_n_03[i+s*4];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_03;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_03[s];
  }
  *suborder_w_d = suborder_w_d_03;

  return;
}
/******************************************************************************/

void tetrahedron_ncc_subrule_04 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE_04 returns a compressed NCC rule 4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.

    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.

    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.

    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
  int i;
  int s;
  int suborder_xyz_n_04[4*3] = {
      0, 0, 0, 3,
      0, 0, 1, 2,
      1, 1, 1, 0
  };
  int suborder_xyz_d_04 = 3;
  int suborder_w_n_04[3] = { 1, 0, 9 };
  int suborder_w_d_04 = 40;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz_n[i+s*4] = suborder_xyz_n_04[i+s*4];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_04;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_04[s];
  }
  *suborder_w_d = suborder_w_d_04;

  return;
}
/******************************************************************************/

void tetrahedron_ncc_subrule_05 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE_05 returns a compressed NCC rule 5.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.

    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.

    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.

    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
  int i;
  int s;
  int suborder_xyz_n_05[4*5] = {
      0, 0, 0, 4,
      0, 0, 3, 1,
      2, 2, 0, 0,
      1, 1, 0, 2,
      1, 1, 1, 1
  };
  int suborder_xyz_d_05 = 4;
  int suborder_w_n_05[5] = { -5, 16, -12, 16, 128 };
  int suborder_w_d_05 = 420;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz_n[i+s*4] = suborder_xyz_n_05[i+s*4];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_05;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_05[s];
  }
  *suborder_w_d = suborder_w_d_05;

  return;
}
/******************************************************************************/

void tetrahedron_ncc_subrule_06 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE_06 returns a compressed NCC rule 6.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.

    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.

    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.

    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
  int i;
  int s;
  int suborder_xyz_n_06[4*6] = {
      0, 0, 0, 5,
      0, 0, 4, 1,
      0, 0, 3, 2,
      1, 1, 0, 3,
      2, 2, 1, 0,
      1, 1, 1, 2
  };
  int suborder_xyz_d_06 = 5;
  int suborder_w_n_06[6] = { 33, -35, 35, 275, -75, 375 };
  int suborder_w_d_06 = 4032;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz_n[i+s*4] = suborder_xyz_n_06[i+s*4];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_06;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_06[s];
  }
  *suborder_w_d = suborder_w_d_06;

  return;
}
/******************************************************************************/

void tetrahedron_ncc_subrule_07 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_NCC_SUBRULE_07 returns a compressed NCC rule 7.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2007

  Author:

    John Burkardt

  Reference:

    Peter Silvester,
    Symmetric Quadrature Formulae for Simplexes,
    Mathematics of Computation,
    Volume 24, Number 109, January 1970, pages 95-100.

  Parameters:

    Input, int SUBORDER_NUM, the number of suborders of the rule.

    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
    the numerators of the barycentric coordinates of the abscissas.

    Output, int *SUBORDER_XYZ_D,
    the denominator of the barycentric coordinates of the abscissas.

    Output, int SUBORDER_W_N[SUBORDER_NUM],
    the numerator of the suborder weights.

    Output, int SUBORDER_W_D,
    the denominator of the suborder weights.
*/
{
  int i;
  int s;
  int suborder_xyz_n_07[4*9] = {
      0, 0, 0, 6,
      0, 0, 5, 1,
      0, 0, 4, 2,
      1, 1, 0, 4,
      3, 3, 0, 0,
      3, 2, 1, 0,
      1, 1, 1, 3,
      2, 2, 2, 0,
      2, 2, 1, 1
  };
  int suborder_xyz_d_07 = 6;
  int suborder_w_n_07[9] = { -7, 24, -30, 0, 40, 30, 180, -45, 0 };
  int suborder_w_d_07 = 1400;

  for ( s = 0; s < suborder_num; s++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      suborder_xyz_n[i+s*4] = suborder_xyz_n_07[i+s*4];
    }
  }
  *suborder_xyz_d = suborder_xyz_d_07;

  for ( s = 0; s < suborder_num; s++ )
  {
    suborder_w_n[s] = suborder_w_n_07[s];
  }
  *suborder_w_d = suborder_w_d_07;

  return;
}
/******************************************************************************/

double tetrahedron_volume ( double tetra[3*4] )

/******************************************************************************/
/*
  Purpose:

    TETRAHEDRON_VOLUME  computes the volume of a tetrahedron.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 August 2005

  Author:

    John Burkardt

  Parameters:

    Input, double TETRA[3*4], the vertices of the tetrahedron.

    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
*/
{
  double a[4*4];
  int i;
  int j;
  double volume;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      a[i+j*4] = tetra[i+j*3];
    }
  }

  i = 3;
  for ( j = 0; j < 4; j++ )
  {
    a[i+j*4] = 1.0;
  }

  volume = fabs ( r8mat_det_4d ( a ) ) / 6.0;

  return volume;
}
/******************************************************************************/

// void timestamp ( )

// /******************************************************************************/
// /*
//   Purpose:

//     TIMESTAMP prints the current YMDHMS date as a time stamp.

//   Example:

//     31 May 2001 09:45:54 AM

//   Licensing:

//     This code is distributed under the GNU LGPL license.

//   Modified:

//     24 September 2003

//   Author:

//     John Burkardt

//   Parameters:

//     None
// */
// {
// # define TIME_SIZE 40

//   static char time_buffer[TIME_SIZE];
//   const struct tm *tm;
//   time_t now;

//   now = time ( NULL );
//   tm = localtime ( &now );

//   strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

//   fprintf ( stdout, "%s\n", time_buffer );

//   return;
// # undef TIME_SIZE
// }

