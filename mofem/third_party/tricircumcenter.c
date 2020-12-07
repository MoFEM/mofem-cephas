// Jonathan R Shewchuk <jrs+@cs.cmu.edu>
/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter()   Find the circumcenter of a triangle.                 */
/*                                                                           */
/*  The result is returned both in terms of x-y coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the x-y coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/

#include <stdlib.h>

void tricircumcenter_tp(a, b, c, circumcenter, xi, eta)
double a[2];
double b[2];
double c[2];
double circumcenter[2];
double *xi;
double *eta;
{
  double xba, yba, xca, yca;
  double balength, calength;
  double denominator;
  double xcirca, ycirca;

  /* Use coordinates relative to point `a' of the triangle. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba;
  calength = xca * xca + yca * yca;

  /* Calculate the denominator of the formulae. */
#ifdef EXACT
  /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
  /*   to ensure a correctly signed (and reasonably accurate) result, */
  /*   avoiding any possibility of division by zero.                  */
  denominator = 0.5 / orient2d(b, c, a);
#else
  /* Take your chances with floating-point roundoff. */
  denominator = 0.5 / (xba * yca - yba * xca);
#endif

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = (yca * balength - yba * calength) * denominator;  
  ycirca = (xba * calength - xca * balength) * denominator;  
  circumcenter[0] = xcirca;
  circumcenter[1] = ycirca;

  if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a     */
    /*   coordinate system with a xi-axis directed from `a' to `b' and    */
    /*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
    /*   are computed by Cramer's Rule for solving systems of linear      */
    /*   equations.                                                       */
    *xi = (xcirca * yca - ycirca * xca) * (2.0 * denominator);
    *eta = (ycirca * xba - xcirca * yba) * (2.0 * denominator);
  }
}

/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter3d()   Find the circumcenter of a triangle in 3D.         */
/*                                                                           */
/*  The result is returned both in terms of xyz coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/

void tricircumcenter3d_tp(a, b, c, circumcenter, xi, eta)
double a[3];
double b[3];
double c[3];
double circumcenter[3];
double *xi;
double *eta;
{
  double xba, yba, zba, xca, yca, zca;
  double balength, calength;
  double xcrossbc, ycrossbc, zcrossbc;
  double denominator;
  double xcirca, ycirca, zcirca;

  /* Use coordinates relative to point `a' of the triangle. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  zba = b[2] - a[2];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  zca = c[2] - a[2];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba + zba * zba;
  calength = xca * xca + yca * yca + zca * zca;
  
  /* Cross product of these edges. */
#ifdef EXACT
  /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
  /*   to ensure a correctly signed (and reasonably accurate) result, */
  /*   avoiding any possibility of division by zero.                  */
  xcrossbc = orient2d(b[1], b[2], c[1], c[2], a[1], a[2]);
  ycrossbc = orient2d(b[2], b[0], c[2], c[0], a[2], a[0]);
  zcrossbc = orient2d(b[0], b[1], c[0], c[1], a[0], a[1]);
#else
  /* Take your chances with floating-point roundoff. */
  xcrossbc = yba * zca - yca * zba;
  ycrossbc = zba * xca - zca * xba;
  zcrossbc = xba * yca - xca * yba;
#endif

  /* Calculate the denominator of the formulae. */
  denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
                       zcrossbc * zcrossbc);

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = ((balength * yca - calength * yba) * zcrossbc -
            (balength * zca - calength * zba) * ycrossbc) * denominator;
  ycirca = ((balength * zca - calength * zba) * xcrossbc -
            (balength * xca - calength * xba) * zcrossbc) * denominator;
  zcirca = ((balength * xca - calength * xba) * ycrossbc -
            (balength * yca - calength * yba) * xcrossbc) * denominator;
  circumcenter[0] = xcirca;
  circumcenter[1] = ycirca;
  circumcenter[2] = zcirca;

  if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a     */
    /*   coordinate system with a xi-axis directed from `a' to `b' and    */
    /*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
    /*   are computed by Cramer's Rule for solving systems of linear      */
    /*   equations.                                                       */

    /* There are three ways to do this calculation - using xcrossbc, */
    /*   ycrossbc, or zcrossbc.  Choose whichever has the largest    */
    /*   magnitude, to improve stability and avoid division by zero. */
    if (((xcrossbc >= ycrossbc) ^ (-xcrossbc > ycrossbc)) &&
        ((xcrossbc >= zcrossbc) ^ (-xcrossbc > zcrossbc))) {
      *xi = (ycirca * zca - zcirca * yca) / xcrossbc;
      *eta = (zcirca * yba - ycirca * zba) / xcrossbc;
    } else if ((ycrossbc >= zcrossbc) ^ (-ycrossbc > zcrossbc)) {
      *xi = (zcirca * xca - xcirca * zca) / ycrossbc;
      *eta = (xcirca * zba - zcirca * xba) / ycrossbc;
    } else {
      *xi = (xcirca * yca - ycirca * xca) / zcrossbc;
      *eta = (ycirca * xba - xcirca * yba) / zcrossbc;
    }
  }
}