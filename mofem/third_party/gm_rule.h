#ifndef GM_RULE_H_
#define GM_RULE_H_

void comp_next ( int n, int k, int a[], int *more, int *h, int *t );
void gm_rule_set ( int rule, int dim_num, int point_num, double w[],
  double x[] );
int gm_rule_size ( int rule, int dim_num );
int i4_choose ( int n, int k );
int i4_huge ( void );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
double *monomial_value ( int dim_num, int point_num, double x[], int expon[] );
double r8_abs ( double x );
double r8_factorial ( int n );
double r8vec_dot ( int n, double a1[], double a2[] );
double *r8vec_uniform_01 ( int n, int *seed );
double simplex_unit_monomial_int ( int dim_num, int expon[] );
double simplex_unit_monomial_quadrature ( int dim_num, int expon[],
  int point_num, double x[], double w[] );
double *simplex_unit_sample ( int dim_num, int n, int *seed );
double *simplex_unit_to_general ( int dim_num, int point_num, double t[],
  double ref[] );
double simplex_unit_volume ( int dim_num );
void timestamp ( void );

#endif
