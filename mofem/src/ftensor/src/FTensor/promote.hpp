/* Traits that allow auto-promotion of int's to double's, double's to
   complex, etc.  Shamelessly stolen from

   http://extreme.indiana.edu/~tveldhui/papers/techniques/

   For now, only int, double, and complex<double> are covered.  If you
   want more, just insert a DECLARE_PROMOTE(A,B,C), where A and B are
   the two types, and C is what they should be coerced to.  */

#pragma once

namespace FTensor
{
  template <class T1, class T2> class promote
  {
  public:
    using V = T1;
  };

#define DECLARE_PROMOTE(A, B, C)                                               \
  template <> class promote<A, B> {                                            \
  public:                                                                      \
    using V = C;                                                               \
  }

  DECLARE_PROMOTE(int, double, double);
  DECLARE_PROMOTE(double, int, double);
  DECLARE_PROMOTE(int, std::complex<double>, std::complex<double>);
  DECLARE_PROMOTE(std::complex<double>, int, std::complex<double>);
  DECLARE_PROMOTE(double, std::complex<double>, std::complex<double>);
  DECLARE_PROMOTE(std::complex<double>, double, std::complex<double>);

  DECLARE_PROMOTE(const double, std::complex<double>, std::complex<double>);
  DECLARE_PROMOTE(const std::complex<double>, double, std::complex<double>);
  DECLARE_PROMOTE(double, const std::complex<double>, std::complex<double>);
  DECLARE_PROMOTE(std::complex<double>, const double, std::complex<double>);
  DECLARE_PROMOTE(const std::complex<double>, const double,
                  std::complex<double>);
  DECLARE_PROMOTE(const double, const std::complex<double>,
                  std::complex<double>);

#ifdef WITH_ADOL_C 
  DECLARE_PROMOTE(adouble, double, adouble);
  DECLARE_PROMOTE(double, adouble, adouble);
  DECLARE_PROMOTE(adouble, int, adouble);
  DECLARE_PROMOTE(int, adouble, adouble);
  DECLARE_PROMOTE(adtl::adouble, double, adtl::adouble);
  DECLARE_PROMOTE(double, adtl::adouble, adtl::adouble);
  DECLARE_PROMOTE(adtl::adouble, int, adtl::adouble);
  DECLARE_PROMOTE(int, adtl::adouble, adtl::adouble);
  DECLARE_PROMOTE(adub, double, adub);
  DECLARE_PROMOTE(double, adub, adub);
  DECLARE_PROMOTE(adub, int, adub);
  DECLARE_PROMOTE(int, adub, adub);
#endif // WITH_ADOL_C

#undef DECLARE_PROMOTE
} // namespace FTensor
