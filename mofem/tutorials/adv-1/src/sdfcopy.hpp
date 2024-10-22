// File: sdf.hpp
#include <petscsnes.h>
#include <boost/math/tools/toms748_solve.hpp>
// double amplitude = 0.0002;
// double indentation = 0.0021745;
double amplitude = 0.0002;
double indentation = 0.002;
const double wave_len = 2.0;
const double w = 2.0 * M_PI / wave_len;

class Wavy_3D {
public:

  double x, y, z, t;

  SmartPetscObj<Mat> A; //< Nested SNES tangent matrix 
  SmartPetscObj<Vec> R; //< Nested SNES residual
  SmartPetscObj<Vec> Chi; //< Nested SNES unknown vector

  ublas::matrix<double, ublas::column_major> dataA;
  SmartPetscObj<SNES> sNes; //< Nested SNES solver

  Wavy_3D() : x(0.0), y(0.0), z(0.0), t(0.0) {}
  void init() {

    PetscInt n = 2;
    dataA.resize(n, n, false);
    Mat A_tmp;
    Vec R_tmp, Chi_tmp;
    CHKERR MatCreateSeqDense(PETSC_COMM_SELF, n, n, &dataA(0, 0), &A_tmp);
    CHKERR VecCreateSeq(PETSC_COMM_SELF, n, &R_tmp);
    CHKERR VecCreateSeq(PETSC_COMM_SELF, n, &Chi_tmp);
    A = SmartPetscObj<Mat>(A_tmp);
    R = SmartPetscObj<Vec>(R_tmp);
    Chi = SmartPetscObj<Vec>(Chi_tmp);
    sNes = createSNES(PETSC_COMM_SELF);

    CHKERR SNESSetFunction(sNes, R, FormFunction, (void *)this);
    CHKERR SNESSetJacobian(sNes, A, A, FormJacobian, (void *)this);

    SNESLineSearch line_search;
    KSP ksp;
    PC pc;

    CHKERR SNESGetKSP(sNes, &ksp);
    CHKERR KSPGetPC(ksp, &pc);
    CHKERR KSPSetType(ksp, KSPPREONLY);
    CHKERR PCSetType(pc, PCLU);


    // CHKERR SNESSetFromOptions(sNes);
    // CHKERR SNESSetType(sNes, SNESNEWTONTR);

    double tol = 1e-10; //1e-10, 1e-12
    CHKERR SNESSetTolerances(sNes, tol, tol, 0, 50, 5000);
    CHKERR SNESGetLineSearch(sNes, &line_search);
    CHKERR SNESLineSearchSetType(line_search, SNESLINESEARCHBT);

    // CHKERR SNESAppendOptionsPrefix(sNes, "wavy_");


  }

  double sDF(double x, double y, double z, double t);
  std::array<double, 3> gradSDF(double x, double y, double z, double t);
  std::array<double, 6> hessSDF(double x, double y, double z, double t);

  static PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *ctx);
  static PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat P, void *ctx);

void destroy() {
    _p_Vec* chi_ptr = Chi.get();
    _p_Vec* r_ptr = R.get();
    _p_Mat* a_ptr = A.get();
    _p_SNES* snes_ptr = sNes.get();
    CHKERR VecDestroy(&chi_ptr);
    CHKERR VecDestroy(&r_ptr);    
    CHKERR MatDestroy(&a_ptr);    
    CHKERR SNESDestroy(&snes_ptr); 
}
};

// Definitions of static member functions
PetscErrorCode Wavy_3D::FormFunction(SNES snes, Vec chi, Vec F, void *ctx) {
    MoFEMFunctionBegin;
    // Wavy_3D *self;
    // self = (Wavy_3D *)ctx;
    Wavy_3D *self = static_cast<Wavy_3D *>(ctx);
    const PetscScalar *array;
    PetscScalar *f;

    CHKERR VecGetArrayRead(chi, &array);
    CHKERR VecGetArray(F, &f);

    double p = array[0];
    double q = array[1];
    double r = amplitude* (1.0 - PetscCosScalar(w * p) * PetscCosScalar(w * q)) - indentation * self->t;
    double lag = 2.0 * (self->z - r);

    f[0] = 2.0 * (p - self->x) - lag * amplitude* w * PetscSinScalar(w * p) * PetscCosScalar(w * q);
    f[1] = 2.0 * (q - self->y) - lag * amplitude* w * PetscCosScalar(w * p) * PetscSinScalar(w * q);

    CHKERR VecRestoreArrayRead(chi, &array);
    CHKERR VecRestoreArray(F, &f);

  MoFEMFunctionReturn(0);
}


PetscErrorCode Wavy_3D::FormJacobian(SNES snes, Vec chi, Mat A, Mat P, void *ctx) {
    MoFEMFunctionBegin;
    Wavy_3D *self = static_cast<Wavy_3D *>(ctx);
    const PetscScalar *array;
    CHKERR VecGetArrayRead(chi, &array);

    double p = array[0];
    double q = array[1];

    CHKERR VecRestoreArrayRead(chi, &array);

    double df1_dp = 2.0 - 2.0 * amplitude * w * w * std::cos(p * w) * std::cos(q * w) *
                    (indentation * self->t + self->z - amplitude * (1.0 - std::cos(p * w) * std::cos(q * w))) +
                    2.0 * amplitude * amplitude * w * w * std::cos(q * w) * std::cos(q * w) * std::sin(p * w) * std::sin(p * w);
    double df1_dq = 2.0 * amplitude * amplitude * w * w * std::cos(p * w) * std::cos(q * w) *
                    std::sin(p * w) * std::sin(q * w) + 2.0 * amplitude * w * w *
                    (indentation * self->t + self->z - amplitude * (1.0 - std::cos(p * w) * std::cos(q * w))) *
                    std::sin(p * w) * std::sin(q * w);
    double df2_dp = df1_dq;
    double df2_dq = 2.0 - 2.0 * amplitude * w * w * std::cos(p * w) * std::cos(q * w) *
                    (indentation * self->t + self->z - amplitude * (1.0 - std::cos(p * w) * std::cos(q * w))) +
                    2.0 * amplitude * amplitude * w * w * std::cos(p * w) * std::cos(p * w) *
                    std::sin(q * w) * std::sin(q * w);

    PetscInt row[2] = {0, 1};
    PetscInt col[2] = {0, 1};
    PetscScalar values[4] = {df1_dp, df1_dq, df2_dp, df2_dq};

    CHKERR MatSetValues(self->A, 2, row, 2, col, values, INSERT_VALUES);
    CHKERR MatAssemblyBegin(self->A, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(self->A, MAT_FINAL_ASSEMBLY);

    MoFEMFunctionReturn(0);
}

double Wavy_3D::sDF(double x, double y, double z, double t) {
    x = fmod(fabs(x), wave_len);
    if (x > wave_len / 2) x = wave_len - x;

    y = fmod(fabs(y), wave_len);
    if (y > wave_len / 2) y = wave_len - y;

    this->x = x;
    this->y = y;
    this->z = z;
    this->t = t;
    double r = amplitude * (1.0 - std::cos(w * x) * std::cos(w * y)) - indentation * t;
    if(std::abs(r - z) < 1e-12) {
        return 0.0;
    }
    // else{
    std::vector<std::pair<double, double>> initial_guesses = {
        {0.0, wave_len / 4.0},
        {wave_len / 2.0, wave_len / 4.0},
        {wave_len / 4.0, 0.0},
        {wave_len / 4.0, wave_len / 2.0}
    };

    double min_sdf_with_sign = std::numeric_limits<double>::max();
    PetscScalar *sol_ptr;
    for (const auto& guess : initial_guesses) {
        double p_initial = guess.first;
        double q_initial = guess.second;


        CHKERR VecGetArray(Chi, &sol_ptr);
        sol_ptr[0] = p_initial;
        sol_ptr[1] = q_initial;
        CHKERR VecRestoreArray(Chi, &sol_ptr);
 
        CHKERR SNESSolve(sNes, NULL, Chi);

       CHKERR  VecGetArray(Chi, &sol_ptr);
        double p = sol_ptr[0];
        double q = sol_ptr[1];
        CHKERR VecRestoreArray(Chi, &sol_ptr);
        double r_result = amplitude * (1.0 - std::cos(w * p) * std::cos(w * q)) - indentation * t;
        double distance = std::sqrt((p - x) * (p - x) + (q - y) * (q - y) + (r_result - z) * (r_result - z));
        double sign = std::copysign(1.0, r_result - z);

        if (distance < std::abs(min_sdf_with_sign)) {
            min_sdf_with_sign = sign * distance;
        }
    }
    return min_sdf_with_sign;
}

std::array<double, 3> Wavy_3D::gradSDF(double x, double y, double z, double t) {
    const double delta = 1e-6; //1e-10
    double df_dx = (sDF(x + delta, y, z, t) - sDF(x - delta, y, z, t)) / (2. * delta);
    double df_dy = (sDF(x, y + delta, z, t) - sDF(x, y - delta, z, t)) / (2. * delta);
    double df_dz = (sDF(x, y, z + delta, t) - sDF(x, y, z - delta, t)) / (2. * delta);
    return {df_dx, df_dy, df_dz};
}

std::array<double, 6> Wavy_3D::hessSDF(double x, double y, double z, double t) {
    const double delta = 1e-4;  
    const double inv_delta2 = 1.0 / (delta * delta);  
    // Second pd
    double d2f_dx2 = (sDF(x + delta, y, z, t) - 2.0 * sDF(x, y, z, t) + sDF(x - delta, y, z, t)) * inv_delta2;
    double d2f_dy2 = (sDF(x, y + delta, z, t) - 2.0 * sDF(x, y, z, t) + sDF(x, y - delta, z, t)) * inv_delta2;
    double d2f_dz2 = (sDF(x, y, z + delta, t) - 2.0 * sDF(x, y, z, t) + sDF(x, y, z - delta, t)) * inv_delta2;

    // Mixed pd
    double d2f_dxdy = (sDF(x + delta, y + delta, z, t) - sDF(x + delta, y - delta, z, t) 
                    - sDF(x - delta, y + delta, z, t) + sDF(x - delta, y - delta, z, t)) / (4.0 * delta * delta);

    double d2f_dxdz = (sDF(x + delta, y, z + delta, t) - sDF(x + delta, y, z - delta, t) 
                    - sDF(x - delta, y, z + delta, t) + sDF(x - delta, y, z - delta, t)) / (4.0 * delta * delta);

    double d2f_dydz = (sDF(x, y + delta, z + delta, t) - sDF(x, y + delta, z - delta, t) 
                    - sDF(x, y - delta, z + delta, t) + sDF(x, y - delta, z - delta, t)) / (4.0 * delta * delta);

    return {d2f_dx2, d2f_dxdy, d2f_dxdz, d2f_dy2, d2f_dydz, d2f_dz2};
}

// std::array<double, 6> Wavy_3D::hessSDF(double x, double y, double z, double t) {
//     const double delta = 1e-4;
//     double d2f_dx2 = (sDF(x + delta, y, z, t) - 2. * sDF(x, y, z, t) + sDF(x - delta, y, z, t)) / (delta * delta);
//     double d2f_dxdy = (sDF(x + delta, y + delta, z, t) - sDF(x + delta, y - delta, z, t) - sDF(x - delta, y + delta, z, t) + sDF(x - delta, y - delta, z, t)) / (4. * delta * delta);
//     double d2f_dxdz = (sDF(x + delta, y, z + delta, t) - sDF(x + delta, y, z - delta, t) - sDF(x - delta, y, z + delta, t) + sDF(x - delta, y, z - delta, t)) / (4. * delta * delta);
//     double d2f_dy2 = (sDF(x, y + delta, z, t) - 2. * sDF(x, y, z, t) + sDF(x, y - delta, z, t)) / (delta * delta);
//     double d2f_dydz = (sDF(x, y + delta, z + delta, t) - sDF(x, y + delta, z - delta, t) - sDF(x, y - delta, z + delta, t) + sDF(x, y - delta, z - delta, t)) / (4. * delta * delta);
//     double d2f_dz2 = (sDF(x, y, z + delta, t) - 2. * sDF(x, y, z, t) + sDF(x, y, z - delta, t)) / (delta * delta);
//     return {d2f_dx2, d2f_dxdy, d2f_dxdz, d2f_dy2, d2f_dydz, d2f_dz2};
// }


class Wavy_2D {

public:
  double F(double p, double x, double y, double t) {
    double q = amplitude * (1.0 - std::cos(w * p)) - indentation * t;
    double lag = 2.0 * (y - q);
    return 2.0 * (p - x) - lag * amplitude * w * std::sin(w * p);
  }

  double sDF(double x, double y, double z, double t) {
    x = std::fmod(std::abs(x), wave_len);
    if (x > wave_len / 2) {
      x = wave_len - x;
    }

    auto F_lambda = [&](double p) { return F(p, x, y, t); };

    try {
    double ytest = amplitude * (1.0 - std::cos(w * x)) - indentation * t;
    if (std::abs(ytest - y) < 1e-14) {
        return 0.0;
      }
      boost::math::tools::eps_tolerance<double> tol(
          std::numeric_limits<double>::digits - 3.);
      double eps = std::numeric_limits<double>::epsilon();
      double a = eps;
      double b = wave_len / 2.0 - eps;
      std::uintmax_t max_iter = 100;
      auto result =
          boost::math::tools::toms748_solve(F_lambda, a, b, tol, max_iter);
      double p0 = (result.first + result.second) / 2.0;

      double q0 = amplitude * (1.0 - std::cos(w * p0)) - indentation * t;


      double distance = std::hypot(p0 - x, q0 - y);
      int sign = std::copysign(1.0, q0 - y);
      return sign * distance;
    } catch (const std::exception &e) {
      double wave_y_d2 = -indentation * t;
      int sgn_d2 = std::copysign(1.0, wave_y_d2 - y);
      double d2 = std::sqrt(x * x + (wave_y_d2 - y) * (wave_y_d2 - y));

      double wave_y_d3 =
          amplitude * (1.0 - std::cos(w * (wave_len / 2.0))) - indentation * t;
      int sgn_d3 = std::copysign(1.0, wave_y_d3 - y);
      double d3 = std::sqrt((wave_len / 2.0 - x) * (wave_len / 2.0 - x) +
                            (wave_y_d3 - y) * (wave_y_d3 - y));

      if (d2 < d3) {
        double sdf = sgn_d2 * d2;
        return sdf;
      } else {
        double sdf = sgn_d3 * d3;

        return sdf;
      }
    }
  }
  // gradsdf
  std::array<double, 3> gradSDF(double x, double y, double z, double t) {
    double delta = 1e-5;
    double df_dx =
        ((sDF(x + delta, y, z, t) - sDF(x - delta, y, z, t)) / (2. * delta));
    double df_dy =
        ((sDF(x, y + delta, z, t) - sDF(x, y - delta, z, t)) / (2. * delta));
    double df_dz = 0.0;

    return {df_dx, df_dy, df_dz};
  }

  // hesssdf
  std::array<double, 6> hessSDF(double x, double y, double z, double t) {
    double delta = 1e-2;
    double d2f_dx2 = ((sDF(x + delta, y, z, t) - 2. * sDF(x, y, z, t) +
                       sDF(x - delta, y, z, t)) /
                      (delta * delta));
    double d2f_dy2 = ((sDF(x, y + delta, z, t) - 2. * sDF(x, y, z, t) +
                       sDF(x, y - delta, z, t)) /
                      (delta * delta));
    double d2f_dxdy =
        ((sDF(x + delta, y + delta, z, t) - sDF(x + delta, y - delta, z, t) -
          sDF(x - delta, y + delta, z, t) + sDF(x - delta, y - delta, z, t)) /
         (4. * delta * delta));
    double d2f_dxdz = 0.0;
    double df2_dydz = 0.0;
    double d2f_dz2 = 0.0;
    return {d2f_dx2, d2f_dxdy, d2f_dxdz, d2f_dy2, df2_dydz, d2f_dz2};
  }
};