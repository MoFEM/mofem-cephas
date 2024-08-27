// File: sdf.hpp
#include <petscsnes.h>
#include <boost/math/tools/toms748_solve.hpp>
double amplitude = 0.0002;
double indentation = 0.00159624;
const double wave_len = 2.0;
const double w = 2.0 * M_PI / wave_len;

// class Wavy_3D {
// public:
//     // Signed Distance Function (SDF)
//     double sDF(double x, double y, double z, double t) {
//         double sdf = amplitude * (1.0 - std::cos(w * x) * std::cos(w * y)) - z -indentation * t;
//         return sdf;
//     }

//     // Gradient of the SDF
//     std::array<double, 3> gradSDF(double x, double y, double z, double t) {
//         double df_dx = amplitude * w * std::sin(w * x) * std::cos(w * y);
//         double df_dy = amplitude * w * std::cos(w * x) * std::sin(w * y);
//         double df_dz = -1.0;

//         return {df_dx, df_dy, df_dz};
//     }

//     // Hessian of the SDF
//     std::array<double, 6> hessSDF(double x, double y, double z, double t) {
//         double w_squared = w * w;

//         double d2f_dx2 = amplitude * w_squared * std::cos(w * x) * std::cos(w * y);
//         double d2f_dy2 = amplitude * w_squared * std::cos(w * x) * std::cos(w * y);
//         double d2f_dz2 = 0.0;
//         double d2f_dxdy = -amplitude * w_squared * std::sin(w * x) * std::sin(w * y);
//         double d2f_dxdz = 0.0;
//         double d2f_dydz = 0.0;

//         return {d2f_dx2, d2f_dxdy, d2f_dxdz, d2f_dy2, d2f_dydz, d2f_dz2};
//     }
// };
// class Wavy_3D {
// public:
//   Wavy_3D() : x(0.0), y(0.0), z(0.0), t(0.0) {}
//   Wavy_3D(double x, double y, double z, double t) : x(x), y(y), z(z), t(t) {}

//   double sDF(double x, double y, double z, double t) {
//     if (x < 0)
//       x = -x;
//     if (x > wave_len)
//       x = x - std::trunc(x / wave_len) * wave_len;
//     if (x > wave_len / 2)
//       x = wave_len - x;

//     if (y < 0)
//       y = -y;
//     if (y > wave_len)
//       y = y - std::trunc(y / wave_len) * wave_len;
//     if (y > wave_len / 2)
//       y = wave_len - y;

//     double r =
//         amplitude * (1.0 - std::cos(w * x) * std::cos(w * y)) - indentation * t;
//     if (std::abs(r - z) < 1e-14) {
//       return 0.0;
//     }

//     std::vector<std::pair<double, double>> initial_guesses = {
//         {0.0, wave_len / 4.0},
//         {wave_len / 2.0, wave_len / 4.0},
//         {wave_len / 4.0, 0.0},
//         {wave_len / 4.0, wave_len / 2.0}
//         // {0.0, 0.0}
//     };

//     double min_sdf_with_sign = std::numeric_limits<double>::max();

//     for (const auto& guess : initial_guesses) {
//         double p_initial = guess.first;
//         double q_initial = guess.second;
//       Wavy_3D wavy_3d(x, y, z, t);
//       Vec sol = wavy_3d.solve(p_initial, q_initial);
//       PetscScalar *sol_ptr;
//       VecGetArray(sol, &sol_ptr);
//       double p = sol_ptr[0];
//       double q = sol_ptr[1];
//       VecRestoreArray(sol, &sol_ptr);

//       double r = amplitude * (1.0 - std::cos(w * p) * std::cos(w * q)) -
//                  indentation * t;
//       double distance =
//           std::sqrt((p - x) * (p - x) + (q - y) * (q - y) + (r - z) * (r - z));
//       double sign = std::copysign(1.0, r - z);

//       if (distance < std::abs(min_sdf_with_sign)) {
//         min_sdf_with_sign = sign * distance;
//       }


//       VecDestroy(&sol);
//     }

//     return min_sdf_with_sign;
//   }

//   std::array<double, 3> gradSDF(double x, double y, double z, double t) {
//     double delta = 1e-7;
//     double df_dx =
//         (sDF(x + delta, y, z, t) - sDF(x - delta, y, z, t)) / (2. * delta);
//     double df_dy =
//         (sDF(x, y + delta, z, t) - sDF(x, y - delta, z, t)) / (2. * delta);
//     double df_dz =
//         (sDF(x, y, z + delta, t) - sDF(x, y, z - delta, t)) / (2. * delta);
//     double magnitude = std::sqrt(df_dx * df_dx + df_dy * df_dy + df_dz * df_dz);
//     return {df_dx, df_dy, df_dz};
//   }

//   std::array<double, 6> hessSDF(double x, double y, double z, double t) {
//     double delta = 1e-4;
//     double d2f_dx2 = (sDF(x + delta, y, z, t) - 2. * sDF(x, y, z, t) +
//                       sDF(x - delta, y, z, t)) /
//                      (delta * delta);
//     double d2f_dxdy =
//         (sDF(x + delta, y + delta, z, t) - sDF(x + delta, y - delta, z, t) -
//          sDF(x - delta, y + delta, z, t) + sDF(x - delta, y - delta, z, t)) /
//         (4. * delta * delta);
//     double d2f_dxdz =
//         (sDF(x + delta, y, z + delta, t) - sDF(x + delta, y, z - delta, t) -
//          sDF(x - delta, y, z + delta, t) + sDF(x - delta, y, z - delta, t)) /
//         (4. * delta * delta);
//     double d2f_dy2 = (sDF(x, y + delta, z, t) - 2. * sDF(x, y, z, t) +
//                       sDF(x, y - delta, z, t)) /
//                      (delta * delta);
//     double d2f_dydz =
//         (sDF(x, y + delta, z + delta, t) - sDF(x, y + delta, z - delta, t) -
//          sDF(x, y - delta, z + delta, t) + sDF(x, y - delta, z - delta, t)) /
//         (4. * delta * delta);
//     double d2f_dz2 = (sDF(x, y, z + delta, t) - 2. * sDF(x, y, z, t) +
//                       sDF(x, y, z - delta, t)) /
//                      (delta * delta);

//     return {d2f_dx2, d2f_dxdy, d2f_dxdz, d2f_dy2, d2f_dydz, d2f_dz2};
//   }

// private:
//   double x, y, z, t;

//   Vec solve(double p_initial, double q_initial) {
//     Vec sol, r;
//     Mat J;
//     SNES snes;

//     CHKERR VecCreate(PETSC_COMM_SELF, &sol);
//     CHKERR VecSetSizes(sol, PETSC_DECIDE, 2);
//     CHKERR VecSetFromOptions(sol);
//     CHKERR VecDuplicate(sol, &r);

//     CHKERR MatCreate(PETSC_COMM_SELF, &J);
//     CHKERR MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
//     CHKERR MatSetFromOptions(J);
//     CHKERR MatSetUp(J);

//     CHKERR SNESCreate(PETSC_COMM_SELF, &snes);
//     CHKERR SNESSetFunction(snes, r, FormFunction, this);
//     CHKERR SNESSetJacobian(snes, J, J, FormJacobian, this);

//     // CHKERR SNESSetType(snes,SNESNEWTONLS); // deafult line serach: Brack tracing
//     CHKERR SNESSetType(snes, SNESNEWTONLS);//SNESQN for Quasi-Newton method with trust-region strategies
//     double tol = boost::math::tools::epsilon<double>();
//     CHKERR SNESSetTolerances(snes, tol, tol, 0.0, PETSC_DEFAULT,
//                              PETSC_DEFAULT);
//     // CHKERR SNESSetTolerances(snes, 1e-8, 1e-8, 1e-8, 30, 1000);
//     KSP ksp;
//     CHKERR SNESGetKSP(snes, &ksp);
//     PC pc;
//     CHKERR KSPGetPC(ksp, &pc);
//     CHKERR PCSetType(pc, PCLU);
//     CHKERR PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

//     CHKERR SNESSetMaxLinearSolveFailures(snes, -1);

//     PetscScalar *sol_ptr;
//     CHKERR VecGetArray(sol, &sol_ptr);
//     sol_ptr[0] = p_initial;
//     sol_ptr[1] = q_initial;
//     CHKERR VecRestoreArray(sol, &sol_ptr);
//     CHKERR SNESAppendOptionsPrefix(snes, "wave_");
//     CHKERR SNESSetFromOptions(snes);
//     CHKERR SNESSolve(snes, NULL, sol);

//     CHKERR VecDestroy(&r);
//     CHKERR MatDestroy(&J);
//     CHKERR SNESDestroy(&snes);

//     return sol;
//   }

//   // fr residual
//   static PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *ctx) {

//     const PetscScalar *sol;
//     PetscScalar *f;
//     VecGetArrayRead(X, &sol);
//     VecGetArray(F, &f);

//     Wavy_3D *self = static_cast<Wavy_3D *>(ctx);

//     double p = sol[0];
//     double q = sol[1];
//     double r = amplitude * (1.0 - std::cos(w * p) * std::cos(w * q)) -
//                indentation * self->t;
//     double lag = 2.0 * (self->z - r);

//     f[0] = 2.0 * (p - self->x) -
//            lag * amplitude * w * std::sin(w * p) * std::cos(w * q);
//     f[1] = 2.0 * (q - self->y) -
//            lag * amplitude * w * std::cos(w * p) * std::sin(w * q);

//     VecRestoreArrayRead(X, &sol);
//     VecRestoreArray(F, &f);

//     return 0;
//   }

//   // Jacobian frm
//   static PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat P,
//                                      void *ctx) { //  p -preconditioner
//     const PetscScalar *sol;
//     PetscScalar v[4];
//     PetscInt idx[2] = {0, 1};
//     VecGetArrayRead(X, &sol);

//     Wavy_3D *self = static_cast<Wavy_3D *>(ctx);

//     double p = sol[0];
//     double q = sol[1];
//     // double df1_dp = 2. -
//     //                 2. * amplitude * w * w * std::cos(p * w) * std::cos(q
//     * w) *
//     //                     (indentation * self->t + self->z -
//     //                      amplitude * (1. - std::cos(p * w) * std::cos(q *
//     w))) +
//     //                 2. * amplitude * amplitude * w * w * std::cos(q * w) *
//     //                     std::cos(q * w) * std::sin(p * w) * std::sin(p *
//     w);
//     // double df1_dq = 2. * amplitude * amplitude * w * w * std::cos(p * w) *
//     //                     std::cos(q * w) * std::sin(p * w) * std::sin(q *
//     w) +
//     //                 2. * amplitude * w * w *
//     //                     (indentation * self->t + self->z -
//     //                      amplitude * (1. - std::cos(p * w) * std::cos(q *
//     w))) *
//     //                     std::sin(p * w) * std::sin(q * w);
//     // // double df2_dp = 2. * amplitude * amplitude * w * w * std::cos(p *
//     w) *
//     // //                     std::cos(q * w) * std::sin(p * w) * std::sin(q
//     * w) +
//     // //                 2. * amplitude * w * w *
//     // //                     (indentation * self->t + self->z -
//     // //                      amplitude * (1. - std::cos(p * w) * std::cos(q
//     * w))) *
//     // //                     std::sin(p * w) * std::sin(q * w);
//     // double df2_dp = df1_dq;
//     // double df2_dq = 2. -
//     //                 2. * amplitude * w * w * std::cos(p * w) * std::cos(q
//     * w) *
//     //                     (indentation * self->t + self->z -
//     //                      amplitude * (1. - std::cos(p * w) * std::cos(q *
//     w))) +
//     //                 2. * amplitude * amplitude * w * w * std::cos(p * w) *
//     //                     std::cos(p * w) * std::sin(q * w) * std::sin(q *
//     w);
//      double df1_dp =
//      -amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y)
//      + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::cos(p*w)*cos(q*w)
//      + 2.0;
//     double df1_dq =
//     amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y)
//     + 1.0) + 2.0*indentation*self->t
//     + 2.0*self->z)*std::sin(p*w)*std::sin(q*w); double df2_dp = df1_dq;
//     double df2_dq =
//     -amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y)
//     + 1.0) + 2.0*indentation*self->t
//     + 2.0*self->z)*std::cos(p*w)*std::cos(q*w) + 2.0; v[0] = df1_dp; v[1] =
//     df1_dq; v[2] = df2_dp; v[3] = df2_dq; MatSetValues(J, 2, idx, 2, idx, v,
//     INSERT_VALUES);

//     VecRestoreArrayRead(X, &sol);
//     MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
//     return 0;
//   }
// };
class Wavy_3D {
public:
  Vec sol, r;
  Mat J;
  SNES snes;
  KSP ksp;
  PC pc;
  double x, y, z, t;

  Wavy_3D() : x(0.0), y(0.0), z(0.0), t(0.0) {}
  void init() {
    VecCreate(PETSC_COMM_SELF, &sol);
    VecSetSizes(sol, PETSC_DECIDE, 2);
    VecSetFromOptions(sol);
    VecDuplicate(sol, &r);

    MatCreate(PETSC_COMM_SELF, &J);
    MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
    MatSetFromOptions(J);
    MatSetUp(J);

    SNESCreate(PETSC_COMM_SELF, &snes);

    SNESSetFunction(snes, r, FormFunction, this);
    SNESSetJacobian(snes, J, J, FormJacobian, this);

    SNESSetType(snes, SNESNEWTONLS);

    SNESGetKSP(snes, &ksp);

    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);
    PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

    // double tol = boost::math::tools::epsilon<double>();
    double tol = 1e-12;
    SNESSetTolerances(snes, tol, tol, 0, PETSC_DEFAULT, PETSC_DEFAULT);
    SNESSetMaxLinearSolveFailures(snes, -1);
    SNESAppendOptionsPrefix(snes, "wave_");
    SNESSetFromOptions(snes);
  }

  double sDF(double x, double y, double z, double t);
  std::array<double, 3> gradSDF(double x, double y, double z, double t);
  std::array<double, 6> hessSDF(double x, double y, double z, double t);

  static PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *ctx);
  static PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat P, void *ctx);

  void destroy() {
    VecDestroy(&sol);
    VecDestroy(&r);
    MatDestroy(&J);
    SNESDestroy(&snes);
  }
};

// Definitions of static member functions
PetscErrorCode Wavy_3D::FormFunction(SNES snes, Vec X, Vec F, void* ctx) {
    Wavy_3D* self = static_cast<Wavy_3D*>(ctx);
    const PetscScalar *sol;
    PetscScalar *f;
    VecGetArrayRead(X, &sol);
    VecGetArray(F, &f);

    double p = sol[0];
    double q = sol[1];
    double r = amplitude* (1.0 - std::cos(w * p) * std::cos(w * q)) - indentation * self->t;
    double lag = 2.0 * (self->z - r);

    f[0] = 2.0 * (p - self->x) - lag * amplitude* w * std::sin(w * p) * std::cos(w * q);
    f[1] = 2.0 * (q - self->y) - lag * amplitude* w * std::cos(w * p) * std::sin(w * q);

    VecRestoreArrayRead(X, &sol);
    VecRestoreArray(F, &f);

    return 0;
}

PetscErrorCode Wavy_3D::FormJacobian(SNES snes, Vec X, Mat J, Mat P, void* ctx) {
    Wavy_3D* self = static_cast<Wavy_3D*>(ctx);
    const PetscScalar *sol;
    PetscScalar v[4];
    PetscInt idx[2] = {0, 1};
    VecGetArrayRead(X, &sol);

    double p = sol[0];
    double q = sol[1];
    double df1_dp = 2.0 - 2.0 * amplitude* w * w * std::cos(p * w) * std::cos(q * w) *
                    (indentation * self->t + self->z - amplitude* (1.0 - std::cos(p * w) * std::cos(q * w))) +
                    2.0 * amplitude* amplitude* w * w * std::cos(q * w) * std::cos(q * w) * std::sin(p * w) * std::sin(p * w);
    double df1_dq = 2.0 * amplitude* amplitude* w * w * std::cos(p * w) * std::cos(q * w) *
                    std::sin(p * w) * std::sin(q * w) + 2.0 * amplitude* w * w *
                    (indentation * self->t + self->z - amplitude* (1.0 - std::cos(p * w) * std::cos(q * w))) *
                    std::sin(p * w) * std::sin(q * w);
    double df2_dp = df1_dq;  // Symmetry
    double df2_dq = 2.0 - 2.0 * amplitude* w * w * std::cos(p * w) * std::cos(q * w) *
                    (indentation * self->t + self->z - amplitude* (1.0 - std::cos(p * w) * std::cos(q * w))) +
                    2.0 * amplitude* amplitude* w * w * std::cos(p * w) * std::cos(p * w) *
                    std::sin(q * w) * std::sin(q * w);
    //  double df1_dp = -amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::cos(p*w)*cos(q*w) + 2.0;
    // double df1_dq = amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::sin(p*w)*std::sin(q*w);
    // double df2_dp = amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::sin(p*w)*std::sin(q*w);
    // double df2_dq = -amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::cos(p*w)*std::cos(q*w) + 2.0;

    v[0] = df1_dp;
    v[1] = df1_dq;
    v[2] = df2_dp;
    v[3] = df2_dq;
    MatSetValues(J, 2, idx, 2, idx, v, INSERT_VALUES);

    VecRestoreArrayRead(X, &sol);
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

    return 0;
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

    std::vector<std::pair<double, double>> initial_guesses = {
        {0.0, wave_len / 4.0},
        {wave_len / 2.0, wave_len / 4.0},
        {wave_len / 4.0, 0.0},
        {wave_len / 4.0, wave_len / 2.0}
                // {0.0, 0.0}
    };

    double min_sdf_with_sign = std::numeric_limits<double>::max();

    for (const auto& guess : initial_guesses) {
        double p_initial = guess.first;
        double q_initial = guess.second;

        // VecCreate(PETSC_COMM_SELF, &sol);
        // VecSetSizes(sol, PETSC_DECIDE, 2);
        // VecSetFromOptions(sol);
        // VecDuplicate(sol, &r);

        // MatCreate(PETSC_COMM_SELF, &J);
        // MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
        // MatSetFromOptions(J);
        // MatSetUp(J);

        // SNESCreate(PETSC_COMM_SELF, &snes);

        // SNESSetFunction(snes, r, FormFunction, this);
        // SNESSetJacobian(snes, J, J, FormJacobian, this);

        // SNESSetType(snes, SNESNEWTONLS);

        // SNESGetKSP(snes, &ksp);

        // KSPGetPC(ksp, &pc);
        // PCSetType(pc, PCLU);
        // PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

        // double tol = boost::math::tools::epsilon<double>();
        // CHKERRQ(SNESSetTolerances(snes, tol, tol, 0, PETSC_DEFAULT, PETSC_DEFAULT));
        // SNESSetMaxLinearSolveFailures(snes, -1);

        PetscScalar *sol_ptr;
        VecGetArray(sol, &sol_ptr);
        sol_ptr[0] = p_initial;
        sol_ptr[1] = q_initial;
        VecRestoreArray(sol, &sol_ptr);
 
        SNESSolve(snes, NULL, sol);

        VecGetArray(sol, &sol_ptr);
        double p = sol_ptr[0];
        double q = sol_ptr[1];
        VecRestoreArray(sol, &sol_ptr);

        double r_result = amplitude * (1.0 - std::cos(w * p) * std::cos(w * q)) - indentation * t;
        double distance = std::sqrt((p - x) * (p - x) + (q - y) * (q - y) + (r_result - z) * (r_result - z));
        double sign = std::copysign(1.0, r_result - z);

        if (distance < std::abs(min_sdf_with_sign)) {
            min_sdf_with_sign = sign * distance;
        }

        // VecDestroy(&sol);
        // VecDestroy(&r);
        // MatDestroy(&J);
        // SNESDestroy(&snes);
    }

    return min_sdf_with_sign;
}

std::array<double, 3> Wavy_3D::gradSDF(double x, double y, double z, double t) {
    const double delta = 1e-4;
    double df_dx = (sDF(x + delta, y, z, t) - sDF(x - delta, y, z, t)) / (2. * delta);
    double df_dy = (sDF(x, y + delta, z, t) - sDF(x, y - delta, z, t)) / (2. * delta);
    double df_dz = (sDF(x, y, z + delta, t) - sDF(x, y, z - delta, t)) / (2. * delta);
    return {df_dx, df_dy, df_dz};
}

std::array<double, 6> Wavy_3D::hessSDF(double x, double y, double z, double t) {
    const double delta = 1e-3;
    double d2f_dx2 = (sDF(x + delta, y, z, t) - 2. * sDF(x, y, z, t) + sDF(x - delta, y, z, t)) / (delta * delta);
    double d2f_dxdy = (sDF(x + delta, y + delta, z, t) - sDF(x + delta, y - delta, z, t) - sDF(x - delta, y + delta, z, t) + sDF(x - delta, y - delta, z, t)) / (4. * delta * delta);
    double d2f_dxdz = (sDF(x + delta, y, z + delta, t) - sDF(x + delta, y, z - delta, t) - sDF(x - delta, y, z + delta, t) + sDF(x - delta, y, z - delta, t)) / (4. * delta * delta);
    double d2f_dy2 = (sDF(x, y + delta, z, t) - 2. * sDF(x, y, z, t) + sDF(x, y - delta, z, t)) / (delta * delta);
    double d2f_dydz = (sDF(x, y + delta, z + delta, t) - sDF(x, y + delta, z - delta, t) - sDF(x, y - delta, z + delta, t) + sDF(x, y - delta, z - delta, t)) / (4. * delta * delta);
    double d2f_dz2 = (sDF(x, y, z + delta, t) - 2. * sDF(x, y, z, t) + sDF(x, y, z - delta, t)) / (delta * delta);
    return {d2f_dx2, d2f_dxdy, d2f_dxdz, d2f_dy2, d2f_dydz, d2f_dz2};
}
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