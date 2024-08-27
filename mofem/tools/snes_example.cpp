#include <petscsnes.h>
#include <MoFEM.hpp>

double amplitude = 0.0002;
double wave_len = 2.0;
double w = 2.0 * M_PI / wave_len;
double indentation = 0.0021745;

class Wavy_3D {
public:
    Vec sol, r;
    Mat J;
    SNES snes;
    KSP ksp;
    PC pc;
    double x, y, z, t;
    

    Wavy_3D() : x(0.0), y(0.0), z(0.0), t(0.0) {}

    double sDF(double x, double y, double z, double t);
    std::array<double, 3> gradSDF(double x, double y, double z, double t);
    std::array<double, 6> hessSDF(double x, double y, double z, double t);

    static PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void* ctx);
    static PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat P, void* ctx);
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
    // double df1_dp = 2.0 - 2.0 * amplitude* w * w * std::cos(p * w) * std::cos(q * w) *
    //                 (indentation * self->t + self->z - amplitude* (1.0 - std::cos(p * w) * std::cos(q * w))) +
    //                 2.0 * amplitude* amplitude* w * w * std::cos(q * w) * std::cos(q * w) * std::sin(p * w) * std::sin(p * w);
    // double df1_dq = 2.0 * amplitude* amplitude* w * w * std::cos(p * w) * std::cos(q * w) *
    //                 std::sin(p * w) * std::sin(q * w) + 2.0 * amplitude* w * w *
    //                 (indentation * self->t + self->z - amplitude* (1.0 - std::cos(p * w) * std::cos(q * w))) *
    //                 std::sin(p * w) * std::sin(q * w);
    // double df2_dp = df1_dq;  // Symmetry
    // double df2_dq = 2.0 - 2.0 * amplitude* w * w * std::cos(p * w) * std::cos(q * w) *
    //                 (indentation * self->t + self->z - amplitude* (1.0 - std::cos(p * w) * std::cos(q * w))) +
    //                 2.0 * amplitude* amplitude* w * w * std::cos(p * w) * std::cos(p * w) *
    //                 std::sin(q * w) * std::sin(q * w);
     double df1_dp = -amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::cos(p*w)*cos(q*w) + 2.0;
    double df1_dq = amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::sin(p*w)*std::sin(q*w);
    double df2_dp = amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::sin(p*w)*std::sin(q*w);
    double df2_dq = -amplitude*w*w*(-2.0*amplitude*(-std::cos(w*self->x)*std::cos(w*self->y) + 1.0) + 2.0*indentation*self->t + 2.0*self->z)*std::cos(p*w)*std::cos(q*w) + 2.0;

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

        double tol = boost::math::tools::epsilon<double>();
        CHKERRQ(SNESSetTolerances(snes, tol, tol, tol, PETSC_DEFAULT, PETSC_DEFAULT));
        SNESSetMaxLinearSolveFailures(snes, -1);

        PetscScalar *sol_ptr;
        VecGetArray(sol, &sol_ptr);
        sol_ptr[0] = p_initial;
        sol_ptr[1] = q_initial;
        VecRestoreArray(sol, &sol_ptr);
        CHKERR SNESAppendOptionsPrefix(snes, "wave_");
        SNESSetFromOptions(snes);
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

        VecDestroy(&sol);
        VecDestroy(&r);
        MatDestroy(&J);
        SNESDestroy(&snes);
    }

    return min_sdf_with_sign;
}

std::array<double, 3> Wavy_3D::gradSDF(double x, double y, double z, double t) {
    const double delta = 1e-8;
    double df_dx = (sDF(x + delta, y, z, t) - sDF(x - delta, y, z, t)) / (2. * delta);
    double df_dy = (sDF(x, y + delta, z, t) - sDF(x, y - delta, z, t)) / (2. * delta);
    double df_dz = (sDF(x, y, z + delta, t) - sDF(x, y, z - delta, t)) / (2. * delta);
    return {df_dx, df_dy, df_dz};
}

std::array<double, 6> Wavy_3D::hessSDF(double x, double y, double z, double t) {
    const double delta = 1e-2;
    double d2f_dx2 = (sDF(x + delta, y, z, t) - 2. * sDF(x, y, z, t) + sDF(x - delta, y, z, t)) / (delta * delta);
    double d2f_dxdy = (sDF(x + delta, y + delta, z, t) - sDF(x + delta, y - delta, z, t) - sDF(x - delta, y + delta, z, t) + sDF(x - delta, y - delta, z, t)) / (4. * delta * delta);
    double d2f_dxdz = (sDF(x + delta, y, z + delta, t) - sDF(x + delta, y, z - delta, t) - sDF(x - delta, y, z + delta, t) + sDF(x - delta, y, z - delta, t)) / (4. * delta * delta);
    double d2f_dy2 = (sDF(x, y + delta, z, t) - 2. * sDF(x, y, z, t) + sDF(x, y - delta, z, t)) / (delta * delta);
    double d2f_dydz = (sDF(x, y + delta, z + delta, t) - sDF(x, y + delta, z - delta, t) - sDF(x, y - delta, z + delta, t) + sDF(x, y - delta, z - delta, t)) / (4. * delta * delta);
    double d2f_dz2 = (sDF(x, y, z + delta, t) - 2. * sDF(x, y, z, t) + sDF(x, y, z - delta, t)) / (delta * delta);
    return {d2f_dx2, d2f_dxdy, d2f_dxdz, d2f_dy2, d2f_dydz, d2f_dz2};
}

// struct SDFCpp {
//     Wavy_3D wavy3d;

//     SDFCpp() = default;
//     virtual ~SDFCpp() = default;

//     void sdfInit(double wave_amplitude, double wave_len, double indentation) {
//         MoFEMFunctionBegin;
//         wavy3d.amplitude = wave_amplitude;
//         wavy3d.wave_len = wave_len;
//         wavy3d.indentation = indentation;
//         wavy3d.w = 2.0 * M_PI / wave_len;

//     }

//     void evalSdf(double x, double y, double z, double t, double& sdf) {
//         MoFEMFunctionBegin;
//         sdf = wavy3d.sDF(x, y, z, t);

//     }

//     void evalGradSdf(double x, double y, double z, double t, std::array<double, 3>& grad_sdf) {
//         MoFEMFunctionBegin;
//         auto grad = wavy3d.gradSDF(x, y, z, t);
//         grad_sdf[0] = grad[0];
//         grad_sdf[1] = grad[1];
//         grad_sdf[2] = grad[2];

//     }

//     void evalHessSdf(double x, double y, double z, double t, std::array<double, 6>& hess_sdf) {
//         MoFEMFunctionBegin;
//         auto hess = wavy3d.hessSDF(x, y, z, t);
//         hess_sdf[0] = hess[0];
//         hess_sdf[1] = hess[1];
//         hess_sdf[2] = hess[2];
//         hess_sdf[3] = hess[3];
//         hess_sdf[4] = hess[4];
//         hess_sdf[5] = hess[5];

//     }
// };
int main(int argc, char **argv) {
  const char help[] = "Solve 3D wavy surface problem.";
  const char param_file[] = "param_file.petsc";
  PetscInitialize(&argc, &argv, NULL, help);

  // Grid size and ranges
  const int grid_size = 20;
  const double x_min = -0.0, x_max = 1.0;
  const double y_min = -0.0, y_max = 1.0;
  const double z_val = -0.5; // Fixed z valuez
  const double t = 1.0;

  // Open output file
  std::ofstream outfile("sdf_gggf.csv");
  if (!outfile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }
  outfile << "x,y,sdf,grad_sdf_x,grad_sdf_y,grad_sdf_z, "
             "hess_sdf_xx,hess_sdf_xy,hess_sdf_xz,hess_sdf_yy,hess_sdf_yz,hess_"
             "sdf_zz\n";

  Wavy_3D wavy_3d;
  double pp = wavy_3d.sDF(1.0, 1.0, 0.4, 1);
  std::cout << "pp: " << pp << std::endl;

  for (int i = 0; i < grid_size; ++i) {
    double x = x_min + i * (x_max - x_min) / (grid_size - 1);
    for (int j = 0; j < grid_size; ++j) {
      double y = y_min + j * (y_max - y_min) / (grid_size - 1);
      double sdf = wavy_3d.sDF(x, y, z_val, t);

      // Calculate the gradient
      auto grad_sdf = wavy_3d.gradSDF(x, y, z_val, t);
      double grad_sdf_x = grad_sdf[0];
      double grad_sdf_y = grad_sdf[1];
      double grad_sdf_z = grad_sdf[2];
      // Calculate the Hessian§
      auto hess_sdf = wavy_3d.hessSDF(x, y, z_val, t);
      double hess_sdf_xx = hess_sdf[0];
      double hess_sdf_xy = hess_sdf[1];
      double hess_sdf_xz = hess_sdf[2];
      double hess_sdf_yy = hess_sdf[3];
      double hess_sdf_yz = hess_sdf[4];
      double hess_sdf_zz = hess_sdf[5];

      // Write the SDF and gradient to the CSV file
      outfile << x << "," << y << "," << sdf << "," << grad_sdf_x << ","
              << grad_sdf_y << "," << grad_sdf_z << "," << hess_sdf_xx << ","
              << hess_sdf_xy << "," << hess_sdf_xz << "," << hess_sdf_yy << ","
              << hess_sdf_yz << "," << hess_sdf_zz << "\n";
      // outfile << x << "," << y << "," << sdf << "\n";
    }
  }

  outfile.close();

  PetscFinalize();
  return 0;
}
