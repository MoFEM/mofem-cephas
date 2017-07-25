/** \file MaterialUnsaturatedFlow.hpp
 * \brief Mix implementation of transport element
 *
 * \ingroup mofem_mix_transport_elem
 *
 */

#ifndef __MATERIALUNSATURATEDFLOW_HPP__
#define __MATERIALUNSATURATEDFLOW_HPP__

namespace MixTransport {

  struct CommonMaterialData: public GenericMaterial {

    CommonMaterialData() {
      kappa = 5.6634e-8;
      mu = 1e-9;
      a = 0.0123;
      thetaS = 0.38;
      thetaR = 0.068;
      m = 0.0826;
      APc = 0.;
      APcZ = 0.;
      APcZZ = 0;
    }

    std::string matName; ///< material name

    double thetaS;      ///< saturated water content
    double thetaR;      ///< residual water contents
    double m;           ///< model parameter

    double mu;          ///< dynamic viscosity [t/mm*s]
    double kappa;       ///< intrinsic conductivity [mm^2]
    double a;           ///< model parameter [MPa]

    // Initial capillary pressure coefficient
    // Pc = APx+APx*Z+APx*Z*Z

    double APc;       ///< Initial capillary pressure coefficient
    double APcZ;      ///< Initial capillary pressure coefficient
    double APcZZ;     ///< Initial capillary pressure coefficient

    void sCaling() {
      a *= (F/1)*(1/(L*L));
      mu *= (M/1)*(1/L)*(1/S)*1e3*1e-3;
      kappa *= ((L*L)/1)*1e-6;
      rho *= (M/1)*(1/(L*L*L));
      g *= (L/1)*(1/(S*S));
      APc *= (F/1)*(1/(L*L));
      APcZ *= (F/1)*(1/(L*L));
      APcZZ *= (F/1)*(1/(L*L));
    }

    double initalPcEval() const {
      return APc+APcZ*z+APcZZ*z*z;
    }

    void addOptions(po::options_description& o,const std::string& prefix) {
      o.add_options()
      ((prefix+".M").c_str(),po::value<double>(&M)->default_value(M))
      ((prefix+".L").c_str(),po::value<double>(&L)->default_value(L))
      ((prefix+".S").c_str(),po::value<double>(&S)->default_value(S))
      ((prefix+".F").c_str(),po::value<double>(&F)->default_value(F))
      ((prefix+".g").c_str(),po::value<double>(&g)->default_value(g))
      ((prefix+".rho").c_str(),po::value<double>(&rho)->default_value(rho))
      ((prefix+".kappa").c_str(),po::value<double>(&kappa)->default_value(kappa))
      ((prefix+".mu").c_str(),po::value<double>(&mu)->default_value(mu))
      ((prefix+".a").c_str(),po::value<double>(&a)->default_value(a))
      ((prefix+".thetaS").c_str(),po::value<double>(&thetaS)->default_value(thetaS))
      ((prefix+".thetaR").c_str(),po::value<double>(&thetaR)->default_value(thetaR))
      ((prefix+".m").c_str(),po::value<double>(&m)->default_value(m))
      ((prefix+".APc").c_str(),po::value<double>(&APc)->default_value(APc))
      ((prefix+".APcZ").c_str(),po::value<double>(&APcZ)->default_value(APcZ))
      ((prefix+".APcZZ").c_str(),po::value<double>(&APcZZ)->default_value(APcZZ));
      // TODO Add more parameters
    }

    void printMatParameters(const int id,const std::string& prefix) {
      PetscPrintf(
        PETSC_COMM_WORLD,"Mat name %s-%s block id %d\n",
        prefix.c_str(),matName.c_str(),id
      );
      PetscPrintf(PETSC_COMM_WORLD,"M=%6.4g [kg]\n",L);
      PetscPrintf(PETSC_COMM_WORLD,"L=%6.4g [m]\n",L);
      PetscPrintf(PETSC_COMM_WORLD,"S=%6.4g [s]\n",S);
      PetscPrintf(PETSC_COMM_WORLD,"F=%6.4g [N]\n",F);
      PetscPrintf(PETSC_COMM_WORLD,"g=%6.4g\n",g);
      PetscPrintf(PETSC_COMM_WORLD,"rho=%6.4g\n",rho);
      PetscPrintf(PETSC_COMM_WORLD,"thetaS=%6.4g\n",thetaS);
      PetscPrintf(PETSC_COMM_WORLD,"thetaR=%6.4g\n",thetaR);
      PetscPrintf(PETSC_COMM_WORLD,"mu=%6.4g\n",mu);
      PetscPrintf(PETSC_COMM_WORLD,"a=%6.4g\n",a);
      PetscPrintf(PETSC_COMM_WORLD,"m=%6.4g\n",m);
      PetscPrintf(PETSC_COMM_WORLD,"kappa=%6.4g\n",kappa);
      PetscPrintf(PETSC_COMM_WORLD,"APc=%6.4g\n",APc);
      PetscPrintf(PETSC_COMM_WORLD,"APcZ=%6.4g\n",APcZ);
      PetscPrintf(PETSC_COMM_WORLD,"APcZZ=%6.4g\n",APcZZ);
    }

  };

  struct SimpleDarcyProblem: public CommonMaterialData {

    SimpleDarcyProblem(const CommonMaterialData &data):
    CommonMaterialData(data) {
    }

    PetscErrorCode calK() {
      PetscFunctionBegin;
      K = kappa*(rho*g)/mu;
      PetscFunctionReturn(0);
    };

    PetscErrorCode calDiffK() {
      PetscFunctionBegin;
      diffK = 0;
      PetscFunctionReturn(0);
    };

    PetscErrorCode calC() {
      PetscFunctionBegin;
      C = rho;
      PetscFunctionReturn(0);
    }

    PetscErrorCode calDiffC() {
      PetscFunctionBegin;
      diffC = 0;
      PetscFunctionReturn(0);
    }

  };

  struct MaterialVanGenuchten: public CommonMaterialData {

    MaterialVanGenuchten(const CommonMaterialData &data):
    CommonMaterialData(data) {
    }

    inline void calSe() {
      if(Pc<=0) {
        Se = 1;
      } else {
        Se = pow(1+pow(Pc/a,1/(1-m)),-m);
      }
    }

    void printSe(const double b,const double e,const double s,const std::string prefix) {
      Pc = b;
      for(;Pc<e;Pc+=s) {
        calSe();
        cout << prefix+": " << Pc << " " << Se << endl;
      }
    }

    inline void calDiffSe() {
      if(Pc<=0) {
        diffSe = 0;
      } else {
        diffSe =
        -pow(0.1e1 + pow(Pc / a, 0.1e1 / (double) (1 - m)), -(double) m) * (double) m *
        pow(Pc / a, 0.1e1 / (double) (1 - m)) / (double) (1 - m) / Pc / (0.1e1 +
          pow(Pc / a, 0.1e1 / (double) (1 - m)));
        }
    }

    inline void calDiff2Se() {
      if(Pc<=0) {
        diff2Se = 0;
      } else {
        diff2Se = -(double) (int) pow((double) (-1 + m), (double) (-2)) * pow(Pc, -0.2e1) * (double) m *
        (-pow(0.1e1 + pow(Pc / a, -0.1e1 / (double) (-1 + m)), (double) (-m - 2)) * (double) m *
        pow(Pc / a, -(double) (2 / (-1 + m))) + pow(0.1e1 +
          pow(Pc / a, -0.1e1 / (double) (-1 + m)), (double) (-m - 1)) * (double) m *
          pow(Pc / a, -0.1e1 / (double) (-1 + m)) - pow(0.1e1 + pow(Pc / a, -0.1e1 / (double) (-1 + m)), (double) (-m - 2)) *
          pow(Pc / a, -(double) (2 / (-1 + m))));
      }
    }

    inline void calKappaR() {
      if(Se==0) {
        kappaR = 0;
      } else if(Se==1) {
        kappaR = 1;
      } else {
        const double a = pow(Se,1/m);
        const double b = pow(1-a,m);
        kappaR = sqrt(Se)*pow(1-b,2);
      }
    }

    inline void  calDiffKappaR() {
      if(Se==0) {
        diffSe = 0;
      } else if(Se==1) {
        diffSe = 0;
      } else {
        diffKappaR =
        pow(Se, -0.1e1 / 0.2e1) * pow(0.1e1 - pow(0.1e1 - pow(Se, 0.1e1 / m), m), 0.2e1) / 0.2e1 + 0.2e1 *
        pow(Se, -0.1e1 / 0.2e1) * (0.1e1 - pow(0.1e1 - pow(Se, 0.1e1 / m), m)) * pow(0.1e1 - pow(Se, 0.1e1 / m), m) *
        pow(Se, 0.1e1 / m) / (0.1e1 - pow(Se, 0.1e1 / m));
      }
    }

    inline void calTheta() {
      theta = Se*(thetaS-thetaR)+thetaR;
    }

    void printTheta(const double b,const double e,const double s,const std::string prefix) {
      Pc = b;
      for(;Pc<e;Pc+=s) {
        calSe();
        calTheta();
        cout << prefix+": " << Pc << " " << theta << endl;
      }
    }

    inline void calDiffTheta() {
      diffTheta = diffSe*(thetaS-thetaR);
    }

    inline void calDiff2Theta() {
      diff2Theta = diff2Se*(thetaS-thetaR);
    }

    PetscErrorCode calK() {
      PetscFunctionBegin;
      calSe();
      calKappaR();
      K = (rho*g/mu)*kappa*kappaR;
      PetscFunctionReturn(0);
    }

    void printKappaR(const double b,const double e,const double s,const std::string prefix) {
      Pc = b;
      for(;Pc<e;Pc+=s) {
        calSe();
        calKappaR();
        cout << prefix+": " << Pc << " " << kappaR << endl;
      }
    }

    PetscErrorCode calDiffK() {
      PetscFunctionBegin;
      calSe();
      calDiffSe();
      calDiffKappaR();
      diffK = (rho*g/mu)*kappa*diffKappaR*diffSe;
      PetscFunctionReturn(0);
    };

    void printDiffKappaR(const double b,const double e,const double s,const std::string prefix) {
      Pc = b;
      for(;Pc<e;Pc+=s) {
        calSe();
        calDiffKappaR();
        calDiffSe();
        cout << prefix+": " << Pc << " " << Se << " " << diffKappaR << " " << diffKappaR*diffSe << endl;
      }
    }

    PetscErrorCode calC() {
      PetscFunctionBegin;
      calDiffSe();
      calDiffTheta();
      C = -rho*diffTheta;
      PetscFunctionReturn(0);
    }

    void printC(const double b,const double e,const double s,const std::string prefix) {
      Pc = b;
      for(;Pc<e;Pc+=s) {
        calC();
        cout << prefix+": " << Pc << " " << C << endl;
      }
    }

    PetscErrorCode calDiffC() {
      PetscFunctionBegin;
      calDiff2Se();
      calDiff2Theta();
      diffC = -rho*diff2Theta;
      PetscFunctionReturn(0);
    }

  private:

    double Se;
    double diffSe;
    double diff2Se;
    double diff2Theta;
    double kappaR;
    double diffKappaR;
    double theta;
    double diffTheta;

  };

}

#endif //__MATERIALUNSATURATEDFLOW_HPP__
