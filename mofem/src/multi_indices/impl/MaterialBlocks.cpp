/** \file MaterialBlocks.cpp
 * \brief Structures for material blocks
 */


namespace MoFEM {

std::ostream &operator<<(std::ostream &os, const BlockSetAttributes &e) {
  os << std::endl << "Block attributes" << std::endl;
  os << "-------------------" << std::endl;
  os << "User attribute 1 = " << e.data.User1 << std::endl;
  os << "User attribute 2 = " << e.data.User2 << std::endl;
  os << "User attribute 3 = " << e.data.User3 << std::endl;
  os << "User attribute 4 = " << e.data.User4 << std::endl;
  os << "User attribute 5 = " << e.data.User5 << std::endl;
  os << "User attribute 6 = " << e.data.User6 << std::endl;
  os << "User attribute 7 = " << e.data.User7 << std::endl;
  os << "User attribute 8 = " << e.data.User7 << std::endl;
  os << "User attribute 9 = " << e.data.User7 << std::endl;
  os << "User attribute 10 = " << e.data.User10 << std::endl << std::endl;

  os << "User attribute 11 = " << e.data.User11 << std::endl << std::endl;
  os << "User attribute 12 = " << e.data.User12 << std::endl << std::endl;
  os << "User attribute 13 = " << e.data.User13 << std::endl << std::endl;
  os << "User attribute 14 = " << e.data.User14 << std::endl << std::endl;
  os << "User attribute 15 = " << e.data.User15 << std::endl << std::endl;
  os << "User attribute 16 = " << e.data.User16 << std::endl << std::endl;
  os << "User attribute 17 = " << e.data.User17 << std::endl << std::endl;
  os << "User attribute 18 = " << e.data.User18 << std::endl << std::endl;
  os << "User attribute 19 = " << e.data.User19 << std::endl << std::endl;
  os << "User attribute 20 = " << e.data.User20 << std::endl << std::endl;

  os << "User attribute 21 = " << e.data.User21 << std::endl << std::endl;
  os << "User attribute 22 = " << e.data.User22 << std::endl << std::endl;
  os << "User attribute 23 = " << e.data.User23 << std::endl << std::endl;
  os << "User attribute 24 = " << e.data.User24 << std::endl << std::endl;
  os << "User attribute 25 = " << e.data.User25 << std::endl << std::endl;
  os << "User attribute 26 = " << e.data.User26 << std::endl << std::endl;
  os << "User attribute 27 = " << e.data.User27 << std::endl << std::endl;
  os << "User attribute 28 = " << e.data.User28 << std::endl << std::endl;
  os << "User attribute 29 = " << e.data.User29 << std::endl << std::endl;
  os << "User attribute 30 = " << e.data.User30 << std::endl << std::endl;

  os << "User attribute 31 = " << e.data.User31 << std::endl << std::endl;
  os << "User attribute 32 = " << e.data.User32 << std::endl << std::endl;
  os << "User attribute 33 = " << e.data.User33 << std::endl << std::endl;
  os << "User attribute 34 = " << e.data.User34 << std::endl << std::endl;
  os << "User attribute 35 = " << e.data.User35 << std::endl << std::endl;
  os << "User attribute 36 = " << e.data.User36 << std::endl << std::endl;
  os << "User attribute 37 = " << e.data.User37 << std::endl << std::endl;
  os << "User attribute 38 = " << e.data.User38 << std::endl << std::endl;
  os << "User attribute 39 = " << e.data.User39 << std::endl << std::endl;
  os << "User attribute 40 = " << e.data.User40 << std::endl << std::endl;

  os << "User attribute 41 = " << e.data.User41 << std::endl << std::endl;
  os << "User attribute 42 = " << e.data.User42 << std::endl << std::endl;
  os << "User attribute 43 = " << e.data.User43 << std::endl << std::endl;
  os << "User attribute 44 = " << e.data.User44 << std::endl << std::endl;
  os << "User attribute 45 = " << e.data.User45 << std::endl << std::endl;
  os << "User attribute 46 = " << e.data.User46 << std::endl << std::endl;
  os << "User attribute 47 = " << e.data.User47 << std::endl << std::endl;
  os << "User attribute 48 = " << e.data.User48 << std::endl << std::endl;
  os << "User attribute 49 = " << e.data.User49 << std::endl << std::endl;
  os << "User attribute 50 = " << e.data.User50 << std::endl << std::endl;

  return os;
}

std::ostream &operator<<(std::ostream &os, const Mat_Elastic &e) {
  os << std::endl << "Material Properties" << std::endl;
  os << "-------------------" << std::endl;
  os << "Young's modulus  = " << e.data.Young << std::endl;
  os << "Poisson's ratio  = " << e.data.Poisson << std::endl;
  os << "Thermal expansion = " << e.data.ThermalExpansion << std::endl;
  os << "User attribute 1 = " << e.data.User1 << std::endl;
  os << "User attribute 2 = " << e.data.User2 << std::endl;
  os << "User attribute 3 = " << e.data.User3 << std::endl;
  os << "User attribute 4 = " << e.data.User4 << std::endl;
  os << "User attribute 5 = " << e.data.User5 << std::endl;
  os << "User attribute 6 = " << e.data.User6 << std::endl;
  os << "User attribute 7 = " << e.data.User7 << std::endl << std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os,
                         const Mat_Elastic_EberleinHolzapfel1 &e) {
  os << std::endl << "Material Properties" << std::endl;
  os << "-------------------" << std::endl;
  os << "Young's modulus  = " << e.data.Young << std::endl;
  os << "Poisson's ratio  = " << e.data.Poisson << std::endl;
  os << "k1 = " << e.data.k1 << std::endl;
  os << "k2 = " << e.data.k2 << std::endl;
  os << "a0_x = " << e.data.a0x << std::endl;
  os << "a0_y = " << e.data.a0y << std::endl;
  os << "a0_z = " << e.data.a0z << std::endl;
  os << "a1_x = " << e.data.a1x << std::endl;
  os << "a1_y = " << e.data.a1y << std::endl;
  os << "a1_Z = " << e.data.a1z << std::endl << std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const Mat_Thermal &e) {
  os << std::endl << "Material Properties" << std::endl;
  os << "-------------------" << std::endl;
  os << "Conductivity  = " << e.data.Conductivity << std::endl;
  os << "User attribute 1 = " << e.data.HeatCapacity << std::endl;
  os << "User attribute 2 = " << e.data.User2 << std::endl;
  os << "User attribute 3 = " << e.data.User3 << std::endl;
  os << "User attribute 4 = " << e.data.User4 << std::endl;
  os << "User attribute 5 = " << e.data.User5 << std::endl;
  os << "User attribute 6 = " << e.data.User6 << std::endl;
  os << "User attribute 7 = " << e.data.User7 << std::endl;
  os << "User attribute 8 = " << e.data.User8 << std::endl << std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const Mat_Moisture &e) {
  os << std::endl << "Material Properties" << std::endl;
  os << "-------------------" << std::endl;
  os << "Diffusivity  = " << e.data.Diffusivity << std::endl;
  os << "Viscosity = " << e.data.Viscosity << std::endl;
  os << "Permeability = " << e.data.Permeability << std::endl;
  os << "User attribute 3 = " << e.data.User3 << std::endl;
  os << "User attribute 4 = " << e.data.User4 << std::endl;
  os << "User attribute 5 = " << e.data.User5 << std::endl;
  os << "User attribute 6 = " << e.data.User6 << std::endl;
  os << "User attribute 7 = " << e.data.User7 << std::endl;
  os << "User attribute 8 = " << e.data.User8 << std::endl << std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const Block_BodyForces &e) {
  os << std::endl << "Block Body Forces" << std::endl;
  os << "-------------------" << std::endl;
  os << "density  = " << e.data.density << std::endl;
  os << "acceleration_x = " << e.data.acceleration_x << std::endl;
  os << "acceleration_y = " << e.data.acceleration_y << std::endl;
  os << "acceleration_z = " << e.data.acceleration_z << std::endl;
  os << "User attribute 4 = " << e.data.User4 << std::endl;
  os << "User attribute 5 = " << e.data.User5 << std::endl;
  os << "User attribute 6 = " << e.data.User6 << std::endl;
  os << "User attribute 7 = " << e.data.User7 << std::endl;
  os << "User attribute 8 = " << e.data.User8 << std::endl << std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const Mat_Elastic_TransIso &e) {
  os << std::endl << "Material Properties" << std::endl;
  os << "-------------------" << std::endl;
  os << "Young's modulus in xy plane (Ep)     = " << e.data.Youngp << std::endl;
  os << "Young's modulus in z-direction (Ez)  = " << e.data.Youngz << std::endl;
  os << "Poisson's ratio in xy plane (vp)     = " << e.data.Poissonp
     << std::endl;
  os << "Poisson's ratio in z-direction (vpz) = " << e.data.Poissonpz
     << std::endl;
  os << "Shear modulus in z-direction (Gzp)   = " << e.data.Shearzp << std::endl
     << std::endl;
  return os;
}

std::ostream &operator<<(std::ostream &os, const Mat_Interf &e) {
  os << std::endl << "Material Properties" << std::endl;
  os << "-------------------" << std::endl;
  os << "Elastic module	= " << e.data.alpha << std::endl << std::endl;
  os << "Damage coupling	= " << e.data.beta << std::endl << std::endl;
  os << "Strengh		= " << e.data.ft << std::endl << std::endl;
  os << "Fracture energy	= " << e.data.Gf << std::endl << std::endl;

  return os;
}

} // namespace MoFEM
