/** \file MaterialBlocks.cpp
 * \brief Structures for material blocks
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

namespace MoFEM {

std::ostream &operator<<(std::ostream &os, const BlockSetAttributes &e) {
  os << std::endl << "Blcok attributes" << std::endl;
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
