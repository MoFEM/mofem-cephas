/** \file BCData.cpp
 * \brief Structures for boundary structures
 * 
 * \note Structures are native for Cubit only
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

std::ostream &operator<<(std::ostream &os, const DisplacementCubitBcData &e) {
  os << "\n";
  os << "D i s p l a c e m e n t \n \n";
  os << "Flag for X-Translation (0/1): " << (int)e.data.flag1 << "\n";
  os << "Flag for Y-Translation (0/1): " << (int)e.data.flag2 << "\n";
  os << "Flag for Z-Translation (0/1): " << (int)e.data.flag3 << "\n";
  os << "Flag for X-Rotation (0/1): " << (int)e.data.flag4 << "\n";
  os << "Flag for Y-Rotation (0/1): " << (int)e.data.flag5 << "\n";
  os << "Flag for Z-Rotation (0/1): " << (int)e.data.flag6 << "\n \n";

  if (e.data.flag1 == 1)
    os << "Displacement magnitude (X-Translation): " << e.data.value1 << "\n";
  else
    os << "Displacement magnitude (X-Translation): N/A"
       << "\n";
  if (e.data.flag2 == 1)
    os << "Displacement magnitude (Y-Translation): " << e.data.value2 << "\n";
  else
    os << "Displacement magnitude (Y-Translation): N/A"
       << "\n";
  if (e.data.flag3 == 1)
    os << "Displacement magnitude (Z-Translation): " << e.data.value3 << "\n";
  else
    os << "Displacement magnitude (Z-Translation): N/A"
       << "\n";
  if (e.data.flag4 == 1)
    os << "Displacement magnitude (X-Rotation): " << e.data.value4 << "\n";
  else
    os << "Displacement magnitude (X-Rotation): N/A"
       << "\n";
  if (e.data.flag5 == 1)
    os << "Displacement magnitude (Y-Rotation): " << e.data.value5 << "\n";
  else
    os << "Displacement magnitude (Y-Rotation): N/A"
       << "\n";
  if (e.data.flag6 == 1)
    os << "Displacement magnitude (Z-Rotation): " << e.data.value6 << "\n \n";
  else
    os << "Displacement magnitude (Z-Rotation): N/A"
       << "\n \n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const ForceCubitBcData &e) {
  os << "\n";
  os << "F o r c e \n \n";
  os << "Force magnitude: " << e.data.value1 << "\n";
  os << "Moment magnitude: " << e.data.value2 << "\n";
  os << "Force direction vector (X-component): " << e.data.value3 << "\n";
  os << "Force direction vector (Y-component): " << e.data.value4 << "\n";
  os << "Force direction vector (Z-component): " << e.data.value5 << "\n";
  os << "Moment direction vector (X-component): " << e.data.value6 << "\n";
  os << "Moment direction vector (Y-component): " << e.data.value7 << "\n";
  os << "Moment direction vector (Z-component): " << e.data.value8 << "\n \n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const VelocityCubitBcData &e) {
  os << "\n";
  os << "V e l o c i t y \n \n";
  if (e.data.flag1 == 1)
    os << "Velocity magnitude (X-Translation): " << e.data.value1 << "\n";
  else
    os << "Velocity magnitude (X-Translation): N/A"
       << "\n";
  if (e.data.flag2 == 1)
    os << "Velocity magnitude (Y-Translation): " << e.data.value2 << "\n";
  else
    os << "Velocity magnitude (Y-Translation): N/A"
       << "\n";
  if (e.data.flag3 == 1)
    os << "Velocity magnitude (Z-Translation): " << e.data.value3 << "\n";
  else
    os << "Velocity magnitude (Z-Translation): N/A"
       << "\n";
  if (e.data.flag4 == 1)
    os << "Velocity magnitude (X-Rotation): " << e.data.value4 << "\n";
  else
    os << "Velocity magnitude (X-Rotation): N/A"
       << "\n";
  if (e.data.flag5 == 1)
    os << "Velocity magnitude (Y-Rotation): " << e.data.value5 << "\n";
  else
    os << "Velocity magnitude (Y-Rotation): N/A"
       << "\n";
  if (e.data.flag6 == 1)
    os << "Velocity magnitude (Z-Rotation): " << e.data.value6 << "\n \n";
  else
    os << "Velocity magnitude (Z-Rotation): N/A"
       << "\n \n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const AccelerationCubitBcData &e) {
  os << "\n";
  os << "A c c e l e r a t i o n \n \n";
  if (e.data.flag1 == 1)
    os << "Acceleration magnitude (X-Translation): " << e.data.value1 << "\n";
  else
    os << "Acceleration magnitude (X-Translation): N/A"
       << "\n";
  if (e.data.flag2 == 1)
    os << "Acceleration magnitude (Y-Translation): " << e.data.value2 << "\n";
  else
    os << "Acceleration magnitude (Y-Translation): N/A"
       << "\n";
  if (e.data.flag3 == 1)
    os << "Acceleration magnitude (Z-Translation): " << e.data.value3 << "\n";
  else
    os << "Acceleration magnitude (Z-Translation): N/A"
       << "\n";
  if (e.data.flag4 == 1)
    os << "Acceleration magnitude (X-Rotation): " << e.data.value4 << "\n";
  else
    os << "Acceleration magnitude (X-Rotation): N/A"
       << "\n";
  if (e.data.flag5 == 1)
    os << "Acceleration magnitude (Y-Rotation): " << e.data.value5 << "\n";
  else
    os << "Acceleration magnitude (Y-Rotation): N/A"
       << "\n";
  if (e.data.flag6 == 1)
    os << "Acceleration magnitude (Z-Rotation): " << e.data.value6 << "\n \n";
  else
    os << "Acceleration magnitude (Z-Rotation): N/A"
       << "\n \n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const TemperatureCubitBcData &e) {
  os << "\n";
  os << "T e m p e r a t u r e \n \n";
  if (e.data.flag1 == 1)
    os << "Temperature: " << e.data.value1 << "\n";
  else
    os << "Temperature (default case): N/A"
       << "\n";
  if (e.data.flag2 == 1)
    os << "Temperature (thin shell middle): " << e.data.value2 << "\n";
  else
    os << "Temperature (thin shell middle): N/A"
       << "\n";
  if (e.data.flag3 == 1)
    os << "Temperature (thin shell gradient): " << e.data.value3 << "\n";
  else
    os << "Temperature (thin shell gradient): N/A"
       << "\n";
  if (e.data.flag4 == 1)
    os << "Temperature (thin shell top): " << e.data.value4 << "\n";
  else
    os << "Temperature (thin shell top): N/A"
       << "\n";
  if (e.data.flag5 == 1)
    os << "Temperature (thin shell bottom): " << e.data.value5 << "\n \n";
  else
    os << "Temperature (thin shell bottom): N/A"
       << "\n \n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const PressureCubitBcData &e) {
  os << "\n";
  os << "P r e s s u r e \n \n";
  os << "Pressure flag2: " << (int)e.data.flag2 << "\n";
  os << "Pressure value: " << e.data.value1 << "\n \n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const HeatFluxCubitBcData &e) {
  os << "\n";
  os << "H e a t  F l u x \n \n";
  if (e.data.flag1 == 1)
    os << "Heat flux value: " << e.data.value1 << "\n";
  else
    os << "Heat flux is applied on thin shells"
       << "\n";
  if (e.data.flag2 == 1)
    os << "Heat flux value (thin shell top): " << e.data.value2 << "\n";
  else
    os << "Heat flux value (thin shell top): N/A"
       << "\n";
  if (e.data.flag3 == 1)
    os << "Heat flux value (thin shell bottom): " << e.data.value3 << "\n \n";
  else
    os << "Heat flux value (thin shell bottom): N/A"
       << "\n \n";
  return os;
}

std::ostream &operator<<(std::ostream &os, const CfgCubitBcData &e) {
  os << "\n";
  os << "CFD BC \n \n";
  return os;
}

} // namespace MoFEM
