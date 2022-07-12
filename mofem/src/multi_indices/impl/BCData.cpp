/** \file BCData.cpp
 * \brief Structures for boundary structures
 * 
 * \note Structures are native for Cubit only
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
