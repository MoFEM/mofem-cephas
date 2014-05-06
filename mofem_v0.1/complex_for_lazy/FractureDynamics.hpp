/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Test for linar elastic dynamics.
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __FRACTURE_DYNAMICS_HPP
#define __FRACTURE_DYNAMICS_HPP

/** \brief FE method for elastic dynamics
  *
  * M*u'' + K*u - F = 0
  *
  * F( t, [ dot_u, u], [ dot_u', u'] ) = [ 0 -1 ][ dot_u' ] + [ 1 0 ][ dot_u ] + [ 0    ] = [ 0 ]
  *                                      [ M  0 ][ u'     ]   [ 0 K ][ u     ]   [ F(t) ]   [ 0 ]
  * 
  * Fu  = [ 1 0 ][ dot_u ]
  *       [ 0 K ][ u     ] 
  * 
  * Fu' = [ 0 -1 ][ dot_u' ]
  *       [ M  0 ][ u'     ]
  *
  * G = Fu + a*Fu'
  *
  * dG/(du,ddot_u) = [ -1 0 ] + a*[ 0 1 ]
  *                  [  0 K ] +   [ M 0 ]
  *
 **/

#endif //__FRACTURE_DYNAMICS_HPP


