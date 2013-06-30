/* Copyright (C) 2009, Lukasz Kaczmarczyk (likask AT civil.gla.ac.uk)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of mofem.
 * mofem is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * mofem is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with mofem. If not, see <http://www.gnu.org/licenses/>. */


//Fibres, by R. Eberlein, G.A. Holzapfel, and C.A.J. Schulze-Bauer. 
//An anisotropic model for annulus tissue and enhanced finite element analyses of intact lumbar disc bodies. 
//Computer Methods in Biomechanics and Biomedical Engineering, 4(3):209–229, 2001. ISSN 1025-5842.
typedef struct {
  enum phisical_equation_volume eq_solid;
  double k1,k2;
  //those are in material space
  double fibre_vector_a1[3];
  double fibre_vector_a2[3];
} ctx_EberleinHolzapfel1;


