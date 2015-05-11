/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#ifndef __BCS_RVELAGRANGE_TRAC_HPP
#define __BCS_RVELAGRANGE_TRAC_HPP

namespace MoFEM {
  
  struct BCs_RVELagrange_Trac: public BCs_RVELagrange_Disp {
    
    BCs_RVELagrange_Trac(FieldInterface &m_field): BCs_RVELagrange_Disp(m_field){}

    
    

    
  };

}

#endif
