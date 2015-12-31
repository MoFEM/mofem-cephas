/** \file CoordSysMultiIndices.hpp
 * \ingroup coordsys_multi_indices
 * \brief Coordinate systems attached to DOFs
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __COORDSYSMULTIINDICES_HPP__
#define __COORDSYSMULTIINDICES_HPP__

namespace MoFEM {

  /** \brief Structure for Coordinate system of two-point tensor

  Because tensors have a geometric interpretation, their underlying physical
  meaning is independent of the coordinate system in which they are defined; the
  numerical value of a tensor depends on the coordinate system, and so any
  numerical instantiation of a tensor field is incomplete without a specification
  of its coordinate system.

  In this data structure a generic two-point tensors are considered, see \cite
  marsden1994mathematical. To each field \ref MoFEM::MoFEMField a CoordSys data
  structure can be attached, carrying information about tensor structure
  approximated by the field.

  */
  struct CoordSys {

    EntityHandle meshSet; 		      ///< keeps entities for this meshset
    const int* tagIdData;
    const int* tagCoordSysDim;
    const char* tagCoordSysName; 		///< tag keeps name of the field
    int tagCoordSysNameSize;
    CoordSys(Interface &moab,const EntityHandle meshset);
    inline int getId() const { return *tagIdData; };

    /** \brief Get tensor dimension

    This is general two point tensor \f$\mathbf{T}\f$ of type
    \f[
    \left(
    \begin{array}{cc}
    q & l \\
    p & m
    \end{array}
    \right)
    \f]
    at point \f$\mathbf{X} \in \mathcal{B}\f$ over mapping
    \f[
    \phi: \mathcal{B} \to \mathcal{S}
    \f]
    is a multilinear mapping
    \f[
    \mathbf{T}:
    (
    T^*_\mathbf{X}\mathcal{B}\times\dots\times T^*_\mathbf{X}\mathcal{B}
    ) \times
    (
    T_\mathbf{X}\mathcal{B}\times\dots\times T_\mathbf{X}\mathcal{B}
    ) \times
    (
    T^*_\mathbf{x}\mathcal{S}\times\dots\times T^*_\mathbf{x}\mathcal{S}
    ) \times
    (
    T_\mathbf{x}\mathcal{B}\times\dots\times T_\mathbf{x}\mathcal{B}
    )
    \to \mathbb{R}
    \f]
    where \f$\mathbf{x} = \phi(\mathbf{X})\f$. See details in \cite marsden1994mathematical

    \param d = 0,1,2,3 is equivalent to q,l,o and m respectively.

    */
    inline int getDim(const int d = 0) const {
      return tagCoordSysDim[d];
    };

    virtual PetscErrorCode get_E_Base(const double m[]) const {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode get_E_DualBase(const double m[]) const {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode get_e_Base(const double m[]) const {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode get_e_DualBase(const double m[]) const {
      PetscFunctionBegin;
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      PetscFunctionReturn(0);
    }

    inline EntityHandle getMeshSet() const { return meshSet; };
    inline boost::string_ref getNameRef() const { return boost::string_ref((char *)tagCoordSysName,tagCoordSysNameSize); };
    inline string getName() const { return string((char *)tagCoordSysName,tagCoordSysNameSize); };
  };

  typedef multi_index_container<
    CoordSys,
    indexed_by<
      ordered_non_unique<
        tag<CoordSysID_mi_tag>, const_mem_fun<CoordSys,int,&CoordSys::getId> >,
      ordered_unique<
        tag<Meshset_mi_tag>, member<CoordSys,EntityHandle,&CoordSys::meshSet>
      >,
      ordered_unique<
        tag<CoordSysName_mi_tag >, const_mem_fun<CoordSys,boost::string_ref,&CoordSys::getNameRef>
      >
  > > CoordSys_multiIndex;

}

#endif //__COORDSYSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup coordsys_multi_indices Coordinate system of general tenor field
 * \ingroup mofem
 ******************************************************************************/
