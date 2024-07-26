/** \file BaseFunction.hpp
\brief General implementation of base function

*/



#ifndef __BASEFUNCTION_HPP__
#define __BASEFUNCTION_HPP__

namespace MoFEM {

struct BaseFunctionUnknownInterface : public UnknownInterface {

  using UnknownInterface::UnknownInterface;

  virtual MoFEMErrorCode
  query_interface(boost::typeindex::type_index type_index,
                  UnknownInterface **iface) const = 0;

  virtual ~BaseFunctionUnknownInterface() = default;
};

/**
 * \brief Base class used to exchange data between element data structures and
 * class calculating base functions \ingroup mofem_base_functions
 */
struct BaseFunctionCtx : public BaseFunctionUnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  using BaseFunctionUnknownInterface::BaseFunctionUnknownInterface;
};

/**
 * \brief Base class if inherited used to calculate base functions
 * \ingroup mofem_base_functions
 */
struct BaseFunction : public BaseFunctionUnknownInterface {

  struct DofsSideMapData {
    EntityType type;
    int side;
    int dof;
  };

  /**
   * @brief Map entity stype and side to element/entity dof index
   *
   * Such map is used to establish connection between dofs in the interior for
   * broken specs. Is assume that trace of interior on given side is not zero.
   *
   */
  using DofsSideMap = multi_index_container<

      DofsSideMapData,

      indexed_by<
          ordered_unique<
              tag<TypeSide_mi_tag>,
              composite_key<

                  DofsSideMapData,
                  member<DofsSideMapData, EntityType, &DofsSideMapData::type>,
                  member<DofsSideMapData, int, &DofsSideMapData::side>>>,

          ordered_unique<tag<EntDofIdx_mi_tag>,
                         member<DofsSideMapData, int, &DofsSideMapData::dof>>

          >>;

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 MoFEM::UnknownInterface **iface) const;

  using BaseFunctionUnknownInterface::BaseFunctionUnknownInterface;

  virtual MoFEMErrorCode getValue(MatrixDouble &pts,
                                  boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

  virtual MoFEMErrorCode getValue(MatrixDouble &pts_x, MatrixDouble &pts_t,
                                  boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

};

} // namespace MoFEM

#endif //__BASEFUNCTION_HPP__

/**
 * \defgroup mofem_base_functions Base functions
 *
 * \brief Calculation of base functions at integration points.
 *
 * \ingroup mofem
 **/
