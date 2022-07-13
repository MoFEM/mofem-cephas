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
