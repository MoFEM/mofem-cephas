/** \file Common.hpp
 * \brief Basic structures and data
 */



#ifndef __COMMON_HPP__
#define __COMMON_HPP__

namespace MoFEM {

const EntityHandle no_handle =
    0; ///< No entity handle is indicated by zero handle, i.e. root meshset

} // namespace MoFEM

#include <Exceptions.hpp>
#include <Types.hpp>
#include <Templates.hpp>
#include <PetscSmartObj.hpp>

#endif //__COMMON_HPP__

/**
 * \defgroup mofem MoFEM
 */
