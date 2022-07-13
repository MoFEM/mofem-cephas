/**
 * \file IntegrationRules.hpp
 *
 * Integration rules
 *
 */



namespace MoFEM {

namespace IntRules {

namespace NCC {

#include <triangle_ncc_rule.h>
#include <tetrahedron_ncc_rule.h>

} // namespace NCC

namespace NCO {

#include <triangle_nco_rule.h>
#include <tetrahedron_nco_rule.h>

} // namespace NCO

} // namespace IntRules

}