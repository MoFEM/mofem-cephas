/**
 * \file IntegrationRules.cpp
 *
 * Integration rules
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

namespace MoFEM {

namespace IntRules {

namespace NCC {

#include <triangle_ncc_rule.c>
#include <tetrahedron_ncc_rule.c>

} // namespace NCC

namespace NCO {

#include <triangle_nco_rule.c>
#include <tetrahedron_nco_rule.c>

} // namespace NCO

} // namespace IntRules

} // namespace MoFEM