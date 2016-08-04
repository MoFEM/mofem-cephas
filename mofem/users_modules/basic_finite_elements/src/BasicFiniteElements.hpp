/** \file BasicFiniteElements.hpp
  \ingroup Header file for basic finite elements implementation
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

#ifndef __BASICFINITEELEMENTS_HPP__
#define __BASICFINITEELEMENTS_HPP__

#include <MoFEM.hpp>
using namespace MoFEM;

#ifdef WITH_ADOL_C
  #include <adolc/adolc.h>
#endif // WITH_ADOL_C

extern "C" {
  void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
    double circumcenter[3],double *xi,double *eta,double *zeta);
  void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
    double circumcenter[3],double *xi,double *eta);
}

#include <cholesky.hpp>

#include <MethodForForceScaling.hpp>
#include <ArcLengthTools.hpp>
#include <BodyForce.hpp>
#include <ConstrainMatrixCtx.hpp>
#include <DirichletBC.hpp>
#include <EdgeForce.hpp>
#include <FieldApproximation.hpp>
#include <FluidPressure.hpp>
#include <KelvinVoigtDamper.hpp>
#include <NodalForce.hpp>
#include <NonLinearElasticElement.hpp>
#include <NitscheMethod.hpp>
#include <PCMGSetUpViaApproxOrders.hpp>
#include <PostProcOnRefMesh.hpp>
#include <PostProcStresses.hpp>
#include <PostProcHookStresses.hpp>
#include <Smoother.hpp>
#include <SurfacePressure.hpp>
#include <SurfaceSlidingConstrains.hpp>
#include <ThermalElement.hpp>
#include <ThermalStressElement.hpp>
#include <TimeForceScale.hpp>
#include <UltraWeakTransportElement.hpp>
#include <VolumeCalculation.hpp>
#include <ConvectiveMassElement.hpp>


#endif // __BASICFINITEELEMENTS_HPP__
