#pragma once

#include<CL/sycl.hpp>

#include "parameters.h"
#include "trace.h"
#include "maths.h"
#include "particles.h"
#include "meshes.h"
#include "FESolver.h"

extern SYCL_EXTERNAL bool d_XtoLtet(Particle &part, Volume &volume);

extern SYCL_EXTERNAL bool d_XtoLtet2(Particle &part, Volume &volume);