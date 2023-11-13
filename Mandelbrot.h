/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// MANDELBROT RESULT STRUCTURE
//
// This is packed because it must be storable inside a pixel array.  Speed of
// access is not important here, but space is.

#pragma pack(push, 4)

typedef struct
{
  float64  dwell;    // Mandelbrot "dwell" value.  (8 bytes)
  uint64   iter;     // Number of iterations computed.  (8 bytes)
  uint64   period;   // Period of orbit, if detected.  (8 bytes)
}
MandelbrotResult;    // (24 bytes)

#pragma pack(pop)


//-----------------------------------------------------------------------------
// MANDELBROT CONFIGURATION STRUCTURE
//
// This contains the parameters governing computation of samples.

typedef struct
{
  uint64   iter_max;
  int      mp_prec;
}
MandelbrotConfiguration;


//-----------------------------------------------------------------------------
// MANDELBROT RESULT STATISTICS STRUCTURE
//
// This keeps track of various statistics that are useful to know when choosing
// parameters or just to satisfy curiosity.

typedef struct
{
  // For tracking both interior and exterior points.
  uint64  total_iter;
  uint64  total_probes;
  uint64  total_probes_uniterated;

  // For tracking interior points.
  uint64  interior_iter;
  uint64  interior_probes;
  uint64  interior_probes_uniterated;
  uint64  interior_probes_aperiodic;
  uint64  interior_probes_by_log2_iter[64];

  // For tracking exterior points.
  uint64  exterior_iter;
  uint64  exterior_probes;
  uint64  exterior_probes_uniterated;
  uint64  exterior_probes_by_log2_iter[64];
}
MandelbrotResultStatistics;


//-----------------------------------------------------------------------------
// MANDELBROT CALCULATION STRUCTURE

typedef struct
{
  MandelbrotConfiguration     conf;
  MandelbrotResultStatistics  stats;
}
Mandelbrot;


//-----------------------------------------------------------------------------
// PUBLIC INLINE FUNCTIONS

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_undefined(const MandelbrotResult mr)
{
  return (mr.dwell == 0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_defined(const MandelbrotResult mr)
{
  return (mr.dwell != 0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_exterior(const MandelbrotResult mr)
{
  return (mr.dwell != 0) && (mr.dwell != INFINITY);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_exterior_uniterated(const MandelbrotResult mr)
{
  return (mr.dwell != 0) && (mr.dwell != INFINITY) && (mr.iter == 0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_interior(const MandelbrotResult mr)
{
  return (mr.dwell == INFINITY);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_interior_periodic(const MandelbrotResult mr)
{
  return (mr.dwell == INFINITY) && (mr.period > 0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_interior_aperiodic(const MandelbrotResult mr)
{
  return (mr.dwell == INFINITY) && (mr.period == 0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_interior_iterated(const MandelbrotResult mr)
{
  return (mr.dwell == INFINITY) && (mr.iter > 0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool mandelbrot_result_is_interior_uniterated(const MandelbrotResult mr)
{
  return (mr.dwell == INFINITY) && (mr.iter == 0);
}

//-----------------------------------------------------------------------------
public_inline_function
MandelbrotResult mandelbrot_result_undefined(void)
{
  return (MandelbrotResult)
  {
    .dwell   = 0,
    .iter    = 0,
    .period  = 0,
  };
}

//-----------------------------------------------------------------------------
public_inline_function
MandelbrotResult mandelbrot_result_interior_iterated_periodic(const uint64 iter,
                                                              const uint64 period)
{
  return (MandelbrotResult)
  {
    .dwell   = INFINITY,
    .iter    = iter,
    .period  = period,
  };
}

//-----------------------------------------------------------------------------
public_inline_function
MandelbrotResult mandelbrot_result_interior_iterated_aperiodic(const uint64 iter)
{
  return (MandelbrotResult)
  {
    .dwell   = INFINITY,
    .iter    = iter,
    .period  = 0,
  };
}

//-----------------------------------------------------------------------------
public_inline_function
MandelbrotResult mandelbrot_result_interior_uniterated(void)
{
  return (MandelbrotResult)
  {
    .dwell   = INFINITY,
    .iter    = 0,
    .period  = 0,
  };
}

//-----------------------------------------------------------------------------
public_inline_function
MandelbrotResult mandelbrot_result_exterior_iterated(const uint64 iter,
                                                     const float64 dwell)
{
  return (MandelbrotResult)
  {
    .dwell   = dwell,
    .iter    = iter,
    .period  = 0,
  };
}

//-----------------------------------------------------------------------------
public_inline_function
MandelbrotResult mandelbrot_result_exterior_uniterated(const float64 dwell)
{
  return (MandelbrotResult)
  {
    .dwell   = dwell,
    .iter    = 0,
    .period  = 0,
  };
}

//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_function
  real mandelbrot_max_scalar_value_during_iteration();

extern_public_constructor
  Mandelbrot *mandelbrot_create(uint64 iter_max,
                                int mp_prec);

extern_public_destructor
  void mandelbrot_destroy(Mandelbrot **this);

extern_public_method
  MandelbrotResult mandelbrot_compute(Mandelbrot *this,
                                      const mp_real cx, const mp_real cy,
                                      const mp_real periodicity_epsilon);

