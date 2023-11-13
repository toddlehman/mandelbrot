/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// MANDELBROT CALCULATION STRUCTURE

typedef struct
{
  uint64  iter_max;
  int     precision_bits;
  mpfr_t  periodicity_epsilon;
}
Mandelbrot;


//-----------------------------------------------------------------------------
// MANDELBROT RESULT STRUCTURE
//
// This is packed because it must be storable inside a pixel array.  Speed of
// access is not important here, but space is.

#pragma pack(push, 4)

typedef struct
{
  float64  dwell;
  uint64   iter;
  uint64   period;
}
MandelbrotResult;

#pragma pack(pop)


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
MandelbrotResult mandelbrot_result_exterior(const uint64 iter,
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
// FUNCTION PROTOTYPES

extern_public_constructor
  Mandelbrot *mandelbrot_create(uint64 iter_max,
                                int precision_bits,
                                mpfr_t periodicity_epsilon);

extern_public_destructor
  void mandelbrot_destroy(Mandelbrot **this);

extern_public_method
  MandelbrotResult mandelbrot_compute(Mandelbrot *this, mpfr_t cx, mpfr_t cy);

