/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "MPReal.h"


//-----------------------------------------------------------------------------
// LINEAR INTERPOLATION

//-----------------------------------------------------------------------------
public_function
void mp_lerp(mp_real *x,
             const mp_real t, const mp_real x0, const mp_real x1)
{
  assert(x);
                         // *x = x0 + (t * (x1 - x0));
  mp_sub(*x, x1, x0);    // *x = (x1 - x0);
  mp_mul(*x, *x, t);     // *x *= t;
  mp_add(*x, x0, *x);    // *x += x0;
}

//-----------------------------------------------------------------------------
public_function
void mp_lerp_d(mp_real *x,
               const real t, const mp_real x0, const mp_real x1)
{
  assert(x);
                         // *x = x0 + (t * (x1 - x0));
  mp_sub(*x, x1, x0);    // *x = (x1 - x0);
  mp_mul_d(*x, *x, t);   // *x *= t;
  mp_add(*x, x0, *x);    // *x += x0;
}


//-----------------------------------------------------------------------------
// GET MINIMUM PRECISION FROM STRING

public_function
int mp_get_min_prec_from_string(const char *str)
{
  int digit_count = 0;

  for (const char *p = str; *p; p++)
    if (isdigit(*p))  // KLUDGE
      digit_count++;

  return (int)ceil(digit_count * log(10) / log(2));
}
