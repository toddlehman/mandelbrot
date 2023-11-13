/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "mp_real.h"


//-----------------------------------------------------------------------------
// LINEAR INTERPOLATION

//-----------------------------------------------------------------------------
public_function
void mp_lerp(mp_real *x, mp_real t, mp_real x0, mp_real x1)
{
  assert(x);
                                   // *x = x0 + (t * (x1 - x0));
  mp_sub(*x, x1, x0, MP_ROUND);    // *x = (x1 - x0);
  mp_mul(*x, *x, t,  MP_ROUND);    // *x *= t;
  mp_add(*x, x0, *x, MP_ROUND);    // *x += x0;
}

//-----------------------------------------------------------------------------
public_function
void mp_lerp_d(mp_real *x, real t, mp_real x0, mp_real x1)
{
  assert(x);
                                   // *x = x0 + (t * (x1 - x0));
  mp_sub(*x, x1, x0, MP_ROUND);    // *x = (x1 - x0);
  mp_mul_d(*x, *x, t,  MP_ROUND);  // *x *= t;
  mp_add(*x, x0, *x, MP_ROUND);    // *x += x0;
}


//-----------------------------------------------------------------------------
