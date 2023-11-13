/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// CAMERA STRUCTURE

typedef struct
{
  // Camera position
  mp_real  camera_x;     // x-coordinate
  mp_real  camera_y;     // y-coordinate
  mp_real  camera_z;     // z-coordinate
  mp_real  camera_d;     // distance to viewport

  // Camera orientation
  real     camera_ox;    // x-axis rotational angle (neg=down, pos=up)
  real     camera_oy;    // y-axis rotational angle (neg=left, pos=right)
  real     camera_oz;    // z-axis rotational angle (neg=right, pos=left)
}
Camera;


//-----------------------------------------------------------------------------
// PUBLIC INLINE METHODS


//-----------------------------------------------------------------------------
// PUBLIC METHODS

extern_public_constructor
  Camera *camera_create(mp_real x, mp_real y, mp_real z, mp_real d,
                        real ox, real oy, real oz);

extern_public_destructor
  void camera_destroy(Camera **p_this);

extern_public_method
  bool camera_get_argand_point(Camera *this,
                               real u, real v,
                               mp_real *x, mp_real *y);

