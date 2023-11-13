/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//----------------------------------------------------------------------------
// MAIN CONTROL DATA STRUCTURE

typedef struct
{
  mp_real  camera_x;   // x-coordinate
  mp_real  camera_y;   // y-coordinate
  mp_real  camera_z;   // z-coordinate
  mp_real  camera_d;   // distance to viewport
  real     camera_oz;  // z-axis rotational angle (neg=right, pos=left)
  real     camera_ox;  // x-axis rotational angle (neg=down, pos=up)
  real     camera_oy;  // y-axis rotational angle (neg=left, pos=right)

  mp_real  x_center;
  mp_real  y_center;
  mp_real  xy_min_size;

  uint64   iter_max;

  int      width_pixels;
  int      height_pixels;
  int      supersample_interior_min_depth;
  int      supersample_interior_max_depth;
  int      supersample_exterior_min_depth;
  int      supersample_exterior_max_depth;
  float32  supersample_solidarity;

  bool     output_statistics;
  bool     output_image_text_format;
}
MainControl;
