/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//----------------------------------------------------------------------------
// MAIN CONTROL DATA STRUCTURE

typedef struct
{
  mp_real  x_center;
  mp_real  y_center;
  mp_real  xy_min_size;
  real     roll;
  real     pitch;
  real     yaw;

  uint64   iter_max;

  int      width_pixels;
  int      height_pixels;
  int      supersample_int_min_depth;
  int      supersample_int_max_depth;
  int      supersample_ext_min_depth;
  int      supersample_ext_max_depth;
  float32  supersample_solidarity;

  bool     output_statistics;
  bool     output_image_text_format;
}
MainControl;
