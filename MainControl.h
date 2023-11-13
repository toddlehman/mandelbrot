/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//----------------------------------------------------------------------------
// MAIN CONTROL DATA STRUCTURE

typedef struct
{
  // TODO:  These are currently unused.  They belong in the Mandelbrot
  // iteration object.
  mp_real  julia_x;
  mp_real  julia_y;

  // FIXME:  These 7 parameters belong in a camera object, or as calls to
  // property setters of a camera object.
  mp_real  target_x;
  mp_real  target_y;
  real     target_camera_rho;
  real     target_camera_theta;
  real     target_camera_phi;
  real     viewport_tilt;
  real     viewport_fov;

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
