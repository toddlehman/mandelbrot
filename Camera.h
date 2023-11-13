/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// CAMERA STRUCTURE

typedef struct
{
  // Target point in Argand plane
  mp_real  target_x;     // Real component
  mp_real  target_y;     // Imaginary component

  // Relative position of camera (relative to target point)
  real     target_camera_rho;    // Distance to camera from target point
  real     target_camera_theta;  // Angle in xy plane
  real     target_camera_phi;    // Angle from positive z-axis

  // Viewport parameters
  real     viewport_tilt;      // Angle from vertical (neg=cc-wise, pos=c-wise)
  //real     viewport_scale;     // Reciprocal of magnification
  real     viewport_fov;       // Center angle of frustum (0 to π/2)
  //real     viewport_distance;  // Distance from camera

  // Absolute camera position (derived from the above)
  mp_real  camera_x;     // x-coordinate
  mp_real  camera_y;     // y-coordinate
  mp_real  camera_z;     // z-coordinate
}
Camera;


//-----------------------------------------------------------------------------
// PUBLIC INLINE METHODS


//-----------------------------------------------------------------------------
// PUBLIC METHODS

extern_public_constructor
  Camera *camera_create(mp_real target_x,
                        mp_real target_y,
                        real target_camera_rho,
                        real target_camera_theta,
                        real target_camera_phi,
                        real viewport_tilt,
                        real viewport_fov);

extern_public_destructor
  void camera_destroy(Camera **p_this);

extern_public_method
  bool camera_get_argand_point(const Camera *this,
                               real u, real v,
                               mp_real *x, mp_real *y);

