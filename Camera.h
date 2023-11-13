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

  // Absolute camera position (derived from other values here)
  mp_real  camera_x;     // x-coordinate
  mp_real  camera_y;     // y-coordinate
  mp_real  camera_z;     // z-coordinate

  // Relative position of camera (relative to target point)
  real     target_camera_rho;    // Distance to camera from target point
  real     target_camera_theta;  // Angle in xy plane
  real     target_camera_phi;    // Angle from positive z-axis

  // Viewing angle of camera (additive to target-relative angles)
  real     camera_theta;     // Additional angle in xy plane
  real     camera_phi;       // Additional angle from positive z-axis
  real     camera_roll;      // Angle from vertical (neg=cc-wise, pos=c-wise)

  // Viewport parameters
  real     camera_fov;       // Center angle of frustum (0 to π/2)
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
                        real camera_theta,
                        real camera_phi,
                        real camera_roll,
                        real camera_fov);

extern_public_destructor
  void camera_destroy(Camera **p_this);

extern_public_method
  bool camera_get_argand_point(const Camera *this,
                               real u, real v,
                               mp_real *x, mp_real *y,
                               real *sky_angle);

#if 0  // OBSOLETE -- Not needed ever?
extern_public_method
  bool camera_get_viewport_point(const Camera *this,
                                 mp_real x, mp_real y, mp_real z,
                                 real *u, real *v);
#endif
