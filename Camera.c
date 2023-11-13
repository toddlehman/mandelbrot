/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Camera.h"


//-----------------------------------------------------------------------------
// CREATE CAMERA OBJECT

public_constructor
Camera *camera_create(mp_real target_x,
                      mp_real target_y,
                      real target_camera_rho,
                      real target_camera_theta,
                      real target_camera_phi,
                      real viewport_tilt,
                      real viewport_fov)
{
  #if 0  // OBSOLETE -- Was only for debugging.
  fprintf(stderr, "target=(%f,%f) rho=%f theta=%f phi=%f tilt=%f fov=%f\n",
          mp_get_d(target_x), mp_get_d(target_y),
          (double)target_camera_rho,
          (double)target_camera_theta,
          (double)target_camera_phi,
          (double)viewport_tilt,
          (double)viewport_fov);
  #endif

  Camera *this = mem_alloc_clear(1, sizeof(Camera));

  // Store non-derived values.
  mp_init2(this->target_x, mp_get_prec(target_x));
  mp_init2(this->target_y, mp_get_prec(target_y));
  mp_set(this->target_x, target_x);
  mp_set(this->target_y, target_y);

  this->target_camera_rho   = target_camera_rho;
  this->target_camera_theta = target_camera_theta;
  this->target_camera_phi   = target_camera_phi;

  this->viewport_tilt = viewport_tilt;
  this->viewport_fov  = viewport_fov;

  // Compute target-relative coordinates of camera.
  real phi = this->target_camera_phi;
  real theta = this->target_camera_theta + (TAU * 3 / 4);
  real trx = target_camera_rho * sin(phi) * cos(theta);
  real try = target_camera_rho * sin(phi) * sin(theta);
  real trz = target_camera_rho * cos(phi);

  // Compute absolute position of camera.
  mp_init2(this->camera_x, mp_get_prec(target_x));  // KLUDGE on precision
  mp_init2(this->camera_y, mp_get_prec(target_x));  // KLUDGE on precision
  mp_init2(this->camera_z, mp_get_prec(target_x));  // KLUDGE on precision
  mp_add_d(this->camera_x, this->target_x, trx);
  mp_add_d(this->camera_y, this->target_y, try);
  mp_set_d(this->camera_z,                 trz);

  #if 0  // OBSOLETE -- Was only for debugging.
  fprintf(stderr, "camera=(%f,%f,%f)\n",
          mp_get_d(this->camera_x),
          mp_get_d(this->camera_y),
          mp_get_d(this->camera_z));
  #endif

  return this;
}


//-----------------------------------------------------------------------------
// DESTROY CAMERA OBJECT

extern_public_destructor
void camera_destroy(Camera **p_this)
{
  assert(p_this);
  Camera *this = *p_this;
  assert(this);

  mp_clear(this->target_x);
  mp_clear(this->target_y);

  mp_clear(this->camera_x);
  mp_clear(this->camera_y);
  mp_clear(this->camera_z);

  mem_dealloc(p_this);
}


//-----------------------------------------------------------------------------
// MAP VIEWPORT COORDINATES TO ARGAND PLANE COORDINATES
//
// ENTRY:  u and v are viewport coordinates, typically in the interval
//          [-1,+1] for each axis.
//
// EXIT:   Returns a success code.
//         If successful, then *x and *y contain the Argand plane coordinates.
//         Otherwise, *sky_angle contains the angle in radians of the cast ray
//           as it hits the sky (0 is straight up; π/2 radians is the horizon).
//
public_method
bool camera_get_argand_point(const Camera *this,
                             real u, real v,
                             mp_real *x, mp_real *y,
                             real *sky_angle)
{
  assert(this); assert(x); assert(y);

  // FIXME:  This routine only works for standard-precision floating-point
  // numbers (e.g., real).  It needs to be upgraded to work with mp_real.


  // --- Transform viewport coordinates for camera tilt rotation.
  {
    Vector3 viewport = { .x = u, .y = v, .z = 0 };
    viewport = vector3_rotated_z(viewport, this->viewport_tilt);
    u = viewport.x;
    v = viewport.y;
  }


  // --- Compute viewport radius from distance and field-of-view.

  real viewport_radius = tan(this->viewport_fov / 2) * this->target_camera_rho;


  // --- Start with camera position.

  Vector3 position =
  {
    .x = mp_get_d(this->camera_x),
    .y = mp_get_d(this->camera_y),
    .z = mp_get_d(this->camera_z)
  };


  // --- Cast a ray from the camera to the viewport, using the camera's
  //     orientation (given by 2 of 3 Euler angles).

  Vector3 ray =
  {
    .x = u * viewport_radius,
    .y = v * viewport_radius,
    .z = -this->target_camera_rho
  };

  ray = vector3_rotated_x(ray, this->target_camera_phi);
  ray = vector3_rotated_z(ray, this->target_camera_theta);


  // --- Extend the ray to meet the Argand plane and note the coordinates of
  //     intersection.

  if ((position.z > 0) && (ray.z < 0))
  {
    real t = unlerp(0, position.z, position.z + ray.z);
    mp_set_d(*x, lerp(t, position.x, position.x + ray.x));
    mp_set_d(*y, lerp(t, position.y, position.y + ray.y));
    //fprintf(stderr, "(u,v)=(%f,%f) (x,y)=(%f,%f)\n",
    //        (double)u, (double)v,
    //        mp_get_d(*x), mp_get_d(*y));
    return true;
  }
  else
  {
    mp_set_d(*x, 0);
    mp_set_d(*y, 0);
    //fprintf(stderr, "(u,v)=(%f,%f) (x,y)=(undef,undef)\n",
    //        (double)u, (double)v);
    #if 0  // OBSOLETE
    return false;
    #endif

    // This is a big fat KLUDGE.
    real xy = sqrt((ray.x * ray.x) + (ray.y * ray.y));
    *sky_angle = (PI / 2) - atan(ray.z / xy);
    return false;
  }
}


//-----------------------------------------------------------------------------
// MAP 3-D COORDINATES TO VIEWPORT COORDINATES

#if 0
public_method
bool camera_get_viewport_point(const Camera *this,
                               mp_real x, mp_real y, mp_real z,
                               real *u, real *v)
{
  assert(this); assert(u); assert(v);

  // FIXME:  This routine only works for standard-precision floating-point
  // numbers (e.g., real).  It needs to be upgraded to work with mp_real.


  // --- Compute viewport radius from distance and field-of-view.

  real viewport_radius = tan(this->viewport_fov / 2) * this->target_camera_rho;


  // --- Start with camera position.

  Vector3 position =
  {
    .x = mp_get_d(this->camera_x),
    .y = mp_get_d(this->camera_y),
    .z = mp_get_d(this->camera_z)
  };


  // --- Cast a ray from the camera to the viewport, using the camera's
  //     orientation (given by 2 of 3 Euler angles).

  Vector3 ray =
  {
    .x = u * viewport_radius,
    .y = v * viewport_radius,
    .z = -this->target_camera_rho
  };

  ray = vector3_rotated_x(ray, this->target_camera_phi);
  ray = vector3_rotated_z(ray, this->target_camera_theta);


  // --- Extend the ray to meet the Argand plane and note the coordinates of
  //     intersection.

  if ((position.z > 0) && (ray.z < 0))
  {
    real t = unlerp(0, position.z, position.z + ray.z);
    mp_set_d(*x, lerp(t, position.x, position.x + ray.x));
    mp_set_d(*y, lerp(t, position.y, position.y + ray.y));
    //fprintf(stderr, "(u,v)=(%f,%f) (x,y)=(%f,%f)\n",
    //        (double)u, (double)v,
    //        mp_get_d(*x), mp_get_d(*y));
    return true;
  }
  else
  {
    mp_set_d(*x, 0);
    mp_set_d(*y, 0);
    //fprintf(stderr, "(u,v)=(%f,%f) (x,y)=(undef,undef)\n",
    //        (double)u, (double)v);
    #if 0  // OBSOLETE
    return false;
    #endif

    // This is a big fat KLUDGE.
    real xy = sqrt((ray.x * ray.x) + (ray.y * ray.y));
    *sky_angle = (PI / 2) - atan(ray.z / xy);
    return false;
  }
}
#endif

