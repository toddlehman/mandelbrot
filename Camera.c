/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Camera.h"


//-----------------------------------------------------------------------------
// CREATE CAMERA OBJECT

public_constructor
Camera *camera_create(mp_real x, mp_real y, mp_real z, mp_real d,
                      real ox, real oy, real oz)
{
  Camera *this = mem_alloc_clear(1, sizeof(Camera));

  mp_init2(this->camera_x, mp_get_prec(x));
  mp_set(this->camera_x, x);

  mp_init2(this->camera_y, mp_get_prec(y));
  mp_set(this->camera_y, y);

  mp_init2(this->camera_z, mp_get_prec(z));
  mp_set(this->camera_z, z);

  mp_init2(this->camera_d, mp_get_prec(d));
  mp_set(this->camera_d, d);

  this->camera_ox = ox;
  this->camera_oy = oy;
  this->camera_oz = oz;

  return this;
}


//-----------------------------------------------------------------------------
// DESTROY CAMERA OBJECT

extern_public_destructor
void camera_destroy(Camera **p_this)
{
  assert(p_this);
  Camera *this = *p_this;

  mp_clear(this->camera_x);
  mp_clear(this->camera_y);
  mp_clear(this->camera_z);
  mp_clear(this->camera_d);

  assert(this);
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
//         Otherwise, *x and *y are set to zero.
//
public_method
bool camera_get_argand_point(Camera *this,
                             real u, real v,
                             mp_real *x, mp_real *y)
{
  assert(this); assert(x); assert(y);

  // FIXME:  This routine only works for standard-precision floating-point
  // numbers (e.g., real).  It needs to be upgraded to work with mp_real.


  // --- Start with camera position.

  Vector3 position = (Vector3)
  {
    .x = mp_get_d(this->camera_x),
    .y = mp_get_d(this->camera_y),
    .z = mp_get_d(this->camera_z)
  };


  // --- Cast a ray from the camera to the viewport, using the camera's
  //     orientation (3 Euler angles).

  Vector3 ray = (Vector3)
  {
    .x = u,
    .y = v,
    .z = -mp_get_d(this->camera_d)
  };

  ray = vector3_rotated_z(ray, this->camera_oz);
  ray = vector3_rotated_x(ray, this->camera_ox);
  ray = vector3_rotated_y(ray, this->camera_oy);


  // --- Extend the ray to meet the Argand plane and note the coordinates of
  //     intersection.

  if ((position.z > 0) && (ray.z < 0))
  {
    real t = unlerp(0, position.z, position.z + ray.z);
    mp_set_d(*x, lerp(t, position.x, position.x + ray.x));
    mp_set_d(*y, lerp(t, position.y, position.y + ray.y));
    return true;
  }
  else
  {
    mp_set_d(*x, 0);
    mp_set_d(*y, 0);
    return false;
  }
}

