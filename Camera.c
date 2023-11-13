/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Camera.h"

// FIXME:  Make this configurable by a runtime switch.
#define  FLAT_WORLD       0
#define  SPHERICAL_WORLD  1


//-----------------------------------------------------------------------------
// FORWARD DECLARATIONS

private_method
Vector3 camera_argand_to_sphere(const Camera *this, real a, real b);


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
  #if FLAT_WORLD
    mp_add_d(this->camera_x, this->target_x, trx);
    mp_add_d(this->camera_y, this->target_y, try);
    mp_set_d(this->camera_z,                 trz);
  #elif SPHERICAL_WORLD
    Vector3 s = camera_argand_to_sphere(
                  this, mp_get_d(target_x), mp_get_d(target_y));
    s = vector3_scaled(s, 1); // + target_camera_rho);
    // TODO:  Take theta and phi into account.
    phi += 0; theta += 0;  // (Temporary) Avoid compiler warnings.
    mp_set_d(this->camera_x, s.x + trx);
    mp_set_d(this->camera_y, s.y + try);
    mp_set_d(this->camera_z, s.z + trz);
  #endif

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

public_destructor
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
// MAP ARGAND PLANE TO UNIT SPHERE SURFACE

private_method
Vector3 camera_argand_to_sphere(const Camera *this, real a, real b)
{
  assert(this);

  // Translate entire Argand plane so that target point is at origin.
  a -= mp_get_d(this->target_x);
  b -= mp_get_d(this->target_y);

  // KLUDGE:  This information really belongs in the Mandelbrot class.
  //real scale = 1.25;  // Scale unit square by this much.
  real scale = 2.00;  // Scale unit square by this much.
  real oy = lerp(unlerp(a, -scale, +scale), -PI/2, +PI/2);
  real ox = lerp(unlerp(b, -scale, +scale), +PI/2, -PI/2);

  Vector3 s = { .x = 0, .y = 0, .z = 1 };
  s = vector3_rotated_x(s, ox);
  s = vector3_rotated_y(s, oy);

  return s;
}


//-----------------------------------------------------------------------------
// MAP UNIT SPHERE SURFACE TO ARGAND PLANE

private_method
void camera_sphere_to_argand(const Camera *this, Vector3 s, real *a, real *b)
{
  assert(this);
  assert(fabs(vector3_magnitude(s) - 1) < 0.000001);
  assert(a); assert(b);

  real oy = atan2(s.x, s.z);

       if (s.y < -1) s.y = -1;
  else if (s.y > +1) s.y = +1;
  real ox = asin(s.y);

  // KLUDGE:  This information really belongs in the Mandelbrot class.
  //real scale = 1.25;  // Scale unit square by this much.
  real scale = 2.00;  // Scale unit square by this much.
  *a = lerp(unlerp(oy, -PI/2, +PI/2), -scale, +scale);
  *b = lerp(unlerp(ox, -PI/2, +PI/2), -scale, +scale);

  // Translate entire Argand plane so that target point is at origin.
  *a += mp_get_d(this->target_x);
  *b += mp_get_d(this->target_y);

  // KLUDGE:  This information really belongs in the Mandelbrot class.
  //fprintf(stderr, "s=(%+f,%+f,%+f) oy=%+f ox=%+f (a,b)=(%+f,%+f)\n",
  //        (double)s.x, (double)s.y, (double)s.z, (double)oy, (double)ox,
  //        (double)*a, (double)*b);
}


//-----------------------------------------------------------------------------
// INTERSECT RAY WITH UNIT SPHERE
//
// ENTRY:  point specifies the starting point of the ray.
//         ray specifies the direction of the ray.
//
// EXIT:   Returns a result code.
//          true   *sphere contains the point on the unit sphere closest to the
//                 starting point of the ray.
//          false  The ray does not intersect the unit sphere, in which case
//                 *sphere is undefined (do not rely on it being set to zero).
//
// Arithmetic derived from:
//   <http://en.wikipedia.org/wiki/Line–sphere_intersection>

private_method
bool camera_intersect_ray_with_sphere(const Camera *this,
                                      Vector3 point, Vector3 ray,
                                      Vector3 *sphere)
{
  assert(this);

  ray = vector3_normalized(ray);  // Convert to unit vector.

  real radicand = pow(vector3_dot_product(ray, point), 2)
                - vector3_dot_product(point, point)
                + 1;

  if (radicand >= 0)
  {
    real distance = - vector3_dot_product(ray, point) - sqrt(radicand);
    //fprintf(stderr, "distance=%+f\n", (double)distance);

    if (distance >= 0)
    {
      *sphere = vector3_sum(point, vector3_scaled(ray, distance));
      // Verify that the point is indeed on unit sphere.
      assert(fabs(vector3_magnitude(*sphere) - 1) < 0.000001);
      return true;
    }
    else
    {
      *sphere = (Vector3) { .x = 0, .y = 0, .z = 0 };
      return false;
    }
  }
  else
  {
    *sphere = (Vector3) { .x = 0, .y = 0, .z = 0 };
    return false;
  }
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
  //fprintf(stderr, "(%f,%f,%f)\n",
  //        (double)position.x, (double)position.y, (double)position.z);


  // --- Cast a ray from the camera to the viewport, using the camera's
  //     orientation (given by 2 of 3 Euler angles).

#if FLAT_WORLD

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

    // This is a big fat KLUDGE.
    real xy = sqrt((ray.x * ray.x) + (ray.y * ray.y));
    *sky_angle = (PI / 2) - atan(ray.z / xy);
    return false;
  }

#elif SPHERICAL_WORLD

  Vector3 ray =
  {
    .x = u * viewport_radius,
    .y = v * viewport_radius,
    .z = -this->target_camera_rho
  };

  ray = vector3_rotated_x(ray, this->target_camera_phi);
  ray = vector3_rotated_z(ray, this->target_camera_theta);

  Vector3 s;
  if (camera_intersect_ray_with_sphere(this, position, ray, &s))
  {
    real a, b;
    camera_sphere_to_argand(this, s, &a, &b);
    mp_set_d(*x, a);
    mp_set_d(*y, b);
    //fprintf(stderr, "(%+.3f,%+.3f) --> (%+.3f,%+.3f,%+.3f)+(%+.3f,%+.3f,%+.3f) --> (%+.3f,%+.3f,%+.3f) --> (%+.3f,%+.3f)\n",
    //        (double)u, (double)v,
    //        (double)position.x, (double)position.y, (double)position.z,
    //        (double)ray.x, (double)ray.y, (double)ray.z,
    //        (double)s.x, (double)s.y, (double)s.z,
    //        (double)a, (double)b);
    return true;
  }
  else
  {
    //fprintf(stderr, "(%+.3f,%+.3f) --> (%+.3f,%+.3f,%+.3f)+(%+.3f,%+.3f,%+.3f) --> space\n",
    //        (double)u, (double)v,
    //        (double)position.x, (double)position.y, (double)position.z,
    //        (double)ray.x, (double)ray.y, (double)ray.z);
    *sky_angle = 0;
    return false;
  }

#endif
}

