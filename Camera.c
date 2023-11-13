/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Camera.h"

// TODO:  Make this configurable by a runtime switch.
// TODO:  Factor out into an orthographic vs. perspective switch so that
//        MAP_WORLD and FLAT_WORLD can be combined.
#define  MAP_WORLD        1
#define  FLAT_WORLD       0
#define  SPHERICAL_WORLD  0


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
                      real camera_theta,
                      real camera_phi,
                      real camera_roll,
                      real camera_fov)
{
  #if 0  // OBSOLETE -- Was only for debugging.
  fprintf(stderr, "target=(%f,%f) cd=%f ct=%f cp=%f vt=%f vp=%f vr=%f fov=%f\n",
          mp_get_d(target_x), mp_get_d(target_y),
          (double)target_camera_rho,
          (double)target_camera_theta,
          (double)target_camera_phi,
          (double)camera_theta,
          (double)camera_phi,
          (double)camera_roll,
          (double)camera_fov);
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
  this->camera_theta        = camera_theta;
  this->camera_phi          = camera_phi;
  this->camera_roll         = camera_roll;
  this->camera_fov          = camera_fov;

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
  #if MAP_WORLD
    mp_add_d(this->camera_x, this->target_x, trx);
    mp_add_d(this->camera_y, this->target_y, try);
    mp_set_d(this->camera_z,                 trz);
  #elif FLAT_WORLD
    mp_add_d(this->camera_x, this->target_x, trx);
    mp_add_d(this->camera_y, this->target_y, try);
    mp_set_d(this->camera_z,                 trz);
  #elif SPHERICAL_WORLD
    Vector3 s = camera_argand_to_sphere(
                  this, mp_get_d(target_x), mp_get_d(target_y));
    s = vector3_scaled(s, 1); // + target_camera_rho);
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
  assert(fabs((double)vector3_magnitude(s) - 1) < 0.000001);
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
//         sphere_radius specifies the radius of the sphere.
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
                                      real sphere_radius,
                                      Vector3 *sphere_point)
{
  assert(this);

  ray = vector3_normalized(ray);  // Convert to unit vector.

  real radicand = pow(vector3_dot_product(ray, point), 2)
                - vector3_dot_product(point, point)
                + pow(sphere_radius, 2);

  if (radicand >= 0)
  {
    real distance = - vector3_dot_product(ray, point) - sqrt(radicand);
    //fprintf(stderr, "distance=%+f\n", (double)distance);

    if (distance >= 0)
    {
      *sphere_point = vector3_sum(point, vector3_scaled(ray, distance));
      // Verify that the point is indeed on the sphere.
      assert(fabs((double)(vector3_magnitude(*sphere_point) - sphere_radius))
             < 0.000001);
      return true;
    }
    else
    {
      *sphere_point = (Vector3) { .x = 0, .y = 0, .z = 0 };
      return false;
    }
  }
  else
  {
    *sphere_point = (Vector3) { .x = 0, .y = 0, .z = 0 };
    return false;
  }
}


//-----------------------------------------------------------------------------
// DISTANCE THROUGH SPHERE
//
// ENTRY:  point specifies the starting point of the ray.
//         ray specifies the direction of the ray.
//         sphere_radius specifies the radius of the sphere.
//
// EXIT:   Returns the distance the ray travels through the sphere.

private_method
real camera_distance_through_sphere(const Camera *this,
                                    Vector3 point, Vector3 ray,
                                    real sphere_radius)
{
  assert(this);

  ray = vector3_normalized(ray);  // Convert to unit vector.

  real radicand = pow(vector3_dot_product(ray, point), 2)
                - vector3_dot_product(point, point)
                + pow(sphere_radius, 2);

  if (radicand >= 0)
  {
    real distance_near = - vector3_dot_product(ray, point) - sqrt(radicand);
    real distance_far  = - vector3_dot_product(ray, point) + sqrt(radicand);

    if (distance_near >= 0)
    {
      return distance_far - distance_near;
    }
    else
    {
      return (distance_far >= 0)? distance_far : 0;
    }
  }
  else
  {
    return 0;
  }
}


#if 0  // OBSOLETE
//-----------------------------------------------------------------------------
// COMPUTE DISTANCE TO LINE FROM ORIGIN
//
// Arithmetic derived from:
//   <http://geomalgorithms.com/a02-_lines.html>

private_method
real camera_distance_to_line_from_origin(const Camera *this,
                                         Vector3 point, Vector3 ray)
{
  assert(this);

  real c1 = vector3_dot_product(vector3_negated(point), ray);
  real c2 = vector3_dot_product(ray, ray);

  Vector3 closest = vector3_sum(point, vector3_scaled(ray, c1/c2));

  real distance = vector3_magnitude(closest);
  return distance;
}
#endif


//-----------------------------------------------------------------------------
// MAP VIEWPORT COORDINATES TO ARGAND PLANE COORDINATES
//
// ENTRY:  u and v are viewport coordinates, typically in the interval
//          [-1,+1] for each axis.
//         x and y point to the desired Argand plane coordinates.
//         atmo_scat points to the desired atmospheric scattering parameter.
//          (Set this to null if this value is unneeded.)
//
// EXIT:   Returns a success code.
//           true   *x and *y contain the Argand plane coordinates.
//           false  *x and *y are undefined; the point is outside the plane.
//         If atmo_scat was non-null on entry, *atmo_scat contains the distance
//          traveled through atmosphere (valid regardless of the return value).
//
public_method
bool camera_get_argand_point(const Camera *this,
                             real u, real v,
                             mp_real *x, mp_real *y,
                             real *atmo_scat)
{
  assert(this); assert(x); assert(y);

  // FIXME:  This routine only works for standard-precision floating-point
  // numbers (e.g., real).  It needs to be upgraded to work with mp_real.


  // --- Transform viewport coordinates for camera roll rotation.
  {
    Vector3 viewport = { .x = u, .y = v, .z = 0 };
    viewport = vector3_rotated_z(viewport, this->camera_roll);
    u = viewport.x;
    v = viewport.y;
  }


  // --- Compute viewport radius from distance and field-of-view.

  real viewport_radius = tan(this->camera_fov / 2) * this->target_camera_rho;


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

#if MAP_WORLD

  assert(mp_get_d(this->camera_z) > 0);
  assert(this->target_camera_rho > 0);
  assert(this->target_camera_theta == 0);
  assert(this->target_camera_phi == 0);
  viewport_radius = this->target_camera_rho / 2;
  mp_set_d(*x, position.x + (u * viewport_radius));
  mp_set_d(*y, position.y + (v * viewport_radius));
  unless (atmo_scat == NULL)
    *atmo_scat = 0;
  return true;

#elif FLAT_WORLD

  Vector3 ray =
  {
    .x = u * viewport_radius,
    .y = v * viewport_radius,
    .z = -this->target_camera_rho
  };

  ray = vector3_rotated_x(ray, this->target_camera_phi   + this->camera_phi);
  ray = vector3_rotated_z(ray, this->target_camera_theta + this->camera_theta);


  // --- Extend the ray to meet the Argand plane and note the coordinates of
  //     intersection.

  if ((position.z > 0) && (ray.z < 0))
  {
    real t = unlerp(0, position.z, position.z + ray.z);
    real x_t = lerp(t, position.x, position.x + ray.x);
    real y_t = lerp(t, position.y, position.y + ray.y);
    mp_set_d(*x, x_t);
    mp_set_d(*y, y_t);
    //fprintf(stderr, "(u,v)=(%f,%f) (x,y)=(%f,%f)\n",
    //        (double)u, (double)v, (double)x_t, (double)y_t);
    unless (atmo_scat == NULL)
    {
      // Note to self:  This calculation isn't very good.  If flat worlds are
      // ever actually desired, this should be revisited.
      real ray_xy = sqrt(pow(position.x - x_t, 2) +
                         pow(position.y - y_t, 2) +
                         pow(position.z - 0,   2));
      real damper = 5;
      *atmo_scat = atan(ray_xy / damper) / (PI / 2);
      //fprintf(stderr, "ray_xy:%f atmo_scat:%f\n",
      //        (double)ray_xy, (double)*atmo_scat);
    }
    return true;
  }
  else
  {
    mp_set_d(*x, 0);
    mp_set_d(*y, 0);
    //fprintf(stderr, "(u,v)=(%f,%f) (x,y)=(undef,undef)\n",
    //        (double)u, (double)v);

    // This is a big fat KLUDGE.
    // FIXME:  This is probably no longer workable for flat world.
    unless (atmo_scat == NULL)
    {
      // Note to self:  Surprisingly, this calculation turns out to look pretty
      // nice -- especially for square images looking almost horizontally (i.e.,
      // where the sky becomes deep blue at the top).
      real ray_xy = sqrt((ray.x * ray.x) + (ray.y * ray.y));
      *atmo_scat = 1 - atan(ray.z / ray_xy) / (PI / 2);
      //fprintf(stderr, "ray_xy:%f atmo_scat:%f\n",
      //        (double)ray_xy, (double)*atmo_scat);
    }

    return false;
  }

#elif SPHERICAL_WORLD

  Vector3 ray =
  {
    .x = u * viewport_radius,
    .y = v * viewport_radius,
    .z = -this->target_camera_rho
  };

  ray = vector3_rotated_x(ray, this->target_camera_phi   + this->camera_phi);
  ray = vector3_rotated_z(ray, this->target_camera_theta + this->camera_theta);

  real globe_radius = 1.00;
  real atmosphere_radius = 1.05;
  real max_atmo_dist = 0.640308;  // KLUDGE from observation

  Vector3 s;
  if (camera_intersect_ray_with_sphere(this, position, ray, globe_radius, &s))
  {
    real a, b;
    camera_sphere_to_argand(this, s, &a, &b);
    mp_set_d(*x, a);
    mp_set_d(*y, b);

    unless (atmo_scat == NULL)
    {
      real camera_radius = vector3_magnitude(position);

      real f;
      if (camera_radius > atmosphere_radius)
      {
        // Outside the atmospheric radius.
        Vector3 sa;
        (void)camera_intersect_ray_with_sphere(this, position, ray,
                                               atmosphere_radius, &sa);
        real d = vector3_magnitude(vector3_difference(sa, s));
        f = d / (max_atmo_dist / 2);
      }
      else
      {
        // Inside the atmospheric radius.
        real d = vector3_magnitude(vector3_difference(position, s));
        real f1 = d / (max_atmo_dist / 2);

        real f2 = atan(d / this->target_camera_rho) / (PI / 2);
        f2 = pow(f2, 4) * 0.5;

        f = lerp(unlerp(camera_radius, globe_radius, atmosphere_radius),
                        f2, f1);
      }
      if (f < 0) f = 0;  // FIXME:  Use a clamp function instead.
      if (f > 1) f = 1;

      // Apply a power curve to the computed atmospheric coefficient.  This is
      // just a tuning value.
      real foo = lerp(unlerp(camera_radius, globe_radius, atmosphere_radius),
                      1.0, 2.0);
      if (foo < 1.0) foo = 1.0;  // FIXME:  Use a clamp function instead.
      if (foo > 2.0) foo = 2.0;
      f = pow(f, foo);

      *atmo_scat = f;
    }

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
    unless (atmo_scat == NULL)
    {
      real d = camera_distance_through_sphere(this, position, ray,
                                              atmosphere_radius);
      real camera_radius = vector3_magnitude(position);

      // Temper the atmospheric calculation by the distance light travels
      // through the atmosphere.
      real f;
      if (camera_radius > atmosphere_radius)
      {
        // Outside the atmospheric radius.
        f = d / max_atmo_dist;
      }
      else
      {
        // Inside the atmospheric radius.  Temper the effective maximum
        // atmospheric distance by the current altitude.
        real mad = lerp(unlerp(camera_radius, globe_radius, atmosphere_radius),
                        max_atmo_dist / 2,
                        max_atmo_dist);
        f = d / mad;
      }
      if (f < 0) f = 0;  // FIXME:  Use a clamp function instead.
      if (f > 1) f = 1;

      // Apply a power curve to the computed atmospheric coefficient.
      real foo = lerp(unlerp(camera_radius, globe_radius, atmosphere_radius),
                      0.5, 4.0);
      if (foo < 0.5) foo = 0.5;  // FIXME:  Use a clamp function instead.
      if (foo > 4.0) foo = 4.0;
      f = pow(f, foo);
      f = pow(f, 2);  // Why?

      *atmo_scat = f;

      //fprintf(stderr, "(%+.3f,%+.3f) --> (%+.3f,%+.3f,%+.3f)+(%+.3f,%+.3f,%+.3f) --> atmo_scat %f\n",
      //        (double)u, (double)v,
      //        (double)position.x, (double)position.y, (double)position.z,
      //        (double)ray.x, (double)ray.y, (double)ray.z,
      //        (double)*atmo_scat);
    }

    return false;
  }

#endif
}

