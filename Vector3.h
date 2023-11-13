/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// VECTOR

typedef struct
{
  real  x;
  real  y;
  real  z;
}
Vector3;


//-----------------------------------------------------------------------------
// PUBLIC INLINE FUNCTIONS

//-----------------------------------------------------------------------------
public_inline_function
real vector3_magnitude(const Vector3 a)
{
  return sqrtl((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
}


//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_normalized(const Vector3 a)
{
  real m = vector3_magnitude(a);

  return (Vector3)
  {
    .x = a.x / m,
    .y = a.y / m,
    .z = a.z / m
  };
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_negated(const Vector3 a)
{
  return (Vector3)
  {
    .x = -a.x,
    .y = -a.y,
    .z = -a.z
  };
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_scaled(const Vector3 a, const real s)
{
  return (Vector3)
  {
    .x = a.x * s,
    .y = a.y * s,
    .z = a.z * s
  };
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_sum(const Vector3 a, const Vector3 b)
{
  return (Vector3)
  {
    .x = a.x + b.x,
    .y = a.y + b.y,
    .z = a.z + b.z
  };
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_difference(const Vector3 a, const Vector3 b)
{
  return (Vector3)
  {
    .x = a.x - b.x,
    .y = a.y - b.y,
    .z = a.z - b.z
  };
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_lerp(const real t, const Vector3 a, const Vector3 b)
{
  return (Vector3)
  {
    .x = lerp(t, a.x, b.x),
    .y = lerp(t, a.y, b.y),
    .z = lerp(t, a.z, b.z)
  };
}

//-----------------------------------------------------------------------------
public_inline_function
real vector3_dot_product(const Vector3 a, const Vector3 b)
{
  return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

//-----------------------------------------------------------------------------
public_inline_function
real vector3_normalized_dot_product(const Vector3 a, const Vector3 b)
{
  return vector3_dot_product(vector3_normalized(a),
                             vector3_normalized(b));
}

//-----------------------------------------------------------------------------
public_inline_function
real vector3_angle(const Vector3 a, const Vector3 b)
{
  real cos_angle = vector3_normalized_dot_product(a, b);

  // Very important:  Due to floating-point rounding errors, it is possible to
  // obtain a value slightly less than -1 or slightly greater than 1.  This line
  // of code clamps the value to the unit interval [-1,1].  If this is omitted,
  // then acos() will return NaN for values of cos_angle outside that interval.
       if (cos_angle < -1) cos_angle = -1;
  else if (cos_angle > +1) cos_angle = +1;

  real angle = acosl(cos_angle);

  return angle;
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_cross_product(const Vector3 a, const Vector3 b)
{
  return (Vector3)
  {
    .x = (a.y * b.z) - (a.z * b.y),
    .y = (a.z * b.x) - (a.x * b.z),
    .z = (a.x * b.y) - (a.y * b.x)
  };
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_rotated_z(const Vector3 a, const real oz)
{
  if (oz != 0)
  {
    real cos_oz = cosl(oz), sin_oz = sinl(oz);

    return (Vector3)
    {
      .z = a.z,
      .x = (a.x * cos_oz) - (a.y * sin_oz),
      .y = (a.x * sin_oz) + (a.y * cos_oz)
    };
  }
  else
  {
    return a;
  }
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_rotated_y(const Vector3 a, const real oy)
{
  if (oy != 0)
  {
    real cos_oy = cosl(oy), sin_oy = sinl(oy);

    return (Vector3)
    {
      .y = a.y,
      .z = (a.z * cos_oy) - (a.x * sin_oy),
      .x = (a.z * sin_oy) + (a.x * cos_oy)
    };
  }
  else
  {
    return a;
  }
}

//-----------------------------------------------------------------------------
public_inline_function
Vector3 vector3_rotated_x(const Vector3 a, const real ox)
{
  if (ox != 0)
  {
    real cos_ox = cosl(ox), sin_ox = sinl(ox);

    return (Vector3)
    {
      .x = a.x,
      .y = (a.y * cos_ox) - (a.z * sin_ox),
      .z = (a.y * sin_ox) + (a.z * cos_ox)
    };
  }
  else
  {
    return a;
  }
}


//-----------------------------------------------------------------------------
// PUBLIC FUNCTION PROTOTYPES

//extern_public_function
//  vector3_foo();

