/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// LINEAR RGB STRUCTURE

#pragma pack(push, 4)

typedef struct
{
  float32  r;  // Red color component, where 0 <= r <= 1.
  float32  g;  // Green color component, where 0 <= g <= 1.
  float32  b;  // Blue color component, where 0 <= b <= 1.
}
LinearRGB;

#pragma pack(pop)


//-----------------------------------------------------------------------------
// DEVICE 24-BIT RGB STRUCTURE

#pragma pack(push, 1)

typedef struct
{
  uint8  r;  // Red color component, where 0 <= r <= 255.
  uint8  g;  // Green color component, where 0 <= g <= 255.
  uint8  b;  // Blue color component, where 0 <= b <= 255.
}
DeviceRGB24;

#pragma pack(pop)


//-----------------------------------------------------------------------------
// DEVICE 48-BIT RGB STRUCTURE

#pragma pack(push, 2)

typedef struct
{
  uint16  r;  // Red color component, where 0 <= r <= 65535.
  uint16  g;  // Green color component, where 0 <= g <= 65535.
  uint16  b;  // Blue color component, where 0 <= b <= 65535.
}
DeviceRGB48;

#pragma pack(pop)


//-----------------------------------------------------------------------------
// INTERPOLATE LINEAR RGB VALUES

public_inline_function
LinearRGB linear_rgb_lerp(const float32 t,
                          const LinearRGB color0,
                          const LinearRGB color1)
{
  return (LinearRGB)
  {
    .r = lerp(t, color0.r, color1.r),
    .g = lerp(t, color0.g, color1.g),
    .b = lerp(t, color0.b, color1.b),
  };
}


//-----------------------------------------------------------------------------
// CONVERT LINEAR RGB TO DEVICE RGB

#define DEVICE_GAMMA 2.2

public_inline_function
DeviceRGB24 linear_rgb_to_device_rgb24(const LinearRGB color)
{
  return (DeviceRGB24)
  {
    .r = (uint8)(pow(color.r, 1/DEVICE_GAMMA) * 255.999999),
    .g = (uint8)(pow(color.g, 1/DEVICE_GAMMA) * 255.999999),
    .b = (uint8)(pow(color.b, 1/DEVICE_GAMMA) * 255.999999),
  };
}

public_inline_function
DeviceRGB48 linear_rgb_to_device_rgb48(const LinearRGB color)
{
  return (DeviceRGB48)
  {
    .r = (uint16)(pow(color.r, 1/DEVICE_GAMMA) * 65535.999999),
    .g = (uint16)(pow(color.g, 1/DEVICE_GAMMA) * 65535.999999),
    .b = (uint16)(pow(color.b, 1/DEVICE_GAMMA) * 65535.999999),
  };
}


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_function
LinearRGB linear_rgb_average4(LinearRGB color1, LinearRGB color2,
                              LinearRGB color3, LinearRGB color4);

extern_public_function
float32 linear_rgb_diff4(LinearRGB color1, LinearRGB color2,
                         LinearRGB color3, LinearRGB color4);
