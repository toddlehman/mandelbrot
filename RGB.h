/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// LINEAR RGB STRUCTURE

#pragma pack(push, 4)

typedef struct
{
  float32  r;  // Red   color component, where 0 <= r <= 1.  (4 bytes)
  float32  g;  // Green color component, where 0 <= g <= 1.  (4 bytes)
  float32  b;  // Blue  color component, where 0 <= b <= 1.  (4 bytes)
}
LinearRGB;     // (12 bytes)

#pragma pack(pop)


//-----------------------------------------------------------------------------
// DEVICE 48-BIT RGB STRUCTURE

#pragma pack(push, 2)

typedef struct
{
  uint16  r;  // Red   color component, where 0 <= r <= 65535.  (2 bytes)
  uint16  g;  // Green color component, where 0 <= g <= 65535.  (2 bytes)
  uint16  b;  // Blue  color component, where 0 <= b <= 65535.  (2 bytes)
}
DeviceRGB48;  // (6 bytes)

#pragma pack(pop)


//-----------------------------------------------------------------------------
// DEVICE 24-BIT RGB STRUCTURE

#pragma pack(push, 1)

typedef struct
{
  uint8  r;   // Red   color component, where 0 <= r <= 255.  (1 byte)
  uint8  g;   // Green color component, where 0 <= g <= 255.  (1 byte)
  uint8  b;   // Blue  color component, where 0 <= b <= 255.  (1 byte)
}
DeviceRGB24;  // (3 bytes)

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
// INTERPOLATE LINEAR RGB VALUES

public_inline_function
LinearRGB linear_rgb_swerp(const float32 t,
                           const LinearRGB color0,
                           const LinearRGB color1)
{
  return (LinearRGB)
  {
    .r = swerp(t, color0.r, color1.r),
    .g = swerp(t, color0.g, color1.g),
    .b = swerp(t, color0.b, color1.b),
  };
}


//-----------------------------------------------------------------------------
// COMPUTE AVERAGE OF 2 RGB VALUES

public_inline_function
LinearRGB linear_rgb_average2(const LinearRGB color1, const LinearRGB color2)
{
  return (LinearRGB)
  {
    .r = (color1.r + color2.r) / 2,
    .g = (color1.g + color2.g) / 2,
    .b = (color1.b + color2.b) / 2,
  };
}


//-----------------------------------------------------------------------------
// COMPUTE AVERAGE OF 4 RGB VALUES

public_inline_function
LinearRGB linear_rgb_average4(const LinearRGB color1, const LinearRGB color2,
                              const LinearRGB color3, const LinearRGB color4)
{
  return (LinearRGB)
  {
    .r = (color1.r + color2.r + color3.r + color4.r) / 4,
    .g = (color1.g + color2.g + color3.g + color4.g) / 4,
    .b = (color1.b + color2.b + color3.b + color4.b) / 4,
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
  LinearRGB linear_rgb_average2(const LinearRGB color1,
                                const LinearRGB color2);

extern_public_function
  LinearRGB linear_rgb_average4(const LinearRGB color1,
                                const LinearRGB color2,
                                const LinearRGB color3,
                                const LinearRGB color4);

extern_public_function
  float32 linear_rgb_solidarity2(const LinearRGB color1,
                                 const LinearRGB color2);

extern_public_function
  float32 linear_rgb_solidarity4(const LinearRGB color1,
                                 const LinearRGB color2,
                                 const LinearRGB color3,
                                 const LinearRGB color4);

