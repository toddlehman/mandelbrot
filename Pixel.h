/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"
#import "RGB.h"
#import "Mandelbrot.h"


//-----------------------------------------------------------------------------
// PIXEL STRUCTURE

#pragma pack(push, 4)

typedef struct
{
  MandelbrotResult  mr;               // 24 bytes
  LinearRGB         color;            // 12 bytes

  bool              supersample:8;    //  1 byte
  bool              supersampled:8;   //  1 byte
  uint16            padding;          //  2 bytes
}
Pixel;                                // 40 bytes

#pragma pack(pop)


//-----------------------------------------------------------------------------
// PUBLIC INLINE FUNCTIONS

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_undefined(const Pixel pixel)
{
  return mandelbrot_result_is_undefined(pixel.mr);
}

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_defined(const Pixel pixel)
{
  return mandelbrot_result_is_defined(pixel.mr);
}

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_interior(const Pixel pixel)
{
  return mandelbrot_result_is_interior(pixel.mr);
}

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_exterior(const Pixel pixel)
{
  return mandelbrot_result_is_exterior(pixel.mr);
}

