/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"
#import "RGB.h"
#import "Mandelbrot.h"


//-----------------------------------------------------------------------------
// PIXEL EDGE MASKS

typedef enum
{
  PixelEdge_none    = 0,
  PixelEdge_right   = (1 << 0),
  PixelEdge_top     = (1 << 1),
  PixelEdge_left    = (1 << 2),
  PixelEdge_bottom  = (1 << 3),
}
PixelEdge;


//-----------------------------------------------------------------------------
// PIXEL STRUCTURE

#pragma pack(push, 0)

typedef struct
{
  MandelbrotResult  mr;                                  // 24 bytes
  LinearRGB         color;                               // 12 bytes
                                              // Subtotal:  36 bytes

  bool              supersample:1;                       //  1 bit
  bool              supersampled:1;                      //  1 bit
  bool              probed_top_edge:1;                   //  1 bit
  bool              probed_left_edge:1;                  //  1 bit
  uint              padding:28;                          // 28 bits
                                              // Subtotal:   4 bytes
}
Pixel;                                           // Total:  40 bytes

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

