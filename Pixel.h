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
  LinearRGB  color;                    // 12 bytes

  float32    displacement_3d;          //  4 bytes

  float32    interior_portion;         //  4 bytes

  bool       is_defined:1;             //  1 bit
  bool       is_interior_periodic:1;   //  1 bit
  bool       supersample:1;            //  1 bit
  bool       supersampled:1;           //  1 bit
  bool       probed_top_edge:1;        //  1 bit
  bool       probed_left_edge:1;       //  1 bit
  uint       padding:26;               // 26 bits
                            // Subtotal:   4 bytes
}
Pixel;                         // Total:  24 bytes

#pragma pack(pop)


//-----------------------------------------------------------------------------
// PUBLIC INLINE FUNCTIONS

#if 0  // OBSOLETE
//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_undefined(const Pixel pixel)
{
  return !pixel.is_defined;
}

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_defined(const Pixel pixel)
{
  return pixel.is_defined;
}
#endif

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_interior(const Pixel pixel)
{
  return (pixel.interior_portion == 1.0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_exterior(const Pixel pixel)
{
  return (pixel.interior_portion == 0.0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_mixed(const Pixel pixel)
{
  return (pixel.interior_portion > 0.0) && (pixel.interior_portion < 1.0);
}

//-----------------------------------------------------------------------------
public_inline_function
bool pixel_is_black(const Pixel pixel)
{
  return (pixel.color.r == 0) &&
         (pixel.color.g == 0) &&
         (pixel.color.b == 0);
}

