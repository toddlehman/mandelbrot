/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"
#import "RGB.h"
#import "Pixel.h"
#import "Mandelbrot.h"


//-----------------------------------------------------------------------------
// CYCLICAL COLOR PALETTE STRUCTURE

typedef struct
{
  uint32     exterior_count;
  float64    *exterior_locations;
  LinearRGB  *exterior_colors;
  LinearRGB  interior_iterated_color;
  LinearRGB  interior_uniterated_color;
  LinearRGB  undefined_color;
}
Palette;


//-----------------------------------------------------------------------------
// INLINE METHODS

//-----------------------------------------------------------------------------
public_inline_method
LinearRGB palette_undefined_color(Palette *this)
{
  return this->undefined_color;
}

//-----------------------------------------------------------------------------
public_inline_method
LinearRGB palette_interior_iterated_color(Palette *this)
{
  return this->interior_iterated_color;
}

//-----------------------------------------------------------------------------
public_inline_method
LinearRGB palette_interior_uniterated_color(Palette *this)
{
  return this->interior_uniterated_color;
}


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_constructor
Palette *palette_create(void);

extern_public_destructor
void palette_destroy(Palette **p_this);

extern_public_method
LinearRGB palette_color_from_mandelbrot_result(Palette *this,
                                               MandelbrotResult mr);
