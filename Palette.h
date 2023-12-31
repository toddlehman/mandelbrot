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

  uint32     sky_count;
  float64    *sky_locations;
  LinearRGB  *sky_colors;

  LinearRGB  undefined_color;
  LinearRGB  interior_iterated_periodic_color;
  LinearRGB  interior_iterated_aperiodic_color;
  LinearRGB  interior_uniterated_color;
  LinearRGB  dead_space_color;

  real       map_zoom_level;
}
Palette;


//-----------------------------------------------------------------------------
// INLINE METHODS

#if 0  // OBSOLETE
//-----------------------------------------------------------------------------
public_inline_method
LinearRGB palette_undefined_color(const Palette *this)
{
  return this->undefined_color;
}

//-----------------------------------------------------------------------------
public_inline_method
LinearRGB palette_interior_iterated_color(const Palette *this)
{
  return this->interior_iterated_color;
}

//-----------------------------------------------------------------------------
public_inline_method
LinearRGB palette_interior_uniterated_color(const Palette *this)
{
  return this->interior_uniterated_color;
}
#endif

//-----------------------------------------------------------------------------
public_inline_method
LinearRGB palette_dead_space_color(const Palette *this)
{
  return this->dead_space_color;
}


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_constructor
  Palette *palette_create(real map_zoom_level);

extern_public_destructor
  void palette_destroy(Palette **p_this);

extern_public_method
  LinearRGB palette_color_from_mandelbrot_result(const Palette *this,
                                                 const MandelbrotResult mr);
extern_public_method
  LinearRGB palette_color_from_sky_coefficient(const Palette *this,
                                               const float64 sky_coefficient);
