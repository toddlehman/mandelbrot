/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"
#import "Pixel.h"
#import "Palette.h"


//-----------------------------------------------------------------------------
typedef struct
{
  Pixel    *_pixels;
  Pixel    **pixels;

  Palette  *palette;

  int      di;  // Height in pixels.
  int      dj;  // Width in pixels.

  real     x_center, y_center;
  real     x_min,    y_min;
  real     x_max,    y_max;

  real     pixel_size;
  real     periodicity_epsilon;

  uint64   iter_max;

  int      subsample_scale;
  real     subsample_tolerance;
}
Image;


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_method
Image *image_alloc(real x_center, real y_center, real x_size,
                   uint64 iter_max,
                   int pixel_width, int pixel_height,
                   int subsample_scale, real subsample_tolerance);

extern_public_method
void image_dealloc(Image **this);

extern_public_method
void image_populate(Image *this);

extern_public_method
void image_output(Image *this, bool text_format);

extern_public_method
void image_output_statistics(Image *this);
