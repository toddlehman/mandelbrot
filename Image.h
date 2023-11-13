/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"
#import "Pixel.h"
#import "Palette.h"
#import "Mandelbrot.h"
#import "Camera.h"


//-----------------------------------------------------------------------------
#pragma pack(push, 8)

typedef struct
{
  Pixel    *_pixels;
  Pixel    **pixels;

  Pixel    interior_filler_pixel;

  Palette  *palette;

  Camera   *camera;

  Mandelbrot  *mandelbrot;

  int      di;   // Effective height in pixels.
  int      dj;   // Effective width in pixels.
  int      _di;  // Actual height in pixels (requested value plus 1 extra).
  int      _dj;  // Actual width in pixels (requested value plus 1 extra).

  mp_real  x_center, y_center;
#if 0  // OBSOLETE
  mp_real  x_size,   y_size;
  mp_real  x_min,    y_min;
  mp_real  x_max,    y_max;

  mp_real  pixel_size;
  mp_real  periodicity_epsilon;
#endif

  bool     supersample;
  int      supersample_interior_min_depth;
  int      supersample_interior_max_depth;
  int      supersample_exterior_min_depth;
  int      supersample_exterior_max_depth;
  float32  supersample_solidarity;
}
Image;

#pragma pack(pop)


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_constructor
  Image *image_create(mp_real camera_x,
                      mp_real camera_y,
                      mp_real camera_z,
                      mp_real camera_d,
                      real camera_ox,
                      real camera_oy,
                      real camera_oz,
                      int pixel_width, int pixel_height,
                      int supersample_interior_min_depth,
                      int supersample_interior_max_depth,
                      int supersample_exterior_min_depth,
                      int supersample_exterior_max_depth,
                      float32 supersample_solidarity,
                      uint64 iter_max);

extern_public_destructor
  void image_destroy(Image **this);

extern_public_method
  void image_populate(Image *this);

extern_public_method
  void image_output(Image *this, FILE *stream, bool text_format);

extern_public_method
  void image_output_statistics(Image *this, FILE *stream);

