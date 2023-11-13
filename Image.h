/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"
#import "Pixel.h"
#import "Palette.h"
#import "Mandelbrot.h"
#import "Camera.h"


//-----------------------------------------------------------------------------
typedef enum
{
  ImageWorld_standard   = 0,
  ImageWorld_planar     = 1,
  ImageWorld_spherical  = 2,
}
ImageWorldType;

//-----------------------------------------------------------------------------
typedef enum
{
  ImageColor_standard        = 0,
  ImageColor_period          = 1,
  ImageColor_sample_density  = 2,
}
ImageColorType;

//-----------------------------------------------------------------------------
typedef enum
{
  ImageFormat_ppm_binary  = 0,
  ImageFormat_ppm_text    = 1,
}
ImageFormatType;

//-----------------------------------------------------------------------------
#pragma pack(push, 8)

typedef struct
{
  Pixel    *_pixels;
  Pixel    **pixels;

  Pixel    interior_filler_pixel;

  const Palette  *palette;

  const Camera   *camera;

  Mandelbrot  *mandelbrot;

  int      di;   // Effective height in pixels.
  int      dj;   // Effective width in pixels.
  int      _di;  // Actual height in pixels (requested value plus 1 extra).
  int      _dj;  // Actual width in pixels (requested value plus 1 extra).

  bool     supersample;
  int      supersample_interior_min_depth;
  int      supersample_interior_max_depth;
  int      supersample_exterior_min_depth;
  int      supersample_exterior_max_depth;
  int      supersample_max_depth;
  float32  supersample_solidarity;

  ImageWorldType   world_type;
  ImageColorType   color_type;
  ImageFormatType  format_type;
}
Image;

#pragma pack(pop)


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_constructor
  Image *image_create(int pixel_width, int pixel_height,
                      int supersample_interior_min_depth,
                      int supersample_interior_max_depth,
                      int supersample_exterior_min_depth,
                      int supersample_exterior_max_depth,
                      float32 supersample_solidarity,
                      uint64 iter_max,
                      const Palette *palette,
                      const Camera *camera,
                      ImageWorldType world_type,
                      ImageColorType color_type,
                      ImageFormatType format_type);

extern_public_destructor
  void image_destroy(Image **this);

extern_public_method
  void image_populate(Image *this);

extern_public_method
  void image_output(const Image *this, FILE *stream);

extern_public_method
  void image_output_statistics(const Image *this, FILE *stream);

