/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Image.h"
#import "Mandelbrot.h"
#import "RGB.h"

#define TRACE 0


//-----------------------------------------------------------------------------
// ALLOCATE IMAGE

public_constructor
Image *image_create(int width_pixels, int height_pixels,
                    int supersample_interior_min_depth,
                    int supersample_interior_max_depth,
                    int supersample_exterior_min_depth,
                    int supersample_exterior_max_depth,
                    float32 supersample_solidarity,
                    uint64 iter_max,
                    const Palette *palette,
                    const Camera *camera)
{
  // FIXME:  Move these two assertions to initialization code for each of these
  // modules.
  assert(sizeof(MandelbrotResult) == 24);  // To catch unexpected changes.
  assert(sizeof(Pixel) == 20);  // To catch unexpected changes.

  assert(width_pixels >= 1);
  assert(height_pixels >= 1);
  assert(supersample_interior_min_depth >= 0);
  assert(supersample_interior_max_depth >= 0);
  assert(supersample_interior_min_depth <= supersample_interior_max_depth);
  assert(supersample_exterior_min_depth >= 0);
  assert(supersample_exterior_max_depth >= 0);
  assert(supersample_exterior_min_depth <= supersample_exterior_max_depth);
  assert(supersample_solidarity >= 0);
  assert(supersample_solidarity <= 1);
  assert(iter_max > 0);
  assert(palette != NULL);
  assert(camera != NULL);


  Image *this = mem_alloc_clear(1, sizeof(Image));

  this->camera = camera;
  this->palette = palette;


  // --- Calculate floating-point precision required.

#if 0  // FIXME!!!!
  real bits;
  {
    mp_real log2_xy_min_size;
    mp_init2(log2_xy_min_size, 64);  // (Overkill; only need about 8 or 10)
    mp_log2(log2_xy_min_size, xy_min_size);  // Typically negative.
    bits = -mp_get_d(log2_xy_min_size);
    mp_clear(log2_xy_min_size);
    if (bits < 0) bits = 0;
  }
  // bits now contains the number of bits to the right of the radix point of
  // xy_min_size.  Next, add enough more bits to support the largest scalar
  // value possible during Mandelbrot Set iteration.
  bits += log2(mandelbrot_max_scalar_value_during_iteration());
  // Now add enough more bits to support the resolution of a single pixel.
  bits += log2(MIN(width_pixels, height_pixels));
  // Now add enough more bits to support supersampling.
  bits += MAX(supersample_interior_max_depth, supersample_exterior_max_depth);
  // Now add a few guard bits.  12 is a bare minimum which just barely is
  // sufficient.  16 is a better value.
  bits += 16;
  // Finally, round up and use this as the required precision.
  int mp_prec = ceil(bits);
#endif
  int mp_prec = 68;


  // --- Store image dimensions.

  this->di = height_pixels;
  this->dj = width_pixels;
  this->_di = this->di + 1;
  this->_dj = this->dj + 1;


  // --- Create Mandelbrot Set iteration object.

  this->mandelbrot = mandelbrot_create(iter_max, mp_prec);


  // --- Store supersampling parameters.

  this->supersample_interior_min_depth = supersample_interior_min_depth;
  this->supersample_interior_max_depth = supersample_interior_max_depth;
  this->supersample_exterior_min_depth = supersample_exterior_min_depth;
  this->supersample_exterior_max_depth = supersample_exterior_max_depth;
  this->supersample_max_depth = MAX(this->supersample_interior_max_depth,
                                    this->supersample_exterior_max_depth);
  this->supersample_solidarity = supersample_solidarity;
  this->supersample = (this->supersample_max_depth > 0)
                   && (this->supersample_solidarity > 0);



#if 0  // FIXME -- Do this elsewhere, probably in the population method.
  // --- Tweak image position to account for non-center sampling of pixels.
  //     A supersampling depth of zero results in a tweak of 1/2 pixel.
  //     A depth of 1 results in a tweak of 1/4 pixel.  A depth of 2 results
  //     in a tweak of 1/8 pixel.  In general, a depth of $n$ results in a
  //     tweak of $1/2^{1+n}$ pixel.  This really works beautifully and results
  //     in pixel-perfect symmetry around the x-axis -- for both even and odd
  //     image heights.
  {
    mp_real xy_tweak;
    mp_init2(xy_tweak, mp_prec);
    mp_div_d(xy_tweak, this->pixel_size, pow(2, 1 + supersample_max_depth));

    mp_add(this->target_x, this->target_x, xy_tweak);
    mp_sub(this->target_y, this->target_y, xy_tweak);
    mp_add(this->x_min, this->x_min, xy_tweak);
    mp_sub(this->x_max, this->x_max, xy_tweak);
    mp_add(this->y_min, this->y_min, xy_tweak);
    mp_sub(this->y_max, this->y_max, xy_tweak);

    mp_clear(xy_tweak);
  }
#endif


  // --- Allocate pixel array.

  this->_pixels = mem_alloc_clear(this->_di * this->_dj, sizeof(Pixel));
  this->pixels = mem_alloc_clear(this->_di, sizeof(Pixel **));
  for (int i = 0; i < this->_di; i++)
    this->pixels[i] = this->_pixels + (i * this->_dj);

  for (int i = 0; i < this->_di; i++)
  for (int j = 0; j < this->_dj; j++)
  {
    this->pixels[i][j].is_defined = 0;
  }


  // --- Allocate palette and intialize related resources.

  this->palette = palette_create();

  this->interior_filler_pixel = (Pixel)
  {
    .color                 = palette_color_from_mandelbrot_result(
                               this->palette,
                               mandelbrot_result_interior_uniterated()),
    .interior_portion      = 1.0,
    .is_defined            = true,
    .is_interior_periodic  = false,
    .supersample           = false,
    .supersampled          = false,
    .probed_top_edge       = false,
    .probed_left_edge      = false,
  };

  return this;
}


//-----------------------------------------------------------------------------
// DEALLOCATE IMAGE

public_destructor
void image_destroy(Image **p_this)
{
  assert(p_this);
  Image *this = *p_this;

  assert(this);

  assert(this->pixels);
  mem_dealloc(&this->pixels);

  assert(this->_pixels);
  mem_dealloc(&this->_pixels);

  assert(this->mandelbrot);
  mandelbrot_destroy(&this->mandelbrot);

  mem_dealloc(p_this);
}


//-----------------------------------------------------------------------------
// DETERMINE WHETHER PIXEL NEEDS SUPERSAMPLING

private_inline_method
bool image_pixel_needs_supersampling(const Image *this,
                                     const Pixel pixel,
                                     const float32 solidarity,
                                     const int ss_depth)
{
  // Decide based on the requirement for solidarity, the nature of the pixel,
  // the minimum and maximum supersampling depths, and the current
  // supersampling depth.
  if (pixel_is_interior(pixel))
  {
    return (ss_depth < this->supersample_interior_min_depth) ||
           ((ss_depth < this->supersample_interior_max_depth) &&
            (solidarity < this->supersample_solidarity));
  }
  else
  {
    return (ss_depth < this->supersample_exterior_min_depth) ||
           ((ss_depth < this->supersample_exterior_max_depth) &&
            (solidarity < this->supersample_solidarity));
  }
}


//-----------------------------------------------------------------------------
// COMPUTE NON-SUPERSAMPLED PIXEL

private_method
Pixel image_compute_pixel(const Image *this, real i, real j)
{
  assert(this);

  Pixel pixel;
  MandelbrotResult mr;

  mp_real x, y;
  mp_init2(x, this->mandelbrot->conf.mp_prec);
  mp_init2(y, this->mandelbrot->conf.mp_prec);

  // Map image coordinates to viewport coordinates.  Attempt to set (du,dv)
  // to (2,2) so that (u,v) covers the unit square, but take non-square image
  // into account --- which means du and dv will not always be the same.
  int dij_max = (real)MAX(this->di, this->dj);
  real du = (this->dj / (real)dij_max) * 2.0;
  real dv = (this->di / (real)dij_max) * 2.0;
  real u = lerp(j / this->dj, -du/2, +du/2);
  real v = lerp(i / this->di, +dv/2, -dv/2);

  // Compute size of smallest possible subpixel in viewport coordinate system.
  real uve = (2.0 / (real)dij_max) / pow(2, this->supersample_max_depth);

  // Map viewport coordinates to Argand plane coordinates.
  if (camera_get_argand_point(this->camera, u, v, &x, &y))
  {
    // Determine subpixel size around point (u,v) in the viewport.
    mp_real e, e2, xe, ye;
    mp_init2(e, this->mandelbrot->conf.mp_prec);
    mp_init2(e2, this->mandelbrot->conf.mp_prec);
    mp_init2(xe, this->mandelbrot->conf.mp_prec);
    mp_init2(ye, this->mandelbrot->conf.mp_prec);

    mp_set_d(e, 1e-6);  // Maximum possible epsilon.
    if (camera_get_argand_point(this->camera, u - uve, v, &xe, &ye))
    {
      mp_sub(xe, xe, x); mp_sqr(xe, xe);  // xe = (xe - x) ** 2;
      mp_sub(ye, ye, y); mp_sqr(ye, ye);  // ye = (ye - y) ** 2;
      mp_add(e2, xe, ye);                 // e2 = xe + ye;
      if (mp_less_p(e2, e))               // if (e2 < e)
        mp_set(e, e2);                    //   e = e2;
    }
    if (camera_get_argand_point(this->camera, u + uve, v, &xe, &ye))
    {
      mp_sub(xe, xe, x); mp_sqr(xe, xe);  // xe = (xe - x) ** 2;
      mp_sub(ye, ye, y); mp_sqr(ye, ye);  // ye = (ye - y) ** 2;
      mp_add(e2, xe, ye);                 // e2 = xe + ye;
      if (mp_less_p(e2, e))               // if (e2 < e)
        mp_set(e, e2);                    //   e = e2;
    }
    if (camera_get_argand_point(this->camera, u, v - uve, &xe, &ye))
    {
      mp_sub(xe, xe, x); mp_sqr(xe, xe);  // xe = (xe - x) ** 2;
      mp_sub(ye, ye, y); mp_sqr(ye, ye);  // ye = (ye - y) ** 2;
      mp_add(e2, xe, ye);                 // e2 = xe + ye;
      if (mp_less_p(e2, e))               // if (e2 < e)
        mp_set(e, e2);                    //   e = e2;
    }
    if (camera_get_argand_point(this->camera, u, v + uve, &xe, &ye))
    {
      mp_sub(xe, xe, x); mp_sqr(xe, xe);  // xe = (xe - x) ** 2;
      mp_sub(ye, ye, y); mp_sqr(ye, ye);  // ye = (ye - y) ** 2;
      mp_add(e2, xe, ye);                 // e2 = xe + ye;
      if (mp_less_p(e2, e))               // if (e2 < e)
        mp_set(e, e2);                    //   e = e2;
    }
    mp_sqrt(e, e);

    // Now compute pixel.
    mr = mandelbrot_compute(this->mandelbrot, x, y, e);
    pixel.color = palette_color_from_mandelbrot_result(this->palette, mr);
    pixel.interior_portion = mandelbrot_result_is_interior(mr)? 1.0 : 0.0;
    pixel.is_defined = true;
    pixel.is_interior_periodic = mandelbrot_result_is_interior_periodic(mr);

    mp_clear(ye);
    mp_clear(xe);
    mp_clear(e2);
    mp_clear(e);
  }
  else
  {
    mr = mandelbrot_result_exterior_uniterated(0);
    pixel.color = palette_dead_space_color(this->palette);
    pixel.interior_portion = 0;
    pixel.is_defined = true;
    pixel.is_interior_periodic = false;
  }

  pixel.supersample = false;
  pixel.supersampled = false;
  pixel.probed_top_edge = false;
  pixel.probed_left_edge = false;

  mp_clear(y);
  mp_clear(x);

  // Output progress (for debugging).
  #if TRACE
    fprintf(stderr, "i=%7.2f  j=%7.2f  iter=%9llu\n",
                    (double)i, (double)j, mr.iter);

    if ((i == floor(i)) && (j == floor(j)))
      fprintf(stderr, "pixel(i=%d, j=%d) -> %llu (%.3f,%.3f,%.3f)\n",
              (int)i, (int)j, mr.iter,
              pixel.color.r, pixel.color.g, pixel.color.b);
    else
      fprintf(stderr, "pixel(i=%f, j=%f) -> %llu (%.3f,%.3f,%.3f)\n",
              (double)i, (double)j, mr.iter,
              pixel.color.r, pixel.color.g, pixel.color.b);
  #endif

  return pixel;
}


//-----------------------------------------------------------------------------
// CALCULATE SUPERSAMPLED PIXEL

private_method
Pixel image_compute_supersampled_pixel(const Image *this,
                                       real i, real j, real di, real dj,
                                       Pixel pixel_00,
                                       int ss_depth)
{
  assert(this);
  assert(i >= 0); assert(i + di <= this->di);
  assert(j >= 0); assert(j + dj <= this->dj);
  assert(pixel_00.is_defined);

  // Calculate intermediate subpixel values.
  Pixel pixel_01 = image_compute_pixel(this, i,        j+(dj/2));
  Pixel pixel_10 = image_compute_pixel(this, i+(di/2), j);
  Pixel pixel_11 = image_compute_pixel(this, i+(di/2), j+(dj/2));

  // Adjust required solidarity based on supersampling depth.  The deeper we go
  // in supersampling, the sloppier we allow the results to be.  This is
  // somewhat of a KLUDGE, but it seems to work okay and it allows the maximum
  // supersample depth to be specified as infinite because it short-circuits
  // that limit (as long as the value is not 1 to start with).
  float32 required_solidarity = this->supersample_solidarity;
  required_solidarity = 1 - ((1 - required_solidarity) * pow(2, ss_depth));

  // Assess color variation.
  float32 solidarity = linear_rgb_solidarity4(pixel_00.color, pixel_01.color,
                                              pixel_10.color, pixel_11.color);

  // If color variation exceeds tolerance, then subdivide each subquadrant of
  // the pixel recursively (unless the supersample depth limit has already been
  // reached).  Note that interior and exterior pixels have separate minimum
  // and maximum depth limits.
  if (image_pixel_needs_supersampling(this, pixel_00, solidarity, ss_depth))
    pixel_00 = image_compute_supersampled_pixel(this, i,j, di/2,dj/2,
                                                pixel_00, ss_depth+1);

  if (image_pixel_needs_supersampling(this, pixel_01, solidarity, ss_depth))
    pixel_01 = image_compute_supersampled_pixel(this, i,j+dj/2, di/2,dj/2,
                                                pixel_01, ss_depth+1);

  if (image_pixel_needs_supersampling(this, pixel_10, solidarity, ss_depth))
    pixel_10 = image_compute_supersampled_pixel(this, i+di/2,j, di/2,dj/2,
                                                pixel_10, ss_depth+1);

  if (image_pixel_needs_supersampling(this, pixel_11, solidarity, ss_depth))
    pixel_11 = image_compute_supersampled_pixel(this, i+di/2,j+dj/2, di/2,dj/2,
                                                pixel_11, ss_depth+1);

  // Now return a blend of the four supersamples.
  Pixel pixel;
  pixel.color = linear_rgb_average4(pixel_00.color, pixel_01.color,
                                    pixel_10.color, pixel_11.color);
  pixel.interior_portion = (pixel_00.interior_portion +
                            pixel_01.interior_portion +
                            pixel_10.interior_portion +
                            pixel_11.interior_portion) / 4;
  #if TRACE
  fprintf(stderr, "i=%f j=%f di=%f dj=%f ss_depth=%d -> interior=%f\n",
          (double)i, (double)j, (double)di, (double)dj, ss_depth,
          (double)pixel.interior_portion);
  #endif

  #if 0
  if (ss_depth > 1)  // KLUDGE test.
  {
    pixel.color = linear_rgb_lerp(pow(solidarity, 0.5),
                                  palette_undefined_color(this->palette),
                                  pixel.color);
  }
  #elif 0
  if (ss_depth >= 5)  // Another KLUDGE test.
  {
    if (pixel_is_exterior(pixel) && (solidarity < 1))
      pixel.color = linear_rgb_lerp(solidarity / pow(2, ss_depth - 2),
                                    (LinearRGB) { 1, 1, 1 },
                                    pixel.color);
  }
  #endif

  pixel.is_defined = true;
  pixel.is_interior_periodic = pixel_00.is_interior_periodic &&
                               pixel_01.is_interior_periodic &&
                               pixel_10.is_interior_periodic &&
                               pixel_11.is_interior_periodic;
  pixel.supersampled = true;   // Supersampling accomplished for this pixel.
  pixel.supersample  = false;  // Supersampling not required for this pixel.
  pixel.probed_top_edge = false;
  pixel.probed_left_edge = false;

  return pixel;
}


//-----------------------------------------------------------------------------
// POPULATE PIXEL

private_method
void image_populate_pixel(Image *this, int i, int j,
                          bool probe_top_edge, bool probe_left_edge)
{
  assert(this);
  assert(i >= 0); assert(i < this->_di);
  assert(j >= 0); assert(j < this->_dj);

  Pixel *pixel = &this->pixels[i][j];

  if (!pixel->is_defined)
  {
    *pixel = image_compute_pixel(this, i, j);
    assert(pixel->is_defined);
  }

  // BUG!  This misses interior uniterated pixels as in those detected by
  // the cardiod/disc early-out test.
  // if (pixel->is_interior_periodic && this->supersample)

  if (pixel_is_interior(*pixel) && this->supersample)
  {
    bool all_interior = true;

    if (probe_top_edge && !pixel->probed_top_edge)
    {
      int dsj = (1 << this->supersample_interior_max_depth);
      for (int sj = 1; all_interior && (sj < dsj); sj++)  // (Not sj = 0)
      {
        Pixel subpixel = image_compute_pixel(this, i, j+(real)sj/(real)dsj);
        //all_interior &= subpixel.is_interior_periodic;
        all_interior &= pixel_is_interior(subpixel);
      }
      pixel->probed_top_edge = true;
    }

    if (probe_left_edge && !pixel->probed_left_edge)
    {
      int dsi = (1 << this->supersample_interior_max_depth);
      for (int si = 1; all_interior && (si < dsi); si++)  // (Not si = 0)
      {
        Pixel subpixel = image_compute_pixel(this, i+(real)si/(real)dsi, j);
        //all_interior &= subpixel.is_interior_periodic;
        all_interior &= pixel_is_interior(subpixel);
      }
      pixel->probed_left_edge = true;
    }

    if (!all_interior)
    {
      pixel->supersample = true;  // Mark for future refinement.
      pixel->supersampled = false;
      pixel->interior_portion = 0.5;  // Signify neither exterior nor interior.
    }
  }
}


//-----------------------------------------------------------------------------
// REFINE PIXEL

private_method
void image_refine_pixel(Image *this, int i, int j)
{
  assert(this);
  assert(this->supersample);
  assert(i >= 0); assert(i < this->di);  // *Not* this->_di here.
  assert(j >= 0); assert(j < this->dj);  // *Not* this->_dj here.

  // Set new pixel value.
  Pixel old_pixel = this->pixels[i][j];
  Pixel new_pixel = image_compute_supersampled_pixel(this, i, j, 1.0, 1.0,
                                                     old_pixel, 1);
  this->pixels[i][j] = new_pixel;

  // If the color change was significant, then mark the pixel's eight neighbors
  // as also now needing refinement.  (TODO: How much more quality does this
  // provide?  What is the cost?  A few quick experiments showed a slight
  // increase in quality in exterior regions for a slight increase in cost.)
  float32 solidarity = linear_rgb_solidarity2(old_pixel.color, new_pixel.color);
  if (solidarity < this->supersample_solidarity)
  {
    //#pragma message("Neighbor supersampling trigger temporarily disabled.")
    //#if 0  // TEMPORARILY DISABLED
    int ni0 = MAX(i-1, 0), ni1 = MIN(i+1, this->di-1);
    int nj0 = MAX(j-1, 0), nj1 = MIN(j+1, this->dj-1);
    for (int ni = ni0; ni <= ni1; ni++)
    for (int nj = nj0; nj <= nj1; nj++)
      this->pixels[ni][nj].supersample = true;
    //#endif
  }
}


//-----------------------------------------------------------------------------
// POPULATE BLOCK (RECURSIVE)

private_method
void image_populate_block(Image *this, int i0, int j0, int di, int dj)
{
  #if TRACE
  fprintf(stderr, "block(i0=%d, j0=%d, di=%d, dj=%d)\n", i0, j0, di, dj);
  #endif

  assert(this);
  assert(i0 >= 0); assert(i0 + di <= this->_di);
  assert(j0 >= 0); assert(j0 + dj <= this->_dj);


  // --- Populate the boundary of the block.  Because the right edge of this
  //     block will share pixels with the left edge of the block to its right,
  //     and the bottom edge of this block will share pixels with the top edge
  //     of the block below it, only the left and top edges of the pixels at
  //     the right and bottom boundaries of the block need be calculated in the
  //     case of supersampling (unless later they require it based on other
  //     needs).  Also note that neither the top nor the left edge of the
  //     pixel in the lower-right boundary of the block need be scanned,
  //     because only its upper-left single point is needed.  (TODO:  Explain
  //     this better someday with a picture.)

  int i1 = i0 + di - 1;
  int j1 = j0 + dj - 1;

  bool fill_block = true;

  for (int i = i0; fill_block && (i <= i1); i++)  // Proceed top to bottom.
  {
    for (int j = j0; fill_block && (j <= j1); j++)  // Proceed left to right.
    {
      if ((i > i0) && (i < i1) && (j > j0) && (j < j1))
        { j = j1 - 1; continue; }  // Efficiently skip non-boundary pixels.

      #if TRACE
      fprintf(stderr, "i=%d, j=%d\n", i, j);
      #endif

      bool probe_top_edge  = ((i == i0) || (i == i1)) && (j != j1);
      bool probe_left_edge = ((j == j0) || (j == j1)) && (i != i1);

      image_populate_pixel(this, i, j, probe_top_edge, probe_left_edge);

      // BUG!  This misses interior uniterated pixels as in those detected by
      // the cardiod/disc early-out test.
      //fill_block &= this->pixels[i][j].is_interior_periodic;

      fill_block &= pixel_is_interior(this->pixels[i][j]);
      #if TRACE
      if (!fill_block)
        { fprintf(stderr, "Pixel is not interior; exiting boundary scan\n"); }
      #endif
    }
  }
  #if TRACE
  fprintf(stderr, "done with boundary\n");
  #endif


  // --- Fill interior of block or subdivide.

  if (fill_block)
  {
    // If the boundary of the block consists entirely of M-Set interior pixels,
    // then fill the entire block with interior pixels.  Note that if di <= 2
    // or dj <= 2, there is actually no interior and the loop below is never
    // entered.

    for (int i = i0 + 1; i <= i1 - 1; i++)
    for (int j = j0 + 1; j <= j1 - 1; j++)
    {
      assert(!this->pixels[i][j].is_defined);
      this->pixels[i][j] = this->interior_filler_pixel;
    }
  }
  else if ((di >= 3) || (dj >= 3))
  {
    // Otherwise, split the block into two roughly equal parts (either
    // vertically or horizontally, depending on the aspect ratio), and recurse
    // to populate each part.  Note that the height and width of the block are
    // not necessarily powers of 2 and that this still works perfectly fine.

    if (di >= dj)
    {
      // This is a tall block; split vertically.
      int half_di = di / 2;
      image_populate_block(this, i0,           j0,      half_di + 1, dj);
      image_populate_block(this, i0 + half_di, j0, di - half_di,     dj);
    }
    else
    {
      // This is a wide block; split horizontally.
      int half_dj = dj / 2;
      image_populate_block(this, i0, j0,           di,      half_dj + 1);
      image_populate_block(this, i0, j0 + half_dj, di, dj - half_dj    );
    }
  }
  else if ((di >= 2) || (dj >= 2))
  {
    // Handle a 1x2, 2x1, or 2x2 block.
    for (int i = i0; i <= i1; i++)
    for (int j = j0; j <= j1; j++)
    {
      image_populate_pixel(this, i, j, false, false);
    }
  }
  else if ((di == 1) && (dj == 1))
  {
    // Nothing to do; already done with this 1x1 block.
  }
  else
  {
    assert(false);
  }
}


//-----------------------------------------------------------------------------
// POPULATE IMAGE

public_method
void image_populate(Image *this)
{
  assert(this);

  // Make a first pass to populate every pixel in the image, including the
  // extra row and column at the bottom and right edges if supersampling is
  // enabled.
  image_populate_block(this,
                       0,
                       0,
                       this->supersample? this->_di : this->di,
                       this->supersample? this->_dj : this->dj);

  // Return now if supersampling is not requested.
  if (!this->supersample)
    return;

  // Mark interior pixels for supersampling which have an exterior neighbor
  // pixel in any direction, and vice-versa.  (It suffices simply to look for
  // exterior pixels neighboring interior pixels if both are marked.)
  for (int i = 0; i < this->_di; i++)
  for (int j = 0; j < this->_dj; j++)
  {
    if (pixel_is_interior(this->pixels[i][j]))
    {
      int ni0 = MAX(0, i-1), ni1 = MIN(i+1, this->_di-1);
      int nj0 = MAX(0, j-1), nj1 = MIN(j+1, this->_dj-1);

      for (int ni = ni0; ni <= ni1; ni++)
      for (int nj = nj0; nj <= nj1; nj++)
      {
        if (pixel_is_exterior(this->pixels[ni][nj]))
        {
          this->pixels[ni][nj].supersample = true;
          this->pixels[i][j].supersample = true;
        }
      }
    } 
  }

  // Now go back and mark pixels for supersampling whose rightward and downward
  // neighbors are significantly different in color.  (The "+1"s in the indexes
  // here are safe because the image is padded by one row and column at the
  // bottom and right edges.)
  // FIXME:  This doesn't catch all cases!!!
  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    if (!this->pixels[i][j].supersample)
    {
      float32 solidarity = linear_rgb_solidarity4(
        this->pixels[i+0][j+0].color,
        this->pixels[i+0][j+1].color,
        this->pixels[i+1][j+0].color,
        this->pixels[i+1][j+1].color
      );

      if (solidarity < this->supersample_solidarity)
      {
        this->pixels[i][j].supersample = true;
      }
    }
  }

  // Now, refine pixels until the entire image is refined.  Note that this is
  // potentially a lot of looping and relooping, but that the looping is
  // insignificant compared to the actual pixel calculation.
  bool supersampled;
  do
  {
    supersampled = false;
    for (int i = 0; i < this->di; i++)
    for (int j = 0; j < this->dj; j++)
    {
      if (this->pixels[i][j].supersample)
      {
        if (this->pixels[i][j].supersampled)
        {
          this->pixels[i][j].supersample = false;  // Ignore; already done.
        }
        else
        {
          image_refine_pixel(this, i, j);
          supersampled = true;
        }
      }
    }
  }
  while (supersampled);

  // Finally, for testing, set to gray any pixels that did not require
  // refinement.
#if 0
  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    if (!this->pixels[i][j].supersampled)
    {
      this->pixels[i][j].color.r = 0.5;
      this->pixels[i][j].color.g = 0.5;
      this->pixels[i][j].color.b = 0.5;
    }
  }
#endif
}


//-----------------------------------------------------------------------------
// OUTPUT IMAGE

public_method
void image_output(const Image *this, FILE *stream, bool text_format)
{
  assert(this);

  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    assert(this->pixels[i][j].is_defined);
  }

  fprintf(stream, "%s\n", text_format? "P3" : "P6");
  fprintf(stream, "%d %d\n", this->dj, this->di);
  fprintf(stream, "%d\n", 255);

  int buffer_byte_count = (text_format? strlen("000 000 000 ") : 3) *
                          this->dj;
  byte *buffer = mem_alloc_clear(buffer_byte_count, 1);

  for (int i = 0; i < this->di; i++)
  {
    byte *p = &buffer[0];

    for (int j = 0; j < this->dj; j++)
    {
      #if 0  // FOO
      LinearRGB c = this->pixels[i][j].color;
      for (int n = 0; n < 1; n++)
      {
        c.r = swerp(c.r, 0, 1);
        c.g = swerp(c.g, 0, 1);
        c.b = swerp(c.b, 0, 1);
      }
      this->pixels[i][j].color = c;
      #endif

      DeviceRGB24 color = linear_rgb_to_device_rgb24(this->pixels[i][j].color);

      if (text_format)
      {
        // This is WAY faster than calling sprintf() zillions of times.
        *p++ = '0' + ((color.r / 100) % 10);
        *p++ = '0' + ((color.r /  10) % 10);
        *p++ = '0' + ((color.r /   1) % 10);
        *p++ = ' ';
        *p++ = '0' + ((color.g / 100) % 10);
        *p++ = '0' + ((color.g /  10) % 10);
        *p++ = '0' + ((color.g /   1) % 10);
        *p++ = ' ';
        *p++ = '0' + ((color.b / 100) % 10);
        *p++ = '0' + ((color.b /  10) % 10);
        *p++ = '0' + ((color.b /   1) % 10);
        *p++ = (j < this->dj - 1)? ' ' : '\n';
      }
      else
      {
        *p++ = color.r;
        *p++ = color.g;
        *p++ = color.b;
      }
    }

    (void)fwrite(buffer, buffer_byte_count, 1, stream);
  }

  free(buffer);

  fflush(stream);
}


//-----------------------------------------------------------------------------
// PRINT CLEAN FLOATING-POINT NUMBER

#if 0  // OBSOLETE
private_method
void image_print_clean_real(real value)
{
  char buf[100];
  sprintf(buf, "%.16f", value);
  char *p = buf + strlen(buf) - 1;
  while (*p == '0')
    *p-- = '\0';
  if (*p == '.')
    *p-- = '\0';
  fprintf("%s", buf);
}
#endif


//-----------------------------------------------------------------------------
// COMPUTE RATIO

private_method
float64 ratio(uint64 numerator, uint64 denominator)
{
  return (float64)numerator / (float64)(denominator? denominator : 1);
}


//-----------------------------------------------------------------------------
// COMPUTE PERCENTAGE

private_function
float64 percentage(uint64 numerator, uint64 denominator)
{
  return ratio(numerator, denominator) * 100;
}



//-----------------------------------------------------------------------------
// SPECIALLY FORMAT HIGH-PRECISION REAL NUMBER FOR READABILITY
//
// This uses just the right number of digits and avoids the ugly exponential
// format.
//
// FIXME:  This isn't quite right yet.  For example, the number 0.001 comes out
// as "0.00099945068" when 39-bit precision is in use.  The MPFR exponential
// format fares no better, producing "9.994506835938E-04".  I suppose the only
// possible correct thing to do is to store the original strings passed by the
// command line.

private_function
char *special_format_mp_real(const mp_real x, const bool sign)
{
  // Calculate bits of precision to the right of the decimal point by
  // subtracting the bits of precision to the left of the decimal point from
  // the total bits of precision.
  double x_trunc = floor(abs(mp_get_d(x)));
  int bits_left = (x_trunc >= 1)? (int)ceil(log2(x_trunc)) : 0;
  int bits_right = mp_get_prec(x) - bits_left;
  int digits_right = floor((double)bits_right / log2(10));

  // Configure mp_printf() format.
  char format[16];
  snprintf(format, ELEMENT_COUNT(format),
           sign? "%%+.%dRNF" : "%%.%dRNF",
           digits_right);

  // Format high-precision number.
  static char str[10000];
  mp_snprintf(str, ELEMENT_COUNT(str), format, x);

  #if 0
  // FIXME:  Overwrite the above with the correct-but-ugly exponential format.
  mp_snprintf(str, ELEMENT_COUNT(str),
              sign? "%+RNE" : "%RNE",
              x);
  #endif

  return str;  // Not thread-safe; be careful!
}


//-----------------------------------------------------------------------------
// SPECIALLY FORMAT LARGE NUMBER FOR READABILITY

private_function
char *special_format_number(const uint64 n)
{
  static char str[100];

  const char *units[] =
  {
    /* 0 */  "ones",
    /* 1 */  "thousand",
    /* 2 */  "million",
    /* 3 */  "billion",
    /* 4 */  "trillion",
    /* 5 */  "quadrillion",
    /* 6 */  "quintrillion",
    /* 7 */  "sextillion",
    /* 8 */  "septillion",
    /* 9 */  "octillion",
  };

  float64 f = n;
  int unit_index;
  for (unit_index = 0;
       f >= 1000;
       unit_index++, f /= 1000)
    ;

  assert(unit_index < ELEMENT_COUNT(units));

  snprintf(str, ELEMENT_COUNT(str),
           f >= 100? "%.0f %s":
           f >= 10?  "%.1f %s":
           f >= 1?   "%.2f %s":
                     "%.0f %s",  // This case will only occur when f == 0.
           f,
           units[unit_index]);

  return str;  // Not thread-safe; be careful!
}


//-----------------------------------------------------------------------------
// FORMAT FRACTIONAL PIXEL COUNT FOR READABILITY

private_method
char *special_format_fractional_pixel_count(const Image *this,
                                            const float64 value)
{
  assert(this);

  int max_supersample_depth = this->supersample
                                ?  MAX(this->supersample_interior_max_depth,
                                    this->supersample_exterior_max_depth)
                                : 0;

  // Each successive supersampling depth results in one additional power of
  // 4^-1, e.g., two extra digits base 10, as the fractions progress
  // 0., 0.25, 0.0625, 0.015625, etc.
  char format[16];
  snprintf(format, ELEMENT_COUNT(format), "%%" "%d" "." "%d" "f",
           20 + (2*max_supersample_depth) - (max_supersample_depth == 0? 1:0),
           2 * max_supersample_depth);

  static char str[100];
  snprintf(str, ELEMENT_COUNT(str), format, value);

  return str;  // Not thread-safe; be careful!
}
 

//-----------------------------------------------------------------------------
// OUTPUT STATISTICS

public_method
void image_output_statistics(const Image *this, FILE *stream)
{
  assert(this);

  MandelbrotConfiguration *conf = &(this->mandelbrot->conf);
  MandelbrotResultStatistics *stats = &(this->mandelbrot->stats);

  const char *hbar_thick = "========================================"
                           "=======================================";

  real cpu_time = (real)clock() / (real)CLOCKS_PER_SEC;
  if (cpu_time == 0) cpu_time = 0.001;


  // --- Configuration parameters

  fprintf(stream, "%s\n", hbar_thick);

  fprintf(stream,
          "         MP precision required:  %d bits\n",
          conf->mp_prec);

  mp_fprintf(stream,
             "                      Center x:  %s\n",
             special_format_mp_real(this->camera->target_x, true));

  mp_fprintf(stream,
             "                      Center y:  %s\n",
             special_format_mp_real(this->camera->target_y, true));

#if 0  // OBSOLETE -- FIXME
  mp_fprintf(stream,
             "                        Size x:  %s\n",
             special_format_mp_real(this->x_size, false));

  mp_fprintf(stream,
             "                        Size y:  %s\n",
             special_format_mp_real(this->y_size, false));
#endif

  fprintf(stream,
          "            Maximum iterations:  %llu\n",
          conf->iter_max);

  fprintf(stream,
          "                    Pixel size:  %d x %d\n",
          this->dj, this->di);

  fprintf(stream,
          "  Supersampling interior depth:  %d to %d\n",
          this->supersample_interior_min_depth,
          this->supersample_interior_max_depth);

  fprintf(stream,
          "  Supersampling exterior depth:  %d to %d\n",
          this->supersample_exterior_min_depth,
          this->supersample_exterior_max_depth);

  fprintf(stream,
          "      Supersampling solidarity:  %f\n",
          (double)this->supersample_solidarity);

  fprintf(stream, "\n");


  // --- CPU time

  fprintf(stream, "%s\n", hbar_thick);

  fprintf(stream,
          "                      CPU time: %23.3f seconds\n",
          (double)cpu_time);

  fprintf(stream, "\n");


  // --- Pixels

  uint64 total_pixels_including_overscan =
    this->supersample? (uint64)this->_di * (uint64)this->_dj:
                       (uint64)this->di  * (uint64)this->dj;

  uint64 total_pixels =
                       (uint64)this->di  * (uint64)this->dj;

  float64 total_interior_pixels = 0;
  float64 total_exterior_pixels = 0;
  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    total_interior_pixels +=      this->pixels[i][j].interior_portion;
    total_exterior_pixels += (1 - this->pixels[i][j].interior_portion);
  }

  if (this->supersample)
  fprintf(stream,
          "        Total pixels processed: %19.0f (%s)\n",
          (double)total_pixels_including_overscan,
          special_format_number(total_pixels_including_overscan));

  fprintf(stream,
          "         Total pixels in image: %19.0f (%s)\n",
          (double)total_pixels,
          special_format_number(total_pixels));

  fprintf(stream,
          "         Total interior pixels: %s (%.1f%% of image)\n",
          special_format_fractional_pixel_count(this, total_interior_pixels),
          percentage(total_interior_pixels, total_pixels));

  fprintf(stream,
          "         Total exterior pixels: %s (%.1f%% of image)\n",
          special_format_fractional_pixel_count(this, total_exterior_pixels),
          percentage(total_exterior_pixels, total_pixels));

  fprintf(stream,
          "             Pixels per second: %19.0f (%s)\n",
          (double)total_pixels / (double)cpu_time,
          special_format_number((uint64)((real)total_pixels/(real)cpu_time)));

  fprintf(stream, "\n");


  // --- Probes

  fprintf(stream,
          "                  Total probes: %19llu (%s)\n",
          stats->total_probes,
          special_format_number(stats->total_probes));

  fprintf(stream,
          "         Total interior probes: %19llu (%.1f%% of probes)\n",
          stats->interior_probes,
          percentage(stats->interior_probes, stats->total_probes));

  fprintf(stream,
          "         Total exterior probes: %19llu (%.1f%% of probes)\n",
          stats->exterior_probes,
          percentage(stats->exterior_probes, stats->total_probes));

  fprintf(stream,
          "             Probes per second: %19.0f (%s)\n",
          (double)stats->total_probes / (double)cpu_time,
          special_format_number((uint64)((real)stats->total_probes/(real)cpu_time)));

  fprintf(stream, "\n");


  // --- Iterations

  fprintf(stream,
          "              Total iterations: %19llu (%s)\n",
          stats->total_iter,
          special_format_number(stats->total_iter));

  fprintf(stream,
          "     Total interior iterations: %19llu (%.1f%% of iterations)\n",
          stats->interior_iter,
          percentage(stats->interior_iter, stats->total_iter));

  fprintf(stream,
          "     Total exterior iterations: %19llu (%.1f%% of iterations)\n",
          stats->exterior_iter,
          percentage(stats->exterior_iter, stats->total_iter));

  fprintf(stream,
          "         Iterations per second: %19.0f (%s)\n",
          (double)stats->total_iter / (double)cpu_time,
          special_format_number((uint64)((real)stats->total_iter/(real)cpu_time)));

  fprintf(stream, "\n");

  fprintf(stream,
          "     Avg. iterations per probe: %23.3f\n",
          ratio(stats->total_iter, stats->total_probes));

  fprintf(stream,
          " Avg. iter. per interior probe: %23.3f\n",
          ratio(stats->interior_iter, stats->interior_probes));

  fprintf(stream,
          " Avg. iter. per exterior probe: %23.3f\n",
          ratio(stats->exterior_iter, stats->exterior_probes));

  fprintf(stream,
          "         Avg. probes per pixel: %23.3f\n",
          ratio(stats->total_probes, total_pixels_including_overscan));

  fprintf(stream, "\n");


  // --- Iterations table

  fprintf(stream, "%s\n", hbar_thick);
  fprintf(stream,
          "%13s  %13s  %7s %6s  %13s  %7s %6s\n",
          "Iterations", "Interior", "%age", "%ile", "Exterior", "%age", "%ile");
  fprintf(stream,
          "%13s  %13s  %7s %6s  %13s  %7s %6s\n",
          "----------", "--------", "-------", "-----", "--------", "-------", "-----");

  uint64 accumulated_interior_probes = 0;
  uint64 accumulated_exterior_probes = 0;

  accumulated_interior_probes += stats->interior_probes_uniterated;
  accumulated_exterior_probes += stats->exterior_probes_uniterated;

  fprintf(stream,
    "%13llu  %13llu %7.3f%% %5.1f%%  %13llu %7.3f%% %5.1f%%\n",
    UINT64_C(0), 
    stats->interior_probes_uniterated,
    percentage(stats->interior_probes_uniterated, stats->interior_probes),
    percentage(accumulated_interior_probes, stats->interior_probes),
    stats->exterior_probes_uniterated,
    percentage(stats->exterior_probes_uniterated, stats->exterior_probes),
    percentage(accumulated_exterior_probes, stats->exterior_probes));

  int max_k = (int)floor(log2(conf->iter_max));
  for (int k = 0; k <= max_k; k++)
  {
    accumulated_interior_probes += stats->interior_probes_by_log2_iter[k];
    accumulated_exterior_probes += stats->exterior_probes_by_log2_iter[k];

    fprintf(stream,
      "%13llu  %13llu %7.3f%% %5.1f%%  %13llu %7.3f%% %5.1f%%\n",
      UINT64_C(1) << k,
      stats->interior_probes_by_log2_iter[k],
      percentage(stats->interior_probes_by_log2_iter[k],
                 stats->interior_probes),
      percentage(accumulated_interior_probes, stats->interior_probes),
                 stats->exterior_probes_by_log2_iter[k],
      percentage(stats->exterior_probes_by_log2_iter[k],
                 stats->exterior_probes),
      percentage(accumulated_exterior_probes,
                 stats->exterior_probes));
  }

  accumulated_interior_probes += stats->interior_probes_aperiodic;
  fprintf(stream,
    "%13llu  %13llu %7.3f%% %5.1f%%\n",
    conf->iter_max,
    stats->interior_probes_aperiodic,
    percentage(stats->interior_probes_aperiodic, stats->interior_probes),
    percentage(accumulated_interior_probes, stats->interior_probes));

  fprintf(stream, "\n");


  fflush(stream);
}


//-----------------------------------------------------------------------------
