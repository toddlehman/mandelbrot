/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Image.h"
#import "Mandelbrot.h"
#import "RGB.h"


//-----------------------------------------------------------------------------
// ALLOCATE IMAGE

public_constructor
Image *image_create(mp_real x_center, mp_real y_center, mp_real xy_min_size,
                    int width_pixels, int height_pixels,
                    int supersample_interior_min_depth,
                    int supersample_interior_max_depth,
                    int supersample_exterior_min_depth,
                    int supersample_exterior_max_depth,
                    float32 supersample_solidarity,
                    uint64 iter_max)
{
  // FIXME:  Move these two assertions to initialization code for each of these
  // modules.
  assert(sizeof(MandelbrotResult) == 24);  // To catch unexpected changes.
  assert(sizeof(Pixel) == 20);  // To catch unexpected changes.

  assert(mp_sgn(xy_min_size) > 0);
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

  Image *this = mem_alloc_clear(1, sizeof(Image));

  // Calculate floating-point precision required.
  real bits;
  {
    mp_real log2_xy_min_size;
    mp_init2(log2_xy_min_size, 64);  // (Overkill; only need about 8 or 10)
    mp_log2(log2_xy_min_size, xy_min_size, MP_ROUND);  // Typically negative.
    bits = -mp_get_d(log2_xy_min_size, MP_ROUND);
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
  
  mp_init2(this->x_center, mp_prec); mp_init2(this->y_center, mp_prec);
  mp_init2(this->x_size,   mp_prec); mp_init2(this->y_size,   mp_prec);
  mp_init2(this->x_min,    mp_prec); mp_init2(this->y_min,    mp_prec);
  mp_init2(this->x_max,    mp_prec); mp_init2(this->y_max,    mp_prec);
  mp_init2(this->pixel_size, mp_prec);
  mp_init2(this->periodicity_epsilon, mp_prec);

  if (width_pixels == height_pixels)
  {
    mp_set(this->x_size, xy_min_size, MP_ROUND);
    mp_set(this->y_size, xy_min_size, MP_ROUND);
  }
  else if (width_pixels > height_pixels)
  {
    real aspect_ratio = (real)width_pixels / (real)height_pixels;
    mp_set(this->y_size, xy_min_size, MP_ROUND);
    mp_mul_d(this->x_size, xy_min_size, aspect_ratio, MP_ROUND);
  }
  else if (width_pixels < height_pixels)
  {
    real aspect_ratio = (real)width_pixels / (real)height_pixels;
    mp_set(this->x_size, xy_min_size, MP_ROUND);
    mp_div_d(this->y_size, xy_min_size, aspect_ratio, MP_ROUND);
  }

  this->di = height_pixels;
  this->dj = width_pixels;
  this->_di = this->di + 1;
  this->_dj = this->dj + 1;

  mp_set(this->x_center, x_center, MP_ROUND);
  mp_set(this->y_center, y_center, MP_ROUND);

  {
    mp_real half_x_size; mp_init2(half_x_size, mp_prec);
    mp_real half_y_size; mp_init2(half_y_size, mp_prec);

    mp_div_ui(half_x_size, this->x_size, 2, MP_ROUND);
    mp_div_ui(half_y_size, this->y_size, 2, MP_ROUND);

    mp_sub(this->x_min, this->x_center, half_x_size, MP_ROUND);
    mp_sub(this->y_min, this->y_center, half_y_size, MP_ROUND);

    mp_add(this->x_max, this->x_center, half_x_size, MP_ROUND);
    mp_add(this->y_max, this->y_center, half_y_size, MP_ROUND);

    mp_clear(half_y_size);
    mp_clear(half_x_size);
  }

  // FIXME:  When doing 3D viewport, this will no longer make sense.
  mp_div_ui(this->pixel_size, this->x_size, this->dj, MP_ROUND);

  {
    mp_real periodicity_epsilon;
    mp_init2(periodicity_epsilon, mp_prec);
    mp_div_d(periodicity_epsilon, this->pixel_size, 1e6, MP_ROUND);
    this->mandelbrot = mandelbrot_create(iter_max, mp_prec,
                                         periodicity_epsilon);
    mp_clear(periodicity_epsilon);
  }

  this->supersample_interior_min_depth = supersample_interior_min_depth;
  this->supersample_interior_max_depth = supersample_interior_max_depth;
  this->supersample_exterior_min_depth = supersample_exterior_min_depth;
  this->supersample_exterior_max_depth = supersample_exterior_max_depth;
  this->supersample_solidarity = supersample_solidarity;
  this->supersample = ((supersample_interior_max_depth > 0) ||
                       (supersample_exterior_max_depth > 0))
                   && (supersample_solidarity > 0);

  this->_pixels = mem_alloc_clear(this->_di * this->_dj, sizeof(Pixel));
  this->pixels = mem_alloc_clear(this->_di, sizeof(Pixel **));
  for (int i = 0; i < this->_di; i++)
    this->pixels[i] = this->_pixels + (i * this->_dj);

  for (int i = 0; i < this->_di; i++)
  for (int j = 0; j < this->_dj; j++)
  {
    this->pixels[i][j].is_defined = 0;
  }

  this->palette = palette_create();

  this->interior_filler_pixel = (Pixel)
  {
    #if 0  // OBSOLETE
    .color                 = palette_interior_uniterated_color(this->palette),
    #endif
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

  assert(this->palette);
  palette_destroy(&this->palette);

  assert(this->mandelbrot);
  mandelbrot_destroy(&this->mandelbrot);

  mp_clear(this->x_center);
  mp_clear(this->y_center);
  mp_clear(this->x_size);
  mp_clear(this->y_size);
  mp_clear(this->x_min);
  mp_clear(this->y_min);
  mp_clear(this->x_max);
  mp_clear(this->y_max);
  mp_clear(this->pixel_size);
  mp_clear(this->periodicity_epsilon);

  mem_dealloc(p_this);
}


//-----------------------------------------------------------------------------
// MAP IMAGE COORDINATES TO M-SET COORDINATES

private_method
bool image_get_argand_point(Image *this,
                            real i, real j,
                            mp_real *x, mp_real *y)
{
  assert(this); assert(x); assert(y);
  assert(i >= 0); assert(i <= (real)this->di);
  assert(j >= 0); assert(j <= (real)this->dj);

  //mp_lerp_d(*x, (real)j / (real)this->dj, this->x_min, this->x_max);
  //mp_lerp_d(*y, (real)i / (real)this->di, this->y_max, this->y_min);

  mp_mul_d(*x, this->x_size, (real)j / (real)this->dj, MP_ROUND);
  mp_add(*x, this->x_min, *x, MP_ROUND);

  mp_mul_d(*y, this->y_size, (real)i / (real)this->di, MP_ROUND);
  mp_sub(*y, this->y_max, *y, MP_ROUND);

  return true;
}


//-----------------------------------------------------------------------------
// DETERMINE WHETHER PIXEL NEEDS SUPERSAMPLING

private_inline_method
bool image_pixel_needs_supersampling(Image *this, Pixel pixel,
                                     float32 solidarity, int ss_depth)
{
  // Start with specified solidarity from configuration.
  float32 required_solidarity = this->supersample_solidarity;

  // Adjust required solarity downward (allow sloppiness) based on current
  // supersampling depth.  For example:
  //
  // Depth         0      1      2      3      4      5      6      7      8   
  // -------------------------------------------------------------------------
  // Required    0.999  0.998  0.996  0.992  0.984  0.968  0.936  0.872  0.744
  // Required    0.990  0.980  0.960  0.920  0.840  0.680  0.360    -      -
  // Required    0.980  0.960  0.920  0.840  0.680  0.360    -      -      -
  // Required    0.970  0.940  0.880  0.760  0.520  0.040    -      -      -
  // Required    0.960  0.920  0.840  0.680  0.360    -      -      -      -
  // Required    0.950  0.900  0.800  0.600  0.200    -      -      -      -
  // Required    0.940  0.880  0.760  0.520  0.040    -      -      -      -
  // Required    0.930  0.860  0.720  0.440    -      -      -      -      -
  // Required    0.920  0.840  0.680  0.360    -      -      -      -      -
  // Required    0.910  0.820  0.640  0.280    -      -      -      -      -
  // Required    0.900  0.800  0.600  0.200    -      -      -      -      -
  // Required    0.850  0.700  0.400    -      -      -      -      -      -
  // Required    0.800  0.600  0.200    -      -      -      -      -      -
  // Required    0.750  0.500    -      -      -      -      -      -      -
  // Required    0.700  0.400    -      -      -      -      -      -      -
  // Required    0.650  0.300    -      -      -      -      -      -      -
  // Required    0.600  0.200    -      -      -      -      -      -      -
  // Required    0.550  0.100    -      -      -      -      -      -      -
  // Required    0.500  0.000    -      -      -      -      -      -      -
  //
  // NOTE:  Is this adjustment really a good thing?  Is it useful?  Or should
  // it be removed?
  //
  required_solidarity = 1 - ((1 - required_solidarity) * (1 << ss_depth));

  // Now decide based on the adjusted requirement for solidarity, the nature of
  // the pixel, the minimum and maximum supersampling depths, and the current
  // supersampling depth.
  if (pixel_is_interior(pixel))
  {
    return (ss_depth < this->supersample_interior_min_depth) ||
           ((ss_depth < this->supersample_interior_max_depth) &&
            (solidarity < required_solidarity));
  }
  else
  {
    return (ss_depth < this->supersample_exterior_min_depth) ||
           ((ss_depth < this->supersample_exterior_max_depth) &&
            (solidarity < required_solidarity));
  }
}


//-----------------------------------------------------------------------------
// COMPUTE NON-SUPERSAMPLED PIXEL

private_method
Pixel image_compute_pixel(Image *this, real i, real j)
{
  assert(this);

  Pixel pixel;

  // Map image coordinates to Argand plane coordinates.
  mp_real x, y;
  mp_init2(x, this->mandelbrot->conf.mp_prec);
  mp_init2(y, this->mandelbrot->conf.mp_prec);
  (void)image_get_argand_point(this, i, j, &x, &y);

  // Calculate iterations.
  MandelbrotResult mr = mandelbrot_compute(this->mandelbrot, x, y);
  mp_clear(y);
  mp_clear(x);

  // Output progress (for debugging).
  if (false) //if (this->show_progress)
    fprintf(stderr, "i=%7.2f  j=%7.2f  iter=%9llu\n",
                    (double)i, (double)j, mr.iter);

  // Map iterations to color.
  pixel.color = palette_color_from_mandelbrot_result(this->palette, mr);

  // Define interior portion.
  pixel.interior_portion = mandelbrot_result_is_interior(mr)? 1.0 : 0.0;

  // Set flags.
  pixel.is_defined = true;
  pixel.is_interior_periodic = mandelbrot_result_is_interior_periodic(mr);
  pixel.supersample = false;
  pixel.supersampled = false;
  pixel.probed_top_edge = false;
  pixel.probed_left_edge = false;

  #if 0
  if ((i == floor(i)) && (j == floor(j)))
    fprintf(stderr, "pixel(i=%d, j=%d) -> %llu (%.3f,%.3f,%.3f)\n",
            (int)i, (int)j, pixel.mr.iter,
            pixel.color.r, pixel.color.g, pixel.color.b);
  else
    fprintf(stderr, "pixel(i=%f, j=%f) -> %llu (%.3f,%.3f,%.3f)\n",
            i, j, pixel.mr.iter,
            pixel.color.r, pixel.color.g, pixel.color.b);
  #endif

  return pixel;
}


//-----------------------------------------------------------------------------
// CALCULATE SUPERSAMPLED PIXEL

private_method
Pixel image_compute_supersampled_pixel(Image *this,
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
  //fprintf(stderr, "i=%f j=%f di=%f dj=%f ss_depth=%d -> interior=%f\n",
  //        (double)i, (double)j, (double)di, (double)dj, ss_depth,
  //        (double)pixel.interior_portion);

  #if 0
  if (ss_depth > 1)  // KLUDGE
  {
    pixel.color = linear_rgb_lerp(pow(solidarity, 0.5),
                                  palette_undefined_color(this->palette),
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

  if (pixel->is_interior_periodic && this->supersample)
  {
    bool all_interior = true;

    if (probe_top_edge && !pixel->probed_top_edge)
    {
      int dsj = (1 << this->supersample_interior_max_depth);
      for (int sj = 1; all_interior && (sj < dsj); sj++)  // (Not sj = 0)
      {
        Pixel subpixel = image_compute_pixel(this, i, j+(real)sj/(real)dsj);
        all_interior &= subpixel.is_interior_periodic;
      }
      pixel->probed_top_edge = true;
    }

    if (probe_left_edge && !pixel->probed_left_edge)
    {
      int dsi = (1 << this->supersample_interior_max_depth);
      for (int si = 1; all_interior && (si < dsi); si++)  // (Not si = 0)
      {
        Pixel subpixel = image_compute_pixel(this, i+(real)si/(real)dsi, j);
        all_interior &= subpixel.is_interior_periodic;
      }
      pixel->probed_left_edge = true;
    }

    if (!all_interior)
    {
      pixel->supersample = true;  // Mark for future refinement.
      pixel->supersampled = false;
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
  //fprintf(stderr, "block(i0=%d, j0=%d, di=%d, dj=%d)\n", i0, j0, di, dj);

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
      //fprintf(stderr, "i=%d, j=%d\n", i, j);

      if ((i > i0) && (i < i1) && (j > j0) && (j < j1))
        { j = j1 - 1; continue; }  // Efficiently skip non-boundary pixels.

      bool probe_top_edge  = ((i == i0) || (i == i1)) && (j != j1);
      bool probe_left_edge = ((j == j0) || (j == j1)) && (i != i1);

      image_populate_pixel(this, i, j, probe_top_edge, probe_left_edge);

      fill_block &= this->pixels[i][j].is_interior_periodic;
    }
  }
  //fprintf(stderr, "done with boundary\n");


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
void image_output(Image *this, FILE *stream, bool text_format)
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
char *special_format_mp_real(mp_real x, bool sign)
{
  // Calculate bits of precision to the right of the decimal point by
  // subtracting the bits of precision to the left of the decimal point from
  // the total bits of precision.
  double x_trunc = floor(abs(mp_get_d(x, MP_ROUND)));
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
char *special_format_number(uint64 n)
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
                       "%f %s",
           f,
           units[unit_index]);

  return str;  // Not thread-safe; be careful!
}


//-----------------------------------------------------------------------------
// FORMAT FRACTIONAL PIXEL COUNT FOR READABILITY

private_method
char *special_format_fractional_pixel_count(Image *this, float64 value)
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
void image_output_statistics(Image *this, FILE *stream)
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
             special_format_mp_real(this->x_center, true));

  mp_fprintf(stream,
             "                      Center y:  %s\n",
             special_format_mp_real(this->y_center, true));

  mp_fprintf(stream,
             "                        Size x:  %s\n",
             special_format_mp_real(this->x_size, false));

  mp_fprintf(stream,
             "                        Size y:  %s\n",
             special_format_mp_real(this->y_size, false));

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
