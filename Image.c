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
                    int supersample_int_min_depth,
                    int supersample_int_max_depth,
                    int supersample_ext_min_depth,
                    int supersample_ext_max_depth,
                    float32 supersample_solidarity,
                    uint64 iter_max)
{
  // FIXME:  Move these two assertions to initialization code for each of these
  // modules.
  assert(sizeof(MandelbrotResult) == 24);  // To catch unexpected changes.
  assert(sizeof(Pixel) == 40);  // To catch unexpected changes.

  assert(mp_sgn(xy_min_size) > 0);
  assert(width_pixels >= 1);
  assert(height_pixels >= 1);
  assert(supersample_int_min_depth >= 0);
  assert(supersample_int_max_depth >= 0);
  assert(supersample_int_min_depth <= supersample_int_max_depth);
  assert(supersample_ext_min_depth >= 0);
  assert(supersample_ext_max_depth >= 0);
  assert(supersample_ext_min_depth <= supersample_ext_max_depth);
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
  bits += MAX(supersample_int_max_depth, supersample_ext_max_depth);
  // Now add a few guard bits.  12 is a bare minimum which just barely is
  // sufficient.  16 is a better value.
  bits += 16;
  // Finally, round up and use this as the required precision.
  int mp_prec = ceil(bits);
  fprintf(stderr, "precision required:  %d bits\n", mp_prec);
  
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

  this->supersample_int_min_depth = supersample_int_min_depth;
  this->supersample_int_max_depth = supersample_int_max_depth;
  this->supersample_ext_min_depth = supersample_ext_min_depth;
  this->supersample_ext_max_depth = supersample_ext_max_depth;
  this->supersample_solidarity = supersample_solidarity;
  this->supersample = ((supersample_int_max_depth > 0) ||
                       (supersample_ext_max_depth > 0))
                   && (supersample_solidarity > 0);

  this->_pixels = mem_alloc_clear(this->_di * this->_dj, sizeof(Pixel));
  this->pixels = mem_alloc_clear(this->_di, sizeof(Pixel **));
  for (int i = 0; i < this->_di; i++)
    this->pixels[i] = this->_pixels + (i * this->_dj);

  for (int i = 0; i < this->_di; i++)
  for (int j = 0; j < this->_dj; j++)
  {
    this->pixels[i][j].mr = mandelbrot_result_undefined();
  }

  this->palette = palette_create();

  this->interior_filler_pixel = (Pixel)
  {
    .mr                = mandelbrot_result_interior_uniterated(),
    .color             = palette_interior_uniterated_color(this->palette),
    .supersample       = false,
    .supersampled      = false,
    .probed_top_edge   = false,
    .probed_left_edge  = false,
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
bool image_compute_argand_point(Image *this,
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
  if (pixel_is_interior)
  {
    return (ss_depth < this->supersample_int_min_depth) ||
           ((ss_depth < this->supersample_int_max_depth) &&
            (solidarity < required_solidarity));
  }
  else
  {
    return (ss_depth < this->supersample_ext_min_depth) ||
           ((ss_depth < this->supersample_ext_max_depth) &&
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
  mp_init2(x, this->mandelbrot->mp_prec);
  mp_init2(y, this->mandelbrot->mp_prec);
  (void)image_compute_argand_point(this, i, j, &x, &y);

  // Calculate iterations.
  pixel.mr = mandelbrot_compute(this->mandelbrot, x, y);
  mp_clear(y);
  mp_clear(x);

  // Output progress (for debugging).
  if (false) //if (this->show_progress)
    fprintf(stderr, "i=%7.2f  j=%7.2f  iter=%9llu\n",
                    (double)i, (double)j, pixel.mr.iter);

  // Map iterations to color.
  pixel.color = palette_color_from_mandelbrot_result(this->palette, pixel.mr);

  // Set flags.
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
  pixel.mr.dwell = -1;
  pixel.mr.period = -1;
  pixel.mr.iter = pixel_00.mr.iter + pixel_01.mr.iter +
                  pixel_10.mr.iter + pixel_11.mr.iter;

  if (0) //ss_depth > 1)  // KLUDGE
  {
    pixel.color = linear_rgb_lerp(solidarity,
        palette_undefined_color(this->palette),
        pixel.color);
  }

  pixel.supersampled = true;   // Supersampling accomplished for this pixel.
  pixel.supersample  = false;  // Supersampling not required for this pixel.

  //pixel.probed_top_edge = false;
  //pixel.probed_left_edge = false;

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

  if (pixel_is_undefined(*pixel))
    *pixel = image_compute_pixel(this, i, j);

  if (pixel_is_interior(*pixel) && this->supersample)
  {
    bool interior = true;

    if (probe_top_edge && !pixel->probed_top_edge)
    {
      int dsj = (1 << this->supersample_int_max_depth);
      for (int sj = 1; interior && (sj < dsj); sj++)  // (Not sj = 0)
      {
        Pixel subpixel = image_compute_pixel(this, i, j+(real)sj/(real)dsj);
        pixel->mr.iter += subpixel.mr.iter;
        interior &= pixel_is_interior(subpixel);
      }
      pixel->probed_top_edge = true;
    }

    if (probe_left_edge && !pixel->probed_left_edge)
    {
      int dsi = (1 << this->supersample_int_max_depth);
      for (int si = 1; interior && (si < dsi); si++)  // (Not si = 0)
      {
        Pixel subpixel = image_compute_pixel(this, i+(real)si/(real)dsi, j);
        pixel->mr.iter += subpixel.mr.iter;
        interior &= pixel_is_interior(subpixel);
      }
      pixel->probed_left_edge = true;
    }

    if (interior)
    {
      //pixel->mr.dwell = INFINITY;  // (Already set)
      pixel->mr.period = -1;
    }
    else
    {
      pixel->mr.dwell = -1;
      pixel->mr.period = -1;
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
  // increase in quality in external regions for a slight increase in cost.)
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

      fill_block &= pixel_is_interior(this->pixels[i][j]);
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
      assert(pixel_is_undefined(this->pixels[i][j]));
      if (pixel_is_undefined(this->pixels[i][j]))
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
    assert(pixel_is_defined(this->pixels[i][j]));
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
// OUTPUT IMAGE STATISTICS

public_method
void image_output_statistics(Image *this, FILE *stream)
{
  assert(this);

  uint64 total_interior_pixels = 0;
  uint64 total_exterior_pixels = 0;
  uint64 total_interior_iter = 0;
  uint64 total_exterior_iter = 0;

  uint64 interior_tally_log2[64];
  uint64 interior_tally_aperiodic = 0;
  uint64 interior_tally_uniterated = 0;
  for (int k = 0; k < ELEMENT_COUNT(interior_tally_log2); k++)
    interior_tally_log2[k] = 0;

  uint64 exterior_tally_log2[64];
  uint64 exterior_tally_aperiodic = 0;
  uint64 exterior_tally_uniterated = 0;
  for (int k = 0; k < ELEMENT_COUNT(exterior_tally_log2); k++)
    exterior_tally_log2[k] = 0;


  // --- Tally pixel categories.

  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    Pixel pixel = this->pixels[i][j];

    if (pixel_is_defined(pixel))
    {
      MandelbrotResult mr = pixel.mr;

      if (mandelbrot_result_is_interior(mr))
      {
        total_interior_pixels++;
        total_interior_iter += mr.iter;

        if (mandelbrot_result_is_interior_iterated(mr))
        {
          if (mandelbrot_result_is_interior_aperiodic(mr))
          {
            interior_tally_aperiodic++;
          }
          else if (mandelbrot_result_is_interior_periodic(mr))
          {
            int k = (int)(log(mr.iter) / log(2));
            interior_tally_log2[k]++;
            // TODO:  Do something about talling the periods.
          }
        }
        else
        {
          interior_tally_uniterated++;
        }
      }
      else if (mandelbrot_result_is_exterior(pixel.mr))
      {
        total_exterior_pixels++;
        total_exterior_iter += mr.iter;

        if (mr.iter == 0)
        {
          exterior_tally_uniterated++;
        }
        #if 0  // OBSOLETE
        else if (mr.iter == this->mandelbrot->iter_max)
        {
          // NOTE: This *can* occur due to loop unrolling, which occasionally
          // may cause overage on the iterations.
          assert(false);  // This should never occur.
          exterior_tally_aperiodic++;
        }
        #endif
        else
        {
          int k = (int)(log(mr.iter) / log(2));
          exterior_tally_log2[k]++;
        }
      }
      else
      {
        assert(false);
      }
    }
    else // if (pixel_is_undefined(pixel))
    {
      assert(false);
    }
  }


  // --- Print table.

  uint64 total_pixels = total_interior_pixels + total_exterior_pixels;
  assert(total_pixels == (uint64)this->di * (uint64)this->dj);

  uint64 total_iter = total_interior_iter + total_exterior_iter;

  mp_fprintf(stream,
             "                     Center:  (%RNf,%RNf)\n",
             this->x_center, this->y_center);

  mp_fprintf(stream,
             "                       Size:  (%RNf,%RNf)\n",
             this->x_size, this->y_size);

  fprintf(stream,
          "         Maximum iterations:  %llu\n",
          this->mandelbrot->iter_max);

  fprintf(stream,
          "                 Pixel size:  %d x %d\n",
          this->dj, this->di);

  fprintf(stream,
          "Supersampling int min depth:  %d\n",
          this->supersample_int_min_depth);

  fprintf(stream,
          "Supersampling int max depth:  %d\n",
          this->supersample_int_max_depth);

  fprintf(stream,
          "Supersampling ext min depth:  %d\n",
          this->supersample_ext_min_depth);

  fprintf(stream,
          "Supersampling ext max depth:  %d\n",
          this->supersample_ext_max_depth);

  fprintf(stream,
          "   Supersampling solidarity:  %f\n",
          (double)this->supersample_solidarity);

  fprintf(stream, "\n");

  fprintf(stream,
          "               Total pixels: %20llu\n",
          total_pixels);

  fprintf(stream,
          "      Total interior pixels: %20llu (%.1f%%)\n",
          total_interior_pixels,
          percentage(total_interior_pixels, total_pixels));

  fprintf(stream,
          "      Total exterior pixels: %20llu (%.1f%%)\n",
          total_exterior_pixels,
          percentage(total_exterior_pixels, total_pixels));

  fprintf(stream, "\n");

  fprintf(stream,
          "           Total iterations: %20llu\n",
          total_iter);

  fprintf(stream,
          "  Total interior iterations: %20llu (%.1f%%)\n",
          total_interior_iter,
          percentage(total_interior_iter, total_iter));

  fprintf(stream,
          "  Total exterior iterations: %20llu (%.1f%%)\n",
          total_exterior_iter,
          percentage(total_exterior_iter, total_iter));

  fprintf(stream, "\n");

  fprintf(stream,
          "         Average iterations: %24.3f\n",
          ratio(total_iter, total_pixels));

  fprintf(stream,
          "Average interior iterations: %24.3f\n",
          ratio(total_interior_iter, total_interior_pixels));

  fprintf(stream,
          "Average exterior iterations: %24.3f\n",
          ratio(total_exterior_iter, total_exterior_pixels));

  fprintf(stream, "\n");


  char buf[100];


  // Print periodic orbit table.

  fprintf(stream, "%25s  %20s  %9s  %9s\n",
          "Period detected after", "Quantity", "Percentage", "Percentile");
  fprintf(stream, "%25s  %20s  %9s  %9s\n",
          "---------------------", "--------", "----------", "----------");
  int max_interior_k = (int)(log(this->mandelbrot->iter_max) / log(2));
  uint64 interior_tally_log2_running_total = 0;

  interior_tally_log2_running_total += interior_tally_uniterated;
  sprintf(buf, "%s %llu", "(early-out)", UINT64_C(0));
  fprintf(stream, "%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    interior_tally_uniterated,
    percentage(interior_tally_uniterated, total_interior_pixels),
    percentage(interior_tally_log2_running_total, total_interior_pixels));
  for (int k = 0; k <= max_interior_k; k++)
  {
    interior_tally_log2_running_total += interior_tally_log2[k];
    fprintf(stream, "%25llu  %20llu  %9.6f%% %10.6f%%\n",
      UINT64_C(1) << k,
      interior_tally_log2[k],
      percentage(interior_tally_log2[k], total_interior_pixels),
      percentage(interior_tally_log2_running_total, total_interior_pixels));
  }
  interior_tally_log2_running_total += interior_tally_aperiodic;
  sprintf(buf, "%s %llu", "(never)", this->mandelbrot->iter_max);
  fprintf(stream, "%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    interior_tally_aperiodic,
    percentage(interior_tally_aperiodic, total_interior_pixels),
    percentage(interior_tally_log2_running_total, total_interior_pixels));
  fprintf(stream, "\n");

  // Print escaped-point table.

  fprintf(stream, "%25s  %20s  %9s  %9s\n",
          "Escaped after", "Quantity", "Percentage", "Percentile");
  fprintf(stream, "%25s  %20s  %9s  %9s\n",
          "-------------", "--------", "----------", "----------");
  int max_exterior_k = (int)(log(this->mandelbrot->iter_max) / log(2));
  uint64 exterior_tally_log2_running_total = 0;

  exterior_tally_log2_running_total += exterior_tally_uniterated;
  sprintf(buf, "%s %llu", "(early-out)", UINT64_C(0));
  fprintf(stream, "%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    exterior_tally_uniterated,
    percentage(exterior_tally_uniterated, total_exterior_pixels),
    percentage(exterior_tally_log2_running_total, total_exterior_pixels));
  for (int k = 0; k <= max_exterior_k; k++)
  {
    exterior_tally_log2_running_total += exterior_tally_log2[k];
    fprintf(stream, "%25llu  %20llu  %9.6f%% %10.6f%%\n",
      UINT64_C(1) << k,
      exterior_tally_log2[k],
      percentage(exterior_tally_log2[k], total_exterior_pixels),
      percentage(exterior_tally_log2_running_total, total_exterior_pixels));
  }
  exterior_tally_log2_running_total += exterior_tally_aperiodic;
  sprintf(buf, "%s %llu", "(never)", this->mandelbrot->iter_max);
  fprintf(stream, "%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    exterior_tally_aperiodic,
      percentage(exterior_tally_aperiodic, total_exterior_pixels),
      percentage(exterior_tally_log2_running_total, total_exterior_pixels));

  fflush(stream);
}

