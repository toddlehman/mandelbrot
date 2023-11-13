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
                    int subsample_min_depth, int subsample_max_depth,
                    real subsample_solidarity,
                    uint64 iter_max)
{
  assert(mp_sgn(xy_min_size) > 0);
  assert(width_pixels >= 1);
  assert(height_pixels >= 1);
  assert(subsample_min_depth >= 0);
  assert(subsample_max_depth >= 0);
  assert(subsample_min_depth <= subsample_max_depth);
  assert(subsample_solidarity >= 0);
  assert(subsample_solidarity <= 1);
  assert(iter_max > 0);

  Image *this = mem_alloc_clear(1, sizeof(Image));

  int mp_prec = mp_get_prec(x_center);  // KLUDGE for now
  mp_prec += ceil(log2(MAX(width_pixels, height_pixels)));
  mp_prec += ceil((real)subsample_max_depth / 2);
  fprintf(stderr, "precision required:  %d bits\n", mp_prec);
  mp_prec = MAX(mp_prec, 64);
  mp_prec = 64;
  fprintf(stderr, "precision using:  %d bits\n", mp_prec);

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

  this->subsample_min_depth = subsample_min_depth;
  this->subsample_max_depth = subsample_max_depth;
  this->subsample_solidarity = subsample_solidarity;

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
// COMPUTE NON-SUBSAMPLED PIXEL

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

  // Map iterations to color.
  pixel.color = palette_color_from_mandelbrot_result(this->palette, pixel.mr);

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
// POPULATE PIXEL

private_method
void image_populate_pixel(Image *this, int i, int j)
{
  assert(this);
  assert(i >= 0); assert(i < this->_di);
  assert(j >= 0); assert(j < this->_dj);

  unless (pixel_is_defined(this->pixels[i][j]))
  {
    this->pixels[i][j] = image_compute_pixel(this, i, j);
  }
}



//-----------------------------------------------------------------------------
// CALCULATE SUBSAMPLED PIXEL

private_method
Pixel image_compute_subsampled_pixel(Image *this,
                                     real i, real j, real di, real dj,
                                     Pixel pixel_00,
                                     int current_subsample_depth)
{
  assert(this);
  assert(i >= 0); assert(i + di <= this->di);
  assert(j >= 0); assert(j + dj <= this->dj);

  // Calculate intermediate subpixel values.
  Pixel pixel_01 = image_compute_pixel(this, i,        j+(dj/2));
  Pixel pixel_10 = image_compute_pixel(this, i+(di/2), j);
  Pixel pixel_11 = image_compute_pixel(this, i+(di/2), j+(dj/2));

  // Adjust required solidarity based on subsampling depth.  The deeper we go
  // in subsampling, the sloppier we allow the results to be.  This is somewhat
  // of a KLUDGE, but it seems to work okay and it allows the maximum subsample
  // depth to be specified as infinite because it short-circuits that limit (as
  // long as the value is not 1 to start with).
  float32 required_solidarity = this->subsample_solidarity;
  required_solidarity = 1 - ((1 - required_solidarity) *
                             pow(2, current_subsample_depth));

  // Assess color variation.
  float32 solidarity = linear_rgb_solidarity4(pixel_00.color, pixel_01.color,
                                              pixel_10.color, pixel_11.color);

  // If color variation exceeds tolerance, then subdivide each subquadrant of
  // the pixel recursively (unless the subsample depth limit has already been
  // reached).
  bool subsample = (current_subsample_depth < this->subsample_min_depth) ||
                   //((solidarity < this->subsample_solidarity) &&
                   ((solidarity < required_solidarity) &&
                    (current_subsample_depth < this->subsample_max_depth));
  if (subsample)
  {
    pixel_00 = image_compute_subsampled_pixel(this, i, j, di/2, dj/2,
                                          pixel_00, current_subsample_depth+1);

    pixel_01 = image_compute_subsampled_pixel(this, i, j+dj/2, di/2, dj/2,
                                          pixel_01, current_subsample_depth+1);

    pixel_10 = image_compute_subsampled_pixel(this, i+di/2, j, di/2, dj/2,
                                          pixel_10, current_subsample_depth+1);

    pixel_11 = image_compute_subsampled_pixel(this, i+di/2, j+dj/2, di/2, dj/2,
                                          pixel_11, current_subsample_depth+1);
  }

  // Now return a blend of the four subsamples.

  Pixel pixel;
  pixel.color = linear_rgb_average4(pixel_00.color, pixel_01.color,
                                    pixel_10.color, pixel_11.color);
  pixel.mr.dwell = -1;
  pixel.mr.period = -1;
  pixel.mr.iter = pixel_00.mr.iter + pixel_01.mr.iter +
                  pixel_10.mr.iter + pixel_11.mr.iter;

  if (0) //current_subsample_depth > 1)  // KLUDGE
  {
    pixel.color = linear_rgb_lerp(solidarity,
        palette_undefined_color(this->palette),
        pixel.color);
  }

  return pixel;
}


//-----------------------------------------------------------------------------
// REPOPULATE PIXEL

private_method
void image_repopulate_pixel(Image *this, int i, int j)
{
  assert(this);
  assert(i >= 0); assert(i < this->di);  // *Not* this->_di here.
  assert(j >= 0); assert(j < this->dj);  // *Not* this->_dj here.

  float32 solidarity = linear_rgb_solidarity4(
    this->pixels[i+0][j+0].color,
    this->pixels[i+0][j+1].color,
    this->pixels[i+1][j+0].color,
    this->pixels[i+1][j+1].color
  );

  if (solidarity < this->subsample_solidarity)
  {
    this->pixels[i][j] = image_compute_subsampled_pixel(
      this, i, j, 1, 1,
      this->pixels[i+0][j+0],
      1);
  }
  else
  {
    //this->pixels[i][j].color = (LinearRGB) { .r = .0, .g = .0, .b = .0 };
  }
}


//-----------------------------------------------------------------------------
// POPULATE BLOCK (RECURSIVE)

private_method
void image_populate_block(Image *this, int i0, int j0, int di, int dj)
{
  //fprintf(stderr,
  //        "block(i0=%d, j0=%d, di=%d, dj=%d)\n",
  //        i0, j0, di, dj);

  assert(this);
  assert(i0 >= 0); assert(i0 + di <= this->_di);
  assert(j0 >= 0); assert(j0 + dj <= this->_dj);


  // --- Populate the boundary of the block.
  //
  // TODO:  Reimplement this with early-out as soon as it's known that filling
  // is not going to happen.  This will result in smaller regions sooner, which
  // will greatly help cache coherency.

  bool fill_block = true;
  int i1 = i0 + di - 1;
  int j1 = j0 + dj - 1;
  if (true)  // Top edge
  {
    for (int j = j0; j <= j1; j++)
    {
      image_populate_pixel(this, i0, j);
      fill_block &= pixel_is_interior(this->pixels[i0][j]);
    }
  }
  if (di > 1)  // Bottom edge
  {
    for (int j = j0; j <= j1; j++)
    {
      image_populate_pixel(this, i1, j);
      fill_block &= pixel_is_interior(this->pixels[i1][j]);
    }
  }
  if (true)  // Left edge
  {
    for (int i = i0 + 1; i <= i1 - 1; i++)
    {
      image_populate_pixel(this, i, j0);
      fill_block &= pixel_is_interior(this->pixels[i][j0]);
    }
  }
  if (dj > 1)  // Right edge
  {
    for (int i = i0 + 1; i <= i1 - 1; i++)
    {
      image_populate_pixel(this, i, j1);
      fill_block &= pixel_is_interior(this->pixels[i][j1]);
    }
  }


  // --- Fill interior of block or subdivide.

  if (fill_block)
  {
    // If the boundary of the block consists entirely of M-Set interior pixels,
    // then fill the entire block with interior pixels.

    // FIXME:  This pixel value shouldn't be recomputed every time here.
    Pixel pixel = (Pixel)
    {
      .mr     = mandelbrot_result_interior_uniterated(),
      .color  = palette_interior_uniterated_color(this->palette)
    };

    for (int i = i0 + 1; i <= i1 - 1; i++)
    for (int j = j0 + 1; j <= j1 - 1; j++)
    {
      assert(pixel_is_undefined(this->pixels[i][j]));
      if (pixel_is_undefined(this->pixels[i][j]))
        this->pixels[i][j] = pixel;
    }
  }
  else if ((di > 2) && (dj > 2))
  {
    // Otherwise, split the block in two parts (either vertically or
    // horizontally, depending on the shape), and recurse for each part.
    // Note that the parts may not always be half the block, although they
    // typically are.

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
}


//-----------------------------------------------------------------------------
// POPULATE IMAGE

public_method
void image_populate(Image *this)
{
  assert(this);

  bool subsampling = (this->subsample_max_depth > 0) &&
                     (this->subsample_solidarity > 0);

  // Make a first pass to populate every pixel in the image, including the
  // extra row and column at the bottom and right edges if subsampling is
  // enabled.
  image_populate_block(this,
                       0,
                       0,
                       subsampling? this->_di : this->di,
                       subsampling? this->_dj : this->dj);

  // Return now if subsampling is not requested.
  if (!subsampling)
    return;

  // Go back and subsample pixels whose rightward and downward neighbors are
  // significantly different in color.  (The "+1"s in the indexes here are safe
  // because the image is padded by one row and column at the bottom and right
  // edges.)
  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    image_repopulate_pixel(this, i, j);

#if 0  // OBSOLETE
    if (solidarity < this->subsample_solidarity)
    {
      image_repopulate_pixel(this, i, j);

      LinearRGB color = linear_rgb_average4(
        this->pixels[i+0][j+0].color,
        this->pixels[i+0][j+1].color,
        this->pixels[i+1][j+0].color,
        this->pixels[i+1][j+1].color
      );

      this->pixels[i][j].color = linear_rgb_lerp(solidarity,
        palette_undefined_color(this->palette), color);
    }
#endif
  }
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

  fprintf(stream, "         Maximum iterations:  %llu\n", this->mandelbrot->iter_max);

  fprintf(stream, "                 Pixel size:  %d x %d\n", this->dj, this->di);

  fprintf(stream, "      Subsampling min depth:  %d\n", this->subsample_min_depth);
  fprintf(stream, "      Subsampling max depth:  %d\n", this->subsample_max_depth);
  fprintf(stream, "     Subsampling solidarity:  %f\n", this->subsample_solidarity);
  fprintf(stream, "\n");

  fprintf(stream, "               Total pixels: %20llu\n", total_pixels);
  fprintf(stream, "      Total interior pixels: %20llu (%.1f%%)\n", total_interior_pixels,
          (double)total_interior_pixels / (double)total_pixels * 100);
  fprintf(stream, "      Total exterior pixels: %20llu (%.1f%%)\n", total_exterior_pixels,
          (double)total_exterior_pixels / (double)total_pixels * 100);
  fprintf(stream, "\n");
  fprintf(stream, "           Total iterations: %20llu\n", total_iter);
  fprintf(stream, "  Total interior iterations: %20llu (%.1f%%)\n", total_interior_iter,
          (double)total_interior_iter / (double)total_iter * 100);
  fprintf(stream, "  Total exterior iterations: %20llu (%.1f%%)\n", total_exterior_iter,
          (double)total_exterior_iter / (double)total_iter * 100);
  fprintf(stream, "\n");

  fprintf(stream, "         Average iterations: %20.3f\n",
          (double)total_iter / (double)total_pixels);
  fprintf(stream, "Average interior iterations: %20.3f\n",
          (double)total_interior_iter / (double)total_interior_pixels);
  fprintf(stream, "Average exterior iterations: %20.3f\n",
          (double)total_exterior_iter / (double)total_exterior_pixels);
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
    (double)interior_tally_uniterated / (double)total_interior_pixels * 100,
    (double)interior_tally_log2_running_total / (double)total_interior_pixels * 100);
  for (int k = 0; k <= max_interior_k; k++)
  {
    interior_tally_log2_running_total += interior_tally_log2[k];
    fprintf(stream, "%25llu  %20llu  %9.6f%% %10.6f%%\n",
      UINT64_C(1) << k,
      interior_tally_log2[k],
      (double)interior_tally_log2[k] / (double)total_interior_pixels * 100,
      (double)interior_tally_log2_running_total / (double)total_interior_pixels * 100);
  }
  interior_tally_log2_running_total += interior_tally_aperiodic;
  sprintf(buf, "%s %llu", "(never)", this->mandelbrot->iter_max);
  fprintf(stream, "%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    interior_tally_aperiodic,
    (double)interior_tally_aperiodic / (double)total_interior_pixels * 100,
    (double)interior_tally_log2_running_total / (double)total_interior_pixels * 100);
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
    (double)exterior_tally_uniterated / (double)total_exterior_pixels * 100,
    (double)exterior_tally_log2_running_total / (double)total_exterior_pixels * 100);
  for (int k = 0; k <= max_exterior_k; k++)
  {
    exterior_tally_log2_running_total += exterior_tally_log2[k];
    fprintf(stream, "%25llu  %20llu  %9.6f%% %10.6f%%\n",
      UINT64_C(1) << k,
      exterior_tally_log2[k],
      (double)exterior_tally_log2[k] / (double)total_exterior_pixels * 100,
      (double)exterior_tally_log2_running_total / (double)total_exterior_pixels * 100);
  }
  exterior_tally_log2_running_total += exterior_tally_aperiodic;
  sprintf(buf, "%s %llu", "(never)", this->mandelbrot->iter_max);
  fprintf(stream, "%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    exterior_tally_aperiodic,
    (double)exterior_tally_aperiodic / (double)total_exterior_pixels * 100,
    (double)exterior_tally_log2_running_total / (double)total_exterior_pixels * 100);
}


