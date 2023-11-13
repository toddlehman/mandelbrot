/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Image.h"
#import "Mandelbrot.h"
#import "RGB.h"


//-----------------------------------------------------------------------------
// ALLOCATE IMAGE

public_method
Image *image_alloc(real x_center, real y_center, real x_size,
                   uint64 iter_max,
                   int pixel_width, int pixel_height,
                   int subsample_scale, real subsample_tolerance)
{
  assert(pixel_width >= 2);
  assert(pixel_height >= 2);
  assert(x_size > 0);
  assert(iter_max > 0);
  assert(subsample_scale >= 1);
  assert(subsample_tolerance >= 0);
  assert(subsample_tolerance <= 1);

  real y_size = x_size * (real)pixel_height / (real)pixel_width;

  Image *this = mem_alloc_clear(1, sizeof(Image));

  this->di = pixel_height;
  this->dj = pixel_width;
  this->x_center = x_center;
  this->y_center = y_center;
  this->x_min = x_center - (x_size / 2);
  this->x_max = x_center + (x_size / 2);
  this->y_min = y_center - (y_size / 2);
  this->y_max = y_center + (y_size / 2);
  this->pixel_size = x_size / (real)pixel_width;
  this->periodicity_epsilon = this->pixel_size / 1e6;
  this->iter_max = iter_max;
  this->subsample_scale = subsample_scale;
  this->subsample_tolerance = subsample_tolerance;

  this->_pixels = mem_alloc_clear(this->di * this->dj, sizeof(Pixel));
  this->pixels = mem_alloc_clear(this->di, sizeof(Pixel **));
  for (int i = 0; i < this->di; i++)
    this->pixels[i] = this->_pixels + (i * this->dj);

  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    this->pixels[i][j].mr = mandelbrot_result_undefined();
  }

  this->palette = palette_alloc();

  return this;
}


//-----------------------------------------------------------------------------
// DEALLOCATE IMAGE

public_method
void image_dealloc(Image **p_this)
{
  assert(p_this);
  Image *this = *p_this;

  assert(this);
  assert(this->pixels);
  assert(this->_pixels);
  assert(this->palette);

  mem_dealloc(&this->pixels);
  mem_dealloc(&this->_pixels);
  palette_dealloc(&this->palette);

  mem_dealloc(p_this);
}


#if 0  // OBSOLETE
//-----------------------------------------------------------------------------
// MAP IMAGE COORDINATES TO M-SET COORDINATES

private_method
ArgandPoint image_argand_point(Image *this, int i, int j)
{
  assert(this);

  return (ArgandPoint)
  {
    .x = lerp((real)j / (real)(this->dj - 1), this->x_min, this->x_max),
    .y = lerp((real)i / (real)(this->di - 1), this->y_max, this->y_min)
  };
}
#endif


//-----------------------------------------------------------------------------
// GET PIXEL

#if 0  // OBSOLETE
private_method
Pixel image_get_pixel(Image *this, int i, int j)
{
  assert(this);

  if ((i >= 0) && (i < this->di) &&
      (j >= 0) && (j < this->dj))
  {
    return this->pixels[i][j];
  }
  else
  {
    return (Pixel)
    {
      .mr     = mandelbrot_result_undefined(),
      .color  = palette_undefined_color(this->palette)
    };
  }
}
#endif


//-----------------------------------------------------------------------------
// POPULATE PIXEL

private_method
void image_populate_pixel(Image *this, int i, int j)
{
  assert(this);
  assert(i >= 0); assert(i < this->di);
  assert(j >= 0); assert(j < this->dj);

  Pixel *pixel = &this->pixels[i][j];

  // Early-out if value is already known.
  if (pixel_is_defined(*pixel))
    return;

  // Map image coordinates to Argand plane coordinates.
  #if 0  // OBSOLETE
  ArgandPoint c = image_argand_point(this, i, j);
  #endif
  real cx = lerp((real)j / (real)(this->dj - 1), this->x_min, this->x_max);
  real cy = lerp((real)i / (real)(this->di - 1), this->y_max, this->y_min);

  // Calculate iterations.
  pixel->mr = mandelbrot_compute(cx, cy,
                                 this->periodicity_epsilon,
                                 this->iter_max);

  // Map iterations to color.
  pixel->color = palette_color_from_mandelbrot_result(this->palette, pixel->mr);
}


//-----------------------------------------------------------------------------
// POPULATE BLOCK

private_method
void image_populate_block(Image *this, int i0, int j0, int di, int dj)
{
  //printf("image_populate_block(i0=%d, j0=%d, di=%d, dj=%d)\n",
  //       j0, j0, di, dj);

  assert(this);
  assert(i0 >= 0); assert(i0 + di <= this->di);
  assert(j0 >= 0); assert(j0 + dj <= this->dj);


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

  // Make a first pass to populate every pixel in the image.
  image_populate_block(this, 0, 0, this->di, this->dj);

  // Return now if subsampling is not requested.
  if ((this->subsample_scale == 1) || (this->subsample_tolerance >= 1))
    return;

  // Go back and subsample pixels whose neighbors are significantly different
  // in color.  First mark all such pixels for recalculation.
  for (int i = 0; i < this->di - 1; i++)
  for (int j = 0; j < this->dj - 1; j++)
  {
    float32 diff = linear_rgb_diff4(
      this->pixels[i+0][j+0].color,
      this->pixels[i+0][j+1].color,
      this->pixels[i+1][j+0].color,
      this->pixels[i+1][j+1].color
    );

    if (diff > this->subsample_tolerance)
    {
      //this->pixels[i][j].color = palette_undefined_color(this->palette);

      LinearRGB color = linear_rgb_average4(
        this->pixels[i+0][j+0].color,
        this->pixels[i+0][j+1].color,
        this->pixels[i+1][j+0].color,
        this->pixels[i+1][j+1].color
      );

      this->pixels[i][j].color = linear_rgb_lerp(diff,
        color, palette_undefined_color(this->palette));
    }
  }
}


//-----------------------------------------------------------------------------
// OUTPUT IMAGE

public_method
void image_output(Image *this, bool text_format)
{
  assert(this);

  for (int i = 0; i < this->di; i++)
  for (int j = 0; j < this->dj; j++)
  {
    assert(pixel_is_defined(this->pixels[i][j]));
  }

  printf("%s\n", text_format? "P3" : "P6");
  printf("%d %d\n", this->dj, this->di);
  printf("%d\n", 255);
  fflush(stdout);  // Flush now so as not to conflict with unbuffered write().

  int buffer_byte_count = (text_format? strlen("000 000 000 ") : 3) *
                          this->dj;
  byte *buffer = mem_alloc_clear(buffer_byte_count, 1);

  for (int i = 0; i < this->di; i++)
  {
    byte *p = &buffer[0];

    for (int j = 0; j < this->dj; j++)
    {
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

    write(1, buffer, buffer_byte_count);
  }

  free(buffer);
}


//-----------------------------------------------------------------------------
// PRINT CLEAN FLOATING-POINT NUMBER

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
  printf("%s", buf);
}


//-----------------------------------------------------------------------------
// OUTPUT IMAGE STATISTICS

public_method
void image_output_statistics(Image *this)
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
        else if (mr.iter == this->iter_max)
        {
          assert(false);  // This should never occur.
          exterior_tally_aperiodic++;
        }
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

  printf("                     Center:  (");
  image_print_clean_real(this->x_center);
  printf(",");
  image_print_clean_real(this->y_center);
  printf(")\n");

  printf("                       Size:  (");
  image_print_clean_real(this->x_max - this->x_min);
  printf(",");
  image_print_clean_real(this->y_max - this->y_min);
  printf(")\n");

  printf("         Maximum iterations:  %llu\n", this->iter_max);

  printf("                 Pixel size:  %d x %d\n", this->dj, this->di);

  printf("          Subsampling scale:  %d\n", this->subsample_scale);
  printf("      Subsampling tolerance:  %f\n", this->subsample_tolerance);
  printf("\n");

  printf("               Total pixels: %20llu\n", total_pixels);
  printf("      Total interior pixels: %20llu (%.1f%%)\n", total_interior_pixels,
         (double)total_interior_pixels / (double)total_pixels * 100);
  printf("      Total exterior pixels: %20llu (%.1f%%)\n", total_exterior_pixels,
         (double)total_exterior_pixels / (double)total_pixels * 100);
  printf("\n");
  printf("           Total iterations: %20llu\n", total_iter);
  printf("  Total interior iterations: %20llu (%.1f%%)\n", total_interior_iter,
         (double)total_interior_iter / (double)total_iter * 100);
  printf("  Total exterior iterations: %20llu (%.1f%%)\n", total_exterior_iter,
         (double)total_exterior_iter / (double)total_iter * 100);
  printf("\n");

  printf("         Average iterations: %20.3f\n",
         (double)total_iter / (double)total_pixels);
  printf("Average interior iterations: %20.3f\n",
         (double)total_interior_iter / (double)total_interior_pixels);
  printf("Average exterior iterations: %20.3f\n",
         (double)total_exterior_iter / (double)total_exterior_pixels);
  printf("\n");


  char buf[100];


  // Print periodic orbit table.

  printf("%25s  %20s  %9s  %9s\n",
         "Period detected after", "Quantity", "Percentage", "Percentile");
  printf("%25s  %20s  %9s  %9s\n",
         "---------------------", "--------", "----------", "----------");
  int max_interior_k = (int)(log(this->iter_max) / log(2));
  uint64 interior_tally_log2_running_total = 0;

  interior_tally_log2_running_total += interior_tally_uniterated;
  sprintf(buf, "%s %llu", "(early-out)", UINT64_C(0));
  printf("%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    interior_tally_uniterated,
    (double)interior_tally_uniterated / (double)total_interior_pixels * 100,
    (double)interior_tally_log2_running_total / (double)total_interior_pixels * 100);
  for (int k = 0; k <= max_interior_k; k++)
  {
    interior_tally_log2_running_total += interior_tally_log2[k];
    printf("%25llu  %20llu  %9.6f%% %10.6f%%\n",
      UINT64_C(1) << k,
      interior_tally_log2[k],
      (double)interior_tally_log2[k] / (double)total_interior_pixels * 100,
      (double)interior_tally_log2_running_total / (double)total_interior_pixels * 100);
  }
  interior_tally_log2_running_total += interior_tally_aperiodic;
  sprintf(buf, "%s %llu", "(never)", this->iter_max);
  printf("%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    interior_tally_aperiodic,
    (double)interior_tally_aperiodic / (double)total_interior_pixels * 100,
    (double)interior_tally_log2_running_total / (double)total_interior_pixels * 100);
  printf("\n");

  // Print escaped-point table.

  printf("%25s  %20s  %9s  %9s\n",
         "Escaped after", "Quantity", "Percentage", "Percentile");
  printf("%25s  %20s  %9s  %9s\n",
         "-------------", "--------", "----------", "----------");
  int max_exterior_k = (int)(log(this->iter_max) / log(2));
  uint64 exterior_tally_log2_running_total = 0;

  exterior_tally_log2_running_total += exterior_tally_uniterated;
  sprintf(buf, "%s %llu", "(early-out)", UINT64_C(0));
  printf("%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    exterior_tally_uniterated,
    (double)exterior_tally_uniterated / (double)total_exterior_pixels * 100,
    (double)exterior_tally_log2_running_total / (double)total_exterior_pixels * 100);
  for (int k = 0; k <= max_exterior_k; k++)
  {
    exterior_tally_log2_running_total += exterior_tally_log2[k];
    printf("%25llu  %20llu  %9.6f%% %10.6f%%\n",
      UINT64_C(1) << k,
      exterior_tally_log2[k],
      (double)exterior_tally_log2[k] / (double)total_exterior_pixels * 100,
      (double)exterior_tally_log2_running_total / (double)total_exterior_pixels * 100);
  }
  exterior_tally_log2_running_total += exterior_tally_aperiodic;
  sprintf(buf, "%s %llu", "(never)", this->iter_max);
  printf("%25s  %20llu  %9.6f%% %10.6f%%\n",
    buf,
    exterior_tally_aperiodic,
    (double)exterior_tally_aperiodic / (double)total_exterior_pixels * 100,
    (double)exterior_tally_log2_running_total / (double)total_exterior_pixels * 100);
}


