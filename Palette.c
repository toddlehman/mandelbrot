/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Palette.h"


//-----------------------------------------------------------------------------
// ALLOCATE COLOR PALETTE

public_constructor
Palette *palette_create(void)
{
  Palette *this = mem_alloc_clear(1, sizeof(*this));

  const int max_exterior_count = 20;

  this->exterior_locations = mem_alloc_clear(max_exterior_count + 1,
    sizeof(*(this->exterior_locations)));

  this->exterior_colors = mem_alloc_clear(max_exterior_count + 1,
    sizeof(*(this->exterior_colors)));

  int i = 0;
#if 0
  // Simple, excellent monochrome palette showing lots and lots of detail.
  // A bit boring, though, since there's no color.
  this->exterior_locations[i] = 0.00;
  this->exterior_colors[i++] = (LinearRGB) { 0.00, 0.00, 0.00 };  // Black
  this->exterior_locations[i] = 0.50;
  this->exterior_colors[i++] = (LinearRGB) { 1.00, 1.00, 1.00 };  // White
#elif 0
  // Very nice, but too fiery.  I mean, it's really cool, but the fire gets
  // tiring after a while.  Not enough variety in the photos.
  this->exterior_locations[i] = 0.00;
  this->exterior_colors[i++] = (LinearRGB) { 0.00, 0.00, 0.00 };  // Black
  this->exterior_locations[i] = 0.17;
  this->exterior_colors[i++] = (LinearRGB) { 0.90, 0.00, 0.00 };  // Red
  this->exterior_locations[i] = 0.33;
  this->exterior_colors[i++] = (LinearRGB) { 0.95, 0.80, 0.00 };  // Yellow
  this->exterior_locations[i] = 0.50;
  this->exterior_colors[i++] = (LinearRGB) { 1.00, 1.00, 1.00 };  // White
  this->exterior_locations[i] = 0.67;
  this->exterior_colors[i++] = (LinearRGB) { 0.95, 0.80, 0.00 };  // Yellow
  this->exterior_locations[i] = 0.83;
  this->exterior_colors[i++] = (LinearRGB) { 0.90, 0.00, 0.00 };  // Red
#elif 1
  // Really great all-purpose palette.  Finely tuned.
  this->exterior_locations[i] = 0.00;
  this->exterior_colors[i++] = (LinearRGB) { 0.00, 0.15, 0.85 };  // Blue
  this->exterior_locations[i] = 0.10;
  this->exterior_colors[i++] = (LinearRGB) { 0.30, 0.00, 0.40 };  // Purple
  this->exterior_locations[i] = 0.20;
  this->exterior_colors[i++] = (LinearRGB) { 0.90, 0.00, 0.00 };  // Red
  this->exterior_locations[i] = 0.50;
  this->exterior_colors[i++] = (LinearRGB) { 0.95, 0.80, 0.00 };  // Yellow
  this->exterior_locations[i] = 0.80;
  this->exterior_colors[i++] = (LinearRGB) { 0.00, 0.70, 0.02 };  // Green
  this->exterior_locations[i] = 0.90;
  this->exterior_colors[i++] = (LinearRGB) { 0.00, 0.40, 0.60 };  // Teal
#elif 1
  // Great all-purpose palette.
  this->exterior_locations[i] = 0.00;
  this->exterior_colors[i++] = (LinearRGB) { 0.00, 0.15, 0.85 };  // Blue
  this->exterior_locations[i] = 0.14;
  this->exterior_colors[i++] = (LinearRGB) { 0.30, 0.00, 0.80 };  // Purple
  this->exterior_locations[i] = 0.28;
  this->exterior_colors[i++] = (LinearRGB) { 0.90, 0.00, 0.00 };  // Red
  this->exterior_locations[i] = 0.43;
  this->exterior_colors[i++] = (LinearRGB) { 1.00, 0.25, 0.00 };  // Orange1
  this->exterior_locations[i] = 0.57;
  this->exterior_colors[i++] = (LinearRGB) { 1.00, 0.60, 0.00 };  // Orange2
  this->exterior_locations[i] = 0.71;
  this->exterior_colors[i++] = (LinearRGB) { 0.95, 0.80, 0.00 };  // Yellow
  this->exterior_locations[i] = 0.85;
  this->exterior_colors[i++] = (LinearRGB) { 0.00, 0.70, 0.02 };  // Green
#elif 0
  // Good all-purpose palette.  What I don't like about it, though, is that it
  // doesn't follow a sunset color pattern (blue->purple->red->orange->yellow)
  // and instead goes the other way (blue->green->yellow->orange->red->purple).
  // Also, note that this hasn't been fine-tuned much.
  this->exterior_locations[i] = 0.000;
  this->exterior_colors[i++] = (LinearRGB) { 0.000, 0.145, 0.855 };  // Blue
  this->exterior_locations[i] = 0.125;
  this->exterior_colors[i++] = (LinearRGB) { 0.000, 0.682, 0.018 };  // Green
  this->exterior_locations[i] = 0.250;
  this->exterior_colors[i++] = (LinearRGB) { 0.950, 0.812, 0.000 };  // Yellow
  this->exterior_locations[i] = 0.375;
  this->exterior_colors[i++] = (LinearRGB) { 1.000, 0.348, 0.000 };  // Orange
  this->exterior_locations[i] = 0.500;
  this->exterior_colors[i++] = (LinearRGB) { 0.900, 0.000, 0.000 };  // Red
  this->exterior_locations[i] = 0.625;
  this->exterior_colors[i++] = (LinearRGB) { 0.655, 0.000, 0.245 };  // Mauve
  this->exterior_locations[i] = 0.750;
  this->exterior_colors[i++] = (LinearRGB) { 0.462, 0.000, 0.538 };  // Purple
  this->exterior_locations[i] = 0.825;
  this->exterior_colors[i++] = (LinearRGB) { 0.147, 0.000, 0.853 };  // Indigo
#else
  // Variation on the original.  It's nice, but it's too damn bright and neon.
  this->exterior_locations[i] = 0.000;
  this->exterior_colors[i++] = (LinearRGB) { 0.000, 0.000, 1.000 };  // Blue
  this->exterior_locations[i] = 0.125;
  this->exterior_colors[i++] = (LinearRGB) { 0.000, 0.700, 0.000 };  // Green
  this->exterior_locations[i] = 0.250;
  this->exterior_colors[i++] = (LinearRGB) { 1.000, 1.000, 0.000 };  // Yellow
  this->exterior_locations[i] = 0.375;
  this->exterior_colors[i++] = (LinearRGB) { 1.000, 0.500, 0.000 };  // Orange
  this->exterior_locations[i] = 0.500;
  this->exterior_colors[i++] = (LinearRGB) { 1.000, 0.000, 0.000 };  // Red
  this->exterior_locations[i] = 0.625;
  this->exterior_colors[i++] = (LinearRGB) { 1.000, 0.000, 0.333 };  // Mauve
  this->exterior_locations[i] = 0.750;
  this->exterior_colors[i++] = (LinearRGB) { 0.666, 0.000, 0.666 };  // Purple
  this->exterior_locations[i] = 0.825;
  this->exterior_colors[i++] = (LinearRGB) { 0.333, 0.000, 1.000 };  // Indigo
#endif
  this->exterior_count = i;
  assert(this->exterior_count < max_exterior_count);

  //for (int i = 0; i < this->exterior_count; i++)
  //  this->exterior_locations[i] = (real)i / (real)this->exterior_count;
  this->exterior_colors[i] = this->exterior_colors[0];
  this->exterior_locations[i] = 1.0;

  this->undefined_color =
    (LinearRGB) { 1.0, 1.0, 1.0 };  // White

  this->interior_iterated_periodic_color =
    (LinearRGB) { 0.0, 0.0, 0.0 };  // Black

  this->interior_iterated_aperiodic_color =
    (LinearRGB) { 0.0, 0.0, 0.0 };  // Black (for production)
    //(LinearRGB) { 0.5, 0.5, 0.5 };  // Gray (for debugging)

  this->interior_uniterated_color =
    (LinearRGB) { 0.0 ,0.0 ,0.0 };  // Black (for production)
    //(LinearRGB) { 0.02,0.02,0.02};  // Dark gray (for debugging)

  #if 0
  printf("Color palette:\n");
  printf("%d colors (+1 extra padding on the end)\n", this->exterior_count);
  for (int i = 0; i < this->exterior_count + 1; i++)
  {
    printf("%2d.  %8.6f  (%5.3f, %5.3f, %5.3f)\n",
           i,
           this->exterior_locations[i],
           this->exterior_colors[i].r,
           this->exterior_colors[i].g,
           this->exterior_colors[i].b);
  }
  #endif

  return this;
}


//-----------------------------------------------------------------------------
// DEALLOCATE PALETTE

public_destructor
void palette_destroy(Palette **p_this)
{
  assert(p_this);
  Palette *this = *p_this;

  assert(this);
  assert(this->exterior_locations);
  assert(this->exterior_colors);

  mem_dealloc(&this->exterior_locations);
  mem_dealloc(&this->exterior_colors);

  mem_dealloc(p_this);
}


//-----------------------------------------------------------------------------
// COMPUTE COLOR

private_method
LinearRGB palette_compute_color(Palette *this, real location)
{
  assert(this);
  //assert((location >= 0) && (location <= 1));

  //     if (location < 0) location = 0;
  //else if (location > 1) location = 1;

  if (location > 1) location -= floor(location);

  if (location < 0)
  {
    if (location < -1) location = -1;
    location += 1;
    location = swerp(location, 0, 1);
    return linear_rgb_lerp(location,
                           (LinearRGB) { 0, 0, 0 },
                           this->exterior_colors[0]); 
  }

  for (int i = 0; i < this->exterior_count; i++)
  {
    if ((location >= this->exterior_locations[i]) &&
        (location <  this->exterior_locations[i+1]))
    {
      location = unlerp(location,
                        this->exterior_locations[i],
                        this->exterior_locations[i+1]);

      return linear_rgb_lerp(location,
                             this->exterior_colors[i],
                             this->exterior_colors[i+1]);
    }
  }

  return this->undefined_color;
}


//-----------------------------------------------------------------------------
// MAP MANDELBROT DWELL TO COLOR LOCATION

private_method
real palette_map_dwell_to_color_location(Palette *this, float64 dwell)
{
  float64 f = dwell;

  // **** TODO **** Try f = atan(dwell) / (PI/2) or something involving that.

  #if 0

    f = fmod(f, 200.0) / 200.0;

  #elif 1

    //if (f < 0) f = 0;
    //assert(f >= 0);
    if (f < 0)
    {
      ;
    }
    else if (f < 200)
    {
      f = f / 200;
    }
    else
    {
      f -= 200;
      f = pow(f, 0.5);
      f = fmod(f, 100.0) / 100.0;
    }

  #elif 0

    f = log(f) / log(1.5); f = fmod(f, 8) / 8;

  #elif 0

    f = pow(f, 0.25); f = fmod(f, 2) / 2;

  #endif

  return f;
}


//-----------------------------------------------------------------------------
// COMPUTE COLOR FROM MANDELBROT RESULT

public_method
LinearRGB palette_color_from_mandelbrot_result(Palette *this,
                                               MandelbrotResult mr)
{
  assert(this);

  if (mandelbrot_result_is_defined(mr))
  {
    if (mandelbrot_result_is_exterior(mr))
    {
      float64 location = palette_map_dwell_to_color_location(this, mr.dwell);
      return palette_compute_color(this, location);
    }
    else if (mandelbrot_result_is_interior(mr))
    {
      if (mandelbrot_result_is_interior_iterated(mr))
      {
        if (mandelbrot_result_is_interior_aperiodic(mr))
        {
          return this->interior_iterated_aperiodic_color;
        }
        else
        {
          return this->interior_iterated_periodic_color;
        }
      }
      else if (mandelbrot_result_is_interior_uniterated(mr))
      {
        return this->interior_uniterated_color;
      }
    }
  }

  return this->undefined_color;
}


//-----------------------------------------------------------------------------
