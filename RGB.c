/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "RGB.h"


//-----------------------------------------------------------------------------
// COMPUTE AVERAGE OF 2 RGB VALUES

public_function
LinearRGB linear_rgb_average2(LinearRGB color1, LinearRGB color2)
{
  return (LinearRGB)
  {
    .r = (color1.r + color2.r) / 2,
    .g = (color1.g + color2.g) / 2,
    .b = (color1.b + color2.b) / 2,
  };
}


//-----------------------------------------------------------------------------
// COMPUTE AVERAGE OF 4 RGB VALUES

public_function
LinearRGB linear_rgb_average4(LinearRGB color1, LinearRGB color2,
                              LinearRGB color3, LinearRGB color4)
{
  return (LinearRGB)
  {
    .r = (color1.r + color2.r + color3.r + color4.r) / 4,
    .g = (color1.g + color2.g + color3.g + color4.g) / 4,
    .b = (color1.b + color2.b + color3.b + color4.b) / 4,
  };
}


//-----------------------------------------------------------------------------
// COMPUTE MAXIMAL DIFFERENCE OF 2 RGB VALUES

public_function
float32 linear_rgb_diff2(LinearRGB color1, LinearRGB color2)
{
  #define SET_MAX3(x, a,b,c)  x = MAX(a,b), x = MAX(x,c);

  float32 min_r = MIN(color1.r, color2.r);
  float32 min_g = MIN(color1.g, color2.g);
  float32 min_b = MIN(color1.b, color2.b);

  float32 max_r = MAX(color1.r, color2.r);
  float32 max_g = MAX(color1.g, color2.g);
  float32 max_b = MAX(color1.b, color2.b);

  float32 diff_r = max_r - min_r;
  float32 diff_g = max_g - min_g;
  float32 diff_b = max_b - min_b;

  float32 diff; SET_MAX3(diff, diff_r, diff_g, diff_b);

  return diff;

  #undef SET_MAX3
}


//-----------------------------------------------------------------------------
// COMPUTE MAXIMAL DIFFERENCE OF 4 RGB VALUES

public_function
float32 linear_rgb_diff4(LinearRGB color1, LinearRGB color2,
                         LinearRGB color3, LinearRGB color4)
{
  #define SET_MIN4(x, a,b,c,d)  x = MIN(a,b), x = MIN(x,c), x = MIN(x,d)
  #define SET_MAX4(x, a,b,c,d)  x = MAX(a,b), x = MAX(x,c), x = MAX(x,d)

  #define SET_MAX3(x, a,b,c)  x = MAX(a,b), x = MAX(x,c);

  float32 min_r; SET_MIN4(min_r, color1.r, color2.r, color3.r, color4.r);
  float32 min_g; SET_MIN4(min_g, color1.g, color2.g, color3.g, color4.g);
  float32 min_b; SET_MIN4(min_b, color1.b, color2.b, color3.b, color4.b);

  float32 max_r; SET_MAX4(max_r, color1.r, color2.r, color3.r, color4.r);
  float32 max_g; SET_MAX4(max_g, color1.g, color2.g, color3.g, color4.g);
  float32 max_b; SET_MAX4(max_b, color1.b, color2.b, color3.b, color4.b);

  float32 diff_r = max_r - min_r;
  float32 diff_g = max_g - min_g;
  float32 diff_b = max_b - min_b;

  float32 diff; SET_MAX3(diff, diff_r, diff_g, diff_b);

  return diff;

  #undef SET_MIN4
  #undef SET_MAX4
  #undef SET_MAX3
}


//-----------------------------------------------------------------------------
