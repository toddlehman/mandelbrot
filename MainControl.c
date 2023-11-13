/*-----------------------------------------------------------------------------

  MANDELBROT SET IMAGE GENERATOR

  Todd S. Lehman
  September 2013

  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.

-----------------------------------------------------------------------------*/

#import "Common.h"
#import "Memory.h"
#import "DeepReal.h"
#import "MPReal.h"
#import "RGB.h"
#import "Palette.h"
#import "Mandelbrot.h"
#import "Pixel.h"
#import "Image.h"
#import "MainControl.h"


//-----------------------------------------------------------------------------
// ABORT WITH USAGE MESSAGE

private_const char *program_usage_format =
   "%s: Usage goes here.\n";

private_function
void usage_exit(const char *program_path)
{
  fprintf(stderr, program_usage_format, program_path);
  exit(1);
}


//-----------------------------------------------------------------------------
// CONVERT DEGREES TO RADIANS

private_function
real degrees_to_radians(const real degrees)
{
  real radians = degrees / 360 * TAU;
  return radians;
}


//-----------------------------------------------------------------------------
// VALIDATE STRING NUMERIC VALUE

private_function
bool validate_numeric_value(const char *str, bool allow_radix)
{
  assert(str);

  // Validate syntax of numeric value.  Value must have zero or one leading
  // plus or minus signs, followed by one or more base-10 digits.  If the type
  // is non-integer, up to one radix point is allowed.
  int sign = 0, radix = 0, digits = 0;
  for (const char *p = str; *p; p++)
  {
         if (*p == '-')        { if (sign++ || digits) return false; }
    else if (*p == '+')        { if (sign++ || digits) return false; }
    else if (*p == '.')        { if (radix++) return false; }
    else if (isdigit(*p))      { digits++; }
    else                       { return false; }
  }

  if (!digits)                 { return false; }
  if (radix && !allow_radix)   { return false; }

  return true;
}


//-----------------------------------------------------------------------------
// PARSE AND STORE HIGH-PRECISION NUMERIC VALUE

private_function
bool fetch_mp_real_value(const char *str, mp_real *value)
{
  assert(str);

  if (validate_numeric_value(str, true))
  {
    mp_init2(*value, mp_get_min_prec_from_string(str));
    mp_set_str(*value, str, 10);
    return true;
  }
  else
  {
    return false;
  }
}


//-----------------------------------------------------------------------------
// PARSE AND STORE LOW-PRECISION NUMERIC VALUE

private_function
bool fetch_real_value(const char *str, real *value)
{
  assert(str);

  if (validate_numeric_value(str, true))
  {
    *value = strtold(str, NULL);
    return true;
  }
  else
  {
    return false;
  }
}


//-----------------------------------------------------------------------------
// PARSE AND STORE UINT64 NUMERIC VALUE

private_function
bool fetch_uint64_value(const char *str, uint64 *value)
{
  assert(str);

  if (validate_numeric_value(str, false))
  {
    *value = strtoll(str, NULL, 10);
    return true;
  }
  else
  {
    return false;
  }
}


//-----------------------------------------------------------------------------
// PARSE AND STORE INT NUMERIC VALUE

private_function
bool fetch_int_value(const char *str, int *value)
{
  assert(str);

  if (validate_numeric_value(str, false))
  {
    *value = atoi(str);
    return true;
  }
  else
  {
    return false;
  }
}


//-----------------------------------------------------------------------------
// COMPARE STRING PREFIX FOR MATCH
private_function
const char *argument_matches(const char *str, const char *prefix)
{
  if (strncmp(str, prefix, strlen(prefix)) == 0)
    return str + strlen(prefix);
  else
    return NULL;
}


//-----------------------------------------------------------------------------
// PROCESS COMMAND-LINE ARGUMENT

private_method
bool process_argument(MainControl *this, const char *str)
{
  assert(this); assert(str);

  const char *key = str;
  const char *value = NULL;

  bool valid = false;
  if ((value = argument_matches(key, "jx=")))
  {
    mp_clear(this->julia_x);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->julia_x);
  }
  else if ((value = argument_matches(key, "jy=")))
  {
    mp_clear(this->julia_y);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->julia_y);
  }
  else if ((value = argument_matches(key, "x=")))
  {
    mp_clear(this->target_x);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->target_x);
  }
  else if ((value = argument_matches(key, "y=")))
  {
    mp_clear(this->target_y);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->target_y);
  }
  else if ((value = argument_matches(key, "cd=")))
  {
    valid = fetch_real_value(value, &this->target_camera_rho);
  }
  else if ((value = argument_matches(key, "ct=")))
  {
    valid = fetch_real_value(value, &this->target_camera_theta);
    this->target_camera_theta = degrees_to_radians(this->target_camera_theta);
  }
  else if ((value = argument_matches(key, "cp=")))
  {
    valid = fetch_real_value(value, &this->target_camera_phi);
    this->target_camera_phi = degrees_to_radians(this->target_camera_phi);
  }
  else if ((value = argument_matches(key, "vt=")))
  {
    valid = fetch_real_value(value, &this->camera_theta);
    this->camera_theta = degrees_to_radians(this->camera_theta);
  }
  else if ((value = argument_matches(key, "vp=")))
  {
    valid = fetch_real_value(value, &this->camera_phi);
    this->camera_phi = degrees_to_radians(this->camera_phi);
  }
  else if ((value = argument_matches(key, "vr=")))
  {
    valid = fetch_real_value(value, &this->camera_roll);
    this->camera_roll = degrees_to_radians(this->camera_roll);
  }
  else if ((value = argument_matches(key, "fov=")))
  {
    valid = fetch_real_value(value, &this->camera_fov);
    this->camera_fov = degrees_to_radians(this->camera_fov);
  }
  else if ((value = argument_matches(key, "n=")))
  {
    valid = fetch_uint64_value(value, &this->iter_max);
  }
  else if ((value = argument_matches(key, "di=")))
  {
    valid = fetch_int_value(value, &this->height_pixels);
  }
  else if ((value = argument_matches(key, "dj=")))
  {
    valid = fetch_int_value(value, &this->width_pixels);
  }
  else if ((value = argument_matches(key, "dij=")))
  {
    valid = fetch_int_value(value, &this->height_pixels);
    this->width_pixels = this->height_pixels;
  }
  else if ((value = argument_matches(key, "ssimin=")))
  {
    valid = fetch_int_value(value, &this->supersample_interior_min_depth);
  }
  else if ((value = argument_matches(key, "ssimax=")))
  {
    valid = fetch_int_value(value, &this->supersample_interior_max_depth);
  }
  else if ((value = argument_matches(key, "ssemin=")))
  {
    valid = fetch_int_value(value, &this->supersample_exterior_min_depth);
  }
  else if ((value = argument_matches(key, "ssemax=")))
  {
    valid = fetch_int_value(value, &this->supersample_exterior_max_depth);
  }
  else if ((value = argument_matches(key, "sss=")))
  {
    real supersample_solidarity;
    valid = fetch_real_value(value, &supersample_solidarity);
    this->supersample_solidarity = (float32)supersample_solidarity;
  }
  else if ((value = argument_matches(key, "ss=")))
  {
    // "ss=" is a so-called "convenience" parameter.
    int supersample_max_depth;
    valid = fetch_int_value(value, &supersample_max_depth);
    this->supersample_interior_min_depth = supersample_max_depth;
    this->supersample_interior_max_depth = supersample_max_depth;
    this->supersample_exterior_min_depth = 0;
    this->supersample_exterior_max_depth = supersample_max_depth;
    if ((supersample_max_depth > 0) && (this->supersample_solidarity == 0))
      this->supersample_solidarity = 0.999;  // Default value.
  }
  else if ((value = argument_matches(key, "--stats")))
  {
    this->output_statistics = true;
    valid = true;
  }
  else if ((value = argument_matches(key, "--text")))
  {
    this->output_image_text_format = true;
    valid = true;
  }

  return valid;
}


//-----------------------------------------------------------------------------
// PRINT STRUCTURE SIZES FOR DEBUGGING

void print_struct_sizes(void)
{
  #define PRINT_STRUCT_SIZE(x) \
    fprintf(stderr, "sizeof(%s) = %d\n", #x, (int)sizeof(x));

  PRINT_STRUCT_SIZE(mp_real);
  PRINT_STRUCT_SIZE(DeepReal);
  PRINT_STRUCT_SIZE(LinearRGB);
  PRINT_STRUCT_SIZE(DeviceRGB48);
  PRINT_STRUCT_SIZE(DeviceRGB24);
  PRINT_STRUCT_SIZE(Palette);
  PRINT_STRUCT_SIZE(Mandelbrot);
  PRINT_STRUCT_SIZE(MandelbrotResult);
  PRINT_STRUCT_SIZE(Pixel);
  PRINT_STRUCT_SIZE(Image);
  PRINT_STRUCT_SIZE(MainControl);

  #undef PRINT_STRUCT_SIZE
}


//-----------------------------------------------------------------------------
// MAIN CONTROL

public_function
int main (int arg_count, const char *args[])
{
  MainControl _this, *this = &_this;


  // --- Set program defaults.

  mp_init2(this->julia_x,  32);
  mp_init2(this->julia_y,  32);
  mp_init2(this->target_x, 32);
  mp_init2(this->target_y, 32);
  mp_set_d(this->julia_x,     0.00);
  mp_set_d(this->julia_y,     0.00);
  mp_set_d(this->target_x,   -0.75);
  mp_set_d(this->target_y,    0.00);
  this->target_camera_rho   = 2.0;
  this->target_camera_theta = degrees_to_radians(0);
  this->target_camera_phi   = degrees_to_radians(0);
  this->camera_theta        = degrees_to_radians(0);
  this->camera_phi          = degrees_to_radians(0);
  this->camera_roll         = degrees_to_radians(0);
  this->camera_fov          = degrees_to_radians(90);
  this->iter_max = 10000;
  this->width_pixels  = 8;
  this->height_pixels = 8;
  this->supersample_interior_min_depth = 0;
  this->supersample_interior_max_depth = 0;
  this->supersample_exterior_min_depth = 0;
  this->supersample_exterior_max_depth = 0;
  this->supersample_solidarity = 0;
  this->output_statistics = false;
  this->output_image_text_format = isatty(fileno(stdout));


  // --- Parse command line arguments.

  bool valid = true;
  for (int arg_index = 1; arg_index < arg_count; arg_index++)
  {
    const char *arg = args[arg_index];
    if (!process_argument(this, arg))
    {
      valid = false;
      fprintf(stderr, "%s: Invalid argument\n", arg);
    }
  }
  if (!valid)
    usage_exit(args[0]);

  if (!(this->supersample_interior_min_depth <=
        this->supersample_interior_max_depth))
    error_exit("Supersampling interior minimum depth must be "
               "less than or equal to maximum depth.");
  if (!(this->supersample_exterior_min_depth <=
        this->supersample_exterior_max_depth))
    error_exit("Supersampling exterior minimum depth must be "
               "less than or equal to maximum depth.");


  // --- Dump initial statistics text.

  if (this->output_statistics)
  {
    fprintf(stderr, "Command line:\n");
    for (int i = 0; i < arg_count; i++)
      fprintf(stderr, "%s%c", args[i], (i < arg_count - 1)? ' ':'\n');
    fprintf(stderr, "\n");
    fflush(stderr);
  }


  // --- Create image.

  Palette *palette = palette_create();

  Camera *camera = camera_create(
    this->target_x,
    this->target_y,
    this->target_camera_rho,
    this->target_camera_theta,
    this->target_camera_phi,
    this->camera_theta,
    this->camera_phi,
    this->camera_roll,
    this->camera_fov);

  Image *image = image_create(
    this->width_pixels,
    this->height_pixels,
    this->supersample_interior_min_depth,
    this->supersample_interior_max_depth,
    this->supersample_exterior_min_depth,
    this->supersample_exterior_max_depth,
    this->supersample_solidarity,
    this->iter_max,
    palette,
    camera);

  image_populate(image);


  // --- Output statistics and output image.

  if (this->output_statistics)
    image_output_statistics(image, stderr);

  image_output(image, stdout, this->output_image_text_format);


  // --- Clean up.

  image_destroy(&image);

  camera_destroy(&camera);

  palette_destroy(&palette);

  mp_clear(this->julia_x);
  mp_clear(this->julia_y);
  mp_clear(this->target_x);
  mp_clear(this->target_y);

  return 0;
}

