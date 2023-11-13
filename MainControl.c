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
    mp_set_str(*value, str, 10, MP_ROUND);
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
  if ((value = argument_matches(key, "cx=")))
  {
    mp_clear(this->x_center);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->x_center);
  }
  else if ((value = argument_matches(key, "cy=")))
  {
    mp_clear(this->y_center);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->y_center);
  }
  else if ((value = argument_matches(key, "cr=")))
  {
    mp_clear(this->xy_min_size);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->xy_min_size);
    mp_mul_d(this->xy_min_size, this->xy_min_size, 2, MP_ROUND);
  }
  else if ((value = argument_matches(key, "cd=")))
  {
    mp_clear(this->xy_min_size);  // Assumes previously initialized to default.
    valid = fetch_mp_real_value(value, &this->xy_min_size);
  }
  else if ((value = argument_matches(key, "ar=")))
  {
    valid = fetch_real_value(value, &this->roll);
  }
  else if ((value = argument_matches(key, "ap=")))
  {
    valid = fetch_real_value(value, &this->pitch);
  }
  else if ((value = argument_matches(key, "ay=")))
  {
    valid = fetch_real_value(value, &this->yaw);
  }
  else if ((value = argument_matches(key, "i=")))
  {
    valid = fetch_uint64_value(value, &this->iter_max);
  }
  else if ((value = argument_matches(key, "px=")))
  {
    valid = fetch_int_value(value, &this->width_pixels);
  }
  else if ((value = argument_matches(key, "py=")))
  {
    valid = fetch_int_value(value, &this->height_pixels);
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
  else if ((value = argument_matches(key, "-s")))
  {
    this->output_statistics = true;
    valid = true;
  }
  else if ((value = argument_matches(key, "-t")))
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

  // Set program defaults.
  mp_init2(this->x_center,    32);
  mp_init2(this->y_center,    32);
  mp_init2(this->xy_min_size, 32);
  mp_set_d(this->x_center,   -0.75, MP_ROUND);
  mp_set_d(this->y_center,    0.00, MP_ROUND);
  mp_set_d(this->xy_min_size, 2.75, MP_ROUND);
  this->roll  = 0.0;
  this->pitch = 0.0;
  this->yaw   = 0.0;
  this->iter_max = 10000;
  this->width_pixels = 8;
  this->height_pixels = 8;
  this->supersample_interior_min_depth = 0;
  this->supersample_interior_max_depth = 0;
  this->supersample_exterior_min_depth = 0;
  this->supersample_exterior_max_depth = 0;
  this->supersample_solidarity = 0;
  this->output_statistics = false;
  this->output_image_text_format = isatty(fileno(stdout));

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

  if (this->output_statistics)
  {
    fprintf(stderr, "Command line:\n");
    for (int i = 0; i < arg_count; i++)
      fprintf(stderr, "%s%c", args[i], (i < arg_count - 1)? ' ':'\n');
    fprintf(stderr, "\n");
    fflush(stderr);
  }

  Image *image = image_create(
    this->x_center,
    this->y_center,
    this->xy_min_size,
    this->width_pixels,
    this->height_pixels,
    this->supersample_interior_min_depth,
    this->supersample_interior_max_depth,
    this->supersample_exterior_min_depth,
    this->supersample_exterior_max_depth,
    this->supersample_solidarity,
    this->iter_max);

  image_populate(image);

  if (this->output_statistics)
    image_output_statistics(image, stderr);

  image_output(image, stdout, this->output_image_text_format);

  image_destroy(&image);

  mp_clear(this->x_center);
  mp_clear(this->y_center);
  mp_clear(this->xy_min_size);

  return 0;
}
