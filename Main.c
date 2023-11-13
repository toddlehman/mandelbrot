/*-----------------------------------------------------------------------------

  MANDELBROT SET IMAGE GENERATOR

  Todd S. Lehman
  September 2013

  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.

-----------------------------------------------------------------------------*/

#import "Main.h"
#import "Image.h"
#import "DeepReal.h"


//-----------------------------------------------------------------------------
// ABORT WITH USAGE MESSAGE

private_const char *program_usage_format =
   "%s: Usage goes here.\n";

private_function
void usage_exit(char *program_path)
{
  fprintf(stderr, program_usage_format, program_path);
  exit(1);
}


//-----------------------------------------------------------------------------
// CALCULATE MINIMUM BIT-PRECISION FOR NUMBER IN STRING

private_function
int mp_get_min_prec_from_string(char *str)
{
  real digit_count = 0;
  for (char *p = str; *p; p++)
    if (isdigit(*p))  // KLUDGE
      digit_count++;

  return (int)ceil(digit_count * log(10) / log(2));
}


//-----------------------------------------------------------------------------
// MAIN CONTROL

public_function
int main (int arg_count, char *args[])
{
  bool output_statistics = false;
  bool output_text_format = false;

  int arg_index;
  for (arg_index = 1; arg_index < arg_count; arg_index++)
  {
    char *arg = args[arg_index];

    if ((arg[0] == '-') && isalpha(arg[1]))
    {
      switch (arg[1])
      {
        case 's':
          output_statistics = true;
          break;

        case 't':
          output_text_format = true;
          break;

        default:
          usage_exit(args[0]);
          break;
      }
    }
    else
    {
      break;
    }
  }

  if (output_statistics)
  {
    printf("Command line:\n");
    for (int i = 0; i < arg_count; i++)
      printf("%s%c", args[i], (i < arg_count - 1)? ' ':'\n');
    printf("\n");
  }

  if (arg_count - arg_index == 8)
  {
    char *str_x_center              = args[arg_index+0];
    char *str_y_center              = args[arg_index+1];
    char *str_xy_min_size           = args[arg_index+2];
    char *str_iter_max              = args[arg_index+3];
    char *str_width_pixels          = args[arg_index+4];
    char *str_height_pixels         = args[arg_index+5];
    char *str_subsample_max_depth   = args[arg_index+6];
    char *str_subsample_solidarity  = args[arg_index+7];

    // Calculate minimum precision needed to represent values.
    int mp_prec = 0;
    mp_prec = MAX(mp_prec, 16 + mp_get_min_prec_from_string(str_x_center));
    mp_prec = MAX(mp_prec, 16 + mp_get_min_prec_from_string(str_y_center));
    mp_prec = MAX(mp_prec, 16 + mp_get_min_prec_from_string(str_xy_min_size));

    mp_real x_center, y_center, xy_min_size;
    mp_init2(x_center, mp_prec);
    mp_init2(y_center, mp_prec);
    mp_init2(xy_min_size, mp_prec);

    mp_set_str(x_center, str_x_center, 10, MP_ROUND);
    mp_set_str(y_center, str_y_center, 10, MP_ROUND);
    mp_set_str(xy_min_size, str_xy_min_size, 10, MP_ROUND);

    uint64  iter_max              = atol(str_iter_max);
    int     width_pixels          = atoi(str_width_pixels);
    int     height_pixels         = atoi(str_height_pixels);
    int     subsample_min_depth   = 0;
    int     subsample_max_depth   = atoi(str_subsample_max_depth);
    real    subsample_solidarity  = atof(str_subsample_solidarity);

    Image *image = image_create(x_center, y_center, xy_min_size,
                                width_pixels, height_pixels,
                                subsample_min_depth, subsample_max_depth,
                                subsample_solidarity,
                                iter_max);

    mp_clear(x_center);
    mp_clear(y_center);
    mp_clear(xy_min_size);

    image_populate(image);

    if (output_statistics)
      image_output_statistics(image);
    else
      image_output(image, output_text_format);

    image_destroy(&image);
  }
  else
  {
    usage_exit(args[0]);
  }

  return 0;
}

