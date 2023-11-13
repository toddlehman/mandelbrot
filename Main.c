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

  if (arg_count - arg_index == 7)
  {
    real    x_center             = atof(args[arg_index+0]);
    real    y_center             = atof(args[arg_index+1]);
    real    x_size               = atof(args[arg_index+2]);
    uint64  iter_max             = atol(args[arg_index+3]);
    int     pixel_width          = atoi(args[arg_index+4]);
    int     pixel_height         = atoi(args[arg_index+5]);
    int     subsample_scale      = atoi(args[arg_index+6]);
    real    subsample_tolerance  = 0.5;

    Image *image = image_alloc(x_center, y_center, x_size, iter_max,
                               pixel_width, pixel_height,
                               subsample_scale, subsample_tolerance);
    image_populate(image);

    if (output_statistics)
      image_output_statistics(image);
    else
      image_output(image, output_text_format);

    image_dealloc(&image);
  }
  else
  {
    usage_exit(args[0]);
  }

  return 0;
}

