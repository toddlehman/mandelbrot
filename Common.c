/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// ABORT WITH FATAL ERROR MESSAGE

public_function
void error_exit(char *message)
{
  fprintf(stderr, "%s\n", message);
  exit(1);
}


