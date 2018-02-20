#include "util.h"
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

void error(const char *msg, ...)
{
  va_list args;
  va_start(args, msg);
  vfprintf(stderr, msg, args);
  abort();
}

void print_usage(const char *prog_name)
{
  printf("Usage:\n");
  printf(prog_name);
  printf(" <path to config file>\n");
}
