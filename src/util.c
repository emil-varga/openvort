#include "util.h"
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

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
  fputs(prog_name, stdout);
  printf(" -c | --conf <path to config file> -o | --output <output directory>\n");
}

int setup_outdir(const char *dirname)
{
  DIR *test = opendir(dirname);

  if(test){
      //the directory exist and we can open it
      return 1;
  } else {
      //either it doesn't exist or we can't open it
      //try creating it
      if(mkdir(dirname, 0777) == -1)
	return 0;
  }

  return 1;
}
