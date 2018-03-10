#ifndef UTIL_H
#define UTIL_H

#include <assert.h>
#include <math.h>

/*
 * Assertion that also prints a useful message.
 */
#ifndef NDEBUG
#define assert_msg(expr,msg)\
  if(!(expr)) {\
      printf(msg);\
      assert((expr));\
  }
#else
#define assert_msg(expr,msg) {}
#endif //ifndef NDEBUG

void error(const char *msg, ...);
void print_usage(const char *prog_name);
int setup_outdir(const char *dirname);

#define DEG2RAD(X) ((X)*M_PI/180.0)

#endif
