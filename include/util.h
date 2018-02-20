#ifndef UTIL_H
#define UTIL_H

#include <assert.h>

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

#endif
