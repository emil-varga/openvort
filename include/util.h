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

struct list_elem
{
  void *data;
  struct list_elem *forward;
  struct list_elem *backward;
};

struct list
{
  struct list_elem *head;
};

struct list *new_list();
struct list *add_elem(struct list *l, void *head_data);
void *rem_elem(struct list *l);

void delete_list(struct list *l);

#endif
