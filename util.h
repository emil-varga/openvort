#ifndef UTIL_H
#define UTIL_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

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
