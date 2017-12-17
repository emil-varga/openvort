#include "util.h"
#include <stdlib.h>

struct list *new_list()
{
  struct list *l = (struct list*)malloc(sizeof(struct list));
  l->head = NULL;

  return l;
}

void delete_list(struct list *l)
{
  while(l->head)
    rem_elem(l);

  free(l->head);
  free(l);
}

void *rem_elem(struct list *l)
{
  if(l->head == NULL)
    return NULL;
  
  void *data = l->head->data;

  struct list_elem *old_head = l->head;

  l->head = old_head->backward;

  if(l->head)
    l->head->forward = NULL;

  free(old_head);

  return data;
}

struct list *add_elem(struct list *l, void *data)
{
  struct list_elem *new_elem = (struct list_elem*)malloc(sizeof(struct list_elem));

  new_elem->forward = NULL;
  new_elem->backward = l->head;
  new_elem->data = data;

  if(l->head)
    l->head->forward = new_elem;

  l->head=new_elem;

  return l;
}
