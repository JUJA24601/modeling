#include <stdlib.h>
#include "oh_part.h"

static int psort_comp(const void* x, const void* y);

void particle_sort_(struct S_particle *p, int *n) {
  qsort(p, *n, sizeof(struct S_particle), psort_comp);
}
static int psort_comp(const void* pa, const void* pb) {
  struct S_particle *a = (struct S_particle*)pa;
  struct S_particle *b = (struct S_particle*)pb;

  if ((int) a->z < (int) b->z) return(-1);
  if ((int) a->z > (int) b->z) return(1);
  if ((int) a->y < (int) b->y) return(-1);
  if ((int) a->y > (int) b->y) return(1);
  if (a->x < b->x) return(-1);
  if (a->x > b->x) return(1);
  return(0);
}
