#ifndef SEARCH_H
# define SEARCH_H /* avoid multiple inclusion */

void qtree_free(QTREE_NODE *node);
void qtree_pop_point(DPOINT *p, DATA *d);
void qtree_push_point(DATA *d, DPOINT *p);
void qtree_rebuild(DATA *d);
int qtree_select(DPOINT *where, DATA *d);
/* 2-norm distances from point to block: */
double pb_norm_3D(const DPOINT *where, BBOX bbox);
double pb_norm_2D(const DPOINT *where, BBOX bbox);
double pb_norm_1D(const DPOINT *where, BBOX bbox);
#endif /* SEARCH_H */
