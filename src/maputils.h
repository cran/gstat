#ifndef MAPUTILS_H

void map_name_nr(GRIDMAP *mask, const char *base, char *name, int nr, int max);
int map_equal(GRIDMAP *a, GRIDMAP *b);
int map_convert(int argc, char *argv[]);
int map_cover(int argc, char *argv[]);
int map_nominal(int argc, char *argv[]);
int map_cut(int argc, char *argv[]);
int map_diff(int argc, char *argv[]);
int map_lnh(int argc, char *argv[]);
int map_q(int argc, char *argv[]);

#endif
