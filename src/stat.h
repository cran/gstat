/* stat.c */
double sample_mean(double *list, int n);
double sample_var(double *list, double mean, int n);
double sample_std(double *list, double mean, int n);
double est_quant(double *list, double p, int n);
void calc_r(double *a, double *b, int n, double *r);
int calc_stats(int argc, char *argv[]);
int d_cmp(const double *a, const double *b);
