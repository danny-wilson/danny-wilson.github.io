typedef double DP;

double icdf(const double x);
void nrerror(const char *s);
DP gammln(const DP xx);
DP gammq(const DP a, const DP x);
void gcf(DP &gammcf, const DP a, const DP x, DP &gln);
void gser(DP &gamser, const DP a, const DP x, DP &gln);
