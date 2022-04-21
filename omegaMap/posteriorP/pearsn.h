#include <myerror.h>
#include <cmath>
#include <limits>
using namespace myutils;
using namespace std;

typedef double DP;

DP betacf(const DP a, const DP b, const DP x)
{
	const int MAXIT=100;
	const DP EPS=numeric_limits<DP>::epsilon();
	const DP FPMIN=numeric_limits<DP>::min()/EPS;
	int m,m2;
	DP aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (m > MAXIT) error("a or b too big, or MAXIT too small in betacf");
	return h;
}

DP betai(const DP a, const DP b, const DP x)
{
	DP bt;

	if (x < 0.0 || x > 1.0) error("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

void pearsn(vector<double> &x, vector<double> &y, DP &r, DP &prob, DP &z)
{
	const DP TINY=1.0e-20;
	int j;
	DP yt,xt,t,df;
	DP syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;

	int n=x.size();
	for (j=0;j<n;j++) {
		ax += x[j];
		ay += y[j];
	}
	ax /= n;
	ay /= n;
	for (j=0;j<n;j++) {
		xt=x[j]-ax;
		yt=y[j]-ay;
		sxx += xt*xt;
		syy += yt*yt;
		sxy += xt*yt;
	}
	r=sxy/(sqrt(sxx*syy)+TINY);
	z=0.5*log((1.0+r+TINY)/(1.0-r+TINY));
	df=n-2;
	t=r*sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)));
	prob=betai(0.5*df,0.5,df/(df+t*t));
	// prob=erfcc(fabs(z*sqrt(n-1.0))/1.4142136);
}

