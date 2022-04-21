/*	combineP.exe
	This program computes a single p value for multiple, correlated, discrepancy
	(or test) statistics whose joint reference (or null) distribution is known
	empirically. Marginal p values are also produced, and there is the option to
	generate output for plotting the reference distribution surface.

	To run this program type
		combineP.exe data.txt
	where data.txt is a control file with the following lines
		m = 2									# Number of dimensions for the reference distribution
		n = 2000								# Number of empirical draws from the reference distribution
		x = distribution.txt					# A file containing the empirical draws
		o = 0.5,18.6							# OPTIONAL: observed values of the discrepancy statistics
	distribution.txt is arranged into m columns each of length n (separated by whitespace). Each row of the
	file corresponds to a vector of the discrepancies drawn jointly from the reference distribution.

	If o is omitted from data.txt then observations can be entered at run-time, and 
	a p value surface can be output for plotting purposes. There are run-time options
	to control the resolution of the surface.
*/

#include "combineP.h"
#include <myutils.h>
#include <controlwizard.h>
#include <vector.h>
#include <algorithm>
#include "svd.h"
extern "C" {
#include "mconf.h"
}
#include "pearsn.h"
#include <fstream>
#include <sstream>

extern "C" double ndtri(double y0);

using namespace myutils;
using std::vector;


int main(const int argc, const char* argv[]) {
	if(argc!=2) error("argument must be filename containing:\n\tm = number of test statistics\n\tn = number of values\n\tx = list of columns (length m x n)\n\to = list of observed (length m)");
	const char* filename = argv[1];

	/* Read in the file */
	int m,n,i,j,k;
	vector<double> x_in,X_in,o;
	string xfilename;
	ControlWizard con;
	con.add_item("m",TP_INT,&m);
	con.add_item("n",TP_INT,&n);
	con.add_item("x",TP_STRING,&xfilename);
	con.add_item("o",TP_VEC_DOUBLE,&o);
	con.read_input(filename);
//	if(x_in.size()!=(m*n)) error("x is not of length m x n");
	if(o.size()!=0 && o.size()!=m) error("o is not of length m");
	if(xfilename=="") error("x must specify the file containing a matrix for the reference distribution");

	/* Convert vector x_in to matrix x */
	vector<double> xrow((o.size()==0) ? n : n+1,0.0);
	vector< vector<double> > x(m,xrow);

	j = 0;
	ifstream in(xfilename.c_str());
	if(!in.is_open()) {
		stringstream eText;
		eText << "Could not open file " << xfilename;
		error(eText.str().c_str());
	}
	else while(true) {
		for(i=0;i<m;i++) {
			if(in.eof()) break;
			if(in.fail()) error("Cannot read in file");
			in >> x[i][j];
		}
		if(in.eof()) break;
		j++;
	}
	if(j>n) error("file size exceeds expected number of rows");
	if(j<n) error("file contains fewer than expected number of rows");
	in.close();

	if(o.size()==m) {
		for(i=0;i<m;i++)				// the last row of x contains the observed values
			x[i][n] = o[i];
		++n;
	}

	/* Convert to marginal p values */
	vector<double> index(n);
	for(i=0;i<n;i++) index[i] = i;
	vector< vector<double> > X(m,index);

	for(i=0;i<m;i++) {
		cout << "x["<<i<<"] = ";
		for(j=0;j<n;j++) cout << x[i][j] << " ";
		cout << endl;
	}	cout << endl;

	cout << "Marginal p values:" << endl;
	for(i=0;i<m;i++) sort(X[i].begin(),X[i].end(),sort_by_vector<double>(x[i]));

	vector< vector<double> > P(m,index);
	for(i=0;i<m;i++) {
		double pvalue;
		for(j=0;j<n;j++) {
			if(X[i][j]==n-1) pvalue = 1.-(double)(j+1)/(double)(n+1);
			P[i][X[i][j]] = 1.-(double)(j+1)/(double)(n+1);
			x[i][X[i][j]] = icdf(1.-(double)(j+1)/(double)(n+1));
		}
		cout << "P(x" << i << " <= " << x[i][n-1] << ") = " << 1.-pvalue << "\tp value = " << 2.0*MIN(pvalue,1.-pvalue) << endl;
	}
	cout << endl;


/*	for(i=0;i<m;i++) {
		cout << "p["<<i<<"] = ";
		for(j=0;j<n;j++) cout << P[i][j] << " ";
		cout << endl;
	}
	for(i=0;i<m;i++) {
		cout << "x["<<i<<"] = ";
		for(j=0;j<n;j++) cout << x[i][j] << " ";
		cout << endl;
	}*/


	/* Compute the correlation matrix */
	Matrix<double> Sigma(m,m,1.0);
	for(i=0;i<m;i++)
		for(j=0;j<i;j++) {
			double A,B,AA,BB,AB;
			A = B = AA = BB = AB = 0.0;
			for(k=0;k<n;k++) {
				A += x[i][k];			B += x[j][k];
				AA += x[i][k]*x[i][k];	BB += x[j][k]*x[j][k];
				AB += x[i][k] * x[j][k];
			}
			Sigma[i][j] = Sigma[j][i] = (AB - A*B/(double)(n))/sqrt((AA-A*A/(double)(n))*(BB-B*B/(double)(n)));
		}
	cout << "Correlation matrix:" << endl;
	for(i=0;i<m;i++)
		for(j=0;j<i;j++)
			cout << "rho[" << j << "," << i << "] = " << Sigma[j][i] << endl;
	cout << endl;
	
	/* Perform the singular value decomposition */
	// svdcmp(Sigma,w,v) decomposes Sigma = u . d . vT
	Vector<double> w(m,0.0);						// vector of singular values (diagonal of d)
	Matrix<double> vT(m,m);							// transpose matrix of v
	svdcmp(Sigma,w,vT);								// Sigma becomes u
	Matrix<double> d(m,m,0.0);
	for(i=0;i<m;i++) d[i][i] = w[i];
	Matrix<double> Lambda = Sigma * d.map(sqrt);
	Lambda = Lambda.invert();

	/* print out lambda(-1) */
	cout << "LambdaInv:=Matrix([";
	for(i=0;i<m;i++) {
		cout << "[";
		for(j=0;j<m;j++) cout << Lambda[i][j] << ",";
		cout << "\b],";
	}	cout << "\b]):\n" << endl;

	/* Obtain Y, the uncorrelated standard normal variates */
	cout << "Transformed uncorrelated standard normal variates:" << endl;
	vector< vector<double> > Y(m,index);
	for(j=0;j<n;j++)
		for(i=0;i<m;i++) {
			Y[i][j] = 0.0;
			for(k=0;k<m;k++) Y[i][j] += Lambda[i][k]*x[k][j];
			if(j==n-1) {
				double pvalue = gammq(1.0,pow(Y[i][j],2.0));
				cout << "Pr(y"<<i<<" <= " << Y[i][j] << ") = " << 1.-pvalue << "\tp value = " << 2.0*MIN(pvalue,1.-pvalue) << endl;
			}
		}
	cout << endl;

	/* Print out the Y's 
	for(i=0;i<m;i++) {
		cout << "Y[" << i << "]:=[";
		for(j=0;j<n;j++) cout << Y[i][j] << ",";
		cout << "\b]:\n";
	}	cout << endl;

	myutils::pause();*/

	/* Check the Y's are uncorrelated */
	Matrix<double> SigmaY(m,m,1.0);
	for(i=0;i<m;i++)
		for(j=0;j<i;j++) {
			double A,B,AA,BB,AB;
			A = B = AA = BB = AB = 0.0;
			for(k=0;k<n;k++) {
				A += Y[i][k];			B += Y[j][k];
				AA += Y[i][k]*Y[i][k];	BB += Y[j][k]*Y[j][k];
				AB += Y[i][k] * Y[j][k];
			}
			SigmaY[i][j] = SigmaY[j][i] = (AB - A*B/(double)(n))/sqrt((AA-A*A/(double)(n))*(BB-B*B/(double)(n)));
		}
	cout << "Correlation matrix for Y:" << endl;
	for(i=0;i<m;i++)
		for(j=0;j<i;j++)
			cout << "rho[" << j << "," << i << "] = " << SigmaY[j][i] << endl;
	cout << endl;

	/* Calculate the combined p value */
	double chiSq = 0.0;
	for(i=0;i<m;i++) chiSq += pow(Y[i][n-1],2.0);
	cout << "Combined p value:" << endl;
	cout << "P(Z >= " << chiSq << ") = " << gammq((double)m/2.0,chiSq/2.0) << endl;
	cout << endl;

	if(o.size()==m) {
		myutils::pause();
		return 0;
	}

	/*cout << "Writing null distribution:" << endl;
	ofstream goat("P:/goat.txt");
	goat << "a:=[";
	for(j=0;j<n;j++) {
		chiSq = 0.0;
		for(i=0;i<m;i++) chiSq += pow(Y[i][j],2.0);
		goat << chiSq << ",";
	}	goat << "NULL]:" << endl;
	goat.close();*/

	/******************************************************/
	/* Calculate p values based on reference distribution */
	/******************************************************/
	j = 0;
	ifstream in2(xfilename.c_str());
	if(!in2.is_open()) {
		stringstream eText;
		eText << "Could not open file " << xfilename;
		error(eText.str().c_str());
	}
	else while(true) {
		for(i=0;i<m;i++) {
			if(in2.eof()) break;
			if(in2.fail()) error("Cannot read in file");
			in2 >> x[i][j];
		}
		if(in2.eof()) break;
		j++;
	}
	if(j>n) error("file size exceeds expected number of rows");
	if(j<n) error("file contains fewer than expected number of rows");
	in2.close();

	cout << "Calculate p values based on reference distribution:" << endl;
	int loop,nLoop;
	nLoop = 0;

	while(true) {
		for(loop=0;loop<nLoop;loop++) {
			cout << "Enter observed test statistics..." << endl;
			vector<double> xx(m,0.0);
			vector<double> pp(m,0.0);
			for(i=0;i<m;i++) {
				cout << "x" << i << " = ";
				cin >> xx[i];
				for(j=0;j<n-1;j++)
					if(x[i][j]<xx[i]) pp[i] += 1.0;
				pp[i] = (pp[i]+1.0)/(double)(n+1);
				xx[i] = icdf(1.-pp[i]);
			}

			cout << "Marginal p values:" << endl;
			for(i=0;i<m;i++) cout << "Pr(x"<<i<<" <= "<<xx[i]<<") = " << pp[i] <<"\tp value = " << 2.0*MIN(pp[i],1.-pp[i]) << endl;

			vector<double> yy(m,0.0);
			for(i=0;i<m;i++)
				for(k=0;k<m;k++)
					yy[i] += Lambda[i][k]*xx[k];
			chiSq = 0.0;
			for(i=0;i<m;i++) chiSq += pow(yy[i],2.0);
			cout << "Combined p value: P(Z >= " << chiSq << ") = " << gammq((double)m/2.0,chiSq/2.0) << endl;
			cout << endl;
		}
		cout << "How many times> ";
		cin >> nLoop;
		if(nLoop<=0) break;
	}	cout << endl;

	/**********************************************************/
	/* Obtain p value contour plot for reference distribution */
	/**********************************************************/

	cout << "Compute a p value surface:" << endl;
	vector<double> xmin(m,0.0);
	vector<double> xmax(m,0.0);
	for(i=0;i<m;i++) {
		double XMIN,XMAX;
		XMIN = x[i][0];
		XMAX = x[i][0];
		for(j=1;j<n;j++) {
			if(x[i][j]<XMIN) XMIN = x[i][j];
			if(x[i][j]>XMAX) XMAX = x[i][j];
		}
		xmin[i] = XMIN;
		xmax[i] = XMAX;
	}
	int res;
	cout << "What resolution> ";
	cin >> res;
	
	vector<int> position(m,0);
	int update = m;
	cout << "[[";
	while(true) {
		for(i=m-1;i>0;i--)
			if(position[i]==res+1) {
				position[i] = 0;
				++position[i-1];
				cout << "\b],[";
			}
		if(position[0]==res+1) break;
		vector<double> xx(m,0.0);
		vector<double> pp(m,0.0);
		vector<double> pp2(m,0.0);
		vector<double> yy(m,0.0);
		chiSq = 0.0;
		for(i=0;i<m;i++) {
			xx[i] = xmin[i] + (double)position[i]/(double)(res)*(xmax[i]-xmin[i]);
			pp[i] = pp2[i] = 0.0;
			for(j=0;j<n;j++) {
				if(x[i][j]<xx[i]) pp[i] += 1.0;
				if(x[i][j]>xx[i]) pp2[i] += 1.0;
			}
			pp[i] = (pp[i]+1.0)/(double)(n+1);
			pp2[i] = (pp2[i]+1.0)/(double)(n+1);
			pp[i] = 1.-0.5*(pp[i] + 1.0-pp2[i]);
			pp[i] = icdf(pp[i]);
			yy[i] = 0.0;
		}
		for(i=0;i<m;i++) {
			for(k=0;k<m;k++)
				yy[i] += Lambda[i][k]*pp[k];
			chiSq += pow(yy[i],2.0);
		}
		cout << "[" << xx[0];
		for(i=1;i<m;i++) cout << "," << xx[i];
		cout << "," << gammq((double)m/2.0,chiSq/2.0) << "],";
		++position[m-1];
	}	cout << "\b\b]   " << endl;
	/*	Maple commands:
		VIEWS:=[stats[describe,range](map(x->op(map(y->y[1],x)),A)),stats[describe,range](map(x->op(map(y->y[2],x)),A)),0..0.1];
		plots[surfdata](A,axes=framed,symbol=circle,symbolsize=20,labels=[x,y,z],orientation=[-90,0],contours=20,shading=ZGRAYSCALE,style=CONTOUR,view=VIEWS);
		plots[pointplot3d]([seq(op(A[i]),i=1..nops(A))],axes=framed,symbol=circle,symbolsize=20,view=VIEWS);
	*/

	myutils::pause();

	/*************************************************************/
	/* Re-obtain p value contour plot for reference distribution */
	/*************************************************************/
	/* This is to check that the program is working correctly.	 */
	/* The reference distribution is resampled to ensure that a	 */
	/* correct p value surface is obtained.						 */
	/*************************************************************/

	/*cout << endl << "[";
	for(j=0;j<n;j++) {
		vector<double> xx(m,0.0);
		vector<double> pp(m,0.0);
		vector<double> pp2(m,0.0);
		vector<double> yy(m,0.0);
		chiSq = 0.0;
		for(i=0;i<m;i++) {
			xx[i] = x[i][j];
			pp[i] = pp2[i] = 0.0;
			for(k=0;k<n;k++) {
				if(x[i][k]<xx[i]) pp[i] += 1.0;
				if(x[i][k]>xx[i]) pp2[i] += 1.0;
			}
			pp[i] = (pp[i]+1.0)/(double)(n+1);
			pp2[i] = (pp2[i]+1.0)/(double)(n+1);
			pp[i] = 1.-0.5*(pp[i] + 1.0-pp2[i]);
			P[i][j] = pp[i];
			pp[i] = icdf(pp[i]);
			yy[i] = 0.0;
		}
		for(i=0;i<m;i++) {
			for(k=0;k<m;k++)
				yy[i] += Lambda[i][k]*pp[k];
			chiSq += pow(yy[i],2.0);
		}
		cout << "[" << xx[0];
		for(i=1;i<m;i++) cout << "," << xx[i];
		cout << "," << gammq((double)m/2.0,chiSq/2.0) << "],";
	}	cout << "\b]  " << endl;

	myutils::pause();*/
	return 0;
}

double icdf(const double x) {
	return ndtri(x);
}

//void nrerror(const char *s) {
//	myutils::error(s);
//}

