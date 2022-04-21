/************************************/
/*	main.cpp 23rd February 2005		*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson. Please don't	*/
/*	re-distribute this copy!		*/
/*	www.danielwilson.me.uk			*/
/************************************/

#include <myutils.h>
#include <decode.h>
#include <fstream>
#include <vector>
#include <controlwizard.h>
#include <coalesce.h>
#include <lissim.h>
#include <coalesce.h>

using namespace myutils;
using namespace std;

class predP {
public:
	//enum STATISTICS {S,H,pi,Varpi,Rm,Rsq,Dprime,G4,UNKNOWN};
	enum STATISTICS {UNKNOWN=-1,S,H,pi,Varpi,Rm,Rsq,Dprime,G4,U,D};

	omegaMapAnalyse *oM;
	int burnin;
	int thinning;
	int nsim;
	ofstream out;
	vector<STATISTICS> statistics;

	static STATISTICS stringToStatistics(const string st) {
		string tr = st;
		int i;
		for(i=0;i<tr.length();i++) tr[i] = tolower(tr[i]);
		if(tr=="s") return S;
		if(tr=="h") return H;
		if(tr=="pi") return pi;
		if(tr=="varpi") return Varpi;
		if(tr=="rm") return Rm;
		if(tr=="rsq") return Rsq;
		if(tr=="dprime") return Dprime;
		if(tr=="g4") return G4;
		if(tr=="u") return U;
		if(tr=="d") return D;
		return UNKNOWN;
	}
	static string statisticsToString(const STATISTICS st) {
		switch(st) {
		case S:		return "S";
		case H:		return "H";
		case pi:	return "pi";
		case Varpi:	return "Varpi";
		case Rm:	return "Rm";
		case Rsq:	return "Rsq";
		case Dprime:return "Dprime";
		case G4:	return "G4";
		case U:		return "U";
		case D:		return "D";
		}
		return "UNKNOWN";
	}
};

vector<double> calcStats(vector<bool> &calculate, vector< vector<int> > &codon);
double S(vector< vector<int> > &H);
double H(vector< vector<int> > &H);
double pi(vector< vector<int> > &H);
double Varpi(vector< vector<int> > &H);
double Rm(vector< vector<int> > &H);
vector<double> RecCorrelations(vector< vector<int> > &H);
double U(vector< vector<int> > &H);
double D(vector< vector<int> > &H);
//double Rsq(vector< vector<int> > &H);
//double Dprime(vector< vector<int> > &H);
//double G4(vector< vector<int> > &H);
void dnaToCodons(DNA &dna, vector< vector<int> > &codon);
int tripletToCodon(string &tri);
int removeStopCodons(const int a);
string codonToTriplet(const int a);

vector<double> siteS;
vector<double> realS;
NY98_61* changeType;
vector<double> siteT;
vector<double> realT;

int main(int argc, char* argv[]) {
	if(argc!=2) error("SYNTAX: control-file");

	Random ran;
//	ran.setseed(-180247090);
	changeType = new NY98_61(&ran);
	changeType->state_freq = vector<double>(61,1.0);
	changeType->build_C(1.0,2.0,4.0,changeType->state_freq);
//	changeType->initialize(61,&(changeType->C));


	string FASTAfile, OUTfile;
	predP argP;
	vector<string> stats,MCMCfile;
	ControlWizard con;
	con.add_ITEM("FASTA",TP_STRING,&FASTAfile);
	con.add_ITEM("MCMC",TP_VEC_STRING,&MCMCfile);
	con.add_ITEM("burnin",TP_INT,&argP.burnin);
	con.add_ITEM("thinning",TP_INT,&argP.thinning);
	con.add_ITEM("nsim",TP_INT,&argP.nsim);
	con.add_ITEM("outfile",TP_STRING,&OUTfile);
	con.add_ITEM("statistics",TP_VEC_STRING,&stats);
	con.read_input(argv[1]);
	if(con.got_required==false) error("Not all fields were found");

	DNA FASTA(FASTAfile.c_str());
	vector< vector<int> > codon;
	dnaToCodons(FASTA,codon);
	siteS = vector<double>(codon[0].size(),0.0);
	siteT = siteS;

	argP.out.open(OUTfile.c_str());
	if(argP.out.is_open()!=true) error("Could not open outfile");

	int i,j;
	argP.statistics.resize(stats.size());
	vector<bool> calculate(10,false);
	for(i=0;i<stats.size();i++) {
		argP.statistics[i] = predP::stringToStatistics(stats[i]);
		if(argP.statistics[i]==predP::UNKNOWN) error("Unknown statistic");
		calculate[(int)(argP.statistics[i])] = true;
	}

	vector<double> observed = calcStats(calculate,codon);

	for(i=0;i<argP.statistics.size();i++)
		cout << predP::statisticsToString(argP.statistics[i]) << "\t=\t" << observed[(int)(argP.statistics[i])] << endl;
	cout << "Sitewise S:" << endl;
	for(i=0;i<siteS.size();i++) cout << siteS[i] << ",";
	cout << "\b \n";
	realS = siteS;
	realT = siteT;

	//	argP.out << "theta" << "\t";
	//for(i=0;i<oM.L;i++) argP.out << "T" << i << "\t";
	for(i=0;i<argP.statistics.size();i++) {
		argP.out << predP::statisticsToString(argP.statistics[i]);
		if(i==argP.statistics.size()-1) argP.out << endl;
		else argP.out << "\t";
	}

	int mcmc;
	for(mcmc=0;mcmc<MCMCfile.size();mcmc++) {
		omegaMapAnalyse oM;
		argP.oM = &oM;
		cout << "Reading " << MCMCfile[mcmc] << "..." << endl;
		oM.open(MCMCfile[mcmc].c_str());
		Vector<int> blankVectorInt(0);
		oM.initialize((const char*)"",blankVectorInt,1);

		LiSsimulator LS(ran);
		Control ctrl;
	//	coalescentVR co;
		vector<double> pi(61);
		for(i=0;i<61;i++) pi[i] = oM.pi[i];
		NY98_61 deft_mut(1.0,1.0,1.0,pi,&ran);
		vector<Mutation_Matrix*> _mut(oM.L);
		vector<Mutation_Matrix*> mut(oM.L);
		for(i=0;i<oM.L;i++) _mut[i] = new NY98_61(deft_mut);

		vector<int> deft_LScodon(oM.L,-1);
		vector< vector<int> > LScodon(oM.n,deft_LScodon);
		vector<double> rho(oM.L-1,0.0);
	/*	ctrl.nsamp = oM.n;
		ctrl.ntimes = vector<double>(oM.n,0.0);
		ctrl.Negens = 1.0;
		ctrl.seq_len = oM.L;
		ctrl.len = vector<int>(1,oM.L);
		ctrl.r = ctrl.lambda = 0.0;
		ctrl.rmap = vector<double>(oM.L-1,0.0);
		vector< vector<int> > COcodon(LScodon);
		co.initialize(&ctrl,&ran);*/

		cout << "Reconstructing Markov chain..." << endl;
		int iter,sim;
		oBlock *oB;
		vector<double> simulated;
		cout << endl;
		

		siteS = vector<double>(oM.L,0.0);
		double mu,kappa;
		for(iter=0;iter<oM.niter;iter++) {
			oM.iterate(iter);

			if(iter>=argP.burnin && (iter-argP.burnin)%argP.thinning==0) {
				for(sim=0;sim<argP.nsim;sim++) {
					siteT = vector<double>(oM.L,0.0);
	//				cout << "Equating mutation rate matrices:" << endl;
					for(i=0;i<oM.L;i++) {
						mut[i] = _mut[oM.block[i]->start];
	//					cout << oM.block[i]->start << " ";
					}
	//				cout << endl;
					oB = oM.block[0];
					double avgPi = 0.0;
	//				cout << "Updating certain mutation rate matrices:" << endl;
					mu = oB->oMat->mu;
					kappa = oB->oMat->kappa;
					while(oB!=0) {
						((NY98_61*)_mut[oB->start])->update(mu,kappa,oB->oMat->omega);
	//					cout << oB->start << " ";
	//					avgPi += 2.0*(double)(oB->end-oB->start+1)*_mut[oB->start]->expected_rate();
						oB = oB->_3prime;
					}
	//				cout << endl;
					for(i=0;i<oM.L-1;i++) rho[i] = oM.rblock[i]->rho;
					LS.go(LScodon,(vector<Mutation_Matrix*>&)mut,rho);

					/*cout << "Simulated sequences:" << endl;
					for(i=0;i<LScodon.size();i++) {
						for(j=0;j<LScodon[i].size();j++) cout << codonToTriplet(LScodon[i][j]);
						cout << endl;
					}	cout << endl;
					myutils::pause();*/
					simulated = calcStats(calculate,LScodon);
	/*				ctrl.rmap = rho;
					co.go();
					for(i=0;i<oM.L;i++) co.mutate(i,mut[i]);
					for(i=0;i<oM.n;i++)
						for(j=0;j<oM.L;j++)
							COcodon[i][j] = (int)co.genotype[i][j];
					simulated = calcStats(calculate,COcodon);*/

	//				argP.out << avgPi << "\t";
	//				for(i=0;i<oM.L;i++) argP.out << siteT[i] << "\t";
					for(i=0;i<argP.statistics.size();i++) {
						argP.out << simulated[(int)argP.statistics[i]];
						if(i==argP.statistics.size()-1) argP.out << endl;
						else argP.out << "\t";
					}

#ifdef _WIN32
					cout << "\rDone simulation " << sim+1 << " of " << argP.nsim << " for iteration " << iter+argP.thinning << " of " << oM.niter << " for chain " << mcmc+1 << " of " << MCMCfile.size();
#endif
				}
			}
		}	cout << endl << endl;
	}

//	argP.out << observed[(int)predP::pi] << "\t";
//	for(i=0;i<oM.L;i++) argP.out << realT[i] << "\t";
	for(i=0;i<argP.statistics.size();i++) {
		argP.out << observed[(int)argP.statistics[i]];
		if(i==argP.statistics.size()-1) argP.out << endl;
		else argP.out << "\t";
	}
	argP.out.close();

/*	cout << "Sitewise S:" << endl;
	for(i=0;i<oM.L;i++) cout << siteS[i] << ",";
	cout << "\b\n";*/
	myutils::pause();
	return 0;
}

/*	control file:
	
	FASTA = fasta.txt				//	FASTA file name
	MCMC = mcmc.txt					//	output from MCMC chain
	burnin = bu						//	burn-in to use
	thinning = th					//	thinning to use in extracting MCMC chain
	nsim = n						//	number of simulations per iteration
									//	(total # sims = floor((niter-burnin)/thinning)
	outfile = out.txt				//	output for draw from joint distribution of statistics
	statistics = list				//	list of the following summary statistics to use:

					S				//	Number of segregating sites
					H				//	Number of unique haplotypes
					pi				//	Average pairwise diversity
					Varpi			//	Variance in pairwise diversity
					Rm				//	Myers's Rm (similar to Hudson's Rh)
					Rsq				//	Correlation between physical distance and R-squared
					Dprime			//	Correlation between physical distance and D-prime
					G4				//	Correlation between physical distance and the 4-gamete test
*/

vector<double> calcStats(vector<bool> &calculate, vector< vector<int> > &codon) {
	int i;
	if(calculate.size()!=10) error("calcStats(): calculate must be of length 10");
	vector<double> observed(10,0.0);
	if(calculate[0]) observed[0] = S(codon);
	if(calculate[1]) observed[1] = H(codon);
	if(calculate[2]) observed[2] = pi(codon);
	if(calculate[3]) observed[3] = Varpi(codon);
	if(calculate[4]) observed[4] = Rm(codon);
	if(calculate[5] || calculate[6] || calculate[7]) {
		vector<double> temp = RecCorrelations(codon);
		for(i=0;i<3;i++) observed[i+5] = temp[i];
	}
	if(calculate[8]) observed[8] = U(codon);
	if(calculate[9]) observed[9] = D(codon);
	return observed;
}

/* Number of segregating sites */
double S(vector< vector<int> > &H) {
	double result = 0.0;
	if(H.size()==0) return 0.0;

	int i,j;
	for(j=0;j<H[0].size();j++) {
		int hap = H[0][j];
		for(i=1;i<H.size();i++)
			if(H[i][j]!=hap) {
				//++result;
				++siteS[j];
				if(realS.size()==H[0].size()) {
					if(realS[j]==0) ++result;
				}
				else {
					++result;
				}
				break;
			}
		if(i==H.size() && realS.size()==H[0].size())
			if(realS[j]==1) ++result;
	}

	return result;
}

/*	U() counts the number of sites 

	U = sum_{i=1}^{L} I(Ts_i>To_i) / sum_{i=1}^{L} I(Ts_i<>To_i)
	0 <= U <= 1. Under null hypothesis, E(U) = 0.5

	where I is the indicator function,
	where T_i = n_i - s_i
	where n_i = # pairwise non-synonymous differences
	      s_i = # pairwise synonymous differences

	and   To  = for observed dataset
	      Ts  = for simulated dataset

	so	  T < 0	indicative of functional constraint
	      T > 0 indicative of diversifying selection

	and U is a measure of how often sites are correctly
	assigned 'positive' or 'negative' selection.
*/
double U(vector< vector<int> > &H) {
	double result = 0.0;
	if(H.size()==0) return 0.0;

	int h,i,j;
	int type;
	int syn,nsyn,neither;
	double total = 0.0;
	for(j=0;j<H[0].size();j++) {
		syn = nsyn = neither = 0;
		for(h=0;h<H.size();h++)
			for(i=0;i<h;i++) {
				if(H[h][j]==61 || H[i][j]==61) type = 0;
				else type = (int)changeType->C[H[h][j]][H[i][j]];
				switch(type) {
				case 0:	++neither;	break;
				case 1:	++syn;		break;
				case 2:	++syn;		break;
				case 4:	++nsyn;		break;
				case 8:	++nsyn;		break;
				default:++neither;	break;
				}
			}
		siteT[j] = nsyn - syn;
		if(realT.size()==H[0].size()) {
			if(siteT[j]>realT[j]) {
				++result;
				++total;
			}
			else if(siteT[j]<realT[j])
				++total;
		}
	}
	if(total!=0.0) result /= total;
	else result = 0.5;
	
	return result;
}

/* Number of unique haplotypes */
double H(vector< vector<int> > &H) {
	int result = 1;

	if(H.size()==0) return 0.0;
	vector< vector<int>* > uniqueHaps(H.size(),NULL);
	uniqueHaps[0] = &(H[0]);

	int i,ii,j;
	bool unique;
	for(i=1;i<H.size();i++) {
		unique = true;
		for(ii=0;ii<result;ii++) {
			for(j=0;j<H[0].size();j++)
				//if(H[i][j]!=H[ii][j]) break;
				if(H[i][j]!=(*uniqueHaps[ii])[j]) break;
			if(j==H[0].size()) unique = false;
		}
		if(unique==true) {
			uniqueHaps[result] = &(H[i]);
			++result;
		}
	}

	return (double)result;
}

/* Average number of pairwise differences */
double pi(vector< vector<int> > &H) {
	double result = 0.0;

	int i,j,k;
	for(i=0;i<H.size();i++)
		for(j=0;j<i;j++)
			for(k=0;k<H[0].size();k++)
				result += (H[i][k]==H[j][k]) ? 0.0 : 1.0;
	result *= 2.0/(double)(H.size())/(double)(H.size()-1);

	return result;
}

/* Variance in number of pairwise differences */
double Varpi(vector< vector<int> > &H) {
	double E,EE,pi;
	int i,j,k;

	E = EE = 0.0;
	for(i=0;i<H.size();i++)
		for(j=0;j<i;j++) {
			pi = 0.0;
			for(k=0;k<H[0].size();k++)
				pi += (H[i][k]==H[j][k]) ? 0.0 : 1.0;
			E += pi;
			EE += pi*pi;
		}
	E *= 2.0/(double)(H.size())/(double)(H.size()-1);
	EE *= 2.0/(double)(H.size())/(double)(H.size()-1);

	double result = EE - E*E;
	return result;
}

/* Hudson and Kaplan's Rm, the minimum # recombinations
   See Myers and Griffiths(2003) or Hein, Schierup and Wiuf (2005) */
double Rm(vector< vector<int> > &H) {
	if(H.size()==0) return 0.0;
	if(H[0].size()==0) return 0.0;

	/* Determine which sites are biallelic segregating */
	vector<int> sites(H[0].size(),0);
	int i,j,k;
	int S = 0;
	int hap0,hap1;
	bool segregating;
	for(j=0;j<H[0].size();j++) {
		segregating = false;
		hap0 = H[0][j];
		for(i=1;i<H.size();i++) {
			if(!segregating && H[i][j]!=hap0) {
				segregating = true;
				hap1 = H[i][j];
			}
			else if(segregating && H[i][j]!=hap0 && H[i][j]!=hap1) {
				segregating = false;	// define segregating only for biallelic sites
				break;
			}
		}
		if(segregating) {
			sites[S] = j;
			++S;
		}
	}
	if(S<2) return 0.0;

	/* Calculate the compatibility matrix */
	LowerTriangularMatrix<int> B(S,0);	// so j>=k always
	// B[j][k] = 0 for compatible, 1 for incompatible
	bool comb[3];
	for(j=0;j<S;j++)
		for(k=0;k<j;k++)
		{
			hap0 = H[0][sites[j]];
			hap1 = H[0][sites[k]];
			comb[0] = false;				// hap0  hap1'
			comb[1] = false;				// hap0' hap1
			comb[2] = false;				// hap0' hap1'
			for(i=1;i<H.size();i++) {
				if(H[i][sites[j]]==hap0 && H[i][sites[k]]!=hap1) comb[0] = true;
				if(H[i][sites[j]]!=hap0 && H[i][sites[k]]==hap1) comb[1] = true;
				if(H[i][sites[j]]!=hap0 && H[i][sites[k]]!=hap1) comb[2] = true;
				if(comb[0] && comb[1] && comb[2]) break;
			}
			B[j][k] = (comb[0] && comb[1] && comb[2]) ? 1 : 0;			
		}

/*	cout << "Pairwise compatibilities: (1=incompatible):" << endl;
	for(i=0;i<S;i++) {
		for(j=0;j<S;j++)
			cout << B[i][j] << " ";
		cout << endl;
	}	cout << endl;*/

	/* Calculate the dynamic programming partition matrix */
	vector<int> M(S,0);
	int maxM = 0;
	M[S-1] = 0;
	M[S-2] = B[S-1][S-2];
	for(i=S-3;i>=0;i--) {
		M[i] = B[i+1][i] + M[i+1];
		for(k=i+2;k<S;k++) if(B[k][i]+M[k]>M[i]) M[i] = B[k][i]+M[k];
	}

/*	vector<int> M(S,0);
	int maxM = 0;
	M[0] = 0;
	M[1] = B[1][0];
	for(k=2;k<S;k++) {
		M[k] = M[0] + B[k][0];
		for(i=1;i<k;i++)
			if((M[i] + B[k][i]) > M[k]) M[k] = M[i] + B[k][i];
	}*/

	return (double)M[0];
}

vector<double> RecCorrelations(vector< vector<int> > &H) {
	vector<double> result(3,0.0);
	if(H.size()==0) return result;
	if(H[0].size()==0) return result;

	/* Determine which sites are biallelic segregating */
	vector<int> sites(H[0].size(),0);
	int i,j,k;
	int S = 0;
	int hap0,hap1;
	bool segregating;
	for(j=0;j<H[0].size();j++) {
		segregating = false;
		hap0 = H[0][j];
		for(i=1;i<H.size();i++) {
			if(!segregating && H[i][j]!=hap0) {
				segregating = true;
				hap1 = H[i][j];
			}
			else if(segregating && H[i][j]!=hap0 && H[i][j]!=hap1) {
				segregating = false;	// define segregating only for biallelic sites
				break;
			}
		}
		if(segregating) {
			sites[S] = j;
			++S;
		}
	}
	if(S<3) return result;
	
	/* Calculate frequency statistics */
	vector<double> F(S,1.0);							/* F is the marginal frequency of hap0 at site j */
	for(j=0;j<S;j++) {
		hap0 = H[0][sites[j]];
		for(i=1;i<H.size();i++)
			if(H[i][sites[j]]==hap0) F[j]++;
		F[j] /= (double)H.size();
	}

	vector<double> four(4,0.0);							/* G[j][k] is the frequency of AB (G[j][k][0]),	*/
	LowerTriangularMatrix< vector<double> > G(S,four);	/* Ab (1), aB (2), ab (3) for sites j and k		*/
	for(j=0;j<S;j++)
	  for(k=0;k<j;k++) {
		  hap0 = H[0][sites[j]];
		  hap1 = H[0][sites[k]];
		  for(i=0;i<H.size();i++) {
			  if(H[i][sites[j]]==hap0 && H[i][sites[k]]==hap1) ++G[j][k][0];
			  else if(H[i][sites[j]]==hap0 && H[i][sites[k]]!=hap1) ++G[j][k][1];
			  else if(H[i][sites[j]]!=hap0 && H[i][sites[k]]==hap1) ++G[j][k][2];
			  else if(H[i][sites[j]]!=hap0 && H[i][sites[k]]!=hap1) ++G[j][k][3];
			  else warning("Unexpected choice");
		  }
		  for(i=0;i<4;i++) G[j][k][i] /= (double)H.size();
	  }
  
	/* Calculate LD statistics for pairs of sites */
	LowerTriangularMatrix<double> A(S,0.0);			//	rsq
	LowerTriangularMatrix<double> B(S,0.0);			//	Dprime
	LowerTriangularMatrix<double> C(S,0.0);			//	G4
	Matrix<double> D(S,S,0.0);

	double temp;
//	ofstream out("__out.txt");
//	char tab = '\t';
//	out << "locusA" << tab << "locusB" << tab << "rsq" << tab << "Dprime" << tab << "G4" << tab << "dist" << endl;
	for(i=0;i<S;i++) {
		for(j=0;j<i;j++) {
			temp = G[i][j][0] - F[i]*F[j];
			A[i][j] = pow(temp,2.0)/(F[i]*(1.-F[i])*F[j]*(1.-F[j]));
			B[i][j] = (temp < 0.0) ? -temp/MIN(F[i]*F[j],(1.-F[i])*(1.-F[j])) : temp/MIN(F[i]*(1.-F[j]),(1.-F[i])*F[j]);
			C[i][j] = (G[i][j][0]>0.0 && G[i][j][1]>0.0 && G[i][j][2]>0.0 && G[i][j][3]>0.0) ? 1.0 : 0.0;
			D[i][j] = D[j][i] = sites[i] - sites[j];
//			out << i << tab << j << tab << A[i][j] << tab << B[i][j] << tab << C[i][j] << tab << D[i][j] << endl;
		}
	}
//	out.close();

	double  E[4] = {0.0,0.0,0.0,0.0};
	double EE[4] = {0.0,0.0,0.0,0.0};
	double ED[3] = {0.0,0.0,0.0};
	int ctr;
	for(i=0,ctr=0;i<S;i++)
		for(j=0;j<i;j++,ctr++) {
			E[0] += A[i][j]; E[1] += B[i][j]; E[2] += C[i][j]; E[3] += D[i][j];
			EE[0] += A[i][j]*A[i][j]; EE[1] += B[i][j]*B[i][j]; EE[2] += C[i][j]*C[i][j]; EE[3] += D[i][j]*D[i][j];
			ED[0] += A[i][j]*D[i][j]; ED[1] += B[i][j]*D[i][j]; ED[2] += C[i][j]*D[i][j];
		}
	for(k=0;k<3;k++)
		result[k] = (ED[k]-E[k]*E[3]/(double)ctr)/sqrt((EE[k]-E[k]*E[k]/(double)ctr)*(EE[3]-E[3]*E[3]/(double)ctr));

	/* Calculate remaining statistics for correlation coefficients */
/*	double Abar,Bbar,Cbar,Dbar;
	double Adev,Bdev,Cdev,Ddev;
	Abar=Bbar=Cbar=Dbar=Adev=Bdev=Cdev=Ddev=0.0;
	double ctr = 0.0;
	for(i=0;i<S;i++)
		for(j=0;j<i;j++,ctr++) {
			Abar += A[i][j]; Bbar += B[i][j]; Cbar += C[i][j]; Dbar += D[i][j];
		}
	Abar /= ctr; Bbar /= ctr; Cbar /= ctr; Dbar /= ctr;
	for(i=0;i<S;i++)
		for(j=0;j<i;j++) {
			Adev += pow(A[i][j]-Abar,2.0);
			Bdev += pow(B[i][j]-Bbar,2.0);
			Cdev += pow(C[i][j]-Cbar,2.0);
			Ddev += pow(D[i][j]-Dbar,2.0);
		}

	double dist;
	for(j=0;j<S;j++)
		for(k=0;k<j;k++) {
			dist = D[j][k];
			result[0] += A[j][k] * dist;
			result[1] += B[j][k] * dist;
			result[2] += C[j][k] * dist;
	}
	double Aadd = ctr*Abar*Dbar; double Adiv = sqrt(Adev*Ddev);
	double Badd = ctr*Bbar*Dbar; double Bdiv = sqrt(Bdev*Ddev);
	double Cadd = ctr*Cbar*Dbar; double Cdiv = sqrt(Cdev*Ddev);
	result[0] -= Aadd; if(Adiv!=0) result[0] /= Adiv;
	result[1] -= Badd; if(Bdiv!=0) result[1] /= Bdiv;
	result[2] -= Cadd; if(Cdiv!=0) result[2] /= Cdiv;*/

/*	if(!(result[0]!=numeric_limits<double>::quiet_NaN())) error("Problem in Rsq");
	if(!(result[1]!=numeric_limits<double>::quiet_NaN())) error("Problem in Dprime");
	if(!(result[2]!=numeric_limits<double>::quiet_NaN())) error("Problem in G4");*/
	return result;
}

/* Tajima's D */
double D(vector< vector<int> > &H) {
	double D = 0.0;
	int i,j,k,l,n,L;
	n = H.size();
	L = H[0].size();
	double a1,a2,b1,b2,c1,c2,e1,e2,khat,S;
	bool segregating;
	khat = S = 0.0;
	vector<string> vec(n,"NNN");
	for(k=0;k<L;k++) {
		for(i=0;i<n;i++) vec[i] = codonToTriplet(H[i][k]);
		for(l=0;l<3;l++) {
			segregating = false;
			for(i=0;i<n;i++)
				for(j=0;j<i;j++) {
					if(vec[i][l]!=vec[j][l]) {
						++khat;
						segregating = true;
					}
				}
			if(segregating) ++S;
		}
	}
	if(S==0) return 0.0;
	khat /= (double)(n*(n-1)/2);
	a1 = a2 = 0.0;
	for(i=1;i<=n-1;i++) {
		a1 += 1./(double)i;
		a2 += 1./(double)(i*i);
	}
	b1 = (double)(n+1)/(double)(3*(n-1));
	b2 = (double)(2*(n*n+n+3))/(double)(9*n*(n-1));
	c1 = b1 - 1./a1;
	c2 = b2 - (double)(n+2)/a1/(double)(n) + a2/a1/a1;
	e1 = c1/a1;
	e2 = c2/(a1*a1+a2);
	D = (khat - S/a1)/sqrt(e1*S+e2*S*(S-1.));
	return D;
}

/* Returns 0-60 for non-STOP codons, 61 for indels */
void dnaToCodons(DNA &dna, vector< vector<int> > &codon) {
	int n = dna.nseq;
	int L = (int)floor((double)dna.lseq/(double)3.0);
	if(n<=0) error("dnaToCodons(): n must be positive");
	if(L<=0) error("dnaToCodons(): L must be positive");
	
	vector<int> default_codon_sequence(L,-1);
	codon = vector< vector<int> >(n,default_codon_sequence);
	int i,j;
	string triplet,tr1,tr2,tr3;
	for(i=0;i<n;i++)
		for(j=0;j<L;j++) {
			tr1 = string(1,dna[i][3*j]);
			tr2 = string(1,dna[i][3*j+1]);
			tr3 = string(1,dna[i][3*j+2]);
			//triplet = dna[i][3*j]+dna[i][3*j+1]+dna[i][3*j+2];
			triplet = tr1 + tr2 + tr3;
			codon[i][j] = removeStopCodons(tripletToCodon(triplet));
		}
}

/* Returns 0-63 for codons and 64 for indels */
int tripletToCodon(string &tri) {
	map<char,int> baseToInt;
	baseToInt['T'] = 1;
	baseToInt['U'] = baseToInt['T'];
	baseToInt['C'] = 2;
	baseToInt['A'] = 3;
	baseToInt['G'] = 4;
	baseToInt['-'] = 5;
	int a = baseToInt[tri[0]];
	int b = baseToInt[tri[1]];
	int c = baseToInt[tri[2]];
	bool indel = false;
	if(a==5) indel = true;
	if(b==5) indel = true;
	if(c==5) indel = true;
	if(indel==true) {
		if(a==5 && b==5 && c==5) return 64;
		else return -1;
	}
	/* return a value from 0 to 63 */
	return (a-1)*16 + (b-1)*4 + c - 1;
}

int removeStopCodons(const int a) {
	if(a<0) error("removeStopCodons(): invalid codon number (<0)");
	if(a<10) return a;
	if(a==10) error("removeStopCodons(): found stop codon TAA/UAA");
	if(a==11) error("removeStopCodons(): found stop codon TAG/UAG");
	if(a<14) return a-2;
	if(a==14) error("removeStopCodons(): found stop codon TGA/UGA");
	/* a==64 is for indels. Change this to 61 */
	if(a>64) error("removeStopCodons(): invalid codon number (>64)");
	return a-3;
}

string codonToTriplet(const int a) {
	if(a==64) return string(3,'-');
	char intToBase[6] = {'N','T','C','A','G','-'};
	string tri(3,'N');
	tri[0] = intToBase[a/16 + 1];
	tri[1] = intToBase[(a%16)/4 + 1];
	tri[2] = intToBase[(a%16)%4 + 1];
	return tri;
}
