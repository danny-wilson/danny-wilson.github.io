#ifndef _LiS_sim_H_
#define _LiS_sim_H_

#include <myutils.h>
#include <DNA.h>
#include <mutation.h>
#include <vector>

using namespace std;
using namespace myutils;

class LiSsimulator {
public:
	Random *ran;
	vector<double> nrec;

public:
	LiSsimulator(Random &ran_in) {
		ran = &ran_in;
	}

	LiSsimulator& go(vector< vector<int> > &hap, Mutation_Matrix &mut, const double rho) {
		int k;
		nrec = vector<double>(hap[0].size(),0.0);
		for(k=0;k<hap.size();k++) conditional(hap,k,mut,rho);
		return *this;
	}
	LiSsimulator& go(vector< vector<int> > &hap, Mutation_Matrix &mut, vector<double> rho) {
		int k;
		nrec = vector<double>(hap[0].size(),0.0);
		if(rho.size()!=hap[0].size()-1) error("LiSsimulator::go(): rho should be of length L-1");
		for(k=0;k<hap.size();k++) conditional(hap,k,mut,rho);
		return *this;
	}
	LiSsimulator& go(vector< vector<int> > &hap, vector<Mutation_Matrix*> &mut, vector<double> rho) {
		int k;
		nrec = vector<double>(hap[0].size(),0.0);
		if(rho.size()!=hap[0].size()-1) error("LiSsimulator::go(): rho should be of length L-1");
		for(k=0;k<hap.size();k++) conditional(hap,k,mut,rho);
		return *this;
	}

	/*	k is the number of haps conditioned on, assumed to be the first k
		DNA sequences in hap. hap is assumed to be compatible with mut		*/
	LiSsimulator& conditional(vector< vector<int> > &hap, const int k, Mutation_Matrix &mut, const double rho) {
		if(k>=hap.size()) error("LiSsimulator::conditional(): k >= n");
		if(k<0) error("LiSsimulator::conditional(): k < 0");
		int lseq = hap[0].size();
		nrec.resize(lseq,0);

		int i;
		if(k==0) {
			for(i=0;i<lseq;i++) hap[0][i] = mut.draw();
		}
		else if(k==1) {
			double time;
			for(i=0;i<lseq;i++) {
				time = 2.0*ran->exponential(1.0/(double)k);
				hap[k][i] = mut.mutate_edge(hap[0][i],time);
			}
		}
		else {
			double time;
			int from = ran->discrete(0,k-1);
			//double nextRec = ran->exponential(2.0/((double)k*rho));
			double nextRec = ran->exponential((double)k/rho);
			for(i=0;i<lseq;i++) {
				if((double)i>nextRec) {
					nrec[i-1]++;
					from = ran->discrete(0,k-1);
					//nextRec = (double)i + ran->exponential(2.0/((double)k*rho[i]));
					nextRec = (double)i + ran->exponential((double)k/rho);
				}
				time = 2.0*ran->exponential(1.0/(double)k);
				hap[k][i] = mut.mutate_edge(hap[from][i],time);
			}
		}
		return *this;
	}

	/*	k is the number of haps conditioned on, assumed to be the first k
		DNA sequences in hap. hap is assumed to be compatible with mut		*/
	LiSsimulator& conditional(vector< vector<int> > &hap, const int k, Mutation_Matrix &mut, vector<double> rho) {
		if(k>=hap.size()) error("LiSsimulator::conditional(): k >= n");
		if(k<0) error("LiSsimulator::conditional(): k < 0");
		int lseq = hap[0].size();
		nrec.resize(lseq,0);

		int i;
		if(k==0) {
			for(i=0;i<lseq;i++) hap[0][i] = mut.draw();
		}
		else if(k==1) {
			double time;
			for(i=0;i<lseq;i++) {
				time = 2.0*ran->exponential(1.0/(double)k);
				hap[k][i] = mut.mutate_edge(hap[0][i],time);
			}
		}
		else {
			double time;
			int from = ran->discrete(0,k-1);
			//double nextRec = ran->exponential(2.0/((double)k*rho[0]));
			double nextRec = ran->exponential((double)k/rho[0]);
			for(i=0;i<lseq;i++) {
				if((double)i>nextRec) {
					nrec[i-1]++;
					from = ran->discrete(0,k-1);
					//if(i<lseq-1) nextRec = (double)i + ran->exponential(2.0/((double)k*rho[i]));
					if(i<lseq-1) nextRec = (double)i + ran->exponential((double)k/rho[i]);
				}
				else if(i>0 && i<(lseq-1)) {
					nextRec -= i;
					nextRec *= rho[i-1]/rho[i];
					nextRec += i;
				}
				time = 2.0*ran->exponential(1.0/(double)k);
				hap[k][i] = mut.mutate_edge(hap[from][i],time);
			}
		}
		return *this;
	}

	/*	k is the number of haps conditioned on, assumed to be the first k
		DNA sequences in hap. hap is assumed to be compatible with mut		*/
	LiSsimulator& conditional(vector< vector<int> > &hap, const int k, vector<Mutation_Matrix*> &mut, vector<double> rho) {
		if(k>=hap.size()) error("LiSsimulator::conditional(): k >= n");
		if(k<0) error("LiSsimulator::conditional(): k < 0");
		int lseq = hap[0].size();
		nrec.resize(lseq,0);

		int i;
		if(k==0) {
			for(i=0;i<lseq;i++) hap[0][i] = mut[i]->draw();
		}
		else if(k==1) {
			double time;
			for(i=0;i<lseq;i++) {
				time = 2.0*ran->exponential(1.0/(double)k);
				hap[k][i] = mut[i]->mutate_edge(hap[0][i],time);
			}
		}
		else {
			double time;
			int from = ran->discrete(0,k-1);
			//double nextRec = ran->exponential(2.0/((double)k*rho[0]));
			double nextRec = ran->exponential((double)k/rho[0]);
			for(i=0;i<lseq;i++) {
				if((double)i>nextRec) {
					nrec[i-1]++;
					from = ran->discrete(0,k-1);
					//if(i<lseq-1) nextRec = (double)i + ran->exponential(2.0/((double)k*rho[i]));
					if(i<lseq-1) nextRec = (double)i + ran->exponential((double)k/rho[i]);
				}
				else if(i>0 && i<(lseq-1)) {
					nextRec -= i;
					nextRec *= rho[i-1]/rho[i];
					nextRec += i;
				}
				time = 2.0*ran->exponential(1.0/(double)k);
				hap[k][i] = mut[i]->mutate_edge(hap[from][i],time);
			}
		}
		return *this;
	}

};


#endif // _LiS_sim_H_

