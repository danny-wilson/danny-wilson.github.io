# README
# ~~~~~~
#
# Daniel Wilson. 15th January 2010
#

Thanks for your interest in isource, which uses an approximation to Wright's island model in the case of asymmetric migration to attribute individuals to putative source populations on the basis of allelic profile.

Running
~~~~~~~
The distribution contains the source code (with makefile), and precompiled executables for Windows XP and Mac OS X. The most up-to-date version is available from http://www.danielwilson.me.uk

The input format is simple. It is a tab-separated file with n+2 columns, where n are the number of loci. The columns are
	First:			Sequence type
	Second-penultimate:	Allelic profile (each locus is a separate column)
	Last:			Group identifier, where 0 corresponds to humans
An example input file is included, called input.txt. For an explanation of the example input file, see below.

To run go to a command line. The syntax is

isource input_file output_file_suffix number_of_iterations thinning_interval parameter_of_the_prior

For example, to run for 100,000 iterations, recording the state of the MCMC every 50 iterations, and utilizing a symmetric Dirichlet prior with parameter 1, you would type the following

isource input.txt output.txt 100000 50 1

Make sure to do multiple runs with different seeds and compare results to ensure convergence. (The seed is set automatically from the system clock, and printed to screen at the beginning of a run.) If convergence is poor or questionable, increase the number of iterations, or reduce the thinning interval (minimum 1). If file sizes are too large, try increasing the thinning interval.

Interpretation
~~~~~~~~~~~~~~
Included in the distribution is a file, isource.R, containing R commands for analysing the output. It is annotated to explain the analyses.


Details on the example data sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are three example training datasets, train1.txt, train2.txt and train3.txt. The group numbers correspond to the following groupings, in the order shown, starting with group 1:

train1.txt CATTLE,CHICKEN,BIRD,ENVIRONMENT,SHEEP,PIG
train2.txt CATTLE,CHICKEN,BIRD,WATER,SHEEP,PIG,SAND
train3.txt CATTLE,CHICKEN,BIRD,WATER,SHEEP,PIG,SAND,RABBIT

To run isource, the training data must be combined with a test data set, i.e. a group of isolates of unknown source. These isolates must be labelled group 0, and appended to the training data in a single text file which is then used as input to isource.

Lancashire.txt contains the human isolates from Wilson et al (2008).
input.txt is a concatenation of Lancashire.txt and train3.txt.

NB:- Frequencies matter! If an ST was sampled multiple times from the same source, then it must appear that many times in the input file. The island model relies on differences in allele frequencies (and combinations of allele frequencies) between source populations in order to do source attribution.

Citing
~~~~~~
If you use the program isource, the training data or the human test data in your publications, please cite the following:

Wilson, D. J., E. Gabriel, A. J. H. Leatherbarrow, J. Cheesbrough, S. Gee, E. Bolton, A. Fox, P. Fearnhead, C. A. Hart and P. J. Diggle (2008)
Tracing the source of campylobacteriosis.
PLoS Genetics 4: e1000203.

More information
~~~~~~~~~~~~~~~~
For more information or contact details, please visit my website http://www.danielwilson.me.uk

