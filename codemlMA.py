#codeml Model Averaging
#Brice A. J. Sarver
#Version 1.2

#CHANGELOG
# 1.0: first stable release
# 1.1: incorporates LRTs and M8a testing
# 1.2: fix to the regexp for the number of codons. should handle any length now. other small changes.

from decimal import *
import os, re, math, argparse
from Bio.Phylo.PAML.chi2 import cdf_chi2

parser = argparse.ArgumentParser(prog="codemlMA", usage="%(prog)s [options]", description="Bayesian Model Averaging among Selection Models from Codeml")
parser.add_argument("filename", help="The path to the principal *.codeml file that contains the results from the models fit.")
parser.add_argument("outfile", help="The name of the output file.", default="results.txt")
parser.add_argument("-m8", "--M8a_test", help="The optional argument to include the test of M8 vs. M8a.  Requires a file with the results.")
args = parser.parse_args()

# 
# parser = OptionParser(usage = "usage: python %prog [-f FILENAME] [-o output_file_name]", version = "%prog 0.1")
# parser.add_option("-f", "--file", dest="filename", help="Path to the *.codeml file to be parsed")
# parser.add_option("-o", "--outfile", dest="outfile", help="Output File")
# (options, args) = parser.parse_args()
# if not options.filename and not options.outfile:
# 	parser.error("You must specify a codeml output file and a file to write results to.")


print "-----------------codemlMA.py-----------------"
print "Model averaging of selection models in codeml"
print "-------------Brice A. J. Sarver--------------"
print "----------------Aug 30th, 2013---------------"
print "------------------Version 1.2----------------\n\n\n"

print "Parsing:", args.filename, "Writing output to:", args.outfile

getcontext().prec = 25
parameters = []
likelihood = []
BIC = []
MA = []
posterior = []

print "Lists initialized; beginning to parse...\n"

resultsf = open(args.filename)
try:
	for line in resultsf:
		if "lnL" in line:
			hold = line
			print hold
			m = re.search('(?<=np:)(\s*)\d+', hold)
			#print m.group()
			t = re.search('(?<=ntime: )\d+', hold)
			#print t.group()
			freep = int(m.group()) - int(t.group()) - 1
			#print freep
			parameters.append(freep)
			l = re.search('-\d+.\d+', hold)
			#print l.group()
			likelihood.append(l.group())
		if "ls = " in line:
			holder = line
			p = re.findall('\d+', holder)
			sites = (p[1])
						
	if args.M8a_test:
		m8a = open(args.M8a_test)
		for line in m8a:
			if "lnL" in line:
				hold = line
				print hold
				m = re.search('(?<=np:)\s*\d+', hold)
				#print m.group()
				t = re.search('(?<=ntime: )\d+', hold)
				#print t.group()
				freep = int(m.group()) - int(t.group()) - 1
				#print freep
				parameters.append(freep)
				l = re.search('-\d+.\d+', hold)
				#print l.group()
				likelihood.append(l.group())
			if "ls = " in line:
				holder = line
				p = re.findall('\d+', holder)
				sites = (p[1])
		m8a.close()
finally:
	resultsf.close()
		
floatparam = [float(i) for i in parameters]
floatlike = [float(i) for i in likelihood]
floats = float(sites)

print "CHECKPOINT: POST-PARSING"
print "The list of parameters:", floatparam
print "The list of log likelihoods:", floatlike
print "The number of sites (codons):", floats

print "\nThis set of codeml runs fit", len(parameters), "models.\n"

for num in range(len(parameters)):
	print "Processing list element:", num
	bictemp = -2 * floatlike[num] + floatparam[num] * math.log(floats)
	print "BIC is:", bictemp
	BIC.append(float(bictemp))
	
print "\nCHECKPOINT: The list of BICs is:", BIC

print"\nBegin calculating the posterior probability...\n"

# for bic in BIC:
# 	num = math.exp(-bic/2)
# 	MA.append(float(num))	
# supersum = sum(MA)

#testing decimal library

holdnum = []
for bic in BIC:
	decnum = Decimal(-0.5) * Decimal(bic)
	#print decnum
	exphold = Decimal(decnum).exp()
	#print exphold
	holdnum.append(exphold)
	

# print holdnum
print "The sum of all the transformed BICs is:", sum(holdnum), "\n"

print "CHECKPOINT: The list of numerators of the posterior calculation is:", holdnum

for value in holdnum:
	post = value/sum(holdnum)
	posterior.append(post)


# for value in MA:
# 	post = value/supersum
# 	posterior.append(post)

print "\nThe sum of all posterior probabilities should equal unity.  The actual sum is:", (sum(posterior))
if sum(posterior) == 1:
	print "\nSummation checks out!\n"
else:
	print "Sum does not equal 1.  Something went wrong...or maybe not.  Remember that these are high-precision libraries, and 0.99999... is acceptable.\n"

print "Traditional likelihood ratio tests:\n"

#teststat = -2 * log ( likelihood1 / likelihood 0)

print "M1a vs. M2a\n"
df = 2
print "The two likelihoods are:", round(floatlike[1], 2), round(floatlike[2], 2)
teststat = abs(2 * (-1 * round(floatlike[2], 2) - -1 * round(floatlike[1], 2)))
print "The test statistic is:", teststat
m12_p = cdf_chi2(df, float(teststat))
print "p-value:", cdf_chi2(df, float(teststat))


print "\nM7 vs. M8\n"
df = 2
print "The two likelihoods are:", round(floatlike[7], 2), round(floatlike[8], 2)
teststat = abs(2 * (-1 * round(floatlike[8], 2) - -1 * round(floatlike[7], 2)))
print "The test statistic is:", teststat
m78_p = cdf_chi2(df, float(teststat))
print "p-value:", cdf_chi2(df, float(teststat))


if args.M8a_test:
	print "\nM8a vs. M8\n"
	df = 1
	print "The two likelihoods are:", round(floatlike[8], 2), round(floatlike[9], 2)
	teststat = abs(2 * (-1 * round(floatlike[9], 2) - -1 * round(floatlike[8], 2)))
	print "The test statistic is:", teststat
	m8a_p = float(cdf_chi2(df, float(teststat))) / 2
	print "p-value:", m8a_p, "(with 1 d.f. and p-value / 2)"



	
print "\nSUMMARY:\n"
print "(The model correspondes to nsites in the codeml control file.)\n"
with open(args.outfile, 'w') as out:
	for y in range(len(posterior)):
		print "The posterior probability for model", y, "is:", posterior[y]
		print >> out, "The posterior probability for model", y, "is:", posterior[y]
	print >> out, "The p-value for the M1a vs. M2a LRT is:", m12_p
	print >> out, "The p-value for the M7 vs. M8 LRT is:", m78_p
	if args.M8a_test:
		print >> out, "The p-value for the M8 vs. M8a LRT is:", m8a_p, "(with 1 d.f. and p-value / 2)"

print "\nExecution complete.\n"