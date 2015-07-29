'''
Created on 10/29/2014
@Ronak Shah

'''
import argparse
import os
import sys
import time
import stat
from pprint import pprint
import re
import pandas as pd

def main():
    parser = argparse.ArgumentParser(prog='CompileHSmetrics.py', description='Compile Results of HS metrics calculations', usage='%(prog)s [options]')
    parser.add_argument("-b", "--bamlist", action="store", dest="bamlist", required=True, metavar='bams.fof', help="Full path bam file of files.") 
    parser.add_argument("-i", "--inputDir", action="store", dest="inputDir", required=True, metavar='/somepath/input', help="Full Path to the input directory.")
    parser.add_argument("-f", "--fractions", action="store", nargs='*', dest="fractions", required=False,default='0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90', metavar='0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90', help="Number of fractions to which the bam file need to be down-sampled; Multiple fraction are separated by space")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
    parser.add_argument("-o", "--outDir", action="store", dest="outDir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
    parser.add_argument("-of", "--outfile", action="store", dest="outfile", required=False,default='AllSample_dfracCoverage.txt', metavar='output.txt', help="Name of the output file.")
    args = parser.parse_args()
    if(args.verbose):
        print "Checking the and processing Arguments.\n"
    (bamIDs,filesToProcess,fractionsToDownSample,fractionsDirs) = ProcessArgs(args)
    MergeResults(bamIDs,filesToProcess,fractionsToDownSample,fractionsDirs,args)
       
def ProcessArgs(args):
    if os.path.isfile(args.bamlist):
        filesToProcess = []
        bamIDs = []
        with open(args.bamlist, 'r') as filecontent:
            for line in filecontent:
                fileName = os.path.basename(line)
                baseName = re.search('(.*)_bc\d+.*', fileName).group(1)
                bamIDs.append(baseName)
                filesToProcess.append(fileName.rstrip('\n'))
    fractionsToDownSample = args.fractions.split(" ")
    fractionsDirs = []
    unAccountedFractions = []
    for fraction in fractionsToDownSample:
        os.chdir(args.outDir)
        fractionOutDir = args.outDir + "/Fraction_" + fraction
        if os.path.isdir(fractionOutDir):
                if(args.verbose):
                    print "The output directory ", fractionOutDir, " exists"
                fractionsDirs.append(fractionOutDir)
        else:    
            if(args.verbose):
                    print"Directory", fractionOutDir," does not exists to merge the results.\n"
                    unAccountedFractions.append(fraction)
                    continue
    for unAccountedFraction in unAccountedFractions:
        fractionsToDownSample.remove(unAccountedFraction)
    return(bamIDs,filesToProcess,fractionsToDownSample,fractionsDirs)

def MergeResults(bamIDs,filesToProcess,fractionsToDownSample,fractionsDirs,args):
    outFileTxt = args.outDir + "/" + args.outfile
    coveragePerFraction = []
    for dir in fractionsDirs:
        fcount = 0
        coveragePerFile = []
        covgIndex = 0
        for file in filesToProcess:
            
            hsMetricsOutFile = file.rstrip('.bam')
            hsMetricsOutFile = dir + "/" + hsMetricsOutFile + ".HSmertics.txt"
            
            with open(hsMetricsOutFile, 'r') as filecontent:
                coverageForSample = 0
                linecount=0
                for line in filecontent:
                    line = line.rstrip()
                    if line:  
                        if line.startswith("#"):
                            linecount = linecount+1
                            continue
                        #print linecount,line
                        data = line.rstrip('\n').split('\t')
                        if(fcount == 0):
                                if(data[0] == "BAIT_SET"):
                                    for content in data:
                                        #print linecount,content ,"C\n"
                                        if(content == 'MEAN_TARGET_COVERAGE'):
                                            covgIndex = covgIndex
                                            #print covgIndex,"CI\n"
                                            break
                                        else:
                                            covgIndex = covgIndex + 1
                                else:
                                    #print fcount,covgIndex,data,"Inside\n"  
                                    coverageForSample = data[covgIndex]  
                        if(fcount > 0):
                            if(data[0] == "BAIT_SET"):
                                continue
                            #print fcount,covgIndex,data,"All\n"  
                            coverageForSample = data[covgIndex]
                        linecount = linecount+1
            fcount = fcount + 1
            coveragePerFile.append(coverageForSample)
           
        coveragePerFraction.append(dict(zip(bamIDs, coveragePerFile)))
    #print coveragePerFraction,"\n"
    mergeDF = pd.DataFrame(coveragePerFraction,index=fractionsToDownSample)
    mergeDF.to_csv(outFileTxt, sep='\t')        
    return

if __name__ == "__main__":
    start_time = time.time() 
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))    