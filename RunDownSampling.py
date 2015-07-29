'''
Created on 10/28/2014
@Ronak Shah

'''
import argparse
import os
import sys
import time
import stat
from subprocess import Popen
import shlex
import shutil
from pprint import pprint
import re
import multiprocessing as mp
from Queue import *
from threading import Thread, Lock

def main():
   parser = argparse.ArgumentParser(prog='RunDownSampling.py', description='Down Samples Each Bam File', usage='%(prog)s [options]')
   parser.add_argument("-i", "--bamlist", action="store", dest="bamlist", required=True, metavar='bams.fof', help="Full path bam file of files.") 
   parser.add_argument("-j", "--java", action="store", dest="java", required=False,default='/dmp/resources/prod/tools/system/java/java-1.7.0-openjdk-1.7.0.9.x86_64/jre/bin/java', metavar='/somepath/java', help="Full path JAVA executables.") 
   parser.add_argument("-g", "--gatk", action="store", dest="gatk", required=False, default='/dmp/resources/prod/tools/bio/gatk/production/GenomeAnalysisTK.jar', metavar='/somepath/gatk.jar', help="Full Path to GATK JAR file")
   parser.add_argument("-p", "--picard", action="store", dest="picard", required=False,default='/dmp/resources/prod/tools/bio/picard/production/', metavar='/somepath/picard.jar', help="Full path to Picard Jar File") 
   parser.add_argument("-r", "--referenceFile", action="store", dest="ref", required=False,default='/dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta', metavar='/somepath/Homo_Sapeins_hg19.fasta', help="Full Path to the reference file with the bwa index.")
   parser.add_argument("-q", "--queue", action="store", dest="queue", required=False,default='clin2.q', metavar='all.q or clin.q', help="Name of the SGE queue")
   parser.add_argument("-qsub", "--qsubPath", action="store", dest="qsub", required=False,default='/common/sge/bin/lx-amd64/qsub', metavar='/somepath/qsub', help="Full Path to the qsub executables of SGE.")
   parser.add_argument("-f", "--fractions", action="store", nargs='*', dest="fractions", required=False,default='0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90', metavar='0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90', help="Number of fractions to which the bam file need to be down-sampled; Multiple fraction are separated by space")
   parser.add_argument("-t", "--threads", action="store", dest="threads", required=False,default='5', metavar='5', help="Number of Threads to be used to run Pindel")
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
   parser.add_argument("-o", "--outDir", action="store", dest="outDir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
   parser.add_argument("-ti", "--targetIntervals", action="store", dest="targetIntervals", required=False, default='/dmp/data/mskdata/interval-lists/production/picard_targets.interval_list', metavar='/somePath/taget.intervals', help="Full Path to Target Interval File")
   parser.add_argument("-bi", "--baitIntervals", action="store", dest="baitIntervals", required=False, default= '/dmp/data/mskdata/interval-lists/production/picard_baits.interval_list', metavar='/somePath/bait.intervals', help="Full Path to Bait Interval File")
   
   args = parser.parse_args()
   jobqueue = Queue()
   #Check how many process the analysis should be ran for
   if(args.verbose):
       print "Checking the number of process given for analysis."
       RunDownSample(args,jobqueue)
       
       
def RunDownSample(args,jobqueue):
    if os.path.isfile(args.bamlist):
        fileToProcess = []
        with open(args.bamlist, 'r') as filecontent:
            for line in filecontent:
                fileToProcess.append(line.rstrip('\n'))
    fractionsToDownSample = args.fractions.split(" ")
    for fraction in fractionsToDownSample:
        os.chdir(args.outDir)
        fractionOutDir = args.outDir + "/Fraction_" + fraction
        if os.path.isdir(fractionOutDir):
                if(args.verbose):
                    print "The output directory ", fractionOutDir, " exists"
        else:    
            if(args.verbose):
                    print"Making directory & changing directory to", fractionOutDir," to run the Downsampling process.\n"
            os.mkdir(fractionOutDir)
            call(['chmod', '755', fractionOutDir])
            
        DownSample(jobqueue,fraction,fileToProcess,fractionOutDir,args)
        RunHSmetrics(jobqueue,fraction,fileToProcess,fractionOutDir,args)

def DownSample(jobqueue,fraction,fileToProcess,fractionOutDir,args):
    os.chdir(fractionOutDir)
    jobs=[]
    rcount = 0
    for file in fileToProcess:
        rcount = rcount + 1 
        OutfileName = fractionOutDir + "/" + os.path.basename(file)
        if(os.path.isfile(OutfileName)):
            continue
        else:
            cmd = args.java + " -Xmx25g -jar " +  args.gatk + " -T PrintReads -R " + args.ref + " -I " + file + " -o " + OutfileName + " -dfrac " + fraction
            #print "cmd:",cmd,"\n"
            qsub_cmd = args.qsub + " -q " + args.queue + " -N " + "DownSample_" + str(rcount)+ "_" +str(fraction) + " -o " + "DownSample_"+ str(rcount)+"_"+str(fraction)+".stdout" + " -e " + "DownSample_"+str(rcount)+"_"+str(fraction)+".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp " + args.threads + " -wd " + fractionOutDir + " -sync y " + " -b y " + cmd 
            print "qsub_cmd:",qsub_cmd,"\n"
            jobs.append(qsub_cmd)
    '''iterate over jobs and put each into the queue in sequence'''
    for job in jobs:
        print "inserting job into the queue: %s"%(job)
        jobqueue.put(job)
    '''start some threads, each one will process one job from the queue'''
    for i in range(mp.cpu_count()-1):
        th = Thread(target=processor,args = (i,jobqueue,))
        th.setDaemon(True)
        th.start()

    '''wait until all jobs are processed before quitting'''
    jobqueue.join() 
    return
def RunHSmetrics(jobqueue,fraction,fileToProcess,fractionOutDir,args):
    os.chdir(fractionOutDir)
    jobs=[]
    rcount = 0
    for file in fileToProcess:
        rcount = rcount + 1 
        OutfileName = fractionOutDir + "/" + os.path.basename(file)
        hsMetricsOutFile = OutfileName.rstrip('.bam')
        hsMetricsOutFile = hsMetricsOutFile + ".HSmetrics.txt"
        PerTargetCovgOutFile = OutfileName.rstrip('.bam')
        PerTargetCovgOutFile = PerTargetCovgOutFile + ".target.covg"
        if(os.path.isfile(hsMetricsOutFile) and os.path.isfile(OutfileName)):
            continue
        else:
            hsmetricsJar = args.picard + "/" + "CalculateHsMetrics.jar"
            hs_cmd = args.java + " -Xmx8g -jar " + hsmetricsJar + " I=" + OutfileName + " O=" + hsMetricsOutFile + " TI=" + args.targetIntervals + " BI=" + args.baitIntervals + " R=" + args.ref + " PER_TARGET_COVERAGE=" + PerTargetCovgOutFile + " VALIDATION_STRINGENCY=LENIENT"
            qsub_hscmd = args.qsub + " -q " + args.queue + " -N " + "HSmetrics_" + str(rcount)+ "_" +str(fraction) + " -o " + "HSmetrics_"+ str(rcount)+"_"+str(fraction)+".stdout" + " -e " + "HSmetrics_"+str(rcount)+"_"+str(fraction)+".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp 2" + " -wd " + fractionOutDir + " -sync y " + " -b y " + hs_cmd 
            print "qsub_cmd:",qsub_hscmd,"\n"
            jobs.append(qsub_hscmd)
    '''iterate over jobs and put each into the queue in sequence'''
    for job in jobs:
        print "inserting job into the queue: %s\n"%(job)
        jobqueue.put(job)
    '''start some threads, each one will process one job from the queue'''
    for i in range(mp.cpu_count()-1):
        th = Thread(target=processor,args = (i,jobqueue,))
        th.setDaemon(True)
        th.start()

    '''wait until all jobs are processed before quitting'''
    jobqueue.join() 
    return
     
def DownSample_depricated(fraction,fileToProcess,fractionOutDir,args):
    os.chdir(fractionOutDir)
    rcount = 0
    for file in fileToProcess:
        rcount = rcount + 1 
        OutfileName = fractionOutDir + "/" + os.path.basename(file)
        if(os.path.isfile(OutfileName)):
            retcode = 1
        else:
            cmd = args.java + " -Xmx25g -jar " +  args.gatk + " -T PrintReads -R " + args.ref + " -I " + file + " -o " + OutfileName + " -dfrac " + fraction
            #print "cmd:",cmd,"\n"
            qsub_cmd = args.qsub + " -q " + args.queue + " -N " + "DownSample_" + str(rcount)+ "_" +str(fraction) + " -o " + "DownSample_"+ str(rcount)+"_"+str(fraction)+".stdout" + " -e " + "DownSample_"+str(rcount)+"_"+str(fraction)+".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp " + args.threads + " -wd " + fractionOutDir + " -sync y " + " -b y " + cmd 
            print "qsub_cmd:",qsub_cmd,"\n"
            qsub_args = shlex.split(qsub_cmd)
            proc = Popen(qsub_args)
            proc.wait()
            retcode = proc.returncode
        if(retcode >= 0):
            hsMetricsOutFile = OutfileName.rstrip('.bam')
            hsMetricsOutFile = hsMetricsOutFile + ".HSmertics.txt"
            PerTargetCovgOutFile = OutfileName.rstrip('.bam')
            PerTargetCovgOutFile = PerTargetCovgOutFile + ".target.covg"
            if(os.path.isfile(hsMetricsOutFile) and os.path.isfile(OutfileName)):
                retcode = 1
            else:
                hsmetricsJar = args.picard + "/" + "CalculateHsMetrics.jar"
                hs_cmd = args.java + " -Xmx8g -jar " + hsmetricsJar + " I= " + OutfileName + " O= " + hsMetricsOutFile + " TI=" + args.targetIntervals + " BI=" + args.baitIntervals + " R=" + args.ref + " PER_TARGET_COVERAGE=" + PerTargetCovgOutFile + "VALIDATION_STRINGENCY=LENIENT"
                qsub_hscmd = args.qsub + " -q " + args.queue + " -N " + "HSmetrics_" + str(rcount)+ "_" +str(fraction) + " -o " + "HSmetrics_"+ str(rcount)+"_"+str(fraction)+".stdout" + " -e " + "HSmetrics_"+str(rcount)+"_"+str(fraction)+".stderr" + " -V -l h_vmem=6G,virtual_free=6G -pe smp 1" + " -wd " + fractionOutDir + " -sync y " + " -b y " + cmd 
                hs_args = shlex.split(qsub_hscmd)
                proc = Popen(hs_args)
                proc.wait()
                retcode = proc.returncode
                if(retcode >= 0):
                    if(args.verbose):
                        print "Finished Running HSmetrics\n"
                else:
                    if(args.verbose):
                        print "HSmetrics is either still running or it errored out with return code", retcode,"\n"        
        else:
            if(args.verbose):
                print "PrintReads is either still running or it errored out with return code", retcode,"\n" 

def processor(i,jobqueue):
    if jobqueue.empty() == True:
        print "the Queue is empty!"
        sys.exit(1)
    try:
        job = jobqueue.get()
        print "I'm operating on job item: %s\n"%i,job
        qsub_args = shlex.split(job)
        Popen(qsub_args)
        jobqueue.task_done()
    except:
        print "Failed to operate on job\n"

          
if __name__ == "__main__":
    start_time = time.time() 
    mp.freeze_support() 
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))