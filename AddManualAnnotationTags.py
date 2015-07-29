''''
Created on 11/03/2014
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
import copy
def main():
    parser = argparse.ArgumentParser(prog='AddManualAnnotationTags.py', description='Add Tag for all different annotations', usage='%(prog)s [options]')
    parser.add_argument("-t", "--titleFile", action="store", dest="titleFile", required=True, metavar='title_file.txt', help="Full path title file.") 
    parser.add_argument("-i", "--inputFile", action="store", dest="inputFile", required=True, metavar='/somepath/input.txt', help="Full Path to the input file.")
    parser.add_argument("-hopl", "--hotspotList", action="store", dest="hotspotList", required=False, default='/dmp/data/mskdata/two-tier-filters/production/hotspot_sqm_val.txt', metavar='/some/path/hotspot.list', help="Full Path to HotSpot List")
    parser.add_argument("-ce", "--ClinicalExons", action="store", dest="ClinicalExons", required=False, default='/dmp/data/mskdata/two-tier-filters/production/clinical_exons.txt', metavar='/some/path/ClinicalExons.list', help="Full Path to Clinical Exons List")
    parser.add_argument("-p", "--panelGeneList", action="store", dest="panelGeneList", required=False, default='/dmp/data/mskdata/interval-lists/production/gene_intervals.list.annotated', metavar='/some/path/genesInThePanel.list', help="Full Path to Genes in the Panel List")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
    parser.add_argument("-o", "--outDir", action="store", dest="outDir", required=True, metavar='/somepath/output', help="Full Path to the output dir.")
    parser.add_argument("-of", "--outfile", action="store", dest="outfile", required=False, default='AllSample_MergedResults.txt', metavar='output.txt', help="Name of the output file.")
    parser.add_argument("-ad", "--alleledepth", action="store", dest="ad", required=False, default='3', metavar='3', help="Tumor allele depth threshold")
    parser.add_argument("-tpvf", "--tnVFProportion", action="store", dest="tpvf", required=False, default='10',metavar='5', help="Tumor-Normal variant frequency proportion threshold ")
    parser.add_argument("-vf", "--variantfrequency", action="store", dest="vf", required=False, default='0.01', metavar='0.01', help="Tumor variant frequency threshold ")
    parser.add_argument("-tnr", "--tnRatio", action="store", dest="tnr", required=False, default='5', metavar='5', help="Tumor-Normal variant frequency ratio threshold ")
    args = parser.parse_args()
    if(args.verbose):
        print "Checking the and processing Arguments.\n"
    (inputFileDF,allIds,normalIds,hotspotSNV,hotspotExons,geneList) = ProcessArgs(args)
    if(args.verbose):
        print "Annotating the Fields\n"
    RunAnnotations(inputFileDF,allIds,normalIds,hotspotSNV,hotspotExons,geneList,args)
    if(args.verbose):
        print "Finished running Annotation\n"
    sys.exit()
def ProcessArgs(args):
    (normalIds,allIds) = ReadTitleFile(args)
    (hotspotSNV)= ReadHotspotFile(args)
    (hotspotExons) = ReadClinicalExonsFile(args)
    (geneList) = ReadGeneListFile(args)
    (inputFileDF) = ReadInputFile(args)
    return(inputFileDF,allIds,normalIds,hotspotSNV,hotspotExons,geneList)

def ReadTitleFile(args):
    titleFileDF = pd.read_csv(args.titleFile, sep='\t', header=0, keep_default_na='True')
    normalIds = []
    allIds = []
    for count, row in titleFileDF.iterrows():
        bcId = row.loc['Barcode']
        poolId = row.loc['Pool']
        sampleId = row.loc['Sample_ID']
        sampleClass = row.loc['Class']
        if(sampleClass == 'Normal'):
            normalIds.append(sampleId)
        allIds.append(sampleId)
    return(normalIds,allIds)            
    
def ReadHotspotFile(args):
    hotspotSNV={}
    if(args.verbose):
        print "Reading Hotspot List File\n"
    with open(args.hotspotList, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            data = line.rstrip('\n').split('\t')
            gene = data[1]
            mutationInfo = data[2].split(' ')
            mutation = mutationInfo[0]
            #Check SNV
            snvMatch = re.compile(r'^([A-Z]\d+)\>*[A-Za-z]*.*')
            snvMatchM = snvMatch.match(mutation)
            if snvMatchM:
                mutation = snvMatchM.group(1)
                #print "Gene:",gene,"Mutation:",mutation,"\n"   
                key = gene + ":" + mutation
                hotspotSNV[key] = mutation  
    return(hotspotSNV)
'''
            #Check FrameShift Mutations
            fsMatch = re.match(r'.*\d+fs',mutation)
            if fsMatch:
                mutation = fsMatch.group(0)
                #print "Gene:",gene,"Mutation:",mutation,"\n"
            #Check Insertions
            insMatch = re.compile(r'^(.*?)\>.*')
            insMatchM = insMatch.match(mutation)
            if insMatchM:
                mutation = insMatchM.group(1)
                #print "Gene:",gene,"Mutation:",mutation,"\n"   
            #Check Deletion
            delMatch =  re.compile(r'^(.*?)del')
            delMatchM = delMatch.match(mutation)
            if delMatchM:
                mutation = delMatchM.group(1)
                #print "Gene:",gene,"Mutation:",mutation,"\n"   
'''

def ReadClinicalExonsFile(args):
    hotspotExons={}
    if(args.verbose):
        print "Reading Clinical Exon's List File\n"
    with open(args.ClinicalExons, 'r') as filecontent:
        header = filecontent.readline()
        for line in filecontent:
            (gene,exon,type) = line.rstrip('\n').split('\t')
            key =  gene + ":" + exon 
            hotspotExons[key] = type
    return(hotspotExons)

def ReadGeneListFile(args):
    geneList={}
    if(args.verbose):
        print "Reading Gene list file\n" 
        with open(args.panelGeneList, 'r') as filecontent:
            for line in filecontent:
                (coord,arm,info) = line.rstrip('\n').split('\t')
                infoSplit = info.split(':')
                gene = infoSplit[0]
                geneList[gene] = 1
    return(geneList)

def ReadInputFile(args):
    if(args.verbose):
        print "Reading Input file\n" 
    InputFileDF = pd.read_csv(args.inputFile, sep='\t', header=0, keep_default_na='True', dtype='object')
    return(InputFileDF)

def RunAnnotations(inputFileDF,allIds,normalIds,hotspotSNV,hotspotExons,geneList,args):
    if(args.verbose):
        print "Adding Annotations to the dataframe\n"
    txt_fh = open(args.outfile, "wb") 
    txt_fh.write("Sample\tChr\tPos\tRef\tAlt\tDP_T\tAD_T\tVF_T\tDP_N\tAD_N\tVF_N\tRatio\tOccuranceIn100prcntBams\tOccuranceInNormals\tVarType\tExonic\tAnnotation\tHotspot\tPanel\n")    
    for count, row in inputFileDF.iterrows():
        SampleVal = str(row.loc['Sample'])
        NormalUsedVal = str(row.loc['NormalUsed'])
        ChromVal = str(row.loc['Chrom'])
        StartVal = str(row.loc['Start'])
        RefVal = str(row.loc['Ref'])
        AltVal = str(row.loc['Alt'])
        VCVal = str(row.loc['VariantClass'])
        ExonVal = str(row.loc['Exon'])
        GeneVal = str(row.loc['Gene'])
        TsIDVal = str(row.loc['TranscriptID'])
        CDNAVal = str(row.loc['cDNAchange'])
        AAcVal =  str(row.loc['AAchange'])
        DP_T = int(row.loc['T_TotalDepth'])
        AD_T = int(row.loc['T_AltCount'])
        VF_T = float(row.loc['T_AltFreq'])
        DP_N = int(row.loc['N_TotalDepth'])
        AD_N = int(row.loc['N_AltCount'])
        VF_N = float(row.loc['N_AltFreq'])
        OccurenceIn100pct = 0
        for ids in allIds:
            newId = ids + "-100"
            val = row.loc[newId]
            #print "VAL:",val
            (DPVal,RDVal,ADVal,VFVal) = val.rstrip().split(';')
            DP = DPVal.split('=')[-1]
            AD = ADVal.split('=')[-1]
            VF = VFVal.split('=')[-1]
            TPVF = int(VF_T)/float(args.tpvf)
            if TPVF < 0.01:
                TPVF = 0.01
            if((int(AD) >= int(args.ad)) and (float(VF) >= float(TPVF))):
                OccurenceIn100pct = OccurenceIn100pct + 1
        OccurenceInNormals = 0
        for ids in normalIds:
            val = row.loc[ids]
            (DPVal,RDVal,ADVal,VFVal) = val.rstrip().split(';')
            DP = DPVal.split('=')[-1]
            AD = ADVal.split('=')[-1]
            VF = VFVal.split('=')[-1]
            #print "DP:",DP,"AD:",AD,"VF:",VF
            TPVF = int(VF_T)/float(args.tpvf)
            if TPVF < 0.01:
                TPVF = 0.01
            if((int(AD) >= int(args.ad)) and (float(VF) >= float(TPVF))):
                OccurenceInNormals = OccurenceInNormals + 1
        Ratio = 0
        if(VF_N > 0.0):
            Ratio = float(VF_N)/float(args.tnr)
        else:
            Ratio = 5
        Hotspot = ''
        Panel = ''
        Exonic = ''
        if ((VCVal == 'frameshift_deletion') or (VCVal == 'frameshift_insertion') or (VCVal == 'nonframeshift_deletion') or (VCVal == 'nonframeshift_insertion') or (VCVal == 'nonsynonymous_SNV') or (VCVal == 'splicing') or (VCVal == 'stopgain_SNV')):
            Exonic = 'exonic-nonsynonymous'
        else:
            Exonic = 'exonic-synonymous'
        if(GeneVal in geneList):
            Panel = 'Panel'
        else:
            Panel = 'NonPanel'
        if((len(RefVal) > len(AltVal)) or (len(RefVal) < len(AltVal))):
            key = GeneVal + ":" + ExonVal
            if key in hotspotExons:
                Hotspot = 'Hotspot'
            else:
                Hotspot = 'NotHotpot'
        if((len(RefVal) == len(AltVal)) and (len(RefVal)==1)):
            AAcMatch = re.compile('^p\.([A-Z]\d+).*')
            AAcMatchM = AAcMatch.match(AAcVal)
            if AAcMatchM:
                key = GeneVal + ":" + AAcMatchM.group(1)
                #print key,"\n"
                if key in hotspotSNV:
                    Hotspot = 'Hotspot'
                else:
                    Hotspot = 'NotHotpot'
            else:
                Hotspot = 'NotHotpot'
        VCValcopy = copy.copy(VCVal)
        VCValcopy.replace('_',' ')
        Annotations = GeneVal + ":" + TsIDVal + ":" + ExonVal + ":" + CDNAVal + ":" + AAcVal
        txt_fh.write(SampleVal + "\t" + ChromVal + "\t" + StartVal+"\t"+RefVal+"\t"+AltVal+"\t"+ str(DP_T) + "\t" + str(AD_T) + "\t" + str(VF_T) + "\t" + str(DP_N) + "\t" + str(AD_N) + "\t" + str(VF_N) + "\t" + str(Ratio) + "\t" + str(OccurenceIn100pct) + "\t" + str(OccurenceInNormals) + "\t"+VCValcopy+"\t"+Exonic+"\t"+Annotations+"\t"+Hotspot+"\t"+Panel+"\n")
    txt_fh.close()
    return

if __name__ == "__main__":
    start_time = time.time() 
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))   