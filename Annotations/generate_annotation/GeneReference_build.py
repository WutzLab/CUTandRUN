#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Apr 24 18:49:14 2015

@author: GDM
"""
import time
import datetime


def getDay():
    return (datetime.date.today()).isoformat()
    
def getCurrTime():
        currTime = time.strftime('%y-%m-%d %H:%M:%S', time.localtime(time.time()))
        return currTime

def computeRunTime(startTime, endTime):
    start = datetime.datetime.strptime(startTime, '%y-%m-%d %H:%M:%S')
    start_sec_float = time.mktime(start.timetuple())
    
    end = datetime.datetime.strptime(endTime, '%y-%m-%d %H:%M:%S')
    end_sec_float = time.mktime(end.timetuple())

    return end_sec_float - start_sec_float

import pandas as pd
import HTSeq
import argparse
import pickle as pickle
# from HaSAPPy.HaSAPPY_time import *


def library_preparation(info):
    ################
    def exonInterval_define(transcript, ref_gene):
        exonInterval = []
        starts = ref_gene.loc[transcript, "exonStarts"].split(",")
        ends = ref_gene.loc[transcript, "exonEnds"].split(",")
        for n in range(ref_gene.loc[transcript, "exonCount"]):
            interval = HTSeq.GenomicInterval(ref_gene.loc[transcript, "chrom"], int(starts[n]), int(ends[n]), ref_gene.loc[transcript, "strand"])
            exonInterval.append(interval)
        return exonInterval
        
    ####
    def define_introns(index, library):
        first = True
        introns = []
        for exon_interval in library.loc[index, "exon_specific"]:
            if first:
                _5 = exon_interval.end
            else:
                _3 = exon_interval.start
                intron_interval = HTSeq.GenomicInterval(exon_interval.chrom, _5, _3, exon_interval.strand)
                introns.append(intron_interval)
                _5 = exon_interval.end
            first = False
        return introns
    ##################
    print('\t- Loading .txt file')
    ref_gene = pd.read_table(info.input, sep="\t")

    for index in ref_gene.index:
        ref_gene.loc[index, "interval"] = HTSeq.GenomicInterval(ref_gene.loc[index, "chrom"], ref_gene.loc[index, "txStart"], ref_gene.loc[index, "txEnd"], ref_gene.loc[index, "strand"])

    gene_collection = ref_gene.groupby(ref_gene["name2"])
    gene_collection = pd.Series(gene_collection.groups)

    columns = ["reference", "genomic_interval", "variants", "exons", "exon_specific", "introns_all"]
    
    gene_data = pd.DataFrame(columns= columns)
    gene_data.index.name = "genes"
    print('\t- Redundant gene intervals management')

    for index in gene_collection.index:
        list_interval = []
        first = True
        for n in gene_collection[index]:
            if first:
                list_interval.append((ref_gene.loc[n, "interval"], [n]))
                first = False
            else:
                count = 0
                finded = False
                for interval in list_interval:    
                    if ref_gene.loc[n, "interval"].overlaps(interval[0]):
                        interval[0].extend_to_include(ref_gene.loc[n, "interval"])
                        list_interval[count][1].append(n)
                        finded = True
                        break
                    count += 1
                if not finded:
                    list_interval.append((ref_gene.loc[n, "interval"], [n]))
                        
        if len(list_interval) == 1:
            gene_data.loc[index, "reference"] = list_interval[0][1]
            gene_data.loc[index, "genomic_interval"] = list_interval[0][0]
        else:
            num = len(list_interval)
            for gene in list_interval:
                index_1 = index + "_" + str(num)
                num -= 1
                gene_data.loc[index_1, "reference"] = gene[1]
                gene_data.loc[index_1, "genomic_interval"] = gene[0]

    print('\t- Define exonic regions')
    
    gene_data["variants"] = [{} for gene in range(len(gene_data))]

    for index in gene_data.index:
        for num in gene_data.loc[index, "reference"]:
            exonInterval = exonInterval_define(num, ref_gene)
            transcript = ref_gene.loc[num, "name"]
            gene_data.loc[index, "variants"][transcript] = exonInterval

    print('\t- Define intronic regions')        

    for index in gene_data.index:
        transcripts = gene_data.loc[index, "variants"]
        gene_exon = HTSeq.GenomicArray("auto", stranded=True)
        gene_all_exon = []
        gene_specific_exon = []
        for transcript in transcripts:
            for interval in transcripts[transcript]:
                gene_exon[interval]+=1
        for interval in gene_exon.steps():
            if interval[1] > 0:
                gene_all_exon.append(interval[0])
                if interval[1] == len(transcripts):
                    gene_specific_exon.append(interval[0])
        gene_data.at[index, "exons"] = gene_all_exon
        gene_data.at[index, "exon_specific"] = gene_specific_exon

    for index in gene_data.index:
        introns = define_introns(index, gene_data)
        gene_data.at[index, "introns_all"] = introns

    print('\t- Saving file in: %s' % info.output)

    with open(info.output, 'wb') as write:
        pickle.dump(gene_data, write)

    print(('Gene models library generated\n\tRunTime: %s' %(computeRunTime(startTime, getCurrTime()))))


parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', help='Provide .txt file containing information on genes model', action='store')
parser.add_argument('-o', '--output', help='Provide .pkl file PATH where library will be stored', action='store')

args = parser.parse_args()

print('\n***Generation of gene models library for HaSAPPy program***')
print('\tInput file: %s' % args.input)

if args.input==None or args.output == None:
    print('\nWARNING: informations provided are not sufficent.\nCheck -h option to have more details on requested parameters')

else:
    startTime = getCurrTime()
    print('Started: %s' % startTime)
    library_preparation(args)
 
