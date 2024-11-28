# Module xCHiPPy
# Analysis of ChIPSeq data sets using genome annotation GTF and HTSeq read count data
# CHANGES:
# 1. print ... TO print(...)
# 2. pd.DataFrame.ix TO .loc
# 3. pd.DataFrame.from_csv( TO pd.read_csv(
# 
# cool names include: 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.colors as colors
import pandas as pd
import seaborn as sns
import HTSeq
import pickle
import json
import itertools
import re
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap

class XchipAnalysis(object):
    """
    This class represents a ChIP analysis project. Gene annotation data in a pandas dataframe to which ChIPseq data can be added.
    Several methods work on these data and include normalization functions as well as import and export tools. For convenience conversion
    functions to and from pandas DataFrames and CSV files are implemented. Heatmaps can be generated for various aspects of the dataset.

    Data is held in the instance variables:

        Annotation (pandas dataframe) :  A HaSAPPy formated pandas dataframe that can be accessed with gene names.
        Experiments (dict) :  A dictionary containing GenomicArrays with read coverage of the experiments. Keys are the experiment names.
        ExperimentStatistics (dict): A dictionary that contains statistics of the experiment. It is accessed via experiment names as keys,
                              and returnes a dictionary for each experiment with several values that contains data for the following keys:
                              "totalReads", "alignedReads", "collectedReads", "readsX", "readsA", "coverageX", "coverageA"
                              Use normalizeExperiment(normFactor) to set a normalization Factor for the experiment that will be used for
                              normatization by multiplying coverage values with normFactor.

    Settings for controling behavious by instance variables:
    
        printVerboseLog (boolean) : Controls error and warning messages printed to the screen (default = True)

    """
    RBGcmap   = LinearSegmentedColormap.from_list("mycmap",['#aa0000','#000000','#00aa00'])
    REDcmap   = LinearSegmentedColormap.from_list("myreds",['#ffffff','#ff0000'])
    GREENcmap = LinearSegmentedColormap.from_list("mygreens",['#ffffff','#00ff00'])

    def __init__(self, Annotation_filename=None):
        """
        Instantiate a new analysis and read in a genome annotation from a pickl file, if Annotation_filename is provided
        must be a valid path to a genome annotation pickl file in HaSAPPY format.
        """
        self.printVerboseLog=True
        if Annotation_filename is None:
            self.Annotation = None
        else:
            self.readAnnotation(Annotation_filename)
        self.Experiments=dict()
        self.ExperimentStatistics=dict()

    def readAnnotation(self, filename):
        """
        Reads genome annotation from a HaSAPPy formatted pickl file specified by filename.
        """
        with open(filename,'rb') as annotation_file:
            self.Annotation = pickle.load(annotation_file)

    def getAnnotation(self):
        """
        Returns a refernce to the current Annotation.
        """
        return self.Annotation

    def get_x_linked_genes(self):
        """
        Returns a list of gene names as string for which chrom=='chrX'.
        Note chrom is checked case insensitive by this function.
        """
        x_linked_genes_list=[]
        for g in self.Annotation.index:
            if self.Annotation.loc[g,'genomic_interval'].chrom.upper()=='CHRX':
                x_linked_genes_list.append(g)
        return(x_linked_genes_list)

    def get_autosomal_genes(self):
        """
        Returns a list of gene names as string for which chrom starts with 'chr' but is not 'chrX'.
        Note chrom is checked case insensitive by this function.
        """
        autosomal_genes_list=[]
        for g in self.Annotation.index:
            chrom=self.Annotation.loc[g,'genomic_interval'].chrom.upper()
            if chrom[:3]=='CHR' and chrom[:4]!='CHRX':
                autosomal_genes_list.append(g)
        return(autosomal_genes_list)

    def get_genes_TSSpos(self, chrom):
        """
        Returns a list of tuples containing gene names as string and the position of the TSS as number
        for genes with matching chrom. Note chrom is checked case insensitive by this function.
        This function is useful for generating gene name lists of genes sorted by chromosomal position.
        """
        genes_list=[]
        for g in self.Annotation.index:
            ch=self.Annotation.loc[g,'genomic_interval'].chrom.upper()
            if chrom.upper()==ch:
                TSS=self.get_TSS(gene)
                genes_list.append(g, TSS.start_d)
        return(genes_list)

    def get_genomic_interval(self, gene):
        """
        Returns the genomic interval for the gene identified by gene, which contains a string with the gene name.
        """
        return self.Annotation.loc[gene,'genomic_interval']

    def get_TSS(self, gene, get_TES=False):
        """
        Returns the HTSeq genomic position for the start site of gene given by its name as a string.
        """
        if get_TES:
            return self.Annotation.loc[gene, 'genomic_interval'].end_d_as_pos
        else:
            return self.Annotation.loc[gene, 'genomic_interval'].start_d_as_pos
        # TSS = HTSeq.GenomicPosition(self.Annotation.loc[gene,'genomic_interval'].chrom, self.Annotation.loc[gene,'genomic_interval'].start_d, self.Annotation.loc[gene,'genomic_interval'].strand)
        # return TSS

    def extend_genomic_interval(self, genomic_interval, delta):
        """
        Returns a new genomic interval that has been extended for a number of bases provided by delta on both sides. The position is floored to 0 to avoid negative genomic
        coordinates, in which case the returned GenomicInterval has a smaller size and is not symmetric.
        """
        if genomic_interval.start > delta:
            st=genomic_interval.start-delta
        else:
            st=0
        gi_ext=HTSeq.GenomicInterval(genomic_interval.chrom, st, genomic_interval.end+delta, genomic_interval.strand)
        return gi_ext

    def promoter_genomic_interval(gene, length):
        """
        Returns a genomic interval either of length symmetrically positioned on the start of a gene, if length is given by an int.
        Else length can be a tupple or a list, whereby index [0] indicates the left length and index [1] the right length of the
        genomic intervall from the TSS. This can be useful for perfoming analyses of different intervals around gene promoters.
        The gene is identified through its gene name supplied as a string.
        """
        if type(length) is int:
            left=length/2
            right=length/2
        else:
            left=length[0]
            right=length[1]
        genomic_interval=self.get_genomic_interval(gene)
        if genomic_interval.start > left:
            st=genomic_interval.start-left
        else:
            st=0
        gi_promoter=HTSeq.GenomicInterval(genomic_interval.chrom, genomic_interval.start_d-st, genomic_interval.start_d+right, genomic_interval.strand)
        return gi_promoter

    def addExperiment(self, experiment_name, filename, aQualValue=0, frgmtLength=None, maxReadsPerPos=None):
        """
        Reads in a dataset from a BAM file identified by filename and adds it to the Experiments dictionary of the current XchipAnalysis.
        Various parameters are added to ExperimentStatistics that can be used for assessing dataset quality and normalization.
        A mapping quality value can be specified by aQualValue to select alignments equal or above a certain treshold quality.
        Returns three read count values of the total number of reads in the BAM file, the number of aligned reads, and the number of reads
        that have alignment scores above aQualValue. For eliminating reads that start at the same position and are likely PCR duplicates
        the maximum number of reads that will be accepted per genomic position (stranded) can be provided via maxReadsPerPos.
        If maxReadsPerPos=None (default) all reads will be used for calculating the coverage.
        """
        aligned_file = HTSeq.BAM_Reader(filename)
        count_aligned = 0
        count_GoodQualityAlignment = 0
        count_total = 0
        count_X = 0
        count_A = 0
        coverage_X = 0
        coverage_A = 0

        # Calculate coverage in a non-stranded integer GenomicArray
        readsLib = HTSeq.GenomicArray('auto', stranded=False, typecode='i') # build a GenomicArray holding the read coverage as overlapping reads per base
        if maxReadsPerPos is None:
            for algnt in aligned_file:
                if algnt.aligned:
                    if count_aligned == 0:
                        if algnt.iv.chrom.startswith('chr'):
                            chromosome_style = ''
                        else:
                            chromosome_style = 'chr'
                    if algnt.aQual >= aQualValue:    # only consider alignments equal or above the selected aQualValue
                        if frgmtLength is not None:
                            if algnt.iv.strand=="-":
                                if algnt.iv.start_d > frgmtLength:
                                    algnt.iv.length=frgmtLength
                            else:
                                algnt.iv.length=frgmtLength
                        read = HTSeq.GenomicInterval('%s%s' %(chromosome_style,str(algnt.iv.chrom)),algnt.iv.start,algnt.iv.end) # we need a non stranded interval
                        readsLib[read] += 1
                        count_GoodQualityAlignment += 1
                        if read.chrom == "chrX":
                            count_X += 1
                            coverage_X += read.length
                        else:
                            count_A += 1
                            coverage_A += read.length
                    count_aligned += 1
                count_total += 1
        else:
            strtLib = HTSeq.GenomicArray('auto', stranded = True, typecode='i') # build a GenomicArray holding the read start positions for eliminating PCR duplicates
            for algnt in aligned_file:
                if algnt.aligned:
                    if count_aligned == 0:
                        if algnt.iv.chrom.startswith('chr'):
                            chromosome_style = ''
                        else:
                            chromosome_style = 'chr'                    
                    if algnt.aQual >= aQualValue:    # only consider alignments equal or above the selected aQualValue
                        if strtLib[algnt.iv.start_d_as_pos] < maxReadsPerPos:
                            strtLib[algnt.iv.start_d_as_pos] += 1
                            if frgmtLength is not None:
                                if algnt.iv.strand=="-":
                                    if algnt.iv.start_d > frgmtLength:
                                        algnt.iv.length=frgmtLength
                                else:
                                    algnt.iv.length=frgmtLength
                            read = HTSeq.GenomicInterval('%s%s' %(chromosome_style,str(algnt.iv.chrom)),algnt.iv.start,algnt.iv.end)
                            readsLib[read] += 1
                            count_GoodQualityAlignment += 1
                            if read.chrom == "chrX":
                                count_X += 1
                                coverage_X += read.length
                            else:
                                count_A += 1
                                coverage_A += read.length
                    count_aligned += 1
                count_total += 1
        self.Experiments[experiment_name]=readsLib
        self.ExperimentStatistics[experiment_name]={"totalReads": count_total, "alignedReads": count_aligned, "collectedReads": count_GoodQualityAlignment, "readsX": count_X, "readsA": count_A, "coverageX": coverage_X, "coverageA": coverage_A}
        if self.printVerboseLog==True:
            string = '\t-Total reads: %i\n\t-Aligned reads: %i\n\t-Aligned Reads trusted: %i\n' %(count_total,count_aligned,count_GoodQualityAlignment)
            print(string)
        return count_total, count_aligned, count_GoodQualityAlignment
    
    
    def normalizeExperiment(self, experiment_name, normFactor):
        """
        Sets a normalization Factor for the experiment in ExperimentStatistics that will be used for normalizing coverage values through multiplication with normFactor.
        """
        self.ExperimentStatistics[experiment_name]["normFactor"] = normFactor
        return

    def getExperiment(self, experiment_name):
        """
        Returns the GenomicArray holding the read coverage for the experiment identified by experiment_name.
        """
        return self.Experiments[experiment_name]

    def getExperimentNormFactor(self, experiment_name):
        """
        Returns the current normalization factor of the experiment if it exists, else None.
        """
        ret = None
        if "normFactor" in self.ExperimentStatistics[experiment_name]:
            ret = self.ExperimentStatistics[experiment_name]["normFactor"]
        return ret

    def getExperimentList(self):
        """
        Returns a list of all experiment names in the current analysis as strings.
        """
        return self.Experiments.keys()

    def generateTSSCoverageProfile(self, experiment_name, gene_list, interval, analyse_TES=False):
        """
        Returns a pandas dataframe containing the promoter read coverage for the genes in gene_list identified by gene names, and
        a numpy ndarray with the sum of the coverage. Columns contain read coverage per base around the start (TSS) of a gene and
        are named with the base position relative to TSS. Rows are named with the gene names. The interval is pecified by a tuple
        the contain the left [0] and right [1] length from TSS. IF "normFactor" is defined in ExperimentStatistics it is used to
        multiply the coverage count values for normalization.
        """
        ilen = interval[1]-interval[0]
        if ilen<=0:
            return None
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        TSScoverage = pd.DataFrame(columns= [x for x in range(interval[0],interval[1])])
        TSSprofile = np.zeros(ilen, dtype=float)
        for gene in gene_list:
            TSS = self.get_TSS(gene, get_TES=analyse_TES) # HTSeq.GenomicPosition(self.Annotation.loc[gene,'genomic_interval'].chrom, self.Annotation.loc[gene,'genomic_interval'].start_d, self.genesAnnotation.loc[gene,'genomic_interval'].strand)
            if TSS.strand == '+':
                st=TSS.start+interval[0]
                en=TSS.start+interval[1]
                if st<0 or en<0:
                    continue
                TSSinterval = HTSeq.GenomicInterval(TSS.chrom, st, en, ".")
                wincvg = np.fromiter(data[TSSinterval], dtype="i", count=ilen)
                if normFactor is not None:
                    wincvg*=normFactor
                TSScoverage.loc[gene] = wincvg
                TSSprofile += wincvg
            else:
                st=TSS.start-interval[1]
                en=TSS.start-interval[0]
                if st<0 or en<0:
                    continue
                TSSinterval = HTSeq.GenomicInterval(TSS.chrom, st, en, ".")
                wincvg = np.fromiter(data[TSSinterval], dtype=float, count=ilen)
                if normFactor is not None:
                    wincvg*=normFactor
                TSScoverage.loc[gene] = wincvg[::-1]
                TSSprofile += wincvg[::-1]
        return TSScoverage, TSSprofile
    
    def generateTSSProfile(self, experiment_name, gene_list, interval, analyse_TES=False):
        """
        Returns a numpy ndarray with the sum of the promoter read coverage for the genes in gene_list identified by gene names.
        Elements contain read coverage per base around the start (TSS) of a gene. If "normFactor" is defined in ExperimentStatistics
        it is used to multiply the coverage count values for normalization.
        Use this to conserve memory or when the coverage for individual genes is not needed, else see generateTSSCoverageProfile().
        """
        ilen = interval[1]-interval[0]
        if ilen<=0:
            return None
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        TSSprofile = np.zeros(ilen, dtype=float)
        for gene in gene_list:
            TSS = self.get_TSS(gene, get_TES=analyse_TES) # HTSeq.GenomicPosition(self.Annotation.loc[gene,'genomic_interval'].chrom, self.Annotation.loc[gene,'genomic_interval'].start_d, self.genesAnnotation.loc[gene,'genomic_interval'].strand)
            if TSS.strand == '+':
                st=TSS.start+interval[0]
                en=TSS.start+interval[1]
                if st<0 or en<0:
                    continue
                TSSinterval = HTSeq.GenomicInterval(TSS.chrom, st, en, ".")
                wincvg = np.fromiter(data[TSSinterval], dtype="i", count=ilen)
                TSSprofile += wincvg
            else:
                st=TSS.start-interval[1]
                en=TSS.start-interval[0]
                if st<0 or en<0:
                    continue
                TSSinterval = HTSeq.GenomicInterval(TSS.chrom, st, en, ".")
                wincvg = np.fromiter(data[TSSinterval], dtype="i", count=ilen)
                TSSprofile += wincvg[::-1]
        if normFactor is not None:
            TSSprofile = TSSprofile * normFactor
        return TSSprofile

    def generateGeneProfile(self, experiment_name, gene_list, numIntervals=1000):
        """
        Returns a numpy ndarray with the sum of the read coverage for the genes in gene_list identified by gene names. The transcription
        unit of the genes is scaled to provide numIntervals points that represent average coverage for equal sized intervals over the gene.
        If the gene is smaller than numIntervals it is skipped from the analysis. The total number of genes used for calculating the sum is returned.
        Columns contain read coverage per base around the start (TSS) of a gene and are named with the base position relative to
        TSS. Rows are named with the gene names. The interval is pecified by a tuple the contain the left [0] and right [1] length
        from TSS. IF "normFactor" is defined in ExperimentStatistics it is used to multiply the coverage count values for normalization.
        Use this to conserve memory or when the coverage for individual genes is not needed, else see generateTSSCoverageProfile().
        """
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        GeneProfile = np.zeros(numIntervals, dtype=float)
        numGenes = 0
        for gene in gene_list:
            gene_iv = self.get_genomic_interval(gene)
            if gene_iv.length > numIntervals:
                interval_length = int(gene_iv.length/numIntervals)
                numGenes += 1
                wincvg = np.fromiter(data[gene_iv], dtype="i", count=gene_iv.length)
                if gene_iv.strand == '-':
                    rcvg=wincvg[::-1]
                    wincvg=rcvg[:]
                wincvg=np.resize(wincvg, [numIntervals, interval_length])
                avg = np.mean(wincvg, axis=1)
                GeneProfile += avg
        if normFactor is not None:
            GeneProfile = GeneProfile * normFactor
        return GeneProfile, numGenes

    def getTSScoverageArray(self, experiment_name, gene, interval, analyse_TES=False):
        """
        Returns a numpy ndarray with the read coverage for the TSS of a gene for an experiment identified by experiment_name.
        Elements contain read coverage per base around the start (TSS) of a gene. If "normFactor" is defined in ExperimentStatistics
        it is used to multiply the coverage count values for normalization. The interval is specified by a tuple containing 5' and
        3' border. If the interval would at negative genomic coordinates, or its size is equal to or less than 0 None is returned.
        """
        ilen = interval[1]-interval[0]
        if ilen<=0:
            return None
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        TSS = self.get_TSS(gene, get_TES=analyse_TES)
        if TSS.strand == '+':
            st=TSS.start+interval[0]
            en=TSS.start+interval[1]
            if st<0 or en<0:
                return None
            TSSinterval = HTSeq.GenomicInterval(TSS.chrom, st, en, ".")
            TSScoverage = np.fromiter(data[TSSinterval], dtype=float)
        else:
            st=TSS.start-interval[1]
            en=TSS.start-interval[0]
            if st<0 or en<0:
                return None
            TSSinterval = HTSeq.GenomicInterval(TSS.chrom, st, en, ".")
            wincvg = np.fromiter(data[TSSinterval], dtype=float)
            TSScoverage = wincvg[::-1]
        if normFactor is not None:
            TSScoverage = TSScoverage * normFactor
        return TSScoverage
    
    def getGeneCoverageArray(self, experiment_name, gene):
        """
        Returns a numpy ndarray with the read coverage for a gene in the experiment identified by experiment_name. Elements
        contain read coverage per base of the transcription unit from the start (TSS) of a gene. If a "normFactor" is defined in
        ExperimentStatistics it is used to multiply the coverage count values for normalization.
        """
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        gene_iv = self.get_genomic_interval(gene)
        wincvg = np.fromiter(data[gene_iv], dtype="i", count=gene_iv.length)
        if gene_iv.strand == '-':
            GeneCvg=wincvg[::-1]
        else:
            GeneCvg=wincvg
        if normFactor is not None:
            GeneCvg = GeneCvg * normFactor
        return GeneCvg

    def getGeneInfo(self, experiment_name, gene, TSSinterval=(-5000,5000)):
        """
        Returns genomic coordinates and size of the gene locus and the average read coverage of the TSS interval and the gene
        transcription unit for the experiment identified by experiment_name. If a "normFactor" is defined in ExperimentStatistics
        it is used to multiply the average coverage count values for normalization.
        """
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        gene_iv = self.Annotation.loc[gene, 'genomic_interval']
        ilen = TSSinterval[1]-TSSinterval[0]                                      # calculate average TSS interval coverage
        if ilen<=0:
            TSScoverage=0.0
        else:
            TSS = gene_iv.start_d_as_pos
            if TSS.strand == '+':
                st=max(0, TSS.start+TSSinterval[0])
                en=max(0, TSS.start+TSSinterval[1])
            else:
                st=max(0, TSS.start-TSSinterval[1])
                en=max(0, TSS.start-TSSinterval[0])
            TSSinterval = HTSeq.GenomicInterval(TSS.chrom, st, en, ".")
            wincvg = np.fromiter(data[TSSinterval], dtype=float)
            TSScoverage = np.mean(wincvg)
        wincvg = np.fromiter(data[gene_iv], dtype=float, count=gene_iv.length)     # calculate average gene coverage
        GeneCvg = np.mean(wincvg)
        if normFactor is not None:
            TSScoverage = TSScoverage * normFactor                                 # normalize average coverage values
            GeneCvg = GeneCvg * normFactor
        return TSScoverage, GeneCvg

    def GeneCoverageAnalysis(self, experimentList, TSSinterval=(-5000,5000), filename=None):
        """
        Returns a pandas DataFrame with the average coverage of a TSS interval and the transcription of all genes in the current Annoatation for the experiments specified
        in experimentList. The TSS interval is specified by a tuple containing the 5' (left) and 3' (right) length from TSS. If "normFactor" is defined in
        ExperimentStatistics it is used to multiply the coverage count values for normalization. If a filename is specified the resulting DataFrame is written to a csv
        file for later retrieval by calling ReadAnalysiFromCSVfile.
        """
        cols = ["chrom","start","end","length"]
        for exp_name in experimentList:
            cols.append(exp_name+"AvgTSSCvg")
            cols.append(exp_name+"AvgGeneCvg")
        infoTable = pd.DataFrame(columns=cols)
        infoTable.index.name="geneID"
        for gene in self.Annotation.index:
            gene_iv = self.Annotation.loc[gene, 'genomic_interval']
            cols = [ str(gene_iv.chrom), str(gene_iv.start), str(gene_iv.end), str(gene_iv.length) ]
            for exp_name in experimentList:
                TSScoverage, GeneCvg = self.getGeneInfo(exp_name, gene, TSSinterval)
                cols.append(str(TSScoverage))
                cols.append(str(GeneCvg))
            infoTable.loc[gene]=cols
        if filename is not None:
            infoTable.to_csv(filename)
        return infoTable

    def ReadAnalysisFromCSVfile(self, filename):
        """
        Reads a csv file from a stored GeneCoverageAnalysis back into a pandas DataFrame. Use this in conjunction with the filename option of GeneCoverageAnalysis.
        """
        return pd.read_csv(filename)

    def getIntervalCoverageArray(self, experiment_name, genomicInterval):
        """
        Returns the average read coverage of the genomic interval for the experiment identified by experiment_name. If a "normFactor" is defined in ExperimentStatistics
        it is used to multiply the average coverage count values for normalization. If the genomic interval has strand="-" specified, the data is inverted to represent values in
        the correct orientation relative to start and end.
        """
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        if genomicInterval.length == 0:
            return None
        else:
            wincvg = np.fromiter(data[genomicInterval], dtype=float)
        if genomicInterval.strand == '-':
            Cvg=wincvg[::-1]
        else:
            Cvg=wincvg
        if normFactor is not None:
            Cvg = Cvg * normFactor                                 # normalize coverage values
        return Cvg

    def FeatureCoverageAnalysis(self, FeatureAnnotationDF, experimentList, filename=None):
        """
        Returns a pandas DataFrame with the average and maximal coverage of features in the FeatureAnnotationDF for the experiments specified in experimentList. If "normFactor" is defined
        in ExperimentStatistics it is used to multiply the coverage count values for normalization. If a filename is specified the resulting DataFrame is written to a csv file for later
        retrieval by calling ReadAnalysiFromCSVfile. Select the features for coverage analysis from an Annotation DataFrame using pandas column selection df1=df[["chr","start","end"]] or
        slicing df1=df.loc[[feature row list], [["chr","start","end"]] for selecting a subset of features within the Annotation for analysis. The genomic interval of the feature is constructed
        from the annotation: chr, start, end. For sorting before or after the coverage analysis refer to pandas documentation df.sort_by_value.
        """
        cols = ["name", "chr","start","end", "length"]
        for exp_name in experimentList:
            cols.append(exp_name+"MaxCvg")
            cols.append(exp_name+"AvgCvg")
        infoTable = pd.DataFrame(columns=cols)
        feature_index=1
        for feature in FeatureAnnotationDF.itertuples():
            if "strand" in feature._asdict():
                strand=feature.strand   # getattr(feature, "strand")
            else:
                strand="."
            giv = HTSeq.GenomicInterval(feature.chr, feature.start, feature.end, strand)
            if "name" in feature._asdict():
                feature_name=feature.name
            else:
                if FeatureAnnotationDF.index.name=="name":
                    feature_name=feature.Index
                else:
                    feature_name="feature"+str(feature_index)
            cols = [ feature_name, str(giv.chrom), str(giv.start), str(giv.end), str(giv.length) ]
            valid_data = True
            for exp_name in experimentList:                    
                Cvg = getIntervalCoverageArray(exp_name, giv)
                if Cvg is not None:
                    cols.append(np.max(Cvg))
                    cols.append(np.mean(Cvg))
                else:
                    valid_data = False
                    break
            if valid_data:
                infoTable=infoTable.append(pd.Series(cols, index=infoTable.columns), ignore_index=True)
                feature_index += 1
        infoTable=infoTable.set_index("name")
        if filename is not None:
            infoTable.to_csv(filename)
        return infoTable
    
    def generateFeatureProfile(self, experiment_name, FeatureAnnotationDF, interval=(-2000,2000)):
        """
        Returns a numpy ndarray with the sum of the read coverage for the features in FeatureAnnotationDF ("chr","start","end"). Elements contain read coverage per base for the interval
        centered on the midpoint of the feature. The 5' and 3' length of the interval is a tuple, default is (-2000,2000). If FeatureAnnotationDF contains strand information it will be
        used and if strand="-" is specified, the data for the corresponding feature is inverted to represent values in the correct orientation relative to start and end. If strand is not
        provided non-stranded features are assumed and no corrective calculation for the orientation is performed.
        If "normFactor" is defined in ExperimentStatistics it is used to multiply the coverage count values for normalization.
        """
        ilen = interval[1]-interval[0]
        if ilen<=0:
            return None
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        FeatureProfile = np.zeros(ilen, dtype=float)
        for feature in FeatureAnnotationDF.itertuples():
            if "strand" in feature._asdict():
                strand=feature.strand   # getattr(feature, "strand")
            else:
                strand="."
            midpoint=int((feature.start + feature.end)/2)
            if strand=="-":
                st=midpoint-interval[1]
                en=midpoint-interval[0]
                if st<0 or en<0:
                    continue
                giv = HTSeq.GenomicInterval(feature.chr, st, en, ".")
                wincvg = np.fromiter(data[giv], dtype=float)
                FeatureProfile += wincvg[::-1]
            else:
                st=midpoint+interval[0]
                en=midpoint+interval[1]
                if st<0 or en<0:
                    continue
                giv = HTSeq.GenomicInterval(feature.chr, st, en, ".")
                wincvg = np.fromiter(data[giv], dtype=float)
                FeatureProfile += wincvg
        if normFactor is not None:
            FeatureProfile = FeatureProfile * normFactor
        return FeatureProfile

    def generateFeatureCoverageProfile(self, experiment_name, FeatureAnnotationDF, interval=(-2000,2000)):
        """
        Returns a pandas dataframe containing the read coverage for the features in FeatureAnnotationDF, and
        a numpy ndarray with the sum of the coverage. Columns contain read coverage per base for interval that is centered at the midpoint of the feature.
        The interval is pecified by a tuple of the 5' [0] and 3' [1] length from midpoint. IF "normFactor" is defined in ExperimentStatistics it is used to
        multiply the coverage count values for normalization. If FeatureAnnotationDF contains strand information it will be used and if strand="-" is specified,
        the data for the corresponding feature is inverted to represent values in the correct orientation relative to start and end. If strand is not
        provided non-stranded features are assumed and no corrective calculation for the orientation is performed.
        """
        ilen = interval[1]-interval[0]
        if ilen<=0:
            return None
        data = self.getExperiment(experiment_name)
        normFactor = self.getExperimentNormFactor(experiment_name)
        FeatureCoverage = pd.DataFrame(columns= [x for x in range(interval[0],interval[1])])
        FeatureCoverage.index.name = "name"
        FeatureProfile = np.zeros(ilen, dtype=float)
        feature_index = 1
        for feature in FeatureAnnotationDF.itertuples():
            if "name" in feature._asdict():
                feature_name=feature.name
            else:
                if FeatureAnnotationDF.index.name=="name":
                    feature_name=feature.Index
                else:
                    feature_name="feature"+str(feature_index)
            if "strand" in feature._asdict():
                strand=feature.strand   # getattr(feature, "strand")
            else:
                strand="."
            midpoint=int((feature.start + feature.end)/2)
            if strand=="-":
                st=midpoint-interval[1]
                en=midpoint-interval[0]
                if st<0 or en<0:
                    continue
                giv = HTSeq.GenomicInterval(feature.chr, st, en, ".")
                wincvg = np.fromiter(data[giv], dtype=float)
                if normFactor is not None:
                    wincvg*=normFactor
                FeatureCoverage.loc[feature_name] = wincvg[::-1]
                FeatureProfile += wincvg[::-1]
                feature_index  += 1
            else:
                st=midpoint+interval[0]
                en=midpoint+interval[1]
                if st<0 or en<0:
                    continue
                giv = HTSeq.GenomicInterval(feature.chr, st, en, ".")
                wincvg = np.fromiter(data[giv], dtype=float)
                if normFactor is not None:
                    wincvg*=normFactor
                FeatureCoverage.loc[feature_name] = wincvg
                FeatureProfile += wincvg
        return FeatureCoverage, FeatureProfile

    def plot_average_coverage(self, data, unit=500, centered=True, title="Average coverage", ax=None, xticklabels=True, xtickrotation=90):
        """
        Generate a plot for an ndarray data of the average coverage over bases in interval. Adjust for x tick labels
        for one per unit of number of x values from 0 (centered = False) or centered around 0 (centered = true).
        returns the axis object of the plot.
        """
        if ax is None:
            ax=plt.gca()
        ax.margins(x=0.0)
        x_max=len(data)
        if centered:
            xlabs = [str(i) for i in range(-x_max/2,x_max/2+1,unit)]
        else:
            xlabs = [str(i) for i in range(0,x_max+1,unit)]
        ax.plot(data)
        if xticklabels:
            ax.set_xticks([i for i in range(0,x_max+1,unit)]) #[0,499,999,1499,1999])
            ax.set_xticklabels(xlabs) # ["-1000","-500","0","500","1000",""])
            for t in ax.get_xticklabels():  # vertical lables for matplotlib 1.5.1
                t.set_rotation(xtickrotation)
        else:
            ax.set_xticks(False)
        ax.set_title(title)
        return ax

    def plot_heatmap_interval(self, data, unit=500, centered=True, title="Genes", ax=None, cbar_ax=None, xticklabels=True, xtickrotation=90, cbar_horizontal=False):
        """ Generate a seaborn heatmap for an ndarray data. Adjust for x tick labels
        for one per unit of number of x values from 0 (centered = False) or centered around 0 (centered = true).
        returns the axis object of the plot.
        """
        if type(data) is np.ndarray:
            x_max=data.shape[1]
        else:
            if type(data) is pd.core.frame.DataFrame:
                x_max=data.columns
            else:
                print("plot_heatmap_interval(data ... a numpy ndarray or pandas DataFrame is needed. Type provided is: ", type(data))
                return None
        if cbar_horizontal:
            cbar_kws={"orientation":"horizontal"}
        else:
            cbar_kws={"orientation":"vertical"}
        if centered:
            xlabs = [str(i) for i in range(-x_max/2,x_max/2+1,unit)]
        else:
            xlabs = [str(i) for i in range(0,x_max+1,unit)]
        if cbar_ax is not None:
            ax=sns.heatmap(data, xticklabels=xticklabels, yticklabels=False, ax=ax, cbar_ax=cbar_ax, cbar_kws=cbar_kws)
            ax.axhline(y = 0, color='k',linewidth = 5)
        else:
            ax=sns.heatmap(data, xticklabels=xticklabels, yticklabels=False, ax=ax)
            ax.axhline(y = 0, color='k',linewidth = 5)
        # ax=sns.heatmap(data, yticklabels=False, ax=ax, cbar_ax=cbar_ax)
        if xticklabels:
            ax.set_xticks([i for i in range(0,x_max+1,unit)]) #[0,499,999,1499,1999])
            ax.set_xticklabels(xlabs) # ["-1000","-500","0","500","1000",""])
            for t in ax.get_xticklabels():  # vertical lables for matplotlib 1.5.1
                 t.set_rotation(xtickrotation)
        ax.set_title(title)
        return ax

    def plot_coverage_heatmap_interval(self, data, unit=500, centered=True, title="Genes", heatmap_title=""):
        fig=plt.figure()
        gs=gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[1,4], width_ratios=[20,1])
        ax=fig.add_subplot(gs[0,0])
        if type(data) is np.ndarray:
            self.plot_average_coverage(data, title=title, centered=True, unit= 500, ax=ax)
        else:
            if type(data) is pd.core.frame.DataFrame:
                self.plot_average_coverage(data.values, title=title, centered=True, unit= 500, ax=ax)
            else:
                print("plot_heatmap_interval(data ... a numpy ndarray or pandas DataFrame is needed. Type provided is: ", type(data))
                return None
        ax=fig.add_subplot(gs[1,0])
        cbar_ax=fig.add_subplot(gs[1,1])
        self.plot_heatmap_interval(data, title=heatmap_title, centered=True, unit= 500, xticklabels=False, ax=ax, cbar_horizontal=False, cbar_ax=cbar_ax)
        fig.tight_layout()
        return fig

    def writeExperimentBEDfile(self, experiment_name, filename=None, folder=None):
        """
        Writes a BED graph file containing the coverage. If filename is not specified the experiment name will be used
        and the extension ".wig" will be added. If folder is not None the path will be prepended to the filename.
        The concatenated string from folder and filename must be a valid path and writable to the user.
        Note that ExperimentStatistics is not saved.
        """
        data = self.getExperiment(experiment_name)
        if filename is None:
            fn = experiment_name + ".wig"
        else:
            fn = filename
        if folder is not None:
            fn = folder + fn
        data.write_bedgraph_file(fn, strand=".")

    def readExperimentBEDfile(self, experiment_name, filename=None, folder=None):
        """
        Reads the coverage of an experiment from a BED graph file into a GenomicArray and adds it to the XchipAnalysis
        object Experiment dictionary using the key experiment_name. Old data associated with experiment_name will be replaced.
        If filename is not specified the experiment name will be used and the extension ".wig" will be added.
        If folder is not None the path will be prepended to the filename. The concatenated string from folder and filename must be
        a valid path and writable to the user. Note that ExperimentStatistics will not be restored.
        """
        if filename is None:
            fn = experiment_name + ".wig"
        else:
            fn = filename
        if folder is not None:
            fn = folder + fn
        data = HTSeq.GenomicArray("auto", stranded=False, typecode='i')
        BEDreader = HTSeq.BED_Reader(fn)                             # read the wiggle file record for record and store in GenomicArray
        for rec in BEDreader:
            data[rec.iv]=int(float(rec.name))
        self.Experiments[experiment_name]=data                       # add the GenomicArray to the dictionary of the XchipAnalysis object

    def saveAnalysis(self, folder, save_annotation=True):
        """
        Saves the data associated with the XchipAnalysis object to a specified folder that must already exist. The coverage of all experiments
        are written to BED graph files with the filename constructed from the experiment name and the extension ".wig". ExperimentStatistics
        and printVerboseLog are written to a pickl file ExperimentStatistics. Old files will be replaced.
        """
        state = dict()
        state["ExperimentStatistics"]=self.ExperimentStatistics
        state["printVerboseLog"]=self.printVerboseLog
        outFile = open(folder + "ExperimentStatistics", "wb")        # save ExperimentStatistics and printVerboseLog as a pickel file
        pickle.dump(state, outFile)
        outFile.close()
        if save_annotation:
            outFile = open(folder + "Annotation", "wb")              # if save_annotation==True (default) then save Annotation pickel file
            pickle.dump(self.Annotation, outFile)
            outFile.close()
        exps = self.getExperimentList()                              # save the coverage for every experiment into a BED graph or wiggle file
        for experiment_name in exps:
            self.writeExperimentBEDfile(experiment_name, folder=folder)
        if self.printVerboseLog:
            print(len(exps),"experiments written to BED graph files in", folder,":")
            print(exps)

    def loadAnalysis(self, folder, load_annotation=True):
        """
        Loads data from a specified folder into an XchipAnalysis object. ExperimentStatistics and printVerboseLog are read from a pickl file
        named ExperimentStatistics.The coverage of all experiments is then red from BED graph files with the filename constructed from the
        experiment name and the extension ".wig". If load_annotation==True the Annotation will be read from a pickl file.
        """
        if load_annotation:
            inFile = open(folder + "Annotation", "rb")               # if load_annotation==True (default) then read Annotation pickel file
            self.Annotation = pickle.load(inFile)
            inFile.close()
        inFile = open(folder + "ExperimentStatistics", "rb")         # read ExperimentStatistics and printVerboseLog from a pickel file
        state = pickle.load(inFile)
        inFile.close()
        self.printVerboseLog=state["printVerboseLog"]
        self.ExperimentStatistics=state["ExperimentStatistics"]
        exps = self.ExperimentStatistics.keys()                      # read the coverage for every experiment from a BED graph or wiggle file
        for experiment_name in exps:
            self.readExperimentBEDfile(experiment_name, filename=None, folder=folder)
        if self.printVerboseLog:
            print(len(exps),"experiments loaded into the analysis from", folder,":")
            print(exps)


# Analysis entry point with parameters below
CHIPSEQ_ANALYSIS = True
LOAD_CURRENT_ANALYSIS = False
ANNOTATION_FILE_PATH = "/Users/.../Annotations/CSV_feature_files_mm10/" 
DATA_FILE_PATH = "/Users/.../Analysis/" 
BAM_FILE_PATH= "/Users/.../Bamfiles/"
HASAPPY_ANNOTATION_PICKLE_FILE_PATH="/Users/.../Annotations/generate_annotation/Annotation"

if CHIPSEQ_ANALYSIS:
    if LOAD_CURRENT_ANALYSIS:
        analysis = XchipAnalysis()
        analysis.loadAnalysis(DATA_FILE_PATH)
        x_genes = analysis.get_x_linked_genes()
        a_genes = analysis.get_autosomal_genes()
        
    else:
        analysis = XchipAnalysis(Annotation_filename=HASAPPY_ANNOTATION_PICKLE_FILE_PATH)
        x_genes = analysis.get_x_linked_genes()
        a_genes = analysis.get_autosomal_genes()
        exps=["N8", "N9", "N10"]

        analysis.addExperiment("Name_nontreated", BAM_FILE_PATH+'name_nontreated_filtered.bam', aQualValue=10, frgmtLength=200, maxReadsPerPos=5 )
        analysis.addExperiment("Name_treated", BAM_FILE_PATH+'name_treated_filtered.bam', aQualValue=10, frgmtLength=200, maxReadsPerPos=5 )
        analysis.addExperiment("Name2_nontreated", BAM_FILE_PATH+'name2_nontreated_filtered.bam', aQualValue=10, frgmtLength=200, maxReadsPerPos=5 )
        analysis.addExperiment("Name2_treated", BAM_FILE_PATH+'name2_treated_filtered.bam', aQualValue=10, frgmtLength=200, maxReadsPerPos=5 )
        


# normalize to autosomal coverage
reads_max = max([ analysis.ExperimentStatistics[en]["coverageA"] for en in analysis.getExperimentList() ])
for en in analysis.getExperimentList():
    nf = float(reads_max) / analysis.ExperimentStatistics[en]["coverageA"]
    analysis.normalizeExperiment(en, nf)s
    print(en,"normalization factor auto=",nf)
    
#saving analysis
analysis.saveAnalysis(DATA_FILE_PATH) # ------------------------ save analysis here---------------------------------->>>
experiment_names=analysis.getExperimentList()
print("Experiments:", experiment_names)


# generated ordered experiment list
Name1 = ["Name_nontreated"]
Name1_treated = ["Name_treated"]
Name2 = ["Name2_nontreated"]
Name2_treated = ["Name2_treated"]
Wildtype_samples = {"ctrl": Name1, "dox": Name1_treated}
Mutant_samples = {"ctrl": Name2 , "dox": Name2_treated}
samples = {"Wildtype": Wildtype_samples, "Mutant": Mutant_samples}

expLst=[]

for genotype in ("Wildtype", "Mutant"):
    for condition in ("ctrl", "dox"):
        for replicate in samples[genotype][condition]:
            expLst.append(replicate)
           

# helper function for rolling mean / moving average if pandas v < 0.17 - use: np.mean(rolling_window(x, 3), axis=-1)
def rolling_window(a, window):
    """
    A helper function for calculating the row-wise rolling mean or moving average of numpy arrays. It works by extending the dimensions and
    using stride tricks. Useful when no recent version of pandas is available as df.rolling is not implemented in pandas v < 0.17

    Use: np.mean(rolling_window(x, 3), axis=-1)
    """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def row_sort_by_mean(a, descending=False):
    """
    Sorts rows of a 2-dimensional numpy array according to their mean ascending (default), for descending set descending=True
    """
    sa=np.argsort(a.mean(axis=1))
    if descending:
        return a[sa[::-1]]
    else:
        return a[sa]
    
def scale_y(ymin, ymax):
    """
    Takes a minimum and maximum y axis value and returns a value pair that is aligned to the lower (min) and higer (max) multiple of 1000. This is useful for generating labels with
    larger values. No adjustment is made for values below 1000. Use this with the format_y_ticklabel function for succintly annotating the y axis.
    """
    scaled_min = int(ymin)
    if ymin >= 1000:
        scaled_min = int(ymin/1000)*1000
    scaled_max = int(ymax+1)
    if ymax >= 1000:
        scaled_max = int((ymax+999)/1000)*1000
    return scaled_min, scaled_max

def format_y_ticklabel(val):
    """ For values greater or equal to 1000 the integer value divided by 1000 followed by a "k" for kilo is returned as a string. No adjustment is made for values below 1000, and
    the value is returned as a string. Use this or succintly annotatin the y axis ticks.
    """
    if val >= 1000:
        return str(int(val)/1000)+"k"
    else:
        return str(val)

def analyseFeatureCoverage(analysis, FeatureAnnotationFile, FeatureName, interval=(-5050,5050), chromosome=None):
    """
    Calculates coverage for features in the FeatureAnnotationFile and stores it in a file, which can be plotted using plotFeatureCoverage. Also cummulative coverage is calculated by
    summation of replicates. FeaturAnnotationFile is a csv formated file that must contain "chr", "start", and "end" for each feature. If the feature annotation contains a strand
    specification and strand="-" the data is inverted to provide the correct orientation relative to start and end. The analysis can be restricted by specifying a chromosome as either
    "chrX" or "autosome".
    """
    Feature_Annotation=pd.read_csv(FeatureAnnotationFile)
    if chromosome == "chrX":
        Feature_Annotation=Feature_Annotation.loc[Feature_Annotation['chr'] == 'chrX']
    else:
        if chromosome.lower()[:7] == "autosom":
            Feature_Annotation=Feature_Annotation.loc[Feature_Annotation['chr'].apply(lambda x: x.upper()[:4] != 'CHRX' and x.upper()[:3] == 'CHR')]
    for en in analysis.getExperimentList():
        FeatureCoverage, FeatureProfile = analysis.generateFeatureCoverageProfile(en, Feature_Annotation, interval=interval)
        FeatureCvg=FeatureCoverage.values
        np.save(DATA_FILE_PATH+FeatureName+"cvg_"+en+".npy", FeatureCvg)
        np.save(DATA_FILE_PATH+FeatureName+"prof_"+en+".npy", FeatureProfile)
    # calculate the cummulative coverage by adding of replicates
    for genotype in ("Wildtype", "Mutant"):
        for condition in ("ctrl", "dox"):
            en=samples[genotype][condition][0]
            FeatureCvg=np.load(DATA_FILE_PATH+FeatureName+"cvg_"+en+".npy")
            for replicate in samples[genotype][condition][1:]:
                cvg=np.load(DATA_FILE_PATH+FeatureName+"cvg_"+replicate+".npy")
                FeatureCvg=np.add(FeatureCvg,cvg)
            np.save(DATA_FILE_PATH+FeatureName+"cvg_"+genotype+"_"+condition+".npy", FeatureCvg)

def plotFeatureCoverage(filename, fig_filename=None, fig_name="Heatmap"):
    """
    Generates a heatmap from coverage data file. Use analyseFeatureCoverage to produce coverage files. Eithe coverage of individual experiments or cummulative coverage summed from
    replicates can be plotted. If fig_filename is specified the figure will be stored in this file.
    Returns a reference to the resulting figure. Note the figure should be closed when no longer used to conserve memory using plt. close("Heatmap").
    """
    for en in analysis.getExperimentList():
        FeatureCvg=np.load(filename)
        FeatureMa=np.mean(rolling_window(FeatureCvg, 100), -1)
        Feature=FeatureMa[:,~np.any(np.isnan(FeatureMa), axis=0)]
        zsFeatureCvg=scipy.stats.zscore(Feature, axis=1)
        fig=plt.figure(name=fig_name)
        sns.heatmap(zsFeatureCvg, cmap=analysis.RBGcmap, xticklabels=False, yticklabels=False) # ax=ax[idx], cbar=True cbar_ax=ax[numChrom]
        if fig_filename is not None:
            fig.savefig(fig_filename)
        return fig

def analyseFeatureProfile(analysis, FeatureAnnotationFile, FeatureName, interval=(-2100,2000), chromosome=None):
    """
    Generates a figure for profiles of coverage of features in FeatureAnnotationFile. FeaturAnnotationFile is a csv formated file that must contain "chr", "start", and "end" for each
    feature. If the feature annotation contains a strand specification and strand="-" the data is inverted to provide the correct orientation relative to start and end. The analysis can
    be restricted by specifying a chromosome as either "chrX" or "autosome". Coverage profiles for Ctrl and Dox or all replicates are plotted in the same figure. Two figures are stored in
    files, whose names are construced from the FeatureName, "_profile_", the genotype, and ".pdf". Useful to assess if features are consistently covered among replicates.
    """
    colorList = ["steelblue", "seagreen", "darkcyan", "darkred", "firebrick", "darkorange"]
    Feature_Prof=dict()
    Feature_Annotation=pd.read_csv(FeatureAnnotationFile)
    if chromosome == "chrX":
        Feature_Annotation=Feature_Annotation.loc[Feature_Annotation['chr'].apply(lambda x: x.upper() == 'CHRX')]
    else:
        if chromosome.lower()[:7] == "autosom":
            Feature_Annotation=Feature_Annotation.loc[Feature_Annotation['chr'].apply(lambda x: x.upper()[:4] != 'CHRX' and x.upper()[:3] == 'CHR')]
    for en in analysis.getExperimentList():
        fp = analysis.generateFeatureProfile(en, Feature_Annotation, (-5000,5000))
        Feature_Prof[en]=fp
    for genotype in ("Wildtype", "Mutant"):             # plot profiles of all replicates and conditions for one genotype into a single figure
        fig=plt.figure("FeatureProfileFigure")
        ax=fig.add_subplot(111)
        col=0
        for condition in ("ctrl", "dox"):
            for replicate in samples[genotype][condition]:
                ax.plot(Feature_Prof[replicate],color=colorList[col])
                col=col+1
        filename = FeatureName + "_profile_"+ genotype
        if chromosome is not None:
            filename = filename + "_" + chromosome
        fig.savefig(filename + ".pdf")
        plt.close("FeatureProfileFigure")

def Heizmappe(analysis, FeatureAnnotationFile, FeatureName, interval=(-5000,5000), window=100, chromosome=None, z_transform=False, sort_rows=True, descending=True, log_color_range=True, xtickrotation=90, cbar=True):
    """
    Generates heatmaps and profiles of cummulative coverage from summing all replicates for genotypes, and ctrl and Dox condition. The figure is formated and stored in a file, whose name is
    constructed from DATA_FILE_PATH, FeatureName, and "_heatmap.png". The interval for calculating the feature read coverage is specified with a tuple of the 5' and 3' extend relative to the
    midpoint of the feature genomic interval from FeatureAnnotationFile (default is an interval -5000 to +5000).
    A window can be specified for calculating the rolling mean of coverage. The calculation is performed in such a was that the window average corresponds to the center of the window. For this
    coverage from an interval is calculated that is further extended by window/2 in both directions. FeaturAnnotationFile is a csv formated file that must contain "chr", "start", and "end" for
    each feature. If the feature annotation contains a strand specification and strand="-" the data is inverted to provide the correct orientation relative to start and end. The analysis can
    be restricted by specifying a chromosome as either "chrX" or "autosome".
    For improved visualization in a heatmap, rows can be sorted according to their average coverage, data can be either z-transformed or a logarithmic color range applied.
    """
    REDcmap = LinearSegmentedColormap.from_list("myreds",['#ffffff','#ff8080','#ff4040','#ff0000','#ff0000'])
    draw_cbar=cbar
    if window is not None:
        analysis_interval=(int(interval[0]-window/2), int(interval[1]+window/2))
    else:
        analysis_interval=interval
    # Generate the coverage data and profiles for plotting
    analyseFeatureCoverage(analysis, FeatureAnnotationFile, FeatureName, interval=analysis_interval, chromosome=chromosome)
    FeatureCvg={}
    FeatureProf={}
    for genotype in ("Wildtype", "Mutant"):
        FeatureCvg[genotype]={}
        for condition in ("ctrl", "dox"):
            FeatureCvg[genotype][condition]=np.load(DATA_FILE_PATH+FeatureName+"cvg_"+genotype+"_"+condition+".npy")
    if sort_rows:
        sa=np.argsort(FeatureCvg["Wildtype"]["ctrl"].mean(axis=1))
        if descending:
            sa=sa[::-1]
        for genotype in ("Wildtype", "Mutant"):
            for condition in ("ctrl", "dox"):
                FeatureCvg[genotype][condition]=FeatureCvg[genotype][condition][sa]  # sort all by wildtype ctrl order
    if window is not None:
        for genotype in ("Wildtype", "Mutant"):
            for condition in ("ctrl", "dox"):
                FeatureCvg[genotype][condition]=np.mean(rolling_window(FeatureCvg[genotype][condition], window), axis=-1)  # calculate the moving average over windows
                #FeatureCvg[genotype][condition]=FeatureCvg[genotype][condition][:,window:]
                FeatureCvg[genotype][condition]=FeatureCvg[genotype][condition][:,~np.any(np.isnan(FeatureCvg[genotype][condition]), axis=0)]
    # Generate cummulative profiles from the combined coverage data of the replicates
    prof_y_min = None
    prof_y_max = None
    for genotype in ("Wildtype", "Mutant"):
        FeatureProf[genotype]={}
        for condition in ("ctrl", "dox"):
            FeatureProf[genotype][condition]=np.sum(FeatureCvg[genotype][condition], axis=0)
            if prof_y_min is None:
                prof_y_min = FeatureProf[genotype][condition].min()
                prof_y_max = FeatureProf[genotype][condition].max()
            else:
                prof_y_min = min(prof_y_min, FeatureProf[genotype][condition].min())
                prof_y_max = max(prof_y_max, FeatureProf[genotype][condition].max())
    # z transform the coverage array row-wise if specified. Note this option is not useful with log_color_range
    if z_transform:
        for genotype in ("Wildtype", "Mutant"):
            for condition in ("ctrl", "dox"):
                FeatureCvg[genotype][condition]=scipy.stats.zscore(FeatureCvg[genotype][condition], axis=1)
    fig=plt.figure(num=FeatureName, figsize=(5, 8), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=3, ncols=3, height_ratios=[1,4,4], width_ratios=[20,20,1])
    vmin=None
    vmax=None
    for genotype in ("Wildtype", "Mutant"):
        for condition in ("ctrl", "dox"):
            if log_color_range and not z_transform:
                FeatureCvg[genotype][condition]=FeatureCvg[genotype][condition] + 1   # add one to avoid problems with log 0, but only if not z transformed data
            if vmin is None:
                vmin=FeatureCvg[genotype][condition].min()
                vmax=FeatureCvg[genotype][condition].max()
            else:
                vmin=min(vmin, FeatureCvg[genotype][condition].min())
                vmax=max(vmax, FeatureCvg[genotype][condition].max())                
    # plot the cummulative ctrl and dox profiles for genotypes
    panel_row=0
    panel_col=0
    colorList = ["darkcyan", "firebrick"]
    for genotype in ("Wildtype", "Mutant"):
        with plt.style.context('seaborn-v0_8-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            x=range(len(FeatureProf[genotype]["ctrl"]))
            y1=FeatureProf[genotype]["ctrl"]
            y2=FeatureProf[genotype]["dox"]
            ax.fill_between(x, y1, y2, where=(y1 >= y2), color=colorList[0], alpha=0.4, interpolate=True)
            ax.fill_between(x, y1, y2, where=(y2 > y1), color=colorList[1], alpha=0.4, interpolate=True)
            ax.fill_between(x, y1, where=(y2 > y1), color=colorList[0],alpha=0.4, interpolate=True)
            ax.fill_between(x, y2, where=(y1 > y2), color=colorList[1],alpha=0.4, interpolate=True)
            cl=0
            for condition in ("ctrl", "dox"):
                ax.plot(FeatureProf[genotype][condition],color=colorList[cl], label= condition, alpha=0.7)
                ax.legend(loc='upper right', fontsize='xx-small')
                cl += 1
            max_col=len(FeatureProf[genotype][condition])
            ax.set_xticks([0, max_col/2, max_col-1])
            ax.set_xticklabels(["-"+str(int(max_col/2000))+"kb", "0", "+"+str(int(max_col/2000))+"kb"])
            ymin, ymax=scale_y(prof_y_min, prof_y_max)
            ax.set_yticks([0, ymax])
            ax.set_yticklabels(["0", format_y_ticklabel(ymax)])
            ax.set_title(genotype)
            panel_col += 1
    # Plot the heatmaps for different genotypes and conditions
    panel_col=0
    for genotype in ("Wildtype", "Mutant"):
        panel_row=1
        for condition in ("ctrl", "dox"):
            ax=fig.add_subplot(gs[panel_row,panel_col])
            if draw_cbar:
                cax=fig.add_subplot(gs[panel_row,2])
                if log_color_range:
                    ax=sns.heatmap(FeatureCvg[genotype][condition], norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap='GnBu', xticklabels=False, yticklabels=False, ax=ax, cbar=True, cbar_ax=cax)
                else:
                    ax=sns.heatmap(FeatureCvg[genotype][condition], vmin=vmin, vmax=vmax, cmap='GnBu', xticklabels=False, yticklabels=False, ax=ax, cbar=True, cbar_ax=cax)
                draw_cbar=False
            else:
                if log_color_range:
                    ax=sns.heatmap(FeatureCvg[genotype][condition], norm=colors.LogNorm(vmin=vmin, vmax=vmax), cmap='GnBu', xticklabels=False, yticklabels=False, ax=ax, cbar=False)
                else:
                    ax=sns.heatmap(FeatureCvg[genotype][condition], vmin=vmin, vmax=vmax, cmap='GnBu', xticklabels=False, yticklabels=False, ax=ax, cbar=False)                    
            ax.set_title(genotype + " " + condition)
            x_max=FeatureCvg[genotype][condition].shape[1]
            xlabs = [str(interval[0]),"0",str(interval[1])]
            ax.set_xticks([0, x_max/2, x_max])
            ax.set_xticklabels([str(interval[0]),"0",str(interval[1])])
            for tk in ax.get_xticklabels():  # vertical lables for matplotlib 1.5.1
                 tk.set_rotation(xtickrotation)
            panel_row+=1
        panel_col+=1
    fig.tight_layout()
    fig.savefig(DATA_FILE_PATH+FeatureName+"_heatmap.png")
    plt.close(FeatureName)  # close the figure to free up resources after saving

FeatureAnnotations = [  (ANNOTATION_FILE_PATH + "ES_H3K27me3_mm10.csv", "ES_H3K27me3"),
                        (ANNOTATION_FILE_PATH + "CpG_islands_mm10.csv", "CpG"),
                        (ANNOTATION_FILE_PATH + "ES_Enhancers_mm10.csv", "ES_Enhcr"),
                        (ANNOTATION_FILE_PATH + "Ezh2_d0_moderate_mm10.csv", "ES_Ezh2_mod"),
                        (ANNOTATION_FILE_PATH + "Ezh2_d0_mm10.csv", "ES_Ezh2_strng"),
                        (ANNOTATION_FILE_PATH + "TS_Xchr_dnase1_mm10.csv", "TS_DHS"),
                        (ANNOTATION_FILE_PATH + "TS_Xchr_faire_mm10.csv", "TS_FAIRE"),
                        (ANNOTATION_FILE_PATH + "TS_Xchr_CTCF_mm10.csv", "CTCF"),
                        (ANNOTATION_FILE_PATH + "TS_Xchr_pol2_mm10.csv", "TS_PolII")  ]

# Functions plotting heatmaps and associated profile curves for features and TSS

def plotFeaturesHeatmaps():
    for FeatureAnnotationFile, FeatureName in FeatureAnnotations:
        Heizmappe(analysis, FeatureAnnotationFile, FeatureName, interval=(-5000,5000), window=100, chromosome="chrX", z_transform=False, sort_rows=True, descending=True, log_color_range=True, xtickrotation=90, cbar=True)

def plotAffectedGeneGroups():
    Heizmappe(analysis, ANNOTATION_FILE_PATH + "topGenes_TSS.csv", "topGenes", interval=(-5000,5000), window=100, chromosome="chrX", z_transform=False, sort_rows=True, descending=True, log_color_range=True, xtickrotation=90, cbar=True)
    Heizmappe(analysis, ANNOTATION_FILE_PATH + "flopGenes_TSS.csv", "flopGenes", interval=(-5000,5000), window=100, chromosome="chrX", z_transform=False, sort_rows=True, descending=True, log_color_range=True, xtickrotation=90, cbar=True)

def plotExpressionLevelGroups():
    Heizmappe(analysis, ANNOTATION_FILE_PATH + "ES_genes_high.csv", "highGenes", interval=(-5000,5000), window=100, chromosome="chrX", z_transform=False, sort_rows=True, descending=True, log_color_range=True, xtickrotation=90, cbar=True)
    Heizmappe(analysis, ANNOTATION_FILE_PATH + "ES_genes_medium.csv", "mediumGenes", interval=(-5000,5000), window=100, chromosome="chrX", z_transform=False, sort_rows=True, descending=True, log_color_range=True, xtickrotation=90, cbar=True)
    Heizmappe(analysis, ANNOTATION_FILE_PATH + "ES_genes_low.csv", "lowGenes", interval=(-5000,5000), window=100, chromosome="chrX", z_transform=False, sort_rows=True, descending=True, log_color_range=True, xtickrotation=90, cbar=True)


# gene analysis


def geneXAsplit(analysis, geneList):
    """ Splits gene names that can be either provided as a list, numpy ndarray, or a pandas DataFrame into a list of x-chromosomal genes and a list of autosomal genes.
    If geneList is a DataFrame the index or a column "name" must exist that contains the gene names. If no "name" column or index is found None will be returned.
    For the list of x-linked genes the annotation of the XchipAnalysis object is queried with the gene name and if chrom=='chrX' holds true the gene is considered x-linked,
    else it will be appended to the autosomal gene list. Two lists of x-chromosomal and aututosomal gene names will be returned.
    """
    genes = None
    # handle different geneList data formats
    if (type(geneList) is list) or (type(geneList) is np.ndarray):
        genes = geneList
    else:
        if type(geneList) is pd.DataFrame:
            if geneList.index.name == "name":
                genes=geneList.index
            else:
                if "name" in geneList.columns:
                    genes = geneList["name"]
    if genes is None:
        return None
    # separate x-linked and autosomal genes
    autosomal_genes_list=[]
    x_linked_genes_list=[]
    for g in genes:
        if g in analysis.Annotation.index:
            chrom=analysis.Annotation.loc[g,'genomic_interval'].chrom.upper()
            if chrom == 'CHRX':
                x_linked_genes_list.append(g)
            else:
                if chrom[:3]=='CHR' and chrom[:4]!='CHRX':    # try to avoid assigning genes that might be on more exotic or unmapped x-linked contigs, if these exist
                    autosomal_genes_list.append(g)            
    return x_linked_genes_list, autosomal_genes_list

def geneXprofiles(analysis, GeneGroupName, x_genes, a_genes, analyse_TES=False):
    """
    Generates cumulative profiles for two lists of genes specified x_genes (chrX), and a_genes (autosomal). Both a TSS profile and a scaled profile of the transcription unit will be produced for each genotype,
    and chrX and autosome. Use geneXAsplit() to split a list of genes into chrX and autosomal gene lists. ChrX and autosomal profiles are plotted into separate panels.
    Replicates for Ctrl and Dox conditions are plotted in the same figure panel. The figure is stored in a file, whose name is constructed form DATA_FILE_PATH, GeneGroupName, and "_profile.png"
    """
    XgeneProf=dict()
    AgeneProf=dict()
    X_TSS_Prof=dict()
    A_TSS_Prof=dict()
    A_TSS_min=None
    A_Gene_min=None
    X_TSS_min=None
    X_Gene_min=None
    A_TSS_max=None
    A_Gene_max=None
    X_TSS_max=None
    X_Gene_max=None
    for en in analysis.getExperimentList():
        XgeneProf[en], X_gene_num = analysis.generateGeneProfile(en, x_genes)
        AgeneProf[en], A_gene_num = analysis.generateGeneProfile(en, a_genes)
        X_TSS_Prof[en] = analysis.generateTSSProfile(en, x_genes, (-5000,5000), analyse_TES=analyse_TES)
        A_TSS_Prof[en] = analysis.generateTSSProfile(en, a_genes, (-5000,5000), analyse_TES=analyse_TES)
        if A_TSS_min is None:
            A_TSS_min = A_TSS_Prof[en].min()
            A_Gene_min = AgeneProf[en].min()
            X_TSS_min = X_TSS_Prof[en].min()
            X_Gene_min = XgeneProf[en].min()
            A_TSS_max = A_TSS_Prof[en].max()
            A_Gene_max = AgeneProf[en].max()
            X_TSS_max = X_TSS_Prof[en].max()
            X_Gene_max = XgeneProf[en].max()
        else:
            A_TSS_min = min(A_TSS_min, A_TSS_Prof[en].min())
            A_Gene_min = min(A_Gene_min, AgeneProf[en].min())
            X_TSS_min = min(X_TSS_min, X_TSS_Prof[en].min())
            X_Gene_min = min(X_Gene_min, XgeneProf[en].min())
            A_TSS_max = max(A_TSS_max, A_TSS_Prof[en].max())
            A_Gene_max = max(A_Gene_max, AgeneProf[en].max())
            X_TSS_max = max(X_TSS_max, X_TSS_Prof[en].max())
            X_Gene_max = max(X_Gene_max, XgeneProf[en].max())
    # plot profiles
    A_TSS_min = int(A_TSS_min/1000)*1000
    A_Gene_min = int(A_Gene_min/1000)*1000
    X_TSS_min = int(X_TSS_min/1000)*1000
    X_Gene_min = int(X_Gene_min/1000)*1000
    A_TSS_max = int((A_TSS_max+999)/1000)*1000
    A_Gene_max =int((A_Gene_max+999)/1000)*1000
    X_TSS_max = int((X_TSS_max+999)/1000)*1000
    X_Gene_max =int((X_Gene_max+999)/1000)*1000
    fig=plt.figure(num=GeneGroupName, figsize=(5, 8), dpi=300, facecolor='w', edgecolor='k')
    gs=gridspec.GridSpec(nrows=4, ncols=2, height_ratios=[1,1,1,1], width_ratios=[1,1])
    colorList = ["steelblue", "seagreen", "darkcyan", "darkred", "firebrick", "darkorange"]
    # plot x gene profiles
    panel_row=0
    panel_col=0
    for genotype in ("Wildtype", "Mutant"):
        with plt.style.context('seaborn-v0_8-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_ylim([X_Gene_min, X_Gene_max])
            ax.set_title(genotype+" chrX genes")
            ax.set_xticks([0, 250, 500, 750, 1000])
            ax.set_xticklabels(["TSS", "", "", "", "TES"])
            ymin, ymax =  scale_y(X_Gene_min, X_Gene_max)
            ax.set_yticks([ymin, ymax])
            ax.set_yticklabels([format_y_ticklabel(ymin), format_y_ticklabel(ymax)])
            col=0
            for condition in ("ctrl", "dox"):
                for replicate in samples[genotype][condition]:
                    ax.plot(XgeneProf[replicate],color=colorList[col], label= condition, alpha=0.7)
                    ax.legend(loc='upper right', fontsize='xx-small')
                    col=col+1
                col=3   # set color for dox to red/orange group
        panel_col += 1
    # plot x gene profiles
    panel_row=1
    panel_col=0
    for genotype in ("Wildtype", "Mutant"):
        with plt.style.context('seaborn-v0_8-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_ylim([A_Gene_min, A_Gene_max])
            ax.set_title(genotype+" autosomal genes")
            ax.set_xticks([0, 250, 500, 750, 1000])
            ax.set_xticklabels(["TSS", "", "", "", "TES"])
            ymin, ymax =  scale_y(A_Gene_min, A_Gene_max)
            ax.set_yticks([ymin, ymax])
            ax.set_yticklabels([format_y_ticklabel(ymin), format_y_ticklabel(ymax)])
            col=0
            for condition in ("ctrl", "dox"):
                for replicate in samples[genotype][condition]:
                    ax.plot(AgeneProf[replicate],color=colorList[col], label= condition, alpha=0.7)
                    ax.legend(loc='upper right', fontsize='xx-small')
                    col=col+1
                col=3   # set color for dox to red/orange group
            ax.set_xticks([0, 250, 500, 750, 1000])
            ax.set_xticklabels(["TSS", "", "", "", "TES"])
            ymin, ymax =  scale_y(A_Gene_min, A_Gene_max)
            ax.set_yticks([ymin, ymax])
            ax.set_yticklabels([format_y_ticklabel(ymin), format_y_ticklabel(ymax)])
        panel_col += 1
    # plot x TSS profiles
    panel_row=2
    panel_col=0
    for genotype in ("Wildtype", "Mutant"):
        with plt.style.context('seaborn-v0_8-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_ylim([X_TSS_min, X_TSS_max])
            if analyse_TES:
                ax.set_title(genotype+" chrX TES")
            else:
                ax.set_title(genotype+" chrX TSS")
            max_col=len(X_TSS_Prof[ samples[genotype]["ctrl"][0] ])
            ax.set_xticks([0, max_col/2, max_col-1])
            if analyse_TES:
                ax.set_xticklabels(["-"+str(int(max_col/2000))+"kb", "TES", "+"+str(int(max_col/2000))+"kb"])
            else:
                ax.set_xticklabels(["-"+str(int(max_col/2000))+"kb", "TSS", "+"+str(int(max_col/2000))+"kb"])
            ymin, ymax =  scale_y(X_TSS_min, X_TSS_max)
            ax.set_yticks([ymin, ymax])
            ax.set_yticklabels([format_y_ticklabel(ymin), format_y_ticklabel(ymax)])
            col=0
            for condition in ("ctrl", "dox"):
                for replicate in samples[genotype][condition]:
                    ax.plot(X_TSS_Prof[replicate],color=colorList[col], label= condition, alpha=0.7)
                    ax.legend(loc='upper right', fontsize='xx-small')
                    col=col+1
                col=3   # set color for dox to red/orange group
        panel_col += 1
    # plot autosomal TSS profiles
    panel_row=3
    panel_col=0
    for genotype in ("Wildtype", "Mutant"):
        with plt.style.context('seaborn-v0_8-white'):
            ax=fig.add_subplot(gs[panel_row, panel_col])
            ax.set_ylim([A_TSS_min, A_TSS_max])
            if analyse_TES:
                ax.set_title(genotype+" autosomal TES")
            else:
                ax.set_title(genotype+" autosomal TSS")
            max_col=len(A_TSS_Prof[ samples[genotype]["ctrl"][0] ])
            ax.set_xticks([0, max_col/2, max_col-1])
            if analyse_TES:
                ax.set_xticklabels(["-"+str(int(max_col/2000))+"kb", "TES", "+"+str(int(max_col/2000))+"kb"])
            else:
                ax.set_xticklabels(["-"+str(int(max_col/2000))+"kb", "TSS", "+"+str(int(max_col/2000))+"kb"])
            ymin, ymax =  scale_y(A_TSS_min, A_TSS_max)
            ax.set_yticks([ymin, ymax])
            ax.set_yticklabels([format_y_ticklabel(ymin), format_y_ticklabel(ymax)])
            col=0
            for condition in ("ctrl", "dox"):
                for replicate in samples[genotype][condition]:
                    ax.plot(A_TSS_Prof[replicate],color=colorList[col], label= condition, alpha=0.7)
                    ax.legend(loc='upper right', fontsize='xx-small')
                    col=col+1
                col=3   # set color for dox to red/orange group
        panel_col += 1
    fig.tight_layout()
    fig.savefig(DATA_FILE_PATH+GeneGroupName+"_profile.png")
    plt.close(GeneGroupName)  # close the figure to free up resources after saving

geneGroupLists = [ (ANNOTATION_FILE_PATH + "ES_K27_Genes.csv", "ES_K27_Genes"),
                   (ANNOTATION_FILE_PATH + "ES_K4K27_Genes.csv", "ES_K4K27_Genes"),
                   (ANNOTATION_FILE_PATH + "ES_K4_Genes.csv", "ES_K4_Genes") ]

# Functions plotting profiles over gene TSS and transcription unit

def analyseGeneGroups():
    for gg, gg_name in geneGroupLists:
        genes = pd.read_csv(gg)
        x_genes, a_genes = geneXAsplit(analysis, genes)
        geneXprofiles(analysis, gg_name, x_genes, a_genes)

def topNflopGeneProfiles():
    topXGenes=pd.read_csv(ANNOTATION_FILE_PATH + "topGenes_TSS.csv")
    top_x_genes, ag = geneXAsplit(analysis, topXGenes)
    geneXprofiles(analysis, "topGenes", top_x_genes, a_genes)
    flopXGenes=pd.read_csv(ANNOTATION_FILE_PATH + "flopGenes_TSS.csv")
    flop_x_genes, ag = geneXAsplit(analysis, flopXGenes)
    geneXprofiles(analysis, "flopGenes", flop_x_genes, a_genes)

def ExpressionGroupProfiles():
    GeneGroupNames=["ES_genes_high", "ES_genes_medium", "ES_genes_low"]
    for GroupName in GeneGroupNames:
        GeneAnnotationDF=pd.read_csv(ANNOTATION_FILE_PATH + GroupName + ".csv")
        Group_x_genes, Group_a_genes = geneXAsplit(analysis, GeneAnnotationDF)
        geneXprofiles(analysis, GroupName, Group_x_genes, Group_a_genes)
    plotExpressionLevelGroups()


