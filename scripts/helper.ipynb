#!/usr/bin/env python

# This script contains helper functions.

# import packages
import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np

# class
class helper_funcs(object):
    '''
    Helper functions.
    '''

    def __init__(self):
        '''
        Constructor.
        '''
        return

    def sample_data(self):
        '''
        Samples with both DNA and RNA information and TP53 status.
        '''
        data_file = "../data/DNA_RNA_sample.csv"
        df = pd.read_csv(data_file, sep=",")
        return df

    def p53_fitness(self):
        '''
        p53 fitness.
        '''
        data_file = "../data/PDAC_autopsy_p53_fitness.tsv"
        df = pd.read_csv(data_file, sep="\t")
        df["Fitness"] = df["Fitness"] / df["Fitness"].max()
        df.rename(columns={"Mutation" : "HGVSp_Short"}, inplace=True)
        df["HGVSp_Short"] = df["HGVSp_Short"].replace("^p.", "", regex=True)
        return df

    def TP53_mut_df(self):
        '''
        TP53 mutation status.
        '''
        # read in file
        data_file = "../data/TP53_mut_ccf_L1_no_0.csv" # TP53_mut_ccf_L1_no_0.csv
        df = pd.read_csv(data_file, sep=",")
        ## rename columns
        df = df[["Sample_ID", "Gene", "HGVSc", "HGVSp_Short","Insertion_z", 
                 "Prob MT p53", "VAF", "ccf_Mcopies","Insertion."]]
        return df

    def fpkm_SINE_dsRNAforce(self):
        '''
        Data on SINE and dsRNA force.
        '''
        data_file = "../data/fpkm_SINE_filtered.csv"
        df = pd.read_csv(data_file, sep=",")
        #df.columns = df.columns.str.replace("T0", "T")
        return df

    def IFN_data(self):
        '''
        IFN data.
        '''
        data_file = "../data/immune_P53_adar1.csv"
        df = pd.read_csv(data_file, sep=",")
        return df

    def ADAR1_editing(self):
        '''
        Locus-specific ADAR1 editing.
        '''
        data_file = "../data/sine_locus_editing_filtered.csv" #sine_locus_editing_filtered.csv
        df = pd.read_csv(data_file, sep=",")
        return df

    def organized_model_data(self):
        '''
        Organized model data.
        '''
        # read in data
        RIPseq_data = pd.read_csv("../data/RIPseq-TE_fraction_ORF1p.tsv", sep="\t")
        ADAR1_data = self.ADAR1_editing()
        fpkmSINE_data = self.fpkm_SINE_dsRNAforce()
        IFN_data = self.IFN_data()
        mutation_data = self.TP53_mut_df()
        sample_data = self.sample_data()

        # prep locus-level data
        locus_data = pd.merge(RIPseq_data, ADAR1_data, left_on="TE_ID", right_on="rep.id").drop(["rep.id"], axis=1)
        locus_data = pd.merge(locus_data, fpkmSINE_data, left_on="TE_ID", right_on="rep.id", suffixes=["_ADAR1_Edit", "_Expression"])

        # organize data
        df = pd.merge(sample_data, mutation_data, on=["Sample_ID"], how='right')
        #df["HGVSc"] = df["HGVSc"].fillna("WT")
        df.loc[df["TP53_call"] == "WT", "Gene"] = "TP53"
        df.loc[df["TP53_call"] == "WT", "HGVSc"] = "WT"
        df.loc[df["HGVSc"] == "WT", "HGVSp_Short"] = "WT"
        df.loc[df["HGVSc"] == "WT", "ccf_Mcopies"] = 0
        df.loc[df["HGVSc"] == "WT", "VAF"] = 0
        df.loc[df["HGVSc"] == "WT", "Prob MT p53"] = 0
        df = df[~df["data"].isna()&~df["Gene"].isna()] 
        samples = df["Sample_ID"].to_list()
        ## function to get data - if column missing return NaN
        #def get_column(dataframe, column):
        #    return dataframe.get(column, pd.Series(index=dataframe.index, 
        #                                           name=column))

        df["ADAR1 Term"] = [( (1 - locus_data[f"{s}_ADAR1_Edit"]) * locus_data[f"{s}_Expression"] ).sum() for s in samples]
        df["ORF1p Term"] = [(locus_data[f"{s}_Expression"] ).sum() for s in samples]
        # add to filter out NA values

####### This is the step where NA values are added 

        df["ADAR1 Term"] = (1 - df["Prob MT p53"]) * df["ADAR1 Term"]
        df["ORF1p Term"] = df["Prob MT p53"] * df["Insertion_z"]*df["ORF1p Term"]

        df["ADAR1 Term"] = (1 - df["ccf_Mcopies"]) * df["ADAR1 Term"]
        df["ORF1p Term"] = df["ccf_Mcopies"] * df["ORF1p Term"]

        df["IFN Term"] = [IFN_data[IFN_data["Sample_ID"] == s]["ifn_median"].sum() for s in samples]


        df = df[~df["data"].isna()&~df["ORF1p Term"].isna()]
        df[["ADAR1 Term (z-score)", "ORF1p Term (z-score)", "IFN Term (z-score)"]] = StandardScaler().fit_transform(df[["ADAR1 Term", "ORF1p Term", "IFN Term"]])
        
        return df
