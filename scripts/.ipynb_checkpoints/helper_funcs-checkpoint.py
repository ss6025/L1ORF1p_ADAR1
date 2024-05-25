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
        df["Sample_ID"] = df["Sample_ID"].str.replace("T0", "T")
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
        data_file = "../data/MAFs/TP53_mut.csv"
        df = pd.read_csv(data_file, sep=",")
        # add HGVSp
        data_file = "../data/MAFs/HGVSc_to_HGVSp_mapping.tsv"
        mapping_df = pd.read_csv(data_file, sep="\t")
        mapping_dict = dict(zip(mapping_df["HGVSc"], mapping_df["HGVSp_Short"]))
        df["HGVSp_Short"] = df["HGVSc"].map(mapping_dict)
        # add fitness
        fit_df = self.p53_fitness()
        fit_dict = dict(zip(fit_df["HGVSp_Short"], fit_df["Fitness"]))
        df["Prob MT p53"] = df["HGVSp_Short"].map(fit_dict).fillna(0)
        # rename columns
        df.rename(columns={"SYMBOL" : "Gene", "Allele.Freq.T" : "Unfiltered VAF", 
            "Allele.Freq.T.filtered" : "VAF"}, inplace=True)
        df = df[["Patient_ID", "Sample_ID", "Gene", "HGVSc", "HGVSp_Short", 
                 "Prob MT p53", "Unfiltered VAF", "VAF"]]
        df["Sample_ID"] = df["Sample_ID"].str.replace("T0", "T")
        ## prob MT p53 is mislabled, correct it mannually
        #df.loc[df['Patient_ID'] == "MPAM05", 'Prob MT p53'] = 1.0
        return df

    def fpkm_SINE_dsRNAforce(self):
        '''
        Data on SINE and dsRNA force.
        '''
        data_file = "../data/fpkm_SINE_filtered.csv"
        df = pd.read_csv(data_file, sep=",")
        df.columns = df.columns.str.replace("T0", "T")
        return df

    def IFN_data(self):
        '''
        IFN data.
        '''
        data_file = "../data/immune_P53_adar1.csv"
        df = pd.read_csv(data_file, sep=",")
        df["Sample_ID"] = df["Sample_ID"].str.replace("T0", "T")
        return df

    def ADAR1_editing(self):
        '''
        Locus-specific ADAR1 editing.
        '''
        data_file = "../data/sine_locus_editing_filtered.csv"
        df = pd.read_csv(data_file, sep=",")
        df.columns = df.columns.str.replace("T0", "T")
        return df

    def organized_model_data(self):
        '''
        Organized model data.
        '''
        # read in data
        RIPseq_data = pd.read_csv("../results/RIPseq-TE_fraction_ORF1p.tsv", sep="\t")
        ADAR1_data = self.ADAR1_editing()
        fpkmSINE_data = self.fpkm_SINE_dsRNAforce()
        IFN_data = self.IFN_data()
        mutation_data = self.TP53_mut_df()
        sample_data = self.sample_data()

        # prep locus-level data
        locus_data = pd.merge(RIPseq_data, ADAR1_data, left_on="TE_ID", right_on="rep.id").drop(["rep.id"], axis=1)
        #fpkmSINE_cols = ["rep.id"] + list(fpkmSINE_data.columns[10:])
        locus_data = pd.merge(locus_data, fpkmSINE_data, left_on="TE_ID", right_on="rep.id", suffixes=["_ADAR1_Edit", "_Expression"])
        locus_data = locus_data[locus_data["TE_ID"].str.contains("Alu")]

        # organize data
        df = pd.merge(sample_data, mutation_data, on=["Sample_ID","Patient_ID"], how='left')
        #df["HGVSc"] = df["HGVSc"].fillna("WT")
        df.loc[df["TP53_call"] == "WT", "Gene"] = "TP53"
        df.loc[df["TP53_call"] == "WT", "HGVSc"] = "WT"
        df.loc[df["HGVSc"] == "WT", "HGVSp_Short"] = "WT"
        df.loc[df["HGVSc"] == "WT", "VAF"] = 0
        df.loc[df["HGVSc"] == "WT", "Prob MT p53"] = 0
        samples = df["Sample_ID"].to_list()

        ## function to get data - if column missing return NaN
        #def get_column(dataframe, column):
        #    return dataframe.get(column, pd.Series(index=dataframe.index, 
        #                                           name=column))

        #df["ADAR1 Term"] = [( (1 - locus_data[f"{s}_ADAR1_Edit"]) * locus_data[f"{s}_Expression"] ).sum() for s in samples]
        df["ADAR1 Term"] = [(1 - IFN_data.loc[IFN_data['Sample_ID'] == s, 'all'].iloc[0]) * locus_data[f"{s}_Expression"].max() for s in samples]
        df["ORF1p Term"] = [( (1 - locus_data["median_frac_ORF1p"]) * locus_data[f"{s}_Expression"] ).sum() for s in samples]
        #df["ADAR1 Term"] = [( (1 - get_column(locus_data, f"{s}_ADAR1_Edit")) * get_column(locus_data, f"{s}_Expression") ).sum() for s in samples]
        #df["ORF1p Term"] = [( (1 - locus_data["median_frac_ORF1p"]) * get_column(locus_data, f"{s}_Expression") ).sum() for s in samples]


####### This is the step where NA values are added 

        df["ADAR1 Term"] = (1 - df["Prob MT p53"]) * df["ADAR1 Term"]
        df["ORF1p Term"] = df["Prob MT p53"] * df["ORF1p Term"]

        df["ADAR1 Term"] = (1 - df["VAF"]) * df["ADAR1 Term"]
        df["ORF1p Term"] = df["VAF"] * df["ORF1p Term"]

        df["IFN Term"] = [IFN_data[IFN_data["Sample_ID"] == s]["ifn_median"].sum() for s in samples]

        #df = df[df["VAF"] != 0]

        #df_train = df.sample(frac=0.80, replace=False, random_state=42)
        #df_test = df.drop(df_train.index)

        df[["ADAR1 Term (z-score)", "ORF1p Term (z-score)", "IFN Term (z-score)"]] = StandardScaler().fit_transform(df[["ADAR1 Term", "ORF1p Term", "IFN Term"]])
        return df
