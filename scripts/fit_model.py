#!/usr/bin/env python

# This script fits the model.

# import packages
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import pearsonr,spearmanr
import helper_funcs
helper_funcs = helper_funcs.helper_funcs()

# read in data
df = helper_funcs.organized_model_data()
df.to_csv("../results/immunity_model.tsv", sep="\t", index=False)
print(df)

# fit model
data_fracs = np.linspace(0.1, 0.9, num=9)
data_iters = np.arange(1, 4)
df_train_concat = []
df_test_concat = []
df_train_results = []
df_test_results = []
for data_frac in data_fracs:
    
    print("Data Fraction:", data_frac)

    #df_train = df.sample(frac=0.20, replace=False, random_state=42)
    #df_test = df.drop(df_train.index)

    for data_iter in data_iters:

        print("Data Iteration:", data_iter)

        df_train = df.sample(frac=data_frac, replace=False, random_state=42)
        df_test = df.drop(df_train.index)

        df_train["Training Data Fraction"] = data_frac
        df_train["Training Data Iteration"] = data_iter

        df_test["Training Data Fraction"] = data_frac
        df_test["Training Data Iteration"] = data_iter

        # gather weight correlations
        weights = np.linspace(0, 1)
        for i,w in enumerate(weights):
            df_train["Total Term"] = w * df_train["ADAR1 Term (z-score)"] + (1 - w) * df_train["ORF1p Term (z-score)"]
            df_train["ADAR1 Weight"] = w
            df_train["ORF1p Weight"] = 1 - w

            df_test["Total Term"] = w * df_test["ADAR1 Term (z-score)"] + (1 - w) * df_test["ORF1p Term (z-score)"]
            df_test["ADAR1 Weight"] = w
            df_test["ORF1p Weight"] = 1 - w

            df_train_concat.append(df_train.copy())
            df_test_concat.append(df_test.copy())

            linear_train_r, linear_train_pvalue = pearsonr(df_train["Total Term"], df_train["IFN Term (z-score)"])
            rank_train_r, rank_train_pvalue = spearmanr(df_train["Total Term"], df_train["IFN Term (z-score)"])

            linear_test_r, linear_test_pvalue = pearsonr(df_test["Total Term"], df_test["IFN Term (z-score)"])
            rank_test_r, rank_test_pvalue = spearmanr(df_test["Total Term"], df_test["IFN Term (z-score)"])

            df_train_results.append([data_frac, data_iter, w, 1 - w, 
                                     linear_train_r, linear_train_pvalue, 
                                     rank_train_r, rank_train_pvalue])
            df_test_results.append([data_frac, data_iter, w, 1 - w, 
                                    linear_test_r, linear_test_pvalue, 
                                    rank_test_r, rank_test_pvalue])

df_train_concat = pd.concat(df_train_concat, axis=0)
df_test_concat = pd.concat(df_test_concat, axis=0)

df_train_results = pd.DataFrame(df_train_results)
df_train_results.columns = ["Training Data Fraction", "Training Data Iteration", "ADAR1 Weight", "ORF1p Weight", 
                           "Linear Correlation", "Linear p-value", "Rank Correlation", "Rank p-value"]
df_test_results = pd.DataFrame(df_test_results)
df_test_results.columns = ["Training Data Fraction", "Training Data Iteration", "ADAR1 Weight", "ORF1p Weight", 
                          "Linear Correlation", "Linear p-value", "Rank Correlation", "Rank p-value"]

#################

# plot example

training_fraction = 0.4

sub_df_train_results = df_train_results[(df_train_results["Training Data Fraction"] == training_fraction) & 
                                        (df_train_results["Training Data Iteration"] == 2)]
sub_df_test_results = df_test_results[(df_test_results["Training Data Fraction"] == training_fraction) & 
                                      (df_test_results["Training Data Iteration"] == 2)]

opt_ADAR1_weight = sub_df_train_results[sub_df_train_results["Rank Correlation"] == sub_df_train_results["Rank Correlation"].max()]["ADAR1 Weight"].to_numpy()[0]

sub_df_train_concat = df_train_concat[(df_train_concat["Training Data Fraction"] == training_fraction) &
                                      (df_train_concat["Training Data Iteration"] == 2) & 
                                      (df_train_concat["ADAR1 Weight"] == opt_ADAR1_weight)]
sub_df_test_concat = df_test_concat[(df_test_concat["Training Data Fraction"] == training_fraction) &
                                    (df_test_concat["Training Data Iteration"] == 2) & 
                                    (df_test_concat["ADAR1 Weight"] == opt_ADAR1_weight)]

opt_ORF1p_weight = sub_df_train_results[sub_df_train_results["Rank Correlation"] == sub_df_train_results["Rank Correlation"].max()]["ORF1p Weight"].to_numpy()[0]

sub_df_train_concat = df_train_concat[(df_train_concat["Training Data Fraction"] == training_fraction) &
                                      (df_train_concat["Training Data Iteration"] == 2) & 
                                      (df_train_concat["ORF1p Weight"] == opt_ORF1p_weight)]
sub_df_test_concat = df_test_concat[(df_test_concat["Training Data Fraction"] == training_fraction) &
                                    (df_test_concat["Training Data Iteration"] == 2) & 
                                    (df_test_concat["ORF1p Weight"] == opt_ORF1p_weight)]

fig, ax = plt.subplots()
sns.lineplot(data=sub_df_train_results, x="ADAR1 Weight", y="Rank Correlation")
plt.savefig("../results/ADAR1-weight_vs_Rank_corr.pdf", bbox_inches="tight")

fig, ax = plt.subplots()
sns.lineplot(data=sub_df_train_results, x="ORF1p Weight", y="Rank Correlation")
plt.savefig("../results/ORF1p-weight_vs_Rank_corr.pdf", bbox_inches="tight")

fig, ax = plt.subplots()
sns.regplot(data=sub_df_train_concat, x="Total Term", y="IFN Term (z-score)")
plt.savefig("../results/Predicted_vs_Observed_IFN-training.pdf", bbox_inches="tight")
print("Training Results:", sub_df_train_results[sub_df_train_results["ADAR1 Weight"] == opt_ADAR1_weight])

fig, ax = plt.subplots()
sns.regplot(data=sub_df_test_concat, x="Total Term", y="IFN Term (z-score)")
plt.savefig("../results/Predicted_vs_Observed_IFN-testing.pdf", bbox_inches="tight")
print("Testing Results:", sub_df_test_results[sub_df_test_results["ADAR1 Weight"] == opt_ADAR1_weight])

fig, ax = plt.subplots()
sns.regplot(data=sub_df_test_concat, x="ORF1", y="IFN Term (z-score)")
plt.savefig("../results/Predicted_vs_Observed_IFN-testing.pdf", bbox_inches="tight")
print("Testing Results:", sub_df_test_results[sub_df_test_results["ADAR1 Weight"] == opt_ADAR1_weight])

plt.show()
