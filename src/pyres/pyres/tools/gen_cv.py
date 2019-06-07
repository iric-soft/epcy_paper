import random
import os

import pandas as pd

from .. utils import other as uo

def format_design(df, args):
    df = df.sort_values(by=[args.SUBGROUP, 'sample'], ascending=[False, True])
    df[args.SUBGROUP] = ["Query" if condition == 1 else "Ref" for condition in df[args.SUBGROUP]]

    return(df)

def main_gen_cv(args, argparser):

    df_design = uo.get_design(args)
    num_sample = len(df_design.index)

    fold = args.FOLD
    if args.LOO:
        fold = num_sample

    path_file_out = os.path.join(args.PATH, "cv", "fold" + str(fold))
    if args.LOO:
        path_file_out = os.path.join(args.PATH, "cv", "foldloo")

    if not os.path.exists(path_file_out):
        os.makedirs(path_file_out)

    num_test = round(num_sample / fold)

    df_design_rnd = df_design.sample(frac=1).reset_index(drop=True)

    cpt_sample = num_sample
    num_rnd = num_test
    cpt_fold = 0
    while cpt_sample > 0:
        path_fold = os.path.join(path_file_out, "dataset_" + str(cpt_fold))
        path_train = os.path.join(path_fold, "train")
        path_test = os.path.join(path_fold, "test")

        if not os.path.exists(path_train):
            os.makedirs(path_train)
        if not os.path.exists(path_test):
            os.makedirs(path_test)

        if num_rnd + 1 >= cpt_sample and not args.LOO:
            num_rnd = cpt_sample

        ids_test = random.sample(range(0, cpt_sample), num_rnd)
        df_test = df_design_rnd.iloc[ids_test,:]
        df_train = df_design[~df_design["sample"].isin(df_test["sample"])]

        df_design_rnd = df_design_rnd.drop(df_design_rnd.index[ids_test])

        df_test = format_design(df_test, args)
        df_train = format_design(df_train, args)

        df_test.to_csv(os.path.join(path_test, "design.tsv"), index=False, sep="\t")
        df_train.to_csv(os.path.join(path_train, "design.tsv"), index=False, sep="\t")

        cpt_sample -= num_rnd
        cpt_fold += 1
