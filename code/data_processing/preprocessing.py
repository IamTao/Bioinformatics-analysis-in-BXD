# -*- coding: utf-8 -*-

import pandas as pd
from os import listdir, makedirs
import shutil
from os.path import join, basename, splitext, exists


def read_csv(path):
    return pd.read_csv(path, delimiter=",", header=-1)


def list_files(path):
    return [join(path, f) for f in listdir(path)]

def build_dir(path):
    if exists(path):
        shutil.rmtree(path)
    makedirs(path)

def detect_two_parts(second_column):
    second_column = second_column.dropna()
    split_line = min(second_column.index)
    return split_line - 1


def get_pheno(end_index, data):
    pheno = data.ix[0: end_index - 1, :]
    pheno.drop([1, 2], inplace=True, axis=1)

    pheno_head = data.ix[end_index: end_index, :]
    pheno_head.drop([1, 2], inplace=True, axis=1)
    pheno.columns = pheno_head.values.tolist()

    # remove invalid information
    remove_invalid(pheno)
    pheno.ix[:, 0] = pheno.ix[:, 0].map(lambda symbol: symbol.replace("/", "-"))
    return pheno


def get_geno(start_index, data):
    geno = data.ix[start_index:, :]
    geno.drop([1, 2], inplace=True, axis=1)
    remove_invalid(geno)
    return geno


def remove_invalid(data):
    index_value = zip(data.ix[:, 0].index, data.ix[:, 0].values)
    invalid_index = [ind for ind, val in index_value
                     if "affy_" in val.lower()]
    data.drop(invalid_index, inplace=True, axis=0)


def write_to_file(path_out_root, file_name, geno, pheno):
    path_out_geno = join(path_out_root, file_name) + "_geno.csv"
    path_out_pheno = join(path_out_root, file_name) + "_pheno.csv"
    geno.to_csv(path_out_geno, sep='\t', index=False, header=False)
    pheno.to_csv(path_out_pheno, sep='\t', index=False, header=True)


def deal_single_geno_pheno(file_path, path_data_out_geno_pheno):
    file_name = basename(splitext(file_path)[0])
    in_data = read_csv(file_path)
    split_index = detect_two_parts(in_data.ix[:, 1])
    pheno = get_pheno(split_index, in_data)
    geno = get_geno(split_index, in_data)
    write_to_file(path_data_out_geno_pheno, file_name, geno, pheno)


def deal_all_geno_pheno(path_data_in_geno_pheno, path_data_out_geno_pheno):
    files = list_files(path_data_in_geno_pheno)
    for file in files:
        print file
        deal_single_geno_pheno(file, path_data_out_geno_pheno)


def main(path_root):
    path_data = join(path_root, "data")
    path_data_in = join(path_data, "input")
    path_data_out = join(path_data, "output")
    path_data_in_geno_pheno = join(path_data_in, "geno_pheno")
    path_data_out_geno_pheno = join(path_data_out, "geno_pheno_splitted_csv")
    build_dir(path_data_out_geno_pheno)
    # Debug
    # files = list_files(path_data_in_geno_pheno)
    # deal_single_geno_pheno(files[21], path_data_out_geno_pheno)

    # Total
    deal_all_geno_pheno(path_data_in_geno_pheno, path_data_out_geno_pheno)

if __name__ == '__main__':
    path_root = "/data/Dropbox/sv/gene_expression/new"
    path_root = "/Users/hali/Dropbox/sv/gene_expression/new"
    path_root = "/home/gene_expression"
    main(path_root)
