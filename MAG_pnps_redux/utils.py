#!/bin/python3
import csv
import pandas as pd
import sys
import numpy as np
sys.path.insert(0, '../particle_association_redux')
from third_get_secretory_CAZenzyme import get_all_signalp

depths2 = ["SRF", "DCM"]
depths3 = ["SRF", "DCM", "MES"]
depths4 = ["SRF", "DCM", "MES", "deep"]
has_deep = ["IN","SAT","NAT","SP","NP"]

data_root="/researchdrive/zhongj2/MAG_pnps_redux"
signal_root="/workspace/data/zhongj/MAG_CAZymes"

def MAG_db_fun():
	MAG_db = {"deep": ["deep"]}
	for ocean in ["RS", "MED", "EAC"]:
		MAG_db[ocean] = depths2
	for ocean in ["IN","SAT","NAT","SP","NP", "ARS", "CPC"]:
		MAG_db[ocean] = depths3
	return MAG_db

def ocean_depths_fun():
	ocean_depth ={"ARS": depths3, "CPC":depths3}
	for ocean in ["RS", "MED", "EAC"]:
		ocean_depth[ocean] = depths2
	for ocean in has_deep:
		ocean_depth[ocean] = depths4
	return ocean_depth

def read_pnps(root, ocean, depth): # ocean and depth is reversed for the Malaspina deep ocean
	cols = ["gene_callers_id", "pnps", "sample_id"]
	depth_pnps = pd.read_csv(f"{root}/{ocean}/PROFILE-{depth}/all_MAGs_pnps/pNpS.txt", sep = "\t").astype(str)
	if ocean == 'deep':
		depth_pnps["sample_id"] = depth_pnps["sample_id"].str.replace(f"{depth}_","")
	depth_pnps = depth_pnps[depth_pnps['pNpS_gene_reference'].astype(float) < 5] # sanity check
	depth_pnps.rename(columns= {"corresponding_gene_call":cols[0],"pNpS_gene_reference":"pnps"}, inplace=True)
	return depth_pnps[cols]

def get_bin_helper(fname):
	tmp = list(csv.reader(open(fname), delimiter=" "))[1:]
	bin_info = [[l[0]] + l[1].split('_', 1) for l in tmp]
	bin_info = pd.DataFrame(bin_info)
	bin_info.columns = ["gene_callers_id", "bin", "contig"]
	return bin_info

# gene = "toxin", "transposase"
def merge_blastp_f(df, root, ocean, gene, merge_method = "left"):
	gene_blastp=pd.read_csv(f"{root}/{ocean}/{gene}_diamond_unique.blastp", sep = "\t", header=None).astype(str)[[1,0]]
	gene_blastp.rename(columns={1:"gene_callers_id", 0:f"{gene}_id"}, inplace=True)
	df = df.merge(gene_blastp, on="gene_callers_id", how = merge_method)
	return df

def signalp_to_df_helper(signalp_pos_f):
	signalp_neg_f = signalp_pos_f.replace("gramPos", "gramNeg")
	pos_signalp_list = get_all_signalp(signalp_pos_f)
	neg_signalp_list = get_all_signalp(signalp_neg_f)
	signalp_set = set(pos_signalp_list + neg_signalp_list)
	signal_CAZyme_df = pd.DataFrame(list(signalp_set))
	signal_CAZyme_df.columns = ["gene_callers_id"]
	return signal_CAZyme_df

