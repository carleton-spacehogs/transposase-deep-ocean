#!/bin/python3
import csv
import pandas as pd
import sys
import numpy as np

depths2 = ["SRF", "DCM"]
depths3 = ["SRF", "DCM", "MES"]
depths4 = ["SRF", "DCM", "MES", "deep"]

data_root="/researchdrive/zhongj2/MAG_pnps_redux"

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
	for ocean in ["IN","SAT","NAT","SP","NP"]:
		ocean_depth[ocean] = depths4
	return ocean_depth
