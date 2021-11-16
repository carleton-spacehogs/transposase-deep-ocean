# cd /workspace/data/zhongj/Transposase_Project/integron_finder_tara_contig
# python3
import pandas as pd
import numpy as np
pNpS_integron_contigs = pd.read_csv(f"./for_pnps_integrons.csv", names = ["pnps", "ocean", "gene_callers_id"]) # 0.054658389874379085,ARS,1,N,N
all_integron = pd.read_csv(f"./all_ocean_merged_integrons.csv")

all_integron_with_pnps = pNpS_integron_contigs.merge(all_integron, on=["ocean", "gene_callers_id"], how="inner")[["ocean", "gene_callers_id","gene_function", "category", "pnps", "from_integron_type"]]

all_integron_with_pnps = all_integron_with_pnps.query("from_integron_type != 'In0'") # getting only the integron protein, exclude the transposases/integerase
all_integron_with_pnps.query("pnps < 10").pnps.describe()
'''
>>> all_integron_with_pnps.pnps.describe()
count    1238.000000
mean             inf
std              NaN
min         0.001616
25%         0.099435
50%         0.227482
75%         0.644687
max              inf
Name: pnps, dtype: float64

>>> all_integron_with_pnps.query("pnps < 10").pnps.describe()
count    1186.000000
mean        0.564597
std         1.066065
min         0.001616
25%         0.096967
50%         0.210317
75%         0.544255
max         9.726007
Name: pnps, dtype: float64
'''

all_integron_with_pnps['log_pnps'] = np.where(all_integron_with_pnps['pnps'] < 0.001, -3, np.log10(all_integron_with_pnps['pnps']))
all_integron_with_pnps['integron_type'] = np.where(all_integron_with_pnps['category'].isnull(), "no_call", 
np.where(all_integron_with_pnps['category'].notnull() & all_integron_with_pnps['category'].str.contains("Defense"), "defense", "non_defense"))

all_integron_with_pnps.to_csv('all_integron_with_pnps.csv', sep=',', index=False)
all_integron_with_pnps.query("category.isnull() & pnps < 10").pnps.describe()
'''
>>> all_integron_with_pnps.query("category.isnull() & pnps < 10").pnps.describe()
count    792.000000
mean       0.669321
std        1.168211
min        0.002900
25%        0.123438
50%        0.267709
75%        0.673285
max        9.726007
Name: pnps, dtype: float64
'''
integron_not_null_category = all_integron_with_pnps.query("category.notnull() & pnps < 10")

import scipy.stats
print(scipy.stats.ttest_ind(np.log2(all_integron_with_pnps.query("category.isnull() & pnps < 10").pnps),
np.log2(integron_not_null_category.pnps)))

integron_not_null_category.pnps.describe()
'''
>>> integron_not_null_category.pnps.describe()
count    394.000000
mean       0.354087
std        0.783266
min        0.001616
25%        0.066882
50%        0.142721
75%        0.306189
max        7.846782
Name: pnps, dtype: float64
'''

defense_mech = integron_not_null_category[integron_not_null_category['category'].str.contains("Defense")]
defense_mech.pnps.describe()
'''
>>> integron_not_null_category[integron_not_null_category['category'].str.contains("Defense")].pnps.describe()
count    57.000000
mean      0.258096
std       0.732204
min       0.003440
25%       0.045668
50%       0.086927
75%       0.218023
max       5.463742
Name: pnps, dtype: float64
'''
non_defense_mech = integron_not_null_category[~integron_not_null_category['category'].str.contains("Defense")]
non_defense_mech.pnps.describe()

defense_mech = integron_not_null_category[integron_not_null_category['category'].str.contains("Defense")]
non_defense_mech = integron_not_null_category[~integron_not_null_category['category'].str.contains("Defense")]
print(scipy.stats.ttest_ind(np.log2(non_defense_mech.pnps), np.log2(defense_mech.pnps)))


tmp = all_integron_with_pnps.query("category.notnull()")
tmp[tmp['category'].str.contains("Defense")]