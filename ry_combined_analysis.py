# %% 
import pandas as pd
import numpy as np
import mosaic.io as mio
from tea.plots import plot_snv_clone
from tea.snv.me import get_exclusivity
import seaborn as sns
import matplotlib.pyplot as plt
from tea.utils import sort_for_var
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# load references
ref = "ref/Rona_CRC_tapestri_design.xlsx"
ref = pd.read_excel(ref)
ref = ref[ref["type"].isin(["snv", "sv"])]
ref['condensed_format'] = "chr" + ref["chromosome"].astype(str) + ":" + ref["start"].astype(int).astype(str)
ref["ann"] = ref["gene"] + " p." + ref["aa_change"]


# %% load H5
h5 = "../pipeline_results/RY_pt.patient_wide.genotyped.h5"
smp = mio.load(h5)
smp.dna.genotype_variants()

# %% get per cell mutational burden
PANLE_SIZE=24110
NUM_CELLS = smp.dna.shape[0]
# filter to mutations in >0.5% of cells
cell_per_mut = smp.dna.get_attribute('mut_filtered', constraint="row").sum(axis=0)
cell_per_mut = cell_per_mut[
    (cell_per_mut > 0.005 * NUM_CELLS) & (cell_per_mut < 0.7 * NUM_CELLS)
    ]
voi = cell_per_mut.index.tolist()

# %% plot single cell heatmap
ATTRIBUTE = "AF_MISSING"
TOPIC = "bulk snvs"
if "_Assay__heatmap" in smp.dna.__dict__:
    del smp.dna._Assay__heatmap
fig = plot_snv_clone(
        smp,
        sample_name="RY_Tapestri",
        story_topic = f'high_conf_mutations-{TOPIC}',
        voi = voi,
        attribute = ATTRIBUTE,
        # vars_to_sort_by = ["chr12:25398285:C/A"],
        # barcode_sort_method = "single_var",
    )
fig

# %%
# get the mutational burden
mut_per_cell = smp.dna.get_attribute('mut_filtered', constraint="row", features=voi).sum(axis=1)
median = mut_per_cell.median()

# scale by total genome coverage


# %% plot mutational burden
fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(mut_burden, total_read_count)
ax.set_xlabel('Mutational Burden')
ax.set_ylabel('Total Read Count')
plt.show()

# %% assign cell labels based on SNVs
# fetch all the putative normal cells
# using voi = ['chr3:178936091:G/A','chr4:55152121:G/A','chr6:117677847:T/C','chr7:92404094:A/G','chr12:25398285:C/A','chr17:7577025:T/TA','chrX:47430853:T/C','chrX:66863172:G/A']
VOI = ['chr3:178936091:G/A','chr4:55152121:G/A','chr12:25398285:C/A','chr17:7577025:T/TA'] # <-- only use KRAS and TP53 to define tumors
tumor_cells_ry1 = smp.dna.barcodes()[
    (smp.dna.get_attribute('alt_read_count', constraint='row', features=VOI
        ) > 0).any(axis=1) & (smp.dna.row_attrs["sample_name"] == "RY")
    ]
tumor_cells_ry3 = smp.dna.barcodes()[
    (
        smp.dna.get_attribute('alt_read_count', constraint='row', features=VOI
        ) > 0
        ).any(axis=1) & (smp.dna.row_attrs["sample_name"] == "RY3")
    ]

# assign labels
smp.dna.row_attrs['label'] = np.array(list(
    map(
        lambda x: "Core 1 tumor" if x in tumor_cells_ry1 else "Core 2 tumor" if x in tumor_cells_ry3 else "Normal", 
        smp.dna.barcodes())
    ))
cn_clone_palette = dict(zip(
    ["Core 1 tumor", "Core 2 tumor", "Normal"], 
    ["#b0583a", "#443ab0", "#679665"]
    ))
smp.dna.set_palette(cn_clone_palette)

# # %% select variants that are present in the design sheet
# voi = [v for v in smp.dna.ids() if v.rsplit(":",1)[0] in ref["condensed_format"].unique()]
# ann_map = {}
# for v in voi:
#     matched_row = ref[ref["condensed_format"] == v.rsplit(":",1)[0]]
#     ann_map[v] = matched_row["ann"].values[0]

# # get the VOI's pseudobulk VAF
# alt_sum = smp.dna.get_attribute("alt_read_count", constraint="row", features=voi).sum(axis=0)
# dp_sum = smp.dna.get_attribute("DP", constraint="row", features=voi).sum(axis=0)
# vaf = alt_sum / dp_sum
# %% look around the KRAS locus
for i in smp.dna.ids():
    if "chr12:2539" in i:
        print(i)

# %%
# del smp.dna.__dict__['_Assay__heatmap']
ann_map = {
    "chr3:178936091:G/A": "PIK3CA p.E545K",
    "chr4:55152121:G/A": "PDGFRA p.S851S",
    "chr12:25398285:C/A": "KRAS p.G12C",
    "chr17:7577025:T/TA": "TP53 p.K305*"
}
ATTRIBUTE = "AF_MISSING"
TOPIC = "bulk snvs"
fig = plot_snv_clone(
        smp,
        sample_name="RY_Tapestri",
        story_topic = f'high_conf_mutations-{TOPIC}',
        voi = VOI,
        attribute = ATTRIBUTE,
        ann_map = ann_map,
        # vars_to_sort_by = ['chr3:178936091:G/A','chr12:25398285:C/A','chr17:7577025:T/TA'],
        # barcode_sort_method = "single_var",
    )
fig


#%%
fig.write_image("RY_combined_sc_heatmap.pdf", scale=4)
os.makedirs("raw_data_for_figures", exist_ok=True)
smp.dna.get_attribute("AF_MISSING", features = VOI, constraint="row").to_csv(
    "raw_data_for_figures/RY_combined_sc_heatmap_af_missing.csv"
)
# %% get mutual exclusivity of these vars
ADO_PRECISION = 15 # precision parameter
FP = 0.001 # false positive rate

g00 = 0
g01 = 1 
g10 = 1
g11 = 0

E = get_exclusivity(
    smp.dna.get_attribute('DP', constraint='row', features=list(ann_map.keys())),
    smp.dna.get_attribute('alt_read_count', constraint='row', features=list(ann_map.keys())),
    ADO_PRECISION, FP,
    gametes = [(g00,g01), (g10,g11)],
    # rm_irrelevant_cells=False
)

E.rename(index=ann_map, columns=ann_map, inplace=True)
E = 1-E
# %% Plot and save
fig, ax = plt.subplots(figsize=(15,15))         # Sample figsize in inches
sns.heatmap(
    E, annot=True, 
    cmap='Blues', fmt='.2f', vmin=0, vmax=1,
    ax=ax
    )
# ax.set_title("Mutual Exclusivity")
fig.savefig("RY_Tapestri_mutual_exclusivity.pdf")
# %% CNV setup
smp.cnv.metadata['genome_version'] = "hg19"
smp.cnv.get_gene_names()

# %% inspect SV breakpoints
"""
TFG-MET: AMPL1006857
MET-TGF: ["AMPL1006934", "AMPL1006935"]
SMAD4: ["AMPL479589", "AMPL479590", "AMPL479591", "AMPL479592", "AMPL164287", "AMPL479594", "AMPL479595", "AMPL479596", "AMPL479597", "AMPL459008", "AMPL479599", "AMPL26083"]
"""
raw_rc = smp.cnv.get_attribute('read_counts', constraint='row')
rc_binary = (raw_rc == 0).astype(int)
smp.cnv.add_layer('read_counts_binary[zero/nonzero]', rc_binary.values)

attribute="AF_MISSING"
sorted_bars = sort_for_var(
    dna = smp.cnv,
    vars = VOI,
    attribute = attribute,
    method = "hier"
    )

# %% ===== CNV analysis =====
# %% ----- fit NB model -----
from tea.cnv.train import fit_NB_model, extract_NB_params
normal_cells = [i for i in smp.dna.barcodes() if (i not in tumor_cells_ry1) and (i not in tumor_cells_ry3)]
normal_cells_rc = smp.cnv.get_attribute('read_counts', constraint='row').loc[normal_cells,:]
sc_total_rc = normal_cells_rc.sum(axis = 1)
# X = np.log(sc_total_rc).to_frame(name='total')
X = pd.DataFrame(
    [1.0] * len(sc_total_rc), 
    index = sc_total_rc.index, 
    columns=['constant']
    )

nb_res = pd.Series(
    [fit_NB_model(normal_cells_rc[amplicon_i], X, sc_total_rc) for amplicon_i in normal_cells_rc.columns],
    index = normal_cells_rc.columns
    ) # N x 1 array storing the training results of each amplicon
print('[SUCCESS] NB training finished successfully')

# %%
nb_res_df = pd.DataFrame(
    index = nb_res.index,
    columns = ['amplicon', 'converged', 'amplicon_factor', 'alpha', 'beta_zero', 'beta_one', 'method', 'mean', 'variance']
    )
for amplicon_i in nb_res.index:
    nb_res_df.loc[amplicon_i, :] = extract_NB_params(nb_res[amplicon_i], amplicon_i, normal_cells_rc[amplicon_i])
    if nb_res_df.loc[amplicon_i, 'converged'] == False:
        print('=' * 8)
        print(f'[WARNING] {amplicon_i} did not converge')

ref_df = pd.read_csv("ref/ry_amp_gene_map.csv", )
ref_df.index = ref_df['amplicon_number']
nb_res_df['gene'] = ref_df.loc[
    nb_res_df['amplicon'],
    'gene_name'
] # get mapped gene names
nb_res_df.to_csv('cn_analysis/ry_panel_train_normal_results.csv', index=False, header=True)

print('[INFO] NB training results saved to:', 'cn_analysis/ry_panel_train_normal_results.tsv')

# %%
import plotly.express as px
from tea.utils import rgb_string_to_hex

smp.cnv.row_attrs['label'] = smp.dna.row_attrs['label']
smp.cnv.set_palette(cn_clone_palette)
smp.cnv.get_gene_names('cn_analysis/ry_panel_train_normal_results.csv', gene_name_col = 'gene')
smp.cnv.var = pd.DataFrame.from_dict(smp.cnv.col_attrs).set_index("id")
raw_rc = smp.cnv.get_attribute('read_counts', constraint='row')

# ###### HZ warning: I don't think this is correct.
# normalized_rc = raw_rc / nb_res_df.loc[raw_rc.columns]['amplicon_factor'].astype(float) / raw_rc.sum(axis=1).values[:, None]
# ###### ----------------------------------------
# use the naive method instead
# smp.cnv.normalize_reads(
#     method = 'hz', 
#     diploid_cells = np.array(normal_cells)
#     )

rc = smp.cnv.get_attribute('read_counts', constraint='row')
normal_rc = rc.loc[normal_cells, :]

rc /= normal_rc.mean(axis=0) + 1
rc /= (np.array(rc.mean(axis=1))[:, np.newaxis]) + 1
rc *= 2
normalized_rc = rc.to_numpy()
smp.cnv.add_layer('normalized_counts', normalized_rc)

normalized_rc_binary = (normalized_rc == 0).astype(int)
smp.cnv.add_layer('normalized_counts_binary[zero/nonzero]', normalized_rc_binary)


# %%
grouped_bars_order = sort_for_var(
    dna = smp.dna,
    vars = VOI,
    attribute = "AF_MISSING",
    method = "hier"
)
sc_cnv_heatmap = smp.cnv.heatmap(
    'normalized_counts',
    features = ["12"],
    bars_order = grouped_bars_order
)
sc_cnv_heatmap
sc_cnv_heatmap.write_image(
    "RY_Tapestri_snv_ordered_amplicon_normalized_rc_heatmap_chr12.pdf", 
    scale=4, width=600, height=1000
    )
# save the raw matrix
smp.cnv.get_attribute(
    'normalized_counts', constraint='row', 
    features=["AMPL117564", "AMPL401249", "AMPL401250", "AMPL401252", "AMPL401253"]
    ).to_csv(
        "raw_data_for_figures/RY_Tapestri_snv_ordered_amplicon_normalized_rc_heatmap_chr12.csv"
    )
# %%
VOI = ['chr3:178936091:G/A','chr4:55152121:G/A','chr12:25398285:C/A','chr17:7577025:T/TA']
sc_dna_heatmap = smp.dna.heatmap(
    'AF_MISSING',
    features = VOI,
    bars_order = grouped_bars_order
)
sc_dna_heatmap
sc_dna_heatmap.write_image(
    "RY_Tapestri_snv_ordered_amplicon_AF_MISSING_heatmap.pdf", 
    scale=4, width=600, height=1000
    )
del smp.cnv.__dict__['_Assay__heatmap']

# %% get KRAS CN
cn1 = smp[
    (smp.dna.row_attrs["sample_name"] == "RY") &
    (smp.dna.get_attribute("alt_read_count", features=["chr12:25398285:C/A"])[0]>0).values
    ].cnv.get_attribute('normalized_counts', features= ["AMPL401253"], constraint="row")
cn2 = smp[
    (smp.dna.row_attrs["sample_name"] == "RY3") &
    (smp.dna.get_attribute("alt_read_count", features=["chr12:25398285:C/A"])[0]>0).values
    ].cnv.get_attribute('normalized_counts', features= ["AMPL401253"], constraint="row")


# %% KRAS mutation VAF
af1 = smp[
    (smp.dna.row_attrs["sample_name"] == "RY") &
    (smp.dna.get_attribute("alt_read_count", features=["chr12:25398285:C/A"])[0]>0).values
    ].dna.get_attribute("AF", features= ["chr12:25398285:C/A"], constraint="row")
print(af1.median())
# write af1 to a csv
af1.to_csv(
    "raw_data_for_figures/RY1_KRAS_mutation_VAF.csv"
    )
af2 = smp[
    (smp.dna.row_attrs["sample_name"] == "RY3") &
    (smp.dna.get_attribute("alt_read_count", features=["chr12:25398285:C/A"])[0]>0).values
    ].dna.get_attribute("AF", features= ["chr12:25398285:C/A"], constraint="row")
print(af2.median())
# write af2 to a csv
af2.to_csv(
    "raw_data_for_figures/RY3_KRAS_mutation_VAF.csv"
    )

# %%
s1 = smp[smp.dna.row_attrs["sample_name"] == "RY"]
s2 = smp[smp.dna.row_attrs["sample_name"] == "RY3"]

# %%
voi = [v for v in smp.dna.ids() if v.rsplit(":",1)[0] in ref["condensed_format"].unique()]
ann_map = {}
for v in voi:
    matched_row = ref[ref["condensed_format"] == v.rsplit(":",1)[0]]
    ann_map[v] = matched_row["ann"].values[0]

# get the VOI's pseudobulk VAF
alt_sum = smp.dna.get_attribute("alt_read_count", constraint="row", features=voi).sum(axis=0)
dp_sum = smp.dna.get_attribute("DP", constraint="row", features=voi).sum(axis=0)
vaf = alt_sum / dp_sum

# %%
