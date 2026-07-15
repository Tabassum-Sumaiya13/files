import numpy as np
import pandas as pd

from sklearn.model_selection import KFold, GroupKFold, StratifiedGroupKFold, LeaveOneGroupOut
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored

import random
import sys
import argparse

from utils import *
from rsf import get_rsf_results_for_feature_df

parser = argparse.ArgumentParser()
parser.add_argument('--channels', help='path to csv listing all biomarkers/channels in one column', type=str, default='../data/marker_names.csv')
parser.add_argument('--sample-info', help='path to csv with sample/region names and labels', type=str, default='../data/sample_metadata.csv')
parser.add_argument('--expr-data', help='path to pkl file with celltype labels and normalized expression values (arcsinh only)', type=str, default='../data/labeled_arcsinh_norm_data.pkl')
parser.add_argument('--samples', help='path to csv listing which sample/region ids to use', type=str, default='../data/qc_acq_ids_labeled.csv')
parser.add_argument('--cv-folds', help='which type of fold to train. Default is patient folds.', choices=['patient', 'coverslip', 'random'], default='patient')
parser.add_argument('--output-dir', help='name of output directory for rsf results', type=str)
parser.add_argument('--features', help='which features to use for training rsf. Can list multiple values here.', choices=['neighbor_mat', 'celltype_prop', 'biomarker_cell', 'biomarker_region', 'ripley', 'DenVar', 'biomarker_threshold', 'biomarker_interaction_counts'], nargs='+') 
parser.add_argument('--permute-labels', help='set this argument if training rsf using permuted labels.', action='store_true')
parser.add_argument('--precomputed', help='set this argument if rsf estimators have already been trained. Models should be saved in the {output_dir}/estimators directory. If --permute-labels flag is set, this instead looks for the concordance index array in {output_dir}/{title}_{num_permutations}permutations_concordance.npy.', action='store_true')
parser.add_argument('--primary-only', help='set this argument to use primary site samples only', action='store_true')
parser.add_argument('--exclude-normal', help='set this argument to exclude normal mucosa samples', action='store_true')
parser.add_argument('--combine-features', help='set this argument to combine all specified feature sets for training and prediction. Currently only supports combining celltype proportions and neighborhood matrix features.', action='store_true')
parser.add_argument('--ripley-df', help='path to csv file with Ripley K function values.', type=str, default="../data/k_fns_norm_by_uw_qc_labeled.csv")
parser.add_argument('--ripley-vals', help='radii values to use from Ripley K function values, in integer microns.', type=int, nargs='+', default=[30, 80])
parser.add_argument('--denvar-dir', help='path to directory with DenVar cluster values', type=str, default='../data/denvar')
parser.add_argument('--biomarker-pos-frac-df', help='path to csv file with the biomarker positivity fraction values for each sample. Method used in Patwa et al. (2021)', type=str, default='../data/comparisons/biomarker_frac_positivity_qc_labeled.csv')
parser.add_argument('--biomarker-interaction-df', help='path to csv file with the biomarker interaction counts for each sample. Method used in Patwa et al. (2021)', type=str, default='../data/comparisons/interaction_biomarker_features.csv')
parser.add_argument('--n-repeats', help='number of random seed repeats to perform and save. default 0', type=int, default=0)
args = parser.parse_args()

if args.primary_only and args.exclude_normal:
    raise AssertionError("Can only set either the --primary-only or --exclude-normal argument, but both were used.")

# Loading channels
marker_names = list(pd.read_csv(args.channels).iloc[:, 0].values)
num_markers = len(marker_names)

sample_df = pd.read_csv(args.sample_info)
sample_list = np.array(pd.read_csv(args.samples, header=None).iloc[:,0].values, dtype=str)
# Only use at samples in the sample_list
sample_df = sample_df[sample_df.acquisition_id.isin(sample_list)]

if args.primary_only:
    sample_df = sample_df[sample_df.tissuetype == 'Primary tumor']
    print('Using primary tumor samples only')
elif args.exclude_normal:
    sample_df = sample_df[sample_df.tissuetype != 'Normal mucosa']
    print('Excluding normal samples')
else:
    print('Using all samples')

expr = pd.read_pickle(args.expr_data)
clusters = expr.cluster.values
cluster_names = sorted(np.unique(expr.cluster_label))
num_clusters = len(np.unique(expr.cluster))

output_dir = args.output_dir

label_cols = ["survival_day", "survival_status"]

if args.combine_features:
    all_dfs_dict = {}
    num_features = 0

if 'biomarker_region' in args.features:
    print('Average biomarker intensity across the region')

    input_df = create_input_dataframe_for_features("biomarker_region", sample_df, label_cols, sample_list=sample_list, marker_names=marker_names)

    if args.combine_features:
        input_df.rename(columns={x: f"{x}_biomarker_region" for x in marker_names}, inplace=True)
        all_dfs_dict['biomarker_region'] = input_df
        num_features += num_markers
    else:
        get_rsf_results_for_feature_df(input_df, 'biomarker_region', np.arange(0,num_markers), 'avg_biomarker_region', output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=marker_names, n_repeats=args.n_repeats)

if 'biomarker_cell' in args.features:
    print('Average biomarker intensity within cells per region')

    input_df = create_input_dataframe_for_features("biomarker_cell", sample_df, label_cols, expr=expr, marker_names=marker_names)

    if args.combine_features:
        input_df.rename(columns={x: f"{x}_biomarker_cell" for x in marker_names}, inplace=True)
        all_dfs_dict['biomarker_cell'] = input_df
        num_features += num_markers
    else:
        get_rsf_results_for_feature_df(input_df, 'biomarker_cell', np.arange(0,num_markers), 'avg_biomarker_cell', output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=marker_names, n_repeats=args.n_repeats)

if 'celltype_prop' in args.features:
    print('Celltype proportion features')

    input_df = create_input_dataframe_for_features("celltype_prop", sample_df, label_cols, expr=expr, num_clusters=num_clusters)

    if args.combine_features:
        input_df.rename(columns={x: f'{x}_celltype' for x in np.arange(0, num_clusters)}, inplace=True)
        all_dfs_dict['celltype_prop'] = input_df
        num_features += num_clusters
    else:
        get_rsf_results_for_feature_df(input_df, 'celltype_prop', np.arange(0,num_clusters), f'celltype_prop', output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=cluster_names, n_repeats=args.n_repeats)

# Neighborhood matrix
if 'neighbor_mat' in args.features:
    k = 10
    print(f'Neighborhood matrix k = {k}')
    flatten_dim = num_clusters**2
    neighbor_ids = list(np.arange(flatten_dim))
    input_df = create_input_dataframe_for_features("neighbor_mat", sample_df, label_cols, expr=expr, num_clusters=num_clusters, sample_list=sample_list, k=k)
    
    survival_labels = input_df.loc[:,['survival_status', 'survival_day']].values
  
    if args.combine_features:
        input_df.rename(columns={x: f'{x}_neighbor_mat' for x in np.arange(0, flatten_dim)}, inplace=True)
        all_dfs_dict['neighbor_mat'] = input_df
        num_features += flatten_dim
    else:
        get_rsf_results_for_feature_df(input_df, 'neighbor_mat', np.arange(0,flatten_dim), f'neighbor_mat_k{k}', output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, neighbor_ids=neighbor_ids, cluster_names=cluster_names, n_repeats=args.n_repeats)

# Ripley's K function as features. Each biomarker for each sample has a K function curve from 1-401 pixels.
if 'ripley' in args.features:
    radii_um = args.ripley_vals
    um_per_pix = 0.377
    radii_pix = [round(x/um_per_pix) for x in radii_um]
    radii_names = [f'X{r}' for r in radii_pix]
    print(f"Ripley's K function. Using features at {radii_um} microns.")

    ripley_df = pd.read_csv(args.ripley_df, index_col=0)
    ripley_df['acquisition_id'] = ["_".join(x.split("_")[:5]) for x in ripley_df.index]
    ripley_df["biomarker"] = ["_".join(x.split("_")[5:-2]) for x in ripley_df.index]

    ripley_df = ripley_df[ripley_df.acquisition_id.isin(sample_list)]
    ripley_df = ripley_df[ripley_df.biomarker.isin(marker_names)]
    ripley_df = ripley_df.loc[:,['acquisition_id', 'biomarker'] + radii_names]
    ripley_df = ripley_df.pivot(index='acquisition_id', columns='biomarker', values=radii_names)
    ripley_df = ripley_df.reindex(marker_names, axis=1, level=1)
    col_names = [f'{i}_{j}' for i, j in ripley_df.columns]
    ripley_df.columns = col_names

    ripley_df = sample_df.merge(ripley_df, how='inner', left_on='acquisition_id', right_index=True)

    input_df = ripley_df[col_names + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]

    if args.combine_features:
        input_df.rename(columns={x: f"{x}_ripley" for x in col_names}, inplace=True)
        all_dfs_dict['ripley'] = input_df
        num_features += len(col_names)
    else:
        get_rsf_results_for_feature_df(input_df, 'ripley', np.arange(0,len(col_names)), f"ripley_{'-'.join(str(x) for x in radii_um)}_um", output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=col_names, n_repeats=args.n_repeats)

# DenVar clusters as features
if 'DenVar' in args.features:
    DenVar_df = aggregate_DenVar_clusters(args.denvar_dir, sample_list, marker_names)
    DenVar_df = sample_df.merge(DenVar_df, how='inner', on='acquisition_id')
    input_df = DenVar_df[marker_names + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]

    get_rsf_results_for_feature_df(input_df, 'DenVar', np.arange(0, len(marker_names)), f'DenVar_clusters', output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=marker_names, n_repeats=args.n_repeats)

if 'biomarker_threshold' in args.features:
    print("Biomarker positivity fraction based on background thresholds")
    biomarker_pos_frac_df = pd.read_csv(args.biomarker_pos_frac_df)
    biomarker_pos_frac_df.rename(columns={'sample_id':'acquisition_id'}, inplace=True)
    biomarker_pos_frac_df = sample_df.merge(biomarker_pos_frac_df, how='inner', on='acquisition_id')
    input_df = biomarker_pos_frac_df[marker_names + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]

    get_rsf_results_for_feature_df(input_df, 'biomarker_threshold', np.arange(0, len(marker_names)), f'biomarker_pos_frac', output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=marker_names, n_repeats=args.n_repeats)

if 'biomarker_interaction_counts' in args.features:
    print("Biomarker pair interaction counts based on cell-cell positivities")
    interaction_df = pd.read_csv(args.biomarker_interaction_df)
    interaction_features = list(interaction_df.columns.values[:-1])

    interaction_df = sample_df.merge(interaction_df, how='inner', on='acquisition_id')

    input_df = interaction_df[interaction_features + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]

    get_rsf_results_for_feature_df(input_df, 'biomarker_interaction_counts', np.arange(0, len(marker_names)), f'biomarker_interaction_counts', output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=interaction_features, n_repeats=args.n_repeats)

# Combining features for prediction
# TODO add support for combining arbitrary feature sets
if args.combine_features:
    base_order = ['biomarker_region', 'biomarker_cell', 'celltype_prop', 'neighbor_mat', 'ripley']
    reordered_features = [x for x in base_order if x in args.features]
    print(f"Combining {reordered_features[0]} and {reordered_features[1]} features")
    combined_title = "_".join(reordered_features)
    output_title = combined_title
    if 'neighbor_mat' in reordered_features:
        k = 10
        output_title = output_title.replace('neighbor_mat', f'neighbor_mat_k{k}')
    if 'ripley' in reordered_features:
        radii_um = args.ripley_vals
        output_title = output_title.replace('ripley', f"ripley_{'-'.join(str(x) for x in radii_um)}_um")

    num_feature_sets = len(reordered_features)
    keys = ['survival_day', 'survival_status', 'patient_id', 'coverslip_label', 'acquisition_id']
    for i, feature in enumerate(reordered_features):
        if i == 0:
            input_df = all_dfs_dict[feature]
        else:
            input_df = input_df.merge(all_dfs_dict[feature], how='inner', on=keys)

    input_df = input_df[[x for x in input_df.columns if x not in keys] + keys]

    if ('neighbor_mat' in args.features and 'ripley' in args.features) or ('neighbor_mat' in args.features and 'markcorr' in args.features):
        get_rsf_results_for_feature_df(input_df, combined_title, np.arange(0,num_features), output_title, output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, neighbor_ids=neighbor_ids, cluster_names=cluster_names+col_names, n_repeats=args.n_repeats)
    elif 'neighbor_mat' in args.features:
        get_rsf_results_for_feature_df(input_df, combined_title, np.arange(0,num_features), output_title, output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, neighbor_ids=neighbor_ids, cluster_names=cluster_names, n_repeats=args.n_repeats)
    else: # TODO fix the cluster name issue
        get_rsf_results_for_feature_df(input_df, combined_title, np.arange(0,num_features), output_title, output_dir, args.cv_folds, permute_labels=args.permute_labels, precomputed=args.precomputed, cluster_names=cluster_names+col_names, n_repeats=args.n_repeats)



