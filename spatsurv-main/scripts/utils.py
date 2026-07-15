import numpy as np
import pandas as pd

import os
import glob

def update_hpvstatus(df):
    """Update unknown HPV status based on primary site. HPV+ patients can only occur when the
    primary site is located in the Base of Tongue, Tonsil, or Oropharynx.
    Here, 0 = HPV-, 1 = HPV+, and 5 = unknown status.

    Args: 
        df (pd.DataFrame): dataframe with columns 'hpvstatus' and 'primarysite'

    Returns:
        int: value corresponding to updated HPV status
    """
    if pd.isnull(df['hpvstatus']):
        return df['hpvstatus']
    if int(df['hpvstatus']) == 5:
        if df['primarysite'] not in ['Base of Tongue', 'Tonsil', 'Oropharynx']:
            return 0
        else:
            return 5
    else:
        return df['hpvstatus']

def aggregate_avg_biomarkers(data_dir, sample_list, marker_names):
    """Aggregates data on average pixel intensity values across a region for a list of regions and returns a dataframe with all of this information.

    Args:
        data_dir (str): path to where the avg biomarker values exist (should be in data/experssion_biomarkers/overall directory)
        sample_list (list, str): list of acquisition_ids to get data from
        marker_names(list, str): list of marker names

    Returns:
        avg_markers_df (pd.Dataframe): dataframe with shape (num_regions x num markers), where each row gives the average biomarker intensity across an entire region for each biomarker
    """
    num_markers = len(marker_names)
    avg_markers = np.zeros((len(sample_list), num_markers))
    for idx, sample in enumerate(sample_list):
        # Read in dataframe that lists the summed and averaged pixel intensity values across the entire image for each channel
        df = pd.read_csv(os.path.join(data_dir, f'{sample}_cell_info.csv'), index_col=0)
        vals = df.loc['average', marker_names].values # Average pixel values across entire image
        avg_markers[idx, :] = np.arcsinh(vals) # Performing arcsinh normalization
    
    avg_markers_df = pd.DataFrame(avg_markers, columns=marker_names)
    avg_markers_df['acquisition_id'] = sample_list

    return avg_markers_df

def aggregate_DenVar_clusters(data_dir, sample_list, marker_names):
    """Aggregate DenVar cluster values from a list of samples and markers into a single dataframe.

    Args:
        data_dir (str): path to where the DenVar cluster values for each biomarker are
        sample_list (list, str): list of acquisition_ids to get data from
        marker_names(list, str): list of marker names

    Returns:
        DenVar_df (pd.DataFrame): dataframe with shape (num_regions x num_markers), where each row gives the DenVar cluster value for each biomarker
    """
    DenVar_data = np.zeros((len(sample_list), len(marker_names)))
    for i, biomarker in enumerate(marker_names):
        cluster_df = pd.read_csv(os.path.join(data_dir, f'{biomarker}_DenVar_clusters.csv'), index_col=0)
        cluster_df.index = cluster_df['sample'].values
        cluster_df = cluster_df.reindex(sample_list)
        DenVar_data[:,i] = cluster_df.cluster.values

    DenVar_df = pd.DataFrame(DenVar_data, columns=marker_names)
    DenVar_df['acquisition_id'] = sample_list

    return DenVar_df


def create_neighbor_mat(num_clusters, sample_list, k=10):
    """Aggregate celltype neighborhood matrices from a list of samples into a single dataframe.

    Args:
        num_clusters (int): number of 'celltypes' present in the data
        sample_list (list, str): list of samples/regions to get data from
        k (int, optional): which k value to use from neighborhood matrix creation

    Returns:
        neighbor_flatten_df (pd.DataFrame): dataframe with shape (num_samples x num_clusters**2), where each row gives the flattened neighborhood matrix for that region
    """
    num_samples = len(sample_list)
    flatten_dim = num_clusters**2
    neighbor_flatten = np.zeros((num_samples, flatten_dim))
    for i, sample in enumerate(sample_list):
        arrpath = glob.glob(f'../data/celltype_neighborhoods/k{k}/{sample}.npy')
        arr = np.load(arrpath[0])
        arr = arr / arr.sum(1)[:,np.newaxis] 
        arr = np.nan_to_num(arr, nan=0)
        neighbor_flatten[i, :] = arr.flatten()

    neighbor_flatten_df = pd.DataFrame(neighbor_flatten, columns = np.arange(0,flatten_dim,1))
    neighbor_flatten_df['acquisition_id'] = np.array(sample_list)

    return neighbor_flatten_df

def create_feature_names_neighborhood(mat_ids, celltype_names):
    """Creates names for each of the neighborhood matrix features based on the celltype labels.

    Args:
        mat_ids (list[int] or arraylike): array or list of the matrix indices after flattening.
        celltype_names (list[str]): list of celltype labels in the same order as the columns of the matrix.

    Returns:
        all_feature_names (list[str]): list of feature names for each index in `mat_ids`.
    """
    all_feature_names = [None]*len(mat_ids)
    for i, mat_id in enumerate(mat_ids):
        num_celltypes = len(celltype_names)
        row_celltype = celltype_names[int(mat_id // num_celltypes)]
        col_celltype = celltype_names[int(mat_id % num_celltypes)]
        all_feature_names[i] = f'{row_celltype}>{col_celltype}'
    return all_feature_names

def create_input_dataframe_for_features(feature_set, sample_df, label_cols, expr=None, sample_list=None, marker_names=None, num_clusters=None, k=10):
    """Creates the input dataframe for prediction tasks.

    Args: 
        feature_set (str): Specifies which feature set to use. Must be one of ["biomarker_region", "biomarker_cell", "celltype_prop", "neighbor_mat"].
        sample_df (pd.DataFrame): dataframe with metadata for each sample.
        label_cols (list[str]): Columns with the prediction label information. For example, for a survival prediction task, the label columns would be ['survival_day', 'survival_status']
        expr (pd.DataFrame, optional): dataframe with the biomarker expression values from each cell. This is required for feature_set == 'biomarker_cell', 'celltype_prop', and 'neighbor_mat'
        sample_list (list[str], optional): list of samples/regions to use. This is required for feature_set == 'biomarker_region' and 'neighbor_mat'
        marker_names (list[str], optional): list of biomarker names. Required for feature_set == 'biomarker_region' and 'biomarker_cell'
        num_clusters (int, optional): number of celltype clusters. Required for feature_set == 'celltype_prop' and 'neighbor_mat'
        k (int, optional): number of neighbors for neighbor_mat computation.

    Returns:
        input_df (pd.DataFrame): dataframe with shape (num_samples, num_features + label_cols + ['patient_id, 'coverslip_label', 'acquisition_id'])
    """
    if feature_set not in ["biomarker_region", "biomarker_cell", "celltype_prop", "neighbor_mat"]:
        raise ValueError(f"{feature_set} is not a valid feature set. Must be one of 'biomarker_region', 'biomarker_cell', 'celltype_prop', 'neighbor_mat'.")

    if feature_set == "biomarker_region":
        avg_markers = aggregate_avg_biomarkers('../data/expression_biomarkers/', sample_list, marker_names)
        avg_markers = sample_df.merge(avg_markers, how='inner', on='acquisition_id')
        input_df = avg_markers[marker_names + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]

    elif feature_set == "biomarker_cell":
        cell_biomarkers = expr.copy()
        del cell_biomarkers['cell_id']
        del cell_biomarkers['cluster']
        cell_biomarkers = cell_biomarkers.rename({'sample_id': 'acquisition_id'}, axis='columns')
        cell_biomarkers = cell_biomarkers.groupby(by=['acquisition_id']).mean()
        cell_biomarkers = sample_df.merge(cell_biomarkers, how='inner', on='acquisition_id')

        input_df = cell_biomarkers[marker_names + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]

    else: # celltype_prop or neighbor_mat
        celltype_prop = expr.loc[:,['sample_id','cluster']]
        celltype_prop = celltype_prop.rename({'sample_id': 'acquisition_id'}, axis='columns')
        celltype_prop = celltype_prop.groupby(by=['acquisition_id','cluster']).size().unstack(fill_value=0)
        celltype_prop = celltype_prop.div(celltype_prop.sum(axis=1), axis=0)
        celltype_prop = sample_df.merge(celltype_prop, how='inner', on='acquisition_id')
        if feature_set == "celltype_prop":
            celltype_ids = list(np.arange(num_clusters))
            input_df = celltype_prop[celltype_ids + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]
        else: 
            flatten_dim = num_clusters**2
            neighbor_flatten = create_neighbor_mat(num_clusters, sample_list, k=k)
            neighbor_flatten = sample_df.merge(neighbor_flatten, how='inner', on='acquisition_id')
            
            neighbor_ids = list(np.arange(flatten_dim))
            input_df = neighbor_flatten[neighbor_ids + label_cols + ['patient_id', 'coverslip_label', 'acquisition_id']]

    return input_df


