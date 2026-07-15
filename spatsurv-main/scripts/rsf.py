import numpy as np
import pandas as pd
from scipy import stats

from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.inspection import permutation_importance
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored

from tqdm import tqdm
import os
import json
import pickle
import warnings

from utils import create_feature_names_neighborhood

def train_rsf(df, data_cols, title, label_cols=['survival_status', 'survival_day'], splits=10, cv_folds='patient', random_state=1029, permute_labels=False, output_dir='surv_forest_results', patient_col='patient_id', coverslip_col='coverslip_label', acquisition_col='acquisition_id', leave_progress_bar=True, n_repeats=0):
    '''Trains a random survival forest for specified cross-validation folds and computes the concordance index.

    Args:
        df (pd.DataFrame): dataframe with sample data, labels, and groups (patient + coverslip)
        data_cols (list[int]): List of indices indicating which columns within the dataframe are the sample features.
        title (str): name indicating which data features are being used to train the RSF. Used for file naming purposes.
        label_cols (list[str], optional): 2 strings indicating columns within the dataframe, first one indicating survival status, second survival time. Survival status is 1 when corresponding to a death, 0 otherwise. Default is ['survival_status', 'survival_day'].
        splits (int, optional): Number of cross validation splits to use. Default 10. If using leave one group out logic, this parameter is not applied.
        cv_folds (str, optional): Which cross-validation folds to use. Choices are 'patient', 'coverslip', 'random'. Default is 'patient'.
        random_state (int, optional): Random state to replicate behavior. Default 1029.
        permute_labels (bool, optional): When set to True, the survival labels will be shuffled. This is used to provide a baseline of 'random' performance. Default is False.
        output_dir (str, optional): name of the output directory to save results. Defaults to 'surv_forest_results'.
        patient_col (str, optional): Name of the column indicating patient id in the input dataframe. Default is 'patient_id'.
        coverslip_col (str, optional): Name of the column indicating coverslip id in the input dataframe. Default is 'coverslip_label'.
        acquisition_col (str, optional): Name of the column indicating acquisition id in the input dataframe. Default is 'acquisition_id'.
        leave_progress_bar (bool, optional): When True, the progress bar indicating progess through the cross-validation folds remain on the screen once finished. When False, the progress bar disappears from the screen once finished. 

    Returns:
        cv_concord (np.array, shape (`splits`,)): Array of concordance indices for each cross validation fold.
        feature_imp (np.array, shape (`splits`, `len(data_cols)`, 2)): Array that stores the mean and standard deviation of computed feature importances (via permutation importance) for each feature for each cross validation fold.
        pred_risks (np.array, shape (`len(df)`,)): Array that stores the predicted risk scores from the rsf model for all samples. These scores are those that are predicted when the corresponding sample is in the test set of a cross validation fold.
        split_info (np.array, shape (`splits`, `len(df)`)): Array indicating training and test set indices for each cross validation fold. Each row corresponds to a fold, and an element in that row is equal to 1 if the corresponding sample was present in the test set for that fold, 0 otherwise.
        
    '''
    if cv_folds not in ['patient', 'coverslip', 'random']:
        raise ValueError(f"`cv_folds` must be one of 'patient', 'coverslip', 'random', but got {cv_folds}")

    if not os.path.exists(os.path.join(output_dir, 'estimators')):
        os.makedirs(os.path.join(output_dir, 'estimators'))
    
    if n_repeats > 0:
        if not os.path.exists(os.path.join(output_dir, 'repeats')):
            os.makedirs(os.path.join(output_dir, 'repeats'))
        output_dir = os.path.join(output_dir, 'repeats')

    data = df.iloc[:,data_cols].values
    labels_df = df.loc[:,label_cols]
    if permute_labels:
        labels_df = labels_df.sample(frac=1, random_state=random_state//3)
    labels = np.array([tuple(x) for x in labels_df.values], dtype=[('cens', '?'), ('time', '<f8')])

    # Create cross validation folds based on cv_folds argument
    if cv_folds == 'patient':
        groups = df.loc[:, patient_col]
        kf = GroupKFold(n_splits=splits)
        fold_splits = list(kf.split(data, labels, groups))
    elif cv_folds == 'coverslip':
        groups = df.loc[:, coverslip_col]
        kf = LeaveOneGroupOut()
        fold_splits = list(kf.split(data, labels, groups))
        splits = kf.get_n_splits(groups=groups)
    else: # random folds
        kf = KFold(n_splits=splits, shuffle=True, random_state=random_state*3)
        fold_splits = list(kf.split(data, labels))

    if not permute_labels:
        if n_repeats == 0:
            param_dict = {
                    'run_name': output_dir,
                    'acquisition_ids': list(df.acquisition_id.values),
                    'splits': splits,
                    'cv_folds': cv_folds,
                    'random_state': random_state,
                    'permute_labels': permute_labels
                    }
            with open(os.path.join(output_dir, 'run_params.json'), 'w') as f:
                json.dump(param_dict, f)
        else:
            for repeat_i in range(n_repeats):
                param_dict = {
                        'run_name': output_dir,
                        'acquisition_ids': list(df.acquisition_id.values),
                        'splits': splits,
                        'cv_folds': cv_folds,
                        'random_state': random_state+repeat_i+1,
                        'permute_labels': permute_labels
                        }
                with open(os.path.join(output_dir, f'run_params_repeat{repeat_i:02d}.json'), 'w') as f:
                    json.dump(param_dict, f)


    if leave_progress_bar:
        print(f"Splitting by {cv_folds}")
        print(f"number of splits: {splits}")

    ### NORMAL CASE, NO REPEATS
    if n_repeats == 0:
        cv_concord = np.zeros(splits) # concordance index on test data for each split
        pred_risks = np.zeros(len(data)) # store predicted risk from rsf for each sample
        feature_imp = np.zeros((splits, len(data_cols), 2)) # store feature importances and standard deviation for each split
        split_info = np.zeros((splits, len(data)), dtype=int) # store which samples are testing samples for each fold
        for i, (train_index, test_index) in tqdm(enumerate(fold_splits), total=splits, desc="CV folds", leave=leave_progress_bar):
            split_info[i, test_index] = 1
            train_data = data[train_index]
            train_labels = labels[train_index]
            test_data = data[test_index]
            test_labels = labels[test_index]
            if permute_labels:
                censored_labels = np.array([x[0] for x in test_labels])
                if censored_labels.sum() == 0:
                    cv_concord[i] = np.nan
                    continue

            if os.path.exists(f'./{output_dir}/estimators/{title}_cv{i}_rsf.pkl'):
                with open(f'./{output_dir}/estimators/{title}_cv{i}_rsf.pkl', 'rb') as f:
                    rsf = pickle.load(f)
            else:
                rsf = RandomSurvivalForest(n_estimators=100, n_jobs=-1, random_state=random_state)
                rsf.fit(train_data, train_labels)
            cv_concord[i] = rsf.score(test_data, test_labels)
            pred_risks[test_index] = rsf.predict(test_data)

            if permute_labels: # don't compute permutation importances for random label runs
                continue

            # Save RSF estimators
            with open(f'./{output_dir}/estimators/{title}_cv{i}_rsf.pkl', 'wb') as f:
                pickle.dump(rsf, f)

            perm_imps = permutation_importance(rsf, test_data, test_labels, n_repeats=15, random_state=random_state*5//7)
            feature_imp[i,:,0] = perm_imps['importances_mean']
            feature_imp[i,:,1] = perm_imps['importances_std']

        if permute_labels:
            return cv_concord, feature_imp, pred_risks, split_info

        np.save(f'./{output_dir}/{title}_concordance.npy', cv_concord)
        np.save(f'./{output_dir}/{title}_feature_imp.npy', feature_imp)
        np.save(f'./{output_dir}/{title}_pred_risks.npy', pred_risks)
        np.save(f'./{output_dir}/{title}_split_info.npy', split_info)

    ### REPEAT CASE
    else:
        for repeat_i in tqdm(range(n_repeats), total=n_repeats, desc='diff random seeds'):
            cv_concord = np.zeros(splits) # concordance index on test data for each split
            pred_risks = np.zeros(len(data)) # store predicted risk from rsf for each sample
            feature_imp = np.zeros((splits, len(data_cols), 2)) # store feature importances and standard deviation for each split
            split_info = np.zeros((splits, len(data)), dtype=int) # store which samples are testing samples for each fold
            for i, (train_index, test_index) in tqdm(enumerate(fold_splits), total=splits, desc="CV folds", leave=leave_progress_bar):
                split_info[i, test_index] = 1
                train_data = data[train_index]
                train_labels = labels[train_index]
                test_data = data[test_index]
                test_labels = labels[test_index]

                rsf = RandomSurvivalForest(n_estimators=100, n_jobs=-1, random_state=random_state+repeat_i+1)
                rsf.fit(train_data, train_labels)
                cv_concord[i] = rsf.score(test_data, test_labels)
                pred_risks[test_index] = rsf.predict(test_data)

                perm_imps = permutation_importance(rsf, test_data, test_labels, n_repeats=15, random_state=random_state*5//7)
                feature_imp[i,:,0] = perm_imps['importances_mean']
                feature_imp[i,:,1] = perm_imps['importances_std']

            np.save(f'./{output_dir}/{title}_repeat{repeat_i:02d}_concordance.npy', cv_concord)
            np.save(f'./{output_dir}/{title}_repeat{repeat_i:02d}_feature_imp.npy', feature_imp)
            np.save(f'./{output_dir}/{title}_repeat{repeat_i:02d}_pred_risks.npy', pred_risks)
            np.save(f'./{output_dir}/{title}_repeat{repeat_i:02d}_split_info.npy', split_info)

    return cv_concord, feature_imp, pred_risks, split_info

def permute_labels_performance(df, title, num_permutations=10, random_state=1029, splits=10, output_dir='surv_forest_results', **kwargs):
    """Get performance of RSF if labels are permuted. This is used as a baseline to evaluate model perfrmance.
    Saves the concordance indices from each label permutation in an npy file. Shape is (num_permutations x num_splits). Function return value is the mean value of this matrix.
    
    Args:
        df (pd.DataFrame): dataframe with sample data, labels, and groups (patient + coverslip)
        num_permutations (int): number of times to permute labels
        random_state (int): random state for replicating experiments
        **kwargs: other arguments to pass into train_rsf function

    Returns:
        (float) mean concordance index over all permutations
    """
    all_concordances = np.zeros((num_permutations, splits))
    for i in tqdm(range(num_permutations), total=num_permutations, desc="Label permutations"):
        concordances, _, _, _ = train_rsf(df, title=title, splits=splits, random_state=random_state+4*i, permute_labels=True, leave_progress_bar=False, **kwargs)
        all_concordances[i,:] = concordances

    np.save(f'./{output_dir}/{title}_{num_permutations}permutations_concordance.npy', all_concordances)

    return np.nanmean(np.nanmean(all_concordances, axis=1))

def compute_top_features(feature_mat, feature_names, n=5, verbose=True):
    '''Averages the permutation importances for each feature over all cross validation folds
    

    Args:
        feature_mat (arraylike): Array with shape (num_splits, num_features, 2). Lists the feature importance mean and standard deviations as calcualted by permutation important for each feature for each split.
        feature_names (list[str]): List of names for each feature.
        n (int, optional): Number of top-ranked features to print.
        verbose (bool, optional): When True, prints the top `n` features to the screen. Default is True.
    '''
    feature_names = np.array(feature_names)
    avg_feature_imps = np.zeros((feature_mat.shape[1], feature_mat.shape[2]))
    avg_feature_imps[:,0] = np.mean(feature_mat, axis=0)[:,0]
    avg_feature_imps[:,1] = np.sqrt(np.mean(np.power(feature_mat[:,:,1], 2), axis=0))

    sorted_idx = np.argsort(avg_feature_imps[:,0])[::-1]
    if verbose:
        print(f'Top {n} features:')
        for i in range(n):
            idx = sorted_idx[i]
            print(f"{str(idx) + ': ' + str(feature_names[idx]):<15} {avg_feature_imps[idx,0]:.5f} +/- {avg_feature_imps[idx,1]:.5f}")
            #for j in range(feature_mat.shape[0]):
            #    print(f"{'crossval ' + str(j):<15} {feature_mat[j,idx,0]:.5f} +/- {feature_mat[j,idx,1]:.5f}")

    return avg_feature_imps[sorted_idx], feature_names[sorted_idx]

def get_rsf_results_for_feature_df(df, feature_type, data_cols, title, output_dir, cv_folds, label_cols=['survival_status', 'survival_day'], splits=10, permute_labels=False, precomputed=False, patient_col='patient_id', coverslip_col='coverslip_label', acquisition_col='acquisition_id', neighbor_ids=None, cluster_names=None, n_repeats=0):
    """Given the input dataframe created for RSF training, run the RSF analysis pipeline
    according to specified parameters. This includes:
        - Training the random survival forest and saving the results
        - Loading results if already precomputed
        - Training the RSF using permuted labels as a baseline result

    Args:
        df (pd.DataFrame): dataframe with sample data, labels, and groups (patient + coverslip)
        feature_type (str): string indicating which type of features is being analyzed here. Should be one of ['neighbor_mat', 'celltype_prop', 'biomarker_region', 'biomarker_cell'].
        data_cols (list[int]): List of indices indicating which columns within the dataframe are the sample features.
        title (str): name indicating which data features are being used to train the RSF. Used for file naming purposes.
        output_dir (str): name of the output directory to save results. 
        cv_folds (str): Which cross-validation folds to use. Choices are 'patient', 'coverslip', 'random'.
        label_cols (list[str], optional): 2 strings indicating columns within the dataframe, first one indicating survival status, second survival time. Survival status is 1 when corresponding to a death, 0 otherwise. Default is ['survival_status', 'survival_day'].
        splits (int, optional): Number of cross validation splits to use. Default 10. If using leave one group out logic, this parameter is not applied.
        permute_labels (bool, optional): When set to True, the survival labels will be shuffled. This is used to provide a baseline of 'random' performance. Default is False.
        patient_col (str, optional): Name of the column indicating patient id in the input dataframe. Default is 'patient_id'.
        coverslip_col (str, optional): Name of the column indicating coverslip id in the input dataframe. Default is 'coverslip_label'.
        acquisition_col (str, optional): Name of the column indicating acquisition id in the input dataframe. Default is 'acquisition_id'.
        neighbor_ids (list, optional): The numeric ids for each element in the neighborhood matrix. When feature_type == 'neighbor_mat', this is required to create the feature names. 
        cluster_names (list, optional): A list of feature names. When feature_type == 'neighbor_mat' or 'celltype_prop', these are the names of the celltype clusters. When feature_type == 'biomarker_region' or 'biomarker_cell', these are the names of the biomarkers.
        n_repeats (int, optional): number of random seeds to repeat. Default is 0.

    Returns:
        None
    """
    if permute_labels:
        num_permutations = 10
        if precomputed:
            print("Permuted label performance already computed")
            mean_permuted_concord = np.nanmean(np.load(f'{output_dir}/{title}_{num_permutations}permutations_concordance.npy'))
        else:
            print(f'Computing baseline performance by permuting survival labels')
            if cv_folds == 'coverslip':
                splits = len(np.unique(df.loc[:, coverslip_col].values))
            mean_permuted_concord = permute_labels_performance(df, title, num_permutations=num_permutations, splits=splits, output_dir=output_dir, data_cols=data_cols, label_cols=label_cols, cv_folds=cv_folds, patient_col=patient_col, coverslip_col=coverslip_col, acquisition_col=acquisition_col)
        
        print(f'permuted labels mean CI ({num_permutations} permutations): {mean_permuted_concord}')
    
    else:
        if precomputed:
            print("RSF already trained, loading pre-computed results")
            concordances = np.load(f'{output_dir}/{title}_concordance.npy')
            feature_imp = np.load(f'{output_dir}/{title}_feature_imp.npy')
            pred_risks = np.load(f'{output_dir}/{title}_pred_risks.npy')
            split_info = np.load(f'{output_dir}/{title}_split_info.npy')
        else:
            print("Training RSF")
            if n_repeats >0:
                leave_progress_bar = False
            else:
                leave_progress_bar = True
            concordances, feature_imp, pred_risks, split_info = train_rsf(df, data_cols, title, label_cols=label_cols, splits=splits, cv_folds=cv_folds, output_dir=output_dir, patient_col=patient_col, coverslip_col=coverslip_col, acquisition_col=acquisition_col, n_repeats=n_repeats, leave_progress_bar = leave_progress_bar) 
       
        if n_repeats==0:
            if feature_type == 'neighbor_mat':
                feature_names = create_feature_names_neighborhood(neighbor_ids, cluster_names)
            elif feature_type == 'celltype_prop' or feature_type == 'biomarker_region' or feature_type == 'biomarker_cell' or feature_type == 'ripley' or feature_type == 'leesL':
                feature_names = cluster_names
            elif feature_type == 'celltype_prop_neighbor_mat':
                neighbor_feature_names = create_feature_names_neighborhood(neighbor_ids, cluster_names)
                feature_names = cluster_names + neighbor_feature_names
            elif feature_type == 'neighbor_mat_ripley' or feature_type == 'neighbor_mat_markcorr':
                neighbor_feature_names = create_feature_names_neighborhood(neighbor_ids, cluster_names[:int(np.sqrt(len(neighbor_ids)))])
                feature_names = neighbor_feature_names + cluster_names[int(np.sqrt(len(neighbor_ids))):]
                print(cluster_names)
            elif feature_type == 'neighbor_mat_kcelltype':
                neighbor_feature_names = create_feature_names_neighborhood(neighbor_ids, cluster_names[:int(np.sqrt(len(neighbor_ids)))])
                feature_names = neighbor_feature_names + cluster_names[int(np.sqrt(len(neighbor_ids))):]
            else:
                feature_names = cluster_names # TODO fix this issue for arbitrary combinations
            compute_top_features(feature_imp, feature_names, n=10)
        
            print(f"Average concordance index over all folds: {concordances.mean()}")
            
            survival_labels = df.loc[:,['survival_status', 'survival_day']]

    return

def get_results_from_repeats(title, output_dir, confidence=0.95, num_repeats=100):
    """Get the mean and confidence intervals for the concordance indices over all repeated runs

    Args:
        title (str): name indicating which data features are being used to train the RSF. Used for file naming/searching purposes.
        output_dir (str): path of output directory where results are stored

    Returns:
        mean (float): mean concordance value across all repeated runs
        std (float): standard deviation across all runs
        ci (float): 95% confidence interval
    """

    concords = np.zeros(num_repeats) 
    for repeat_i in range(num_repeats):
        concord_arr = np.load(os.path.join(output_dir, f'{title}_repeat{repeat_i:02d}_concordance.npy'))
        concords[repeat_i] = np.mean(concord_arr)

    std = np.std(concords)
    mean = np.mean(concords)
    ci = stats.norm.interval(confidence, loc=mean, scale=std/np.sqrt(num_repeats))

    return mean, std, ci

