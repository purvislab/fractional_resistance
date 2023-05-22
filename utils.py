import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import scoreatpercentile
from sklearn.preprocessing import StandardScaler
from scipy.stats import norm
import seaborn as sns
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, confusion_matrix, roc_curve, r2_score
from transact.TRANSACT import TRANSACT
import phate 
import anndata
from sklearn.mixture import GaussianMixture

def standardize(x =  None):
    """Standardizes data by removing the mean and scaling to unit variance.

    Parameters
    x: pd.DataFrame (default = None)
        data matrix (dimensions = cells x features)
    ----------

    Returns
    X: pd.DataFrame
        standardized data matrix (dimensions = cells x features)
    ----------
    """
    scaler = StandardScaler(with_mean = True, with_std = True)
    X = pd.DataFrame(scaler.fit_transform(x), index = x.index, columns = x.columns)
    return X

def pRB_cutoff(x_vec = None,
            pct = 50,
            filename_save = None):
    """Finds the cuttoff between arrested from cycling cells using the distribution of pRB expression
    Parameters
    x_vec: np.ndarray (default = None)
        array containing pRB expression values
    pct: int (default = 50)
        percentile for the one-sided distribution
    filename_save: str (default = None)
        str specifying filename for saving
    ----------

    Returns
    cut: float 
        cutpoint from the pRB/RB distribution separating pRB high (proliferative) from pRB low (arrested) cells
    ----------
    """
    med = np.median(x_vec)
    x_one_sided = x_vec.copy()              
    # fold over to the right all values left of the median
    x_one_sided[x_vec < med] = 2*med - x_vec[x_vec < med]

    # find the percentile of the one-sided distribution
    a = scoreatpercentile(x_one_sided, pct) - med

    #check how many elements fall into the range
    sel = (x_vec > (med - a)) & (x_vec < (med + a))

    cut = (np.sum(sel) / float(len(x_vec))) #cutpoint

    #look at histograms
    pd.DataFrame(x_vec).plot(kind='hist', bins=50, range=[-2, 5], legend=None) #un normalized data
    plt.tick_params(labelsize=14)
    plt.axvline(x=cut, ymin=0, ymax=3, color='red')
    plt.xlabel('standardized pRB/RB expression', fontsize=16)
    plt.ylabel('frequency', fontsize=16)
    if filename_save is not None:
        plt.savefig(filename_save+'.pdf', bbox_inches = 'tight')

    return cut

def phase_annotations(data = None,
                    n_components: int = 3,
                    features = ['DNA', 'cycA', 'cycB1'],
                    random_state: int = 0):
    """Performs clustering using a Gaussian  Mixture Model: https://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html
        Parameters
        data: pd.DataFrame
            dataframe containing single-cell data (dimensions = cells x features)
        n_components: int (default = 3)
            number of clusters for the GMM
        features: list (default = ['DNA', 'cycA', 'cycB1'])
            list of features to perform GMM clustering on for phase assignment
        random_state: int (default = 0)
            random seed
        ----------

        Returns
        labels: np.array
            array containing cell-cluster assignment
        ----------
    """
    #subset data according to features of interest
    cluster_df = np.log(data[features])
    
    # train gaussian mixture model 
    gmm = GaussianMixture(n_components = n_components, covariance_type = 'tied', n_init = 10, random_state = random_state)
    gmm.fit(cluster_df)
    labels = gmm.predict(cluster_df)
    frame = pd.DataFrame(cluster_df)
    frame['cluster'] = labels

    #plot for assignment
    color=['blue','green','cyan', 'black']

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    for k in range(0, len(color)):
        cluster_df = frame[frame["cluster"] == k]
        ax.scatter(cluster_df["DNA"], cluster_df["cycA"], cluster_df['cycB1'],
                    c = color[k], alpha=0.1)

    ax.set_xlabel(features[0])
    ax.set_ylabel(features[1])
    ax.set_zlabel(features[2])
    plt.show()

    return labels

def preprocess_t47(directory: str = None,
                    filename: str = None,
                    percentile: float = 20,
                    n_components: int = 3,
                    prolif_features = ['DNA', 'cycA', 'cycB1'],
                    random_state: int = 0):
    """Preprocesses T47D data by:
            1. subsetting feature set (uses the intersection of features profiled across tumor and T47D samples)
            2. standardization
            3. phase annotations
                - computes threshold between cycling (pRB high) and arrested cells (pRB low)
                - uses GMM on proliferative feature subset to assign G1, S, G2, M cells

        Parameters
        directory: str (default = None)
            string specifying path of data
        filename: str or list (default = None)
            list of csv file names
        percentile: float (default = 20)
            float specifying the percentile cutoff for the pRB/RB distribution
        n_components: int (default = 3)
            number of clusters for the GMM
        prolif_features: list (default = ['DNA', 'cycA', 'cycB1'])
            list of features to perform GMM clustering on for phase assignment
        random_state: int (default = 0)
            integer referring to the seed for GMM
        ----------

        Returns
        adata: anndata.AnnData
            annotated data object containing preprocessed 4i data (adata.X, -> dimensions: cells x features) and metadata (adata.obs)
        cut: float
            cutpoint from the pRB/RB distribution separating pRB high (proliferative) from pRB low (arrested) cells
        ----------
    """
    feature_names = ['Intensity_IntegratedIntensity_DNA_nuc',
                    '00_CDK2_nuc_median',
                    '05_CDK4_nuc_median',
                    '04_Cdt1_nuc_median',
                    '02_E2F1_nuc_median',
                    '02_Ki67_nuc_median',
                    '00_Rb_nuc_median',
                    '05_cycA2_nuc_median',
                    '01_cycB1_nuc_median',
                    '01_cycD1_nuc_median',
                    '04_cycE_nuc_median',
                    '00_pRB_nuc_median',
                    '02_p21_nuc_median',
                    'pRB_over_RB',
                    'well']

    feature_dict = {'Intensity_IntegratedIntensity_DNA_nuc': 'DNA', '00_CDK2_nuc_median': 'CDK2', '05_CDK4_nuc_median':'CDK4',
                    '04_Cdt1_nuc_median': 'Cdt1', '02_E2F1_nuc_median':'E2F1', '02_Ki67_nuc_median': 'Ki67', '00_Rb_nuc_median':'RB',
                    '05_cycA2_nuc_median':'cycA', '01_cycB1_nuc_median':'cycB1', '01_cycD1_nuc_median':'cycD1', '04_cycE_nuc_median':'cycE',
                    '02_p21_nuc_median':'p21', '00_pRB_nuc_median':'pRB', 'pRB_over_RB':'pRB_over_RB','well':'well'}

    #read in data
    data = pd.DataFrame()
    for i in filename:
        data_ = pd.read_csv(os.path.join(directory, i + '.csv'), index_col = 0)
        data = pd.concat([data, data_], axis = 0)

    data.index = ['t47_' + str(i) for i in range(0, data.shape[0])]
    data['Intensity_IntegratedIntensity_DNA_nuc'] = data['area']*data['00_DNA_nuc_median']
    data['pRB_over_RB'] = data['00_pRB_nuc_median'] / data['00_Rb_nuc_median']

    data = data.loc[:, feature_names]
    data = data.rename(columns = feature_dict)
    di = {b'B02':'0', b'B03':'10', b'B05':'100'}
    data['well'] = data['well'].astype('|S').replace(di)

    data = data[(data['well'] == '0') | (data['well'] == '10') | (data['well'] == '100')]

    metadata = data.loc[:, ['well']]
    data.drop(columns = ['well'], inplace = True)

    #standardize 
    data_std = standardize(data)

    #compute cutoff for proliferative vs arrested cells
    cut = pRB_cutoff(np.asarray(data_std.loc[:, 'pRB_over_RB']), percentile, 't47_pRB_RB')

    #create preprocessed annotated data object
    adata = anndata.AnnData(data_std)
    adata.obs = metadata.copy()
    adata.obs_names = adata.obs_names.astype('str')

    adata.obs['prb_ratio'] = 'nan'
    adata.obs['prb_ratio'].values[data_std['pRB_over_RB'] >= cut] = '1'
    adata.obs['prb_ratio'].values[data_std['pRB_over_RB'] < cut] = '0'

    adata.obs['RB_status'] = 'nan'
    adata.obs['RB_status'].values[data_std['pRB_over_RB'] >= cut] = 'cycling'
    adata.obs['RB_status'].values[data_std['pRB_over_RB'] < cut] = 'arrested'

    #perform GMM clustering and annotate phases
    labels = phase_annotations(data, n_components = n_components, features = prolif_features, random_state = random_state)

    adata.obs['phase'] = np.nan
    adata.obs['phase'][labels == 0] = 'G1'
    adata.obs['phase'][labels == 1] =  'S'
    adata.obs['phase'][labels == 2] = 'G2M'
    adata.obs['phase'][adata.obs['RB_status'] == 'arrested'] = 'G0'

    adata.obs['phase'] = pd.Categorical(adata.obs['phase'])
    return adata, cut

def preprocess_tumor(directory: str = None,
                    filename: str = None,
                    percentile: float = 70,
                    n_components: int = 3,
                    prolif_features = ['DNA', 'cycA', 'cycB1'],
                    random_state: int = 0):
    """Preprocesses tumor data by:
            1. subsetting feature set (uses the intersection of features profiled across tumor and T47D samples)
            2. standardization
            3. phase annotations
                - computes threshold between cycling (pRB high) and arrested cells (pRB low)
                - uses GMM on proliferative feature subset to assign G1, S, G2, M cells

        Parameters
        directory: str (default = None)
            string specifying path of data
        filename: str or list (default = None)
            list of csv file names
        percentile: float (default = 20)
            float specifying the percentile cutoff for the pRB/RB distribution
        n_components: int (default = 3)
            number of clusters for the GMM
        prolif_features: list (default = ['DNA', 'cycA', 'cycB1'])
            list of features to perform GMM clustering on for phase assignment
        random_state: int (default = 0)
            integer referring to the seed for GMM
        ----------

        Returns
        adata: anndata.AnnData
            annotated data object containing preprocessed 4i data (adata.X, -> dimensions: cells x features) and metadata (adata.obs)
        cut: float
            cutpoint from the pRB/RB distribution separating pRB high (proliferative) from pRB low (arrested) cells
        ----------
    """
    feature_names = ['Intensity_IntegratedIntensity_DNA_nuc',
                    'Intensity_MedianIntensity_CDK2_0_nuc',
                    'Intensity_MedianIntensity_CDK4_nuc',
                    'Intensity_MedianIntensity_Cdt1_nuc',
                    'Intensity_MedianIntensity_E2F1_nuc',
                    'Intensity_MedianIntensity_Ki67_nuc',
                    'Intensity_MedianIntensity_Rb_0_nuc',
                    'Intensity_MedianIntensity_cycA_nuc',
                    'Intensity_MedianIntensity_cycB1_nuc',
                    'Intensity_MedianIntensity_cycD_nuc',
                    'Intensity_MedianIntensity_cycE_nuc',
                    'Intensity_MedianIntensity_p21_nuc',
                    'Intensity_MedianIntensity_pRb_0_nuc',
                    'pRB_over_RB',
                    'palbociclib']

    feature_dict = {'Intensity_IntegratedIntensity_DNA_nuc': 'DNA', 'Intensity_MedianIntensity_CDK2_0_nuc': 'CDK2', 'Intensity_MedianIntensity_CDK4_nuc':'CDK4',
                    'Intensity_MedianIntensity_Cdt1_nuc': 'Cdt1', 'Intensity_MedianIntensity_E2F1_nuc':'E2F1',
                    'Intensity_MedianIntensity_Ki67_nuc': 'Ki67', 'Intensity_MedianIntensity_Rb_0_nuc':'RB', 'Intensity_MedianIntensity_cycA_nuc':'cycA',
                    'Intensity_MedianIntensity_cycB1_nuc':'cycB1', 'Intensity_MedianIntensity_cycD_nuc':'cycD1', 'Intensity_MedianIntensity_cycE_nuc':'cycE',
                    'Intensity_MedianIntensity_p21_nuc':'p21', 'Intensity_MedianIntensity_pRb_0_nuc':'pRB', 'pRB_over_RB':'pRB_over_RB', 'palbociclib':'well'}

    #read in data
    data = pd.read_csv(os.path.join(directory, filename+'.csv'))
    data.index = ['tumor_' + str(i) for i in range(0, data.shape[0])]

    data['pRB_over_RB'] = data['Intensity_MedianIntensity_pRb_0_nuc'] / data['Intensity_MedianIntensity_Rb_0_nuc']
    data['Intensity_IntegratedIntensity_DNA_nuc'] = data['AreaShape_Area_nuc']*data['Intensity_MedianIntensity_DNA_00_nuc']

    data = data[(data['Intensity_MedianIntensity_ER_nuc']) >= (data['Intensity_MedianIntensity_ER_nuc']).quantile(0.50)]
    data = data[(data['Intensity_MedianIntensity_PR_nuc']) >= (data['Intensity_MedianIntensity_PR_nuc']).quantile(0.50)]
    data = data.loc[data['tamoxifen'] == 0]

    data = data.loc[:, feature_names]
    data = data.rename(columns = feature_dict)
    data['well'] = data['well'].astype('str')

    data = data[(data['well'] == '0') | (data['well'] == '10') | (data['well'] == '100')]

    metadata = data.loc[:, 'well']
    data.drop(columns  = ['well'], inplace = True)

    #standardize
    data_std = standardize(data)

    #compute cutoff for proliferative vs arrested cells
    cut = pRB_cutoff(np.asarray(data_std.loc[:, 'pRB_over_RB']), percentile, 'tumor_pRB_RB')

    #create preprocessed annotated data object
    adata = anndata.AnnData(data_std)
    adata.obs = pd.DataFrame(metadata).copy()
    adata.obs_names = adata.obs_names.astype('str')

    adata.obs['prb_ratio'] = 'nan'
    adata.obs['prb_ratio'].values[data_std['pRB_over_RB'] >= cut] = '1'
    adata.obs['prb_ratio'].values[data_std['pRB_over_RB'] < cut] = '0'

    adata.obs['RB_status'] = 'nan'
    adata.obs['RB_status'].values[data_std['pRB_over_RB'] >= cut] = 'cycling'
    adata.obs['RB_status'].values[data_std['pRB_over_RB'] < cut] = 'arrested'

    #perform GMM clustering and annotate phases
    labels = phase_annotations(data, n_components = n_components, features = prolif_features, random_state = random_state)

    adata.obs['phase'] = np.nan
    adata.obs['phase'][labels == 0] = 'S'
    adata.obs['phase'][labels == 2] =  'G2M'
    adata.obs['phase'][labels == 1] = 'G1'
    adata.obs['phase'][adata.obs['RB_status'] == 'arrested'] = 'G0'

    adata.obs['phase'] = pd.Categorical(adata.obs['phase'])
    return adata, cut

def preprocess_t47d_rep(directory: str = None,
                    filename: str = None,
                    percentile: float = 70):
    """Preprocesses T47D replicate data by:
            1. subsetting feature set (uses the intersection of features profiled across tumor and T47D samples)
            2. standardization
            3. computes threshold between cycling (pRB high) and arrested cells (pRB low)

        Parameters
        directory: str (default = None)
            string specifying path of data
        filename: str or list (default = None)
            list of csv file names
        percentile: float (default = 20)
            float specifying the percentile cutoff for the pRB/RB distribution
        ----------

        Returns
        adata: anndata.AnnData
            annotated data object containing preprocessed 4i data (adata.X, -> dimensions: cells x features) and metadata (adata.obs)
        cut: float
            cutpoint from the pRB/RB distribution separating pRB high (proliferative) from pRB low (arrested) cells
        ----------
    """
    feature_names = ['Intensity_IntegratedIntensity_DNA',
                    'Intensity_MedianIntensity_CDK2',
                    'Intensity_MedianIntensity_CDK4',
                    'Intensity_MedianIntensity_CDK6',
                    'Intensity_MedianIntensity_Cdt1',
                    'Intensity_MedianIntensity_E2F1',
                    'Intensity_MedianIntensity_Ki67',
                    'Intensity_MedianIntensity_RB',
                    'Intensity_MedianIntensity_cycA',
                    'Intensity_MedianIntensity_cycB1',
                    'Intensity_MedianIntensity_cycD1',
                    'Intensity_MedianIntensity_cycE',
                    'Intensity_MedianIntensity_p21',
                    'Intensity_MedianIntensity_pRB',
                    'pRB_over_RB',
                    'Metadata_well',
                    'ploidy',
                    'RB_status',
                    'phase']

    feature_dict = {'Intensity_IntegratedIntensity_DNA': 'DNA', 'Intensity_MedianIntensity_CDK2': 'CDK2', 'Intensity_MedianIntensity_CDK4':'CDK4',
                    'Intensity_MedianIntensity_CDK6': 'CDK6', 'Intensity_MedianIntensity_Cdt1': 'Cdt1', 'Intensity_MedianIntensity_E2F1':'E2F1',
                    'Intensity_MedianIntensity_Ki67': 'Ki67', 'Intensity_MedianIntensity_RB':'RB', 'Intensity_MedianIntensity_cycA':'cycA',
                    'Intensity_MedianIntensity_cycB1':'cycB1', 'Intensity_MedianIntensity_cycD1':'cycD1', 'Intensity_MedianIntensity_cycE':'cycE',
                    'Intensity_MedianIntensity_p21':'p21', 'Intensity_MedianIntensity_pRB':'pRB', 'pRB_over_RB':'pRB_over_RB',
                    'Metadata_well':'well', 'ploidy':'ploidy', 'RB_status':'RB_status', 'phase':'phase'}

    #read in data
    data = pd.read_csv(os.path.join(directory, filename + '.csv'), index_col = 0)
    data.index = ['t47d_rep_' + str(i) for i in range(0, data.shape[0])]

    data = data.loc[:, feature_names]
    data = data.rename(columns = feature_dict)
    di = {b'B05':'0', b'C05':'1', b'D05':'10', b'E05':'100', b'F05':'1000', b'G05':'10000'}
    data['well'] = data['well'].astype('|S').replace(di)

    data = data[data['phase'].notna()]
    data = data[(data['well'] == '0') | (data['well'] == '10') | (data['well'] == '100')]

    metadata = data.loc[:, ['well', 'ploidy', 'RB_status', 'phase']]
    data.drop(columns = ['well', 'ploidy', 'RB_status', 'phase'], inplace = True)

    #standardize
    data_std = standardize(data)

    #compute cutoff for proliferative vs arrested cells
    cut = pRB_cutoff(np.asarray(data_std.loc[:, 'pRB_over_RB']), percentile, 't47_rep_pRB_RB')

    #create preprocessed annotated data object
    adata = anndata.AnnData(data_std)
    adata.obs = metadata.copy()
    adata.obs_names = adata.obs_names.astype('str')

    adata.obs['phase'] = pd.Categorical(adata.obs['phase'])

    adata.obs['prb_ratio'] = 'nan'
    adata.obs['prb_ratio'].values[data_std['pRB_over_RB'] >= cut] = '1'
    adata.obs['prb_ratio'].values[data_std['pRB_over_RB'] < cut] = '0'

    return adata, cut

def transact_integrate(df1 = None,
                        df2 = None,
                        kernel: str = 'rbf',
                        n_components: int = 14,
                        n_jobs: int = -1):
    """Performs integration with TRANSACT: https://www.pnas.org/doi/10.1073/pnas.2106682118
    Parameters
    df1: pd.DataFrame
        dataframe containing normalized and preprocessed data from first sample (dimensions = cells x features)
    df1: pd.DataFrame
        dataframe containing normalized and preprocessed data from second sample (dimensions = cells x features)
    kernel: str (default = 'rbf')
        kernel function for computing similarity
    n_components: int (default = 15)
        number of principal vectors
    n_jobs: int (default = -1)
        number of tasks
    ----------
    Returns
    integrated_df: pd.DataFrame
        dataframe containing merged data following integration
    ----------
    """
    clf = TRANSACT(kernel = kernel, kernel_params = {'gamma':1/np.sqrt(500)}, n_components = n_components, n_jobs = n_jobs)
    clf.fit(df1, df2, step = 100, with_interpolation = True)

    trans_df1 = pd.DataFrame(clf.transform(df1), index = df1.index)
    trans_df2 = pd.DataFrame(clf.transform(df2), index = df2.index)

    integrated_df = pd.concat([trans_df1, trans_df2])

    return integrated_df

def compute_phate(df = None,
                n_components = 2,
                knn = 100,
                random_state = 0,
                n_jobs = -1):
    """Performs nonlinear dimensionality reduction with PHATE: https://pubmed.ncbi.nlm.nih.gov/31796933/
    Parameters
    df: pd.DataFrame (default = None)
        dataframe to perform nonlinear dimensionality reduction
    n_components: int (default = 5)
        number of components for MDS
    knn: int (default = 100)
        number of nearest neighbors
    n_jobs: int (default = -1)
        number of tasks
    ----------
    Returns
    X_phate: np.ndarray
        PHATE embedding (dimensions = cells x n_components)
    ----------
    """
    phate_op = phate.PHATE(n_components = n_components, knn = knn, n_jobs = n_jobs, random_state = random_state) 
    X_phate = phate_op.fit_transform(df)
    return X_phate

def subset_cells(adata = None,
                labels_key: str = 'phase',
                label_id: str = None):
    """Subsets cells according to cell type (i.e. label_id in adata.obs[labels_key])

    Parameters
    adata: anndata.AnnData (default = None)
        annotated preprocessed data matrix (dimensions = cells x features)
    labels_key: str (default = 'phase)
        string in adata.obs specifying cell type labels
    label_id: str (default = None)
        string specifying 
    ----------

    Returns
    adata_subset: anndata.AnnData
        annotated data matrix containing cells from cell type of interest (dimensions = cells x features)
    ----------
    """
    adata_subset = adata[np.isin(adata.obs[labels_key], label_id), :].copy()
    return adata_subset

def logit_pvalue(model, x):
    """ Calculate z-scores for scikit-learn LogisticRegression.
    parameters:
        model: fitted sklearn.linear_model.LogisticRegression with intercept and large C
        x:     matrix on which the model was fit
    This function uses asymtptics for maximum likelihood estimates.
    """
    p = model.predict_proba(x)
    n = len(p)
    m = len(model.coef_[0]) + 1
    coefs = np.concatenate([model.intercept_, model.coef_[0]])
    x_full = np.matrix(np.insert(np.array(x), 0, 1, axis = 1))
    ans = np.zeros((m, m))
    for i in range(n):
        ans = ans + np.dot(np.transpose(x_full[i, :]), x_full[i, :]) * p[i,1] * p[i, 0]
    vcov = np.linalg.inv(np.matrix(ans))
    se = np.sqrt(np.diag(vcov))
    t =  coefs/se  
    p = (1 - norm.cdf(abs(t))) * 2
    return p

def logistic_reg(X = None,
                y = None,
                feature_names = None,
                n_splits: int = 5,
                random_state: int = 0):
    """Performs logistic regression: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html
        Parameters
        X: np.array
            array containing preprocessed data (dimensions = cells x features)
        y: np.ndarray
            array containing labels for classification
        feature_names: list (default = None)
            list containing feature names (i.e. column names of X)
        n_splits: int (default = 5)
            number of cross validation folds
        random_state: int (default = 0)
            random seed parameter
        ----------

        Returns
        prediction_coef: pd.DataFrame
            dataframe containing coefficients from logistic regression model across folds
        scores: pd.DataFrame
            dataframe containing AUC scores across folds
        tpr_fpr: pd.DataFrame
            dataframe containing fpr and tpr across folds
        pvals: pd.DataFrame
            dataframe containing p-values across folds
        ----------
    """
    le = LabelEncoder()
    y = le.fit_transform(y)
    n_classes = len(np.unique(y))
    y_bin = label_binarize(y, classes=list(np.arange(0, n_classes)))

    scores = pd.DataFrame()
    prediction_coef = pd.DataFrame()
    tpr_fpr = pd.DataFrame()
    pvals = pd.DataFrame()
    
    param_grid = {'C': [0.001, 0.01, 0.1, 1, 10, 100, 1000], 'l1_ratio': np.arange(0, 1, 0.01)}
    cv = StratifiedKFold(n_splits = n_splits, shuffle = True, random_state = random_state)
    for train_ix, test_ix in cv.split(X, y):
        X_train, X_test = X[train_ix, :], X[test_ix, :]
        y_train, y_test = y[train_ix], y[test_ix]

        cv_hyperparam = StratifiedKFold(n_splits = 3, shuffle = True, random_state = random_state)
        model = LogisticRegression(penalty = 'elasticnet', solver = 'saga', max_iter = 500)
        gsearch = GridSearchCV(model, param_grid, scoring = 'accuracy', n_jobs = -1, cv = cv_hyperparam, refit = True)
    
        result = gsearch.fit(X_train, y_train)

        best_model = result.best_estimator_

        y_prob = best_model.predict_proba(X_test)

        if n_classes == 2:
            auc = roc_auc_score(y_test, y_prob[:, 1], multi_class = 'ovr') #y_prob is already on test
            fpr, tpr, _ = roc_curve(y_test,  y_prob[::,1])
            tpr_fpr_ = pd.DataFrame({'fpr': fpr,
                                    'tpr': tpr})

            tpr_fpr = pd.concat([tpr_fpr, tpr_fpr_], axis = 1)
        else:
            auc = roc_auc_score(y_bin[test_ix, :], y_prob, multi_class = 'ovr')
            tpr_fpr = None
            
        scores_ = pd.DataFrame({'auc': [auc]})
        scores = pd.concat([scores_, scores], axis = 0)

        prediction_coef_ = pd.DataFrame(best_model.coef_.transpose(), index = feature_names)
        prediction_coef = pd.concat([prediction_coef_, prediction_coef], axis= 1)

        pvals_ = pd.DataFrame(logit_pvalue(best_model, X_test), index = ['intercept'] + list(feature_names))
        pvals = pd.concat([pvals_, pvals], axis = 1)

    return prediction_coef, scores, tpr_fpr, pvals

def plot_roc_auc(scores = None,
                tpr_fpr = None,
                n_splits: int = 5,
                filename_save: str = None):
    """Plots roc auc curve 
        Parameters
        scores: pd.DataFrame
            dataframe containing auc scores across folds
        tpr_fpr: pd.DataFrame
            dataframe containing fpr tpr across folds
        n_splits: int (default = 5)
            number of cross validation folds
        filename_save: str (default = None)
            string referring to the filename for saving
        ----------

        Returns
        ----------
    """     
    _, axes = plt.subplots(1,1, figsize = (6, 5.5), gridspec_kw={'bottom':0.15})

    mean_auc = np.asarray(scores.mean())
    std_auc = np.asarray(scores.std())

    for i in range(0, n_splits):
        sns.lineplot(x = tpr_fpr['fpr'].iloc[:, i], y = tpr_fpr['tpr'].iloc[:, i], label='ROC fold %d (AUC = %0.2f)' % (i+1, scores.iloc[i]), lw=0.8, ax = axes)

    sns.lineplot(x = tpr_fpr['fpr'].mean(1), y = tpr_fpr['tpr'].mean(1),
            label = r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2.5, ax = axes, color = 'blue')

    axes.plot([0, 1], [0, 1], linestyle='--', color='black',
            label='chance', alpha=.8, lw=2.5)

    axes.legend(loc = 'lower right',prop={'size':12})
    axes.tick_params(labelsize=16)

    axes.set_xlim([0, 1])
    axes.set_ylim([0, 1])

    axes.set_ylabel('True Positive Rate', fontsize = 16)
    axes.set_xlabel('False Positive Rate', fontsize = 16)
    plt.savefig(filename_save+'.pdf', bbox_inches="tight")

def add_value_labels(label, ax, spacing=1):
    """Add labels to the end of each bar in a bar chart.

    Arguments:
        ax (matplotlib.axes.Axes): The matplotlib object containing the axes
            of the plot to annotate.
        spacing (int): The distance between the labels and the bars.
    """

    # For each bar: Place a label
    for rect, l in zip(ax.patches, label):
        # Get X and Y placement of label from rect.
        y_value = rect.get_height()
        x_value = rect.get_x() + rect.get_width() / 2

        # Number of points between bar and label. Change to your liking.
        space = spacing
        # Vertical alignment for positive values
        va = 'bottom'

        # If value of bar is negative: Place label below bar
        if y_value < 0:
            # Invert space to place label below
            space *= -1
            # Vertically align label at top
            va = 'top'

        # Use Y value as label and format number with one decimal place
        #label = "{:.1f}".format(label)

        # Create annotation
        ax.annotate(
            l,                      # Use `label` as label
            (x_value, y_value),         # Place label at end of the bar
            xytext=(0, space),          # Vertically shift label by `space`
            textcoords="offset points", # Interpret `xytext` as offset in points
            ha='center',                # Horizontally center label
            va=va,
            fontsize = 13)                      # Vertically align label differently for
                                        # positive and negative values.

def run_logistic_regression(adata = None, 
                            groups = ['G0', 'G1', 'S', 'G2M'],
                            origin: str = None,
                            condition_key: str = 'condition',
                            labels_key: str = 'phase',
                            n_splits: int = 5,
                            save_directory: str = 'figures',
                            ylim = [-10.5, 2.5]):
    """Runs logistic regression and plots coefficients, AUCs: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html
        Parameters
        adata: anndata
            annotated data object containing single-cell data (dimensions = cells x features)
        groups: list (default = ['G0', 'G1', 'S', 'G2M'])
            list of cell types to fit logistic regression  model on
        origin: str (default = None)
            string referring to the data id in adata.obs['Origin']
        condition_key: str (default = 'condition')
            string referring to the condition id in adata.obs (i.e. untreated and treated ids)
        labels_key: str (default = None)
            string referring to the cell type labels key in adata.obs
        n_splits: int (default = 5)
            number of cross validation folds
        save_directory: str (default = 'figures')
            string referring to path for saving figures and data
        y_lim: list 
            ylim values for plotting the logistic regression coefficients
        ----------

        Returns
        ----------
    """
    make_directory(save_directory)
    for label_id in groups:
        file_id = f'{origin}_{label_id}'
        print(file_id)    
        adata_subset = subset_cells(adata = adata, labels_key = labels_key, label_id = label_id)

        coef, auc, tpr_fpr, pvals = logistic_reg(X = adata_subset.X.copy(), y = np.asarray(adata_subset.obs[condition_key]), feature_names = adata_subset.var_names, n_splits = n_splits)

        # save data
        coef.to_csv(os.path.join(save_directory, file_id+'_coef.csv'))
        tpr_fpr.to_csv(os.path.join(save_directory,file_id+'_tpr_fpr.csv'))
        pvals.to_csv(os.path.join(save_directory,file_id+'_pvals.csv'))
        np.save(os.path.join(save_directory, file_id+'_auc.npy'), auc)

        print(label_id, np.round(np.mean(auc), 3), '+/-', np.round(np.std(auc),3))

        # plot coefficients
        pvals_mean = pvals.mean(1)
        bar_df = pd.DataFrame(coef.mean(1))
        bar_df['pval'] = pvals_mean.drop('intercept')
        pval_symb = np.full(len(bar_df['pval']), '', dtype = "S10")
        pval_symb[bar_df['pval'] < 0.05] = '*'
        pval_symb[bar_df['pval'] < 0.01] = '**'
        pval_symb[bar_df['pval'] < 0.001] = '***'
        bar_df['pval_symb'] = pval_symb.astype('str')
        bar_df = bar_df.reset_index()
        bar_df = bar_df.sort_values(0, ascending = False)

        plot_coef(df = bar_df, ylim = ylim, filename_save = os.path.join(save_directory,file_id+'_coef'))

        # plot roc auc curves       
        plot_roc_auc(scores = auc, tpr_fpr = tpr_fpr, n_splits = n_splits, filename_save = os.path.join(save_directory,file_id+'_roc_auc'))

def plot_coef(df = None,
            ylim = [-10.5, 2.5],
            palette = 'coolwarm_r',
            figsize = (6, 5.5), 
            filename_save: str = None):
    """Plots coefficients from logistic regression analysis 
        Parameters
        df: pd.DataFrame
            dataframe containing coefficients and p-value symbols
        palette: str (default = 'coolwarm')
            cmap for plotting
        figsize: tuple (default = (6, 5.5))
            size for the plot
        filename_save: str (default = None)
            string referring to the filename for saving
        y_lim: list 
            ylim values for plotting the logistic regression coefficients
        ----------

        Returns
        ----------
    """     
    sns.set_style('ticks')
    fig, axes = plt.subplots(1,1, figsize = figsize, gridspec_kw = {'bottom':0.15})
    g = sns.barplot(x = 'index', y = 0, data = df, palette = palette, ax = axes)

    add_value_labels(df['pval_symb'].values, axes)

    axes.set_xticklabels(df['index'], rotation=45, horizontalalignment='right')
    axes.set_ylim(ylim) #-2 1.5 #-1.75 1.25
    axes.tick_params(labelsize=16)
    axes.set_xlabel('feature', fontsize = 16)
    axes.set_ylabel('coefficient', fontsize = 16)
    plt.savefig(filename_save+'.pdf', bbox_inches="tight")

def mean_barplots(adata = None,
                feature_order = ['pRB_over_RB', 'Ki67', 'pRB', 'RB', 'CDK2', 'CDK4', 'cycD1', 'cycE', 'Cdt1', 'E2F1', 'DNA', 'cycA', 'cycB1', 'p21'],
                colors_dict = {'G0':'#5CAD92', 'G1':'#594997', 'G2M':'#E7739A','S':'#0099CC'},
                ylim = [-1, 2.5],
                hue: str = 'phase', 
                hue_order = ['G0', 'G1', 'S', 'G2M'],
                save_directory: str = None,
                filename_save: str = None):
    """Plots mean expression for each marker
        Parameters
        adata: anndata.AnnData
            annotated data object containing preprocessed data (dimensions = cells x features)
        feature_order: list
            list containing markers of interest
        hue: str (default  = 'phase')
            string referring to group on x
        hue_order: list (default = ['G0', 'G1', 'S', 'G2M'])
            list containing order of groups for barplots
        colors_dict: dict (default = {'G0':'#5CAD92', 'G1':'#594997', 'G2M':'#E7739A','S':'#0099CC'})
           dictionary for color ids for groups in barplots
        y_lim: list 
            ylim values for plotting the expression values
        save_directory: str (default = None)
            string referring to directory for saving
        filename_save: str (default = None)
            string referring to the filename for saving
        ----------

        Returns
        ----------
    """
    make_directory(save_directory)
    df = pd.DataFrame(adata.X, columns = adata.var_names)
    df = df.loc[:, feature_order]
    df['condition'] = adata.obs['condition'].values #'well'
    df['phase'] = adata.obs['phase'].values
    df_phase_condition = df.groupby(['phase', 'condition']).median()

    fig, axes = plt.subplots(4, 4, figsize = (15, 15), gridspec_kw={'hspace': 0.45, 'wspace': 0.45, 'bottom':0.15})
    for i, ax in zip(range(0, len(df_phase_condition.columns)), axes.flat):
        g = sns.barplot(x = 'condition', y = df.columns[i], data = df.reset_index(), hue = hue,ax = ax, palette = colors_dict, hue_order = hue_order, errcolor = 'black', errwidth = 1, capsize = 0.05, estimator = np.mean)
        ax.tick_params(labelsize=12)
        g.set_xlabel('condition', fontsize = 14)
        g.set_ylabel('expression', fontsize = 14)
        g.set_title(df_phase_condition.columns[i], fontsize = 14)
        if i != len(df_phase_condition.columns) -1:
            ax.get_legend().remove()
        else:
            ax.legend(bbox_to_anchor=(1.50,0.95))
        g.set_ylim(ylim)

    axes[3,3].set_visible(False)
    axes[3,2].set_visible(False)
    plt.savefig(os.path.join(save_directory, f'{filename_save}.pdf'), bbox_inches = 'tight')

def make_directory(directory: str = None):
    """Creates a directory at the specified path if one doesn't exist.

    Parameters
    ----------
    directory : str
        A string specifying the directory path.

    Returns
    -------
    """
    if not os.path.exists(directory):
        os.makedirs(directory)