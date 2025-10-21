# the codes are modified from SOAPy

import copy
import anndata
import numpy as np
import pandas as pd
from typing import Optional, Union, Literal, Any
# from ..utils import _add_info_from_sample, _get_info_from_sample, _check_adata_type
# from .utils import (
#     _count_edge,
#     _randomize_helper,
#     Iterators,
#     _best_k,
# )
import numba as nb
from joblib import Parallel, delayed
import logging as logg
from anndata import AnnData

def _get_info_from_sample(
        adata: AnnData,
        sample_id: Optional[str] = None,
        key: Optional[str] = None,
        printf: bool = True
):
    """
    Getting information
    """

    if 'SOAPy' not in adata.uns.keys():
        if printf:
            logg.error('SOAPy was not found in adata.uns', exc_info=True)
        raise KeyError()

    if sample_id is None:
        if key not in adata.uns['SOAPy'].keys():
            if printf:
                logg.error(f'{key} analysis of {sample_id} was not present, please conduct {key} analysis first',
                           exc_info=True)
            raise KeyError()
        data = adata.uns['SOAPy'][key]
    else:
        if key not in adata.uns['SOAPy'][sample_id].keys():
            if printf:
                logg.error(f'{key} analysis of {sample_id} was not present, please conduct {key} analysis first',
                           exc_info=True)
            raise KeyError()
        data = adata.uns['SOAPy'][sample_id][key]

    return data


def _scale(
        adata: AnnData,
        library_id: Optional[str] = None,
        res: Optional[str] = None,
) -> float:
    """
    Get the scale of the coordinates to the drawing map
    """
    try:
        if library_id is None:
            library_id = list(adata.uns['spatial'].keys())[0]
        scale = adata.uns['spatial'][library_id]['scalefactors']['tissue_' + res + '_scalef']
    except KeyError:
        scale = 1.0

    return scale

def _count_edge(edge: pd.DataFrame,
                species_of_clusters: int,
                cell_type_dict: dict,
                ) -> np.ndarray:
    """

    Parameters
    ----------
    edge
    species_of_clusters
    cell_type_dict

    Returns
    -------

    """
    enrich_metrics = np.zeros((species_of_clusters, species_of_clusters), dtype=np.float32)

    for row in edge.itertuples():
        enrich_metrics[cell_type_dict[getattr(row, 'cluster_1')], cell_type_dict[getattr(row, 'cluster_2')]] += 1
        enrich_metrics[cell_type_dict[getattr(row, 'cluster_2')], cell_type_dict[getattr(row, 'cluster_1')]] += 1

    return enrich_metrics

def _randomize_helper(adata,
                      cluster_label,
                      species_of_clusters: int,
                      cell_type_dict: dict,
                      ):
    """

    Parameters
    ----------
    adata
    cluster_label
    species_of_clusters
    cell_type_dict

    Returns
    -------

    """
    Series_cluster = adata.obs[cluster_label]
    iter_cluster = Series_cluster.sample(frac=1.0,
                                         ).reset_index(drop=True)
    adata.obs[cluster_label] = iter_cluster.tolist()
    iter_edge = _preprocessing_of_graph(adata,
                                        cluster_label=cluster_label
                                        )

    enrichment = _count_edge(iter_edge, species_of_clusters, cell_type_dict)
    return enrichment

def _preprocessing_of_graph(adata: AnnData,
                            cluster_label: str,
                            ) -> pd.DataFrame:
    """

    Parameters
    ----------
    adata
    cluster_label

    Returns
    -------

    """

    distances, neighborhoods = adata.uns['SOAPy']['distance'], adata.uns['SOAPy']['indices']
    edges = []

    obs = adata.obs[cluster_label]
    obs_value = obs.values
    for neigh in neighborhoods:
        if len(neigh) == 0:
            continue
        point_1 = neigh[0]
        for point_2 in neigh:
            if point_2 == point_1 | point_2 == -1:
                continue
            elif point_2 < point_1:
                edge = [point_2, point_1,
                        obs_value[point_2],
                        obs_value[point_1]]
                edges.append(edge)

            elif point_2 > point_1:
                edge = [point_1, point_2,
                        obs_value[point_1],
                        obs_value[point_2]]
                edges.append(edge)

    df_edge = pd.DataFrame(data=np.array(edges), columns=['point_1', 'point_2',
                                                          'cluster_1', 'cluster_2'])

    df_edge.drop_duplicates(subset=['point_1', 'point_2'], inplace=True)
    return df_edge

def _add_info_from_sample(
        adata: AnnData,
        sample_id: Optional[str] = None,
        keys: Union[str, list, None] = None,
        add: Any = None,
) -> AnnData:
    """
    Storing information
    """

    if isinstance(keys, str):
        keys = [keys]
        add = [add]

    try:
        adata.uns['SOAPy']
    except KeyError:
        adata.uns['SOAPy'] = {}
        logg.warning('adata has not been initialized, and adata.uns[\'SOAPy\'] has been established')

    if sample_id is None:
        for index, key in enumerate(keys):
            adata.uns['SOAPy'][key] = add[index]
            logg.info(f'{key} information has been incorporated into adata.uns[\'SOAPy\'][{key}]')
    else:
        try:
            adata.uns['SOAPy'][sample_id]
        except KeyError:
            adata.uns['SOAPy'][sample_id] = {}
        for index, key in enumerate(keys):
            adata.uns['SOAPy'][sample_id][key] = add[index]
            logg.info(f'{key} information from {sample_id} has been incorporated into adata.uns[\'SOAPy\'][{sample_id}][{key}]')

    return adata

class Iterators:
    def __init__(self,
                 n_iter: int,
                 *args
                 ):
        self.stop = n_iter
        self.params = []
        for param in args:
            self.params.append(param)

    def __iter__(self):
        self.k = 0
        return self

    def __next__(self):
        if self.k < self.stop:
            self.k += 1
            return self.params
        else:
            raise StopIteration

class cell_network(object):
    """
    Collate the most basic neighborhood information.
    """

    def __init__(self,
                 adata,
                 cluster_key,
                 sample_key: Optional[str] = None,
                 sample=None,
                 ):
        if sample_key is not None:
            self.adata = adata[adata.obs[sample_key] == sample, :]
        else:
            self.adata = adata
        self._edge = _get_info_from_sample(self.adata, sample_id=sample, key='edges')
        self.cluster_key = cluster_key
        self.cluster = adata.obs[cluster_key]
        list_cluster = self.cluster.unique().tolist()
        new_list = [str(elem) for elem in list_cluster]
        self._cell_type_map = {v: i for i, v in enumerate(sorted(new_list))}
        self._species_of_clusters = len(self._cell_type_map)

    @property
    def cell_type(self):
        return self._cell_type_map.keys()

    @property
    def species_of_cell_type(self):
        return self._species_of_clusters

class cell_cell_interaction(cell_network):
    """
    Cell-cell interaction class
    """

    def __init__(self,
                 adata: anndata.AnnData,
                 cluster_key: str = 'clusters',
                 exclude: Union[str, dict] = None,
                 sample_key: Optional[str] = None,
                 sample=None,
                 ):

        super().__init__(adata, cluster_key, sample_key, sample)
        self.exclude = exclude

    def _randomize(self,
                   n_jobs: int,
                   num: int):
        adata = copy.deepcopy(self.adata)

        Iterator = Iterators(num,
                             adata,
                             self.cluster_key,
                             self._species_of_clusters,
                             self._cell_type_map)
        perms = Parallel(n_jobs=n_jobs)(delayed(_randomize_helper)(param[0], param[1], param[2], param[3])
                                        for param in Iterator)

        return perms

    def neighborhood_analysis(self,
                              n_jobs: int,
                              num: int,
                              method: str):

        @nb.njit()
        def enhance(matrix: np.ndarray, species: int):
            matrix_sum = np.sum(matrix, axis=0)
            for i in range(species):
                matrix[i, i] = matrix_sum[i] - matrix[i, i]
            for i in range(species):
                for j in range(i):
                    matrix[i, j] = matrix[i, j] / (matrix[i, i] + matrix[j, j] + 1)
                    matrix[j, i] = matrix[j, i] / (matrix[i, i] + matrix[j, j] + 1)
            for i in range(species):
                matrix[i, i] = 0

            return matrix

        mat_edge = _count_edge(self._edge, self._species_of_clusters, self._cell_type_map)
        perms = self._randomize(n_jobs=n_jobs, num=num)

        if method == 'included':
            pass
        elif method == 'excluded':
            mat_edge = enhance(mat_edge, self._species_of_clusters)
            species_list = [self._species_of_clusters] * num
            perms = Parallel(n_jobs=n_jobs)(delayed(enhance)(perm, spec) for perm, spec in zip(perms, species_list))
        perms = np.array(perms)
        # zscore_array = (mat_edge - perms.mean(axis=0, dtype=np.float32)) / perms.std(axis=0, dtype=np.float32, ddof=1)
        # zscore_df = pd.DataFrame(data=zscore_array, index=list(self._cell_type_map.keys()),
        #                          columns=list(self._cell_type_map.keys()))

        # return zscore_df

        # avoidance sample number
        p=np.sum(mat_edge<perms,axis=0)/num
        # interaction sample number
        q=np.sum(mat_edge>=perms,axis=0)/num
        times=np.where(q>=p,1-q,-(1-p))
        p_df=pd.DataFrame(times)
        return p_df

def neighborhood_analysis(
        adata: anndata.AnnData,
        method: Literal['excluded', 'include'] = 'excluded',
        cluster_key: str = 'clusters',
        sample_key: Optional[str] = None,
        sample: Union[int, str, list, None] = None,
        n_jobs: Optional[int] = None,
        n_iters: int = 1000,
        inplace: bool = True,
) -> anndata.AnnData:
    """
    Compute neighborhood enrichment Z-score by permutation test.

    Parameters
    ----------
    adata : anndata.AnnData
        An AnnData object containing spatial omics data and spatial information.
    method : str, optional
        'included': Z-scores of edges between two cell types were counted directly after randomization.
        'excluded': After randomization, remove self-connected edges between cells of the same type and calculate the
            z-score of edges between two cell types.
    cluster_key : str, optional
        The label of cluster in adata.obs.
    sample_key : str, optional
        The keyword of sample id in adata.obs.columns.
    sample: Union[int, str, list, None],
        The samples involved in calculating infiltration score.
    n_jobs : int, optional
        The maximum number of concurrently running jobs.
    n_iters : int, optional
        Number of rounds of random grouping for permutation tests.
    inplace : bool, optional
        Whether to change the original adata.

    Returns
    -------
    - :attr:`anndata.AnnData.uns` ``['SOAPy']['include_method' or 'exclude_method']['dic_crd_poly']`` - neighborhood score

    """
    if sample is not None and sample_key is None:
        logg.error(f'Mult-sample niche analysis cannot be specified without a given sample key')
        raise ValueError
    if sample is None and sample_key is not None:
        sample = adata.obs[sample_key].unique().tolist()
        logg.info(f'Use all samples in the niche calculation')
    if not isinstance(sample, list):
        sample = [sample]
    
    # adata = _check_adata_type(adata, 'spatial', inplace)
    
    for sample_id in sample:
        if sample_id is not None:
            bdata = adata[adata.obs[sample_key] == sample_id, :].copy()
            bdata.uns['SOAPy'] = {}
            bdata.uns['SOAPy'] = copy.deepcopy(adata.uns['SOAPy'][sample_id])
        else:
            bdata = adata
        new_cci = cell_cell_interaction(bdata, cluster_key)
        p_df = new_cci.neighborhood_analysis(n_jobs=n_jobs, num=n_iters, method=method)
        _add_info_from_sample(adata, sample_id=sample_id, keys=method + '_score', add=p_df)

    return adata

