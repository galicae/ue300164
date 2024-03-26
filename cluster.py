# first mute future warnings and only then import pandas
import datetime
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

from tqdm import tqdm
import numpy as np
import pandas as pd
import scanpy as sc


def grouped_obs_mean(adata, group_key, layer=None):
    """
    Helper function to calculate average expression per group in an `AnnData` object. From
    ivirshup's answer to a scanpy issue:
    https://github.com/theislab/scanpy/issues/181#issuecomment-534867254

    Parameters
    ----------
    adata : `AnnData` object
        The object to analyse
    group_key : string
        The name of the `.obs` entry to group by.
    layer : string (default: None)
        The name of the `.layer` to calculate the averages on. Will default to `.X` if none is
        specified.

    Returns
    -------
    out : `pandas.DataFrame`
        A `DataFrame` where the columns are groups (the categories of `group_key`) and the rows are
        genes of the selected layer.
    """
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    out.index = adata.var_names
    return out


def grouped_obs_present(adata, group_key, layer=None):
    """
    Helper function to calculate how many cells express each gene per group in an `AnnData` object.
    Based on ivirshup's answer to a similar scanpy issue:
    https://github.com/theislab/scanpy/issues/181#issuecomment-534867254

    Parameters
    ----------
    adata : `AnnData` object The object to analyse group_key : string The name of the `.obs` entry
        to group by. layer : string (default: None) The name of the `.layer` to calculate the
        averages on. Will default to `.X` if none is specified.

    Returns
    -------
    out : `pandas.DataFrame` A `DataFrame` where the columns are groups (the categories of
        `group_key`) and the rows are genes of the selected layer.
    """
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names,
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel((X > 0).sum(axis=0, dtype=np.float64))
    return out


def find_closest_candidate(averages, candidates):
    """
    Given a list of eligible clusters (where #cells is smaller than a cutoff), this function will
    search for their closest neighbors and return the candidate and its closest neighbor that have
    the least distance.

    Parameters
    ----------
    averages : `pandas.DataFrame`
        A data frame of the average expression of the clusters. Output of `grouped_obs_mean()`.
    candidates : `numpy.array`
        An array of candidate names. The names must be contained in `averages`.

    Returns
    -------
    to_merge : `int`, `str`, or `obj`
        The name of the first cluster to be merged.
    merge_with : `int`, `str`, or `obj`
        The name of the second cluster to be merged.
    """
    corr = np.corrcoef(averages.T)
    clusters = np.array(averages.columns)
    cand_indices = np.array([np.where(clusters == c) for c in candidates]).flatten()
    ordered_neighbs = np.argsort(-np.abs(corr[cand_indices]), axis=1)
    distances = np.sort(-np.abs(corr[cand_indices]), axis=1)
    nearest = ordered_neighbs[:, 1]
    nearest_dists = np.abs(distances[:, 1])
    closest = np.argmax(nearest_dists)
    to_merge = candidates[closest]
    merge_with = nearest[closest]
    return to_merge, merge_with


def differentially_expressed_genes(
    adata,
    groupby,
    group_1,
    group_2,
    significance=0.05,
    fold_cutoff=2,
    method="wilcoxon",
    layer=None,
):
    """
    Find the number of differentially expressed genes between two clusters using the Wilcoxon test.
    Genes with a statistically significant over- or underexpression beyond a fold-change cutoff are
    counted.

    Parameters
    ----------
    adata : `AnnData` object
        The object to analyse
    groupby : str
        The `.obs` category to use.
    group_1 : str
        The first group.
    group_2 : str
        The second group.
    significance : float (default: 0.05)
        The adjusted p-value cutoff for significance of differential expression (Wilcoxon test) or
        the score cut-off (logistic regression)
    fold_cutoff : float (default: 2)
        The fold change to require for differential expression (only relevant for Wilcoxon test).
    method: str
        One of "wilcoxon" and "logreg".

    Returns
    -------
    int
        The number of over- or underexpressed genes between the two clusters.
    """
    use_raw = False
    if layer is None:
        use_raw = True
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[group_1],
        reference=group_2,
        method=method,
        layer=layer,
        use_raw=use_raw,
    )

    if method == "wilcoxon":
        logfoldchanges = adata.uns["rank_genes_groups"]["logfoldchanges"].astype(float)
        pvals_adj = adata.uns["rank_genes_groups"]["pvals_adj"].astype(float)

        overexpressed = (logfoldchanges > np.log(fold_cutoff)) & (
            pvals_adj < significance
        )
        underexpressed = (logfoldchanges < -np.log(fold_cutoff)) & (
            pvals_adj < significance
        )
        return int(np.sum(overexpressed | underexpressed))
    elif method == "logreg":
        return np.sum(
            np.array(adata.uns["rank_genes_groups"]["scores"].astype(float)) > 0
        )


def find_closest_neighbours(adata, clustering, layer=None):
    averages = grouped_obs_mean(adata, clustering, layer=layer)
    corr = np.corrcoef(averages.T)
    clusters = np.array(averages.columns)

    ordered_neighbs = np.argsort(-np.abs(corr), axis=1)
    distances = np.sort(-np.abs(corr), axis=1)
    nearest = clusters[ordered_neighbs[:, 1]]
    nearest_dists = np.abs(distances[:, 1])

    highest_cor = np.argsort(-nearest_dists)

    candidates = pd.DataFrame(
        {
            "from": clusters[highest_cor],
            "to": nearest[highest_cor],
            "dist": nearest_dists[highest_cor],
        }
    )
    return candidates


def find_next_merge(
    adata,
    clustering,
    candidates,
    tested,
    num_genes=20,
    fold_cutoff=2,
    significance=0.05,
    method="wilcoxon",
    layer=None,
    verbose=False,
):
    for _index, row in candidates.iterrows():
        i, j, _d = row
        if verbose:
            print("testing %s, %s" % (i, j))
        if set([i, j]) in tested:
            continue

        deg = differentially_expressed_genes(
            adata,
            clustering,
            i,
            j,
            significance=significance,
            fold_cutoff=fold_cutoff,
            method=method,
            layer=layer,
        )
        if deg < num_genes:
            return i, j
        else:
            tested.append(set([i, j]))
    return None


def _merge_clusters(adata, clustering, to_merge, merge_with, verbose=True):
    numeric = np.array(adata.obs[clustering], dtype=int)
    new_cluster = np.max(numeric) + 1
    numeric[adata.obs[clustering] == to_merge] = new_cluster
    numeric[adata.obs[clustering] == merge_with] = new_cluster
    if verbose:
        print("Merged %s + %s -> %d" % (to_merge, merge_with, new_cluster))
    adata.obs[clustering] = np.array(numeric, dtype=str)
    adata.obs[clustering] = adata.obs[clustering].astype("category")


def merge_underpopulated(adata, merged_clustering, min_cells=20, verbose=False):
    # first lnum_geness clusters with <cutoff> or fewer cells into nearest neighbor clusters (only if necessary)
    vals, counts = np.unique(adata.obs[merged_clustering], return_counts=True)
    candidates = np.array(vals[counts < min_cells], dtype=int)
    if verbose:
        print(
            "Seeking to merge %d cluster(s) with %d or fewer cells into nearest neighbor cluster."
            % (len(candidates), min_cells)
        )
    merged = np.array(adata.obs[merged_clustering], dtype=int)
    while any(counts < min_cells):
        averages = grouped_obs_mean(adata, merged_clustering)
        merged = np.array(adata.obs[merged_clustering], dtype=int)
        to_merge, merge_with = find_closest_candidate(averages, candidates)
        if verbose:
            print("Merging cluster %d with cluster %d..." % (to_merge, merge_with))
        new_cluster = np.max(merged) + 1
        merged[adata.obs[merged_clustering] == to_merge] = new_cluster
        merged[adata.obs[merged_clustering] == merge_with] = new_cluster
        adata.obs[merged_clustering] = merged
        vals, counts = np.unique(adata.obs[merged_clustering], return_counts=True)
        candidates = np.array(vals[counts < min_cells], dtype=int)

    adata.obs[merged_clustering] = np.array(merged, dtype=str)
    adata.obs[merged_clustering] = adata.obs[merged_clustering].astype("category")
    return adata


def merge_underdifferentiated(
    adata,
    merged_clustering,
    num_genes=20,
    fold_difference=2,
    significance=0.05,
    method="wilcoxon",
    layer=None,
    verbose=False,
):
    tested = []
    candidates = find_closest_neighbours(adata, merged_clustering, layer=layer)
    merge_pair = find_next_merge(
        adata,
        merged_clustering,
        candidates,
        tested,
        num_genes=num_genes,
        fold_cutoff=fold_difference,
        significance=significance,
        method=method,
        layer=layer,
        verbose=verbose,
    )

    while merge_pair is not None:
        to_merge, merge_with = merge_pair
        tested.append(set(merge_pair))
        _merge_clusters(adata, merged_clustering, to_merge, merge_with, verbose=verbose)
        candidates = find_closest_neighbours(adata, merged_clustering, layer=layer)
        merge_pair = find_next_merge(
            adata,
            merged_clustering,
            candidates,
            tested,
            num_genes=num_genes,
            fold_cutoff=fold_difference,
            significance=significance,
            layer=layer,
            verbose=verbose,
        )

    return adata


def merge_clusters(
    adata,
    original_clustering,
    num_genes=20,
    min_cells=10,
    fold_difference=2,
    significance=0.05,
    method="wilcoxon",
    layer=None,
    verbose=False,
):
    merged_clustering = original_clustering + "_merged"
    merged = np.array(adata.obs[original_clustering], dtype=int)
    adata.obs[merged_clustering] = merged

    if verbose:
        now = datetime.datetime.now()
        print("## START_MERGE1 %02d:%02d" % (now.hour, now.minute))
    adata = merge_underpopulated(
        adata, merged_clustering, min_cells=min_cells, verbose=verbose
    )
    if verbose:
        now = datetime.datetime.now()
        print("## END_MERGE1 %02d:%02d" % (now.hour, now.minute))

    # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.settings.verbosity = 1

    # in the second loop we will merge clusters that don't satisfy threshold criteria
    if verbose:
        now = datetime.datetime.now()
        print("## START_MERGE2 %02d:%02d" % (now.hour, now.minute))
    adata = merge_underdifferentiated(
        adata,
        merged_clustering,
        num_genes=num_genes,
        fold_difference=fold_difference,
        significance=significance,
        method=method,
        layer=layer,
        verbose=verbose,
    )
    if verbose:
        now = datetime.datetime.now()
        print("## END_MERGE2 %02d:%02d" % (now.hour, now.minute))


def map_clusterings(adata, cluster_from, cluster_to):
    original_clusters = np.array(adata.obs[cluster_from].cat.categories)
    merged_clusters = np.array(adata.obs[cluster_to].cat.categories)
    overlap = np.zeros((len(original_clusters), len(merged_clusters)))
    for i, original in enumerate(tqdm(original_clusters)):
        for j, merged in enumerate(merged_clusters):
            original_ids = adata.obs.index[adata.obs[cluster_from] == original]
            merged_ids = adata.obs.index[adata.obs[cluster_to] == merged]
            overlap[i][j] = len(np.intersect1d(original_ids, merged_ids)) / len(original_ids)
    overlap = pd.DataFrame(overlap, index=original_clusters, columns=merged_clusters)
    return overlap