"""Sample NMF, robust programs, metaprograms, QC, plots and GSEApy ORA.

No expression normalization, centering, identifier conversion or imputation is
performed. Sample NMF resumes validated sample-by-K checkpoints; downstream
stages recreate their output directories and overwrite existing results.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import shutil
import warnings
from collections import Counter
from dataclasses import dataclass, field
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Union, cast

import numpy as np
from joblib import Parallel, delayed
from scipy import sparse
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

PathLike = Union[str, Path]
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "results"


@dataclass
class Program:
    """One component: full gene loadings plus a ranked top-gene signature."""

    program_id: str
    sample_id: str
    k: int
    component: int
    genes: tuple[str, ...]
    weights: np.ndarray
    top_genes: tuple[str, ...]
    study_id: Optional[str] = None


@dataclass
class NMFResult:
    """All selected fits, sample metadata and tested gene universe."""

    programs: dict[str, Program]
    background: tuple[str, ...]
    sample_info: list[dict[str, Any]]
    output_dir: Path


@dataclass
class MPResult:
    """MP gene lists, member program IDs and complete robust-program audit."""

    genes: dict[str, tuple[str, ...]]
    members: dict[str, tuple[str, ...]]
    nmf: NMFResult
    robust_ids: tuple[str, ...]
    unassigned_ids: tuple[str, ...]
    output_dir: Path
    qc: list[dict[str, Any]] = field(default_factory=list)


def _stage(
    root: PathLike, name: str, parameters: Mapping[str, Any], *,
    resume: bool = False, overwrite: bool = False,
    resume_ignore_keys: Sequence[str] = (),
) -> tuple[Path, logging.Logger]:
    # Resume NMF checkpoints or recreate downstream stages, as requested.
    path = Path(root) / name
    normalized_parameters = json.loads(json.dumps(dict(parameters), default=str))
    root_path = Path(root).resolve()
    resolved = path.resolve()
    if resolved.parent != root_path:
        raise ValueError(f"Unsafe stage path: {resolved}")
    existed = path.exists()
    if path.exists():
        if resume:
            if not path.is_dir():
                raise NotADirectoryError(path)
            parameter_path = path / "parameters.json"
            if not parameter_path.is_file():
                raise ValueError(f"Cannot resume {path}: parameters.json is missing")
            saved_parameters = json.loads(parameter_path.read_text(encoding="utf-8"))
            ignored = set(resume_ignore_keys)
            saved_comparable = {key: value for key, value in saved_parameters.items() if key not in ignored}
            current_comparable = {key: value for key, value in normalized_parameters.items() if key not in ignored}
            if saved_comparable != current_comparable:
                missing = "<missing>"
                differences = {
                    key: {
                        "existing": saved_comparable.get(key, missing),
                        "current": current_comparable.get(key, missing),
                    }
                    for key in sorted(set(saved_comparable) | set(current_comparable))
                    if saved_comparable.get(key, missing) != current_comparable.get(key, missing)
                }
                raise ValueError(
                    f"Cannot resume {path}: numerical/input parameters differ:\n"
                    f"{json.dumps(differences, ensure_ascii=False, indent=2)}"
                )
        elif overwrite:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
        else:
            raise FileExistsError(path)
    if not resume or not existed:
        path.mkdir(parents=True)
        _json(path / "parameters.json", normalized_parameters)
    logger = logging.getLogger(str(path.resolve()))
    logger.setLevel(logging.INFO)
    logger.propagate = False
    handler = logging.FileHandler(path / "run.log", mode="a", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    action = "Resumed" if resume and existed else "Started"
    if overwrite and existed:
        action = "Overwrote"
    logger.info("%s %s", action, name)
    return path, logger


def _finish(logger: logging.Logger) -> None:
    logger.info("Stage closed")
    for handler in list(logger.handlers):
        handler.close()
        logger.removeHandler(handler)


def _json(path: Path, value: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, default=str)


def _tsv(path: Path, rows: Iterable[Mapping[str, Any]], columns: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _replace_tsv(path: Path, rows: Iterable[Mapping[str, Any]], columns: Sequence[str]) -> None:
    # Aggregate tables are rebuilt atomically after all sample-by-K files exist.
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _positive_int(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _mp_sort_key(mp_id: str) -> tuple[int, int, str]:
    """Sort conventional MP numeric suffixes naturally, then other labels."""
    suffix = mp_id.rsplit("_", 1)[-1]
    return (0, int(suffix), mp_id) if suffix.isdigit() else (1, 0, mp_id)


def _all_nan_rows(matrix: Any) -> np.ndarray:
    """Return rows whose every value is NaN, including sparse matrices."""
    if sparse.issparse(matrix):
        csr = sparse.csr_matrix(matrix, dtype=np.float64, copy=True)
        csr.sum_duplicates()
        nan_values = np.isnan(csr.data).astype(np.int32)
        nan_matrix = sparse.csr_matrix(
            (nan_values, csr.indices.copy(), csr.indptr.copy()), shape=csr.shape,
        )
        return np.asarray(nan_matrix.sum(axis=1)).ravel() == csr.shape[1]
    return np.isnan(np.asarray(matrix)).all(axis=1)


def _fit_sample(
    matrix: Any, genes: tuple[str, ...], cell_ids: tuple[str, ...], sample_id: str,
    study_id: Optional[str], index: int, ranks: tuple[int, ...], n_runs: int,
    n_genes: int, max_iter: int, tol: float, solver: str, beta_loss: str,
    basis_scale: str, seed: int, threads: int, path: Path,
) -> tuple[list[Program], dict[str, Any], list[dict[str, Any]]]:
    # Fit genes x cells, matching the orientation of the supplied R basis.
    active = np.asarray(matrix.sum(axis=0)).ravel() > 0
    local_genes = tuple(np.asarray(genes)[active].tolist())
    x = matrix[:, active].T
    if min(x.shape) < max(ranks) or len(local_genes) < n_genes:
        raise ValueError(f"Sample {sample_id}: insufficient cells or active genes: {x.shape}")
    if np.any(np.asarray(x.sum(axis=0)).ravel() == 0):
        raise ValueError(f"Sample {sample_id}: zero-expression cells")
    programs: list[Program] = []
    diagnostics: list[dict[str, Any]] = []
    with threadpool_limits(limits=threads):
        for k in ranks:
            best: Any = None
            best_error = math.inf
            best_run = -1
            rank_diagnostics: list[dict[str, Any]] = []
            for run in range(n_runs):
                random_seed = int(np.random.SeedSequence([seed, index, k, run]).generate_state(1)[0])
                model = cast(Any, NMF)(n_components=k, init="random", solver=solver,
                            beta_loss=beta_loss, max_iter=max_iter, tol=tol,
                            random_state=random_seed)
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always", ConvergenceWarning)
                    w = model.fit_transform(x)
                if any(issubclass(item.category, ConvergenceWarning) for item in caught):
                    raise RuntimeError(f"NMF not converged: sample={sample_id}, k={k}, run={run}; increase max_iter or review tol")
                for item in caught:
                    warnings.warn(str(item.message), item.category)
                error = float(model.reconstruction_err_)
                if not np.isfinite(error):
                    raise ValueError(f"Nonfinite NMF objective: {sample_id}, k={k}")
                diagnostic = dict(sample_id=sample_id, k=k, run=run,
                                  seed=random_seed, objective=error, n_iter=int(model.n_iter_))
                diagnostics.append(diagnostic)
                rank_diagnostics.append(diagnostic)
                if error < best_error:
                    best_error, best_run = error, run
                    best = (w.copy(), model.components_.copy())
            w, h = best
            scales = w.sum(axis=0)
            if np.any(scales <= 0):
                raise ValueError(f"Degenerate NMF component in {sample_id}, k={k}")
            if basis_scale == "l1":
                w = w / scales
                h = h * scales[:, None]
            ids = []
            for component in range(k):
                weights = w[:, component]
                if np.count_nonzero(weights > 0) < n_genes:
                    raise ValueError(f"Fewer than {n_genes} positive loadings in {sample_id}, k={k}")
                order = np.lexsort((np.asarray(local_genes), -weights))[:n_genes]
                pid = f"S{index + 1:05d}_K{k}_P{component + 1}"
                ids.append(pid)
                programs.append(Program(pid, sample_id, k, component + 1, local_genes,
                                        weights, tuple(np.asarray(local_genes)[order].tolist()), study_id))
            target = path / f"S{index + 1:05d}_K{k}.npz"
            temporary = target.with_name(f".{target.name}.{os.getpid()}.tmp")
            with temporary.open("wb") as handle:
                np.savez_compressed(
                    handle, basis=w, activity=h.T, genes=np.asarray(local_genes),
                    cell_ids=np.asarray(cell_ids), program_ids=np.asarray(ids),
                    best_run=best_run, objective=best_error,
                    diagnostic_run=np.asarray([row["run"] for row in rank_diagnostics]),
                    diagnostic_seed=np.asarray([row["seed"] for row in rank_diagnostics]),
                    diagnostic_objective=np.asarray([row["objective"] for row in rank_diagnostics]),
                    diagnostic_n_iter=np.asarray([row["n_iter"] for row in rank_diagnostics]),
                )
            os.replace(temporary, target)
    return programs, dict(sample_id=sample_id, study_id=study_id, n_cells=len(cell_ids),
                          n_active_genes=len(local_genes), n_programs=len(programs)), diagnostics


def _existing_diagnostics(path: Path) -> dict[tuple[str, int], list[dict[str, Any]]]:
    diagnostics: dict[tuple[str, int], list[dict[str, Any]]] = {}
    source = path / "fit_diagnostics.tsv"
    if not source.is_file():
        return diagnostics
    with source.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            parsed = dict(
                sample_id=row["sample_id"], k=int(row["k"]), run=int(row["run"]),
                seed=int(row["seed"]), objective=float(row["objective"]),
                n_iter=int(row["n_iter"]),
            )
            diagnostics.setdefault((row["sample_id"], int(row["k"])), []).append(parsed)
    return diagnostics


def _load_completed_fit(
    fit_path: Path, *, genes: tuple[str, ...], cell_ids: tuple[str, ...],
    sample_id: str, study_id: Optional[str], index: int, k: int, n_genes: int,
    previous_diagnostics: Mapping[tuple[str, int], list[dict[str, Any]]],
) -> tuple[list[Program], list[dict[str, Any]]]:
    """Load one complete sample-by-K checkpoint and reject stale/corrupt files."""
    expected_ids = tuple(f"S{index + 1:05d}_K{k}_P{component + 1}" for component in range(k))
    with np.load(fit_path, allow_pickle=False) as fit:
        required = {"basis", "activity", "genes", "cell_ids", "program_ids", "best_run", "objective"}
        if not required.issubset(fit.files):
            raise ValueError("missing required arrays")
        saved_genes = tuple(fit["genes"].astype(str).tolist())
        saved_cells = tuple(fit["cell_ids"].astype(str).tolist())
        saved_ids = tuple(fit["program_ids"].astype(str).tolist())
        basis = np.asarray(fit["basis"], dtype=np.float64)
        activity = np.asarray(fit["activity"], dtype=np.float64)
        if saved_genes != genes or saved_cells != cell_ids or saved_ids != expected_ids:
            raise ValueError("gene, cell or program identifiers differ")
        if basis.shape != (len(genes), k) or activity.shape != (len(cell_ids), k):
            raise ValueError("basis or activity dimensions differ")
        if not np.isfinite(basis).all() or not np.isfinite(activity).all() or np.any(basis < 0) or np.any(activity < 0):
            raise ValueError("basis or activity contains invalid values")
        diagnostic_keys = {"diagnostic_run", "diagnostic_seed", "diagnostic_objective", "diagnostic_n_iter"}
        if diagnostic_keys.issubset(fit.files):
            diagnostic_rows = [
                dict(sample_id=sample_id, k=k, run=int(run), seed=int(seed),
                     objective=float(objective), n_iter=int(n_iter))
                for run, seed, objective, n_iter in zip(
                    fit["diagnostic_run"], fit["diagnostic_seed"],
                    fit["diagnostic_objective"], fit["diagnostic_n_iter"],
                )
            ]
        else:
            diagnostic_rows = list(previous_diagnostics.get((sample_id, k), []))
            if not diagnostic_rows:
                diagnostic_rows = [dict(
                    sample_id=sample_id, k=k, run=int(fit["best_run"]), seed=-1,
                    objective=float(fit["objective"]), n_iter=-1,
                )]
    programs = []
    for component, program_id in enumerate(expected_ids):
        weights = basis[:, component]
        if np.count_nonzero(weights > 0) < n_genes:
            raise ValueError("fewer positive loadings than n_genes")
        order = np.lexsort((np.asarray(genes), -weights))[:n_genes]
        programs.append(Program(
            program_id, sample_id, k, component + 1, genes, weights,
            tuple(np.asarray(genes)[order].tolist()), study_id,
        ))
    return programs, diagnostic_rows

def scale_by_batch(
    adata,
    batch_key,
    target_layer=None,
    max_value=10,
    min_value=None,
):
    """Scale genes within each batch while retaining original cell order.

    Batches containing only one cell are assigned NaN for all genes.
    Genes with zero variance in batches containing at least two cells
    are assigned zero.
    """
    if batch_key not in adata.obs:
        raise KeyError(f"Missing adata.obs column: {batch_key}")

    if adata.obs[batch_key].isna().any():
        raise ValueError(f"adata.obs[{batch_key!r}] contains missing values")

    if target_layer is None:
        source = adata.X
    else:
        if target_layer not in adata.layers:
            raise KeyError(f"Missing adata layer: {target_layer}")
        source = adata.layers[target_layer]

    # Initialize with NaN so singleton batches naturally remain NaN.
    layer_nmf = np.full(adata.shape, np.nan, dtype=np.float32)

    for batch, indices in adata.obs.groupby(
        batch_key,
        observed=True,
        sort=False,
    ).indices.items():
        indices = np.asarray(indices)

        # A single cell has no within-batch variance.
        if len(indices) == 1:
            layer_nmf[indices, :] = np.nan
            continue

        block = source[indices, :]

        if sparse.issparse(block):
            block = block.toarray()
        else:
            block = np.asarray(block).copy()

        block = block.astype(np.float64, copy=False)

        if not np.isfinite(block).all():
            raise ValueError(
                f"Batch {batch!r} contains NaN or infinite input values"
            )

        means = block.mean(axis=0)
        standard_deviations = block.std(axis=0, ddof=1)
        variable_genes = standard_deviations > 0

        block[:, variable_genes] = (
            block[:, variable_genes] - means[variable_genes]
        ) / standard_deviations[variable_genes]

        # Constant genes have no standardized variation.
        block[:, ~variable_genes] = 0

        if max_value is not None:
            block[block > max_value] = max_value

        if min_value is not None:
            block[block < min_value] = min_value

        layer_nmf[indices, :] = block.astype(np.float32)

    return layer_nmf

def run_sample_nmf(
    adata: Any, sample_key: str, *, layer: Optional[str] = None,
    ranks: Sequence[int] = tuple(range(4, 10)), n_runs: int = 10,
    n_genes: int = 50, study_key: Optional[str] = None,
    min_cells_per_sample: int = 50,
    clip_negative: bool = False, solver: str = "mu",
    beta_loss: str = "kullback-leibler", basis_scale: str = "l1",
    max_iter: int = 3000, tol: float = 1e-4, random_state: int = 0,
    n_jobs: int = 1, threads_per_job: Optional[int] = None,
    output_dir: PathLike = DEFAULT_OUTPUT,
) -> NMFResult:
    """Factorize each sample at every requested K and save its best restart.

    Parameters
    ----------
    adata : AnnData
        Cells x genes; var_names must be unique gene identifiers. Subset the
        desired cell population before calling. The object is never modified.
    sample_key, study_key : str, optional
        obs columns defining sample and study; each sample must have one study.
    layer : str or None
        Read adata.layers[layer], or adata.X when None. No normalization or
        centering is done. Cells whose entire expression row is NaN are removed
        internally and audited; any remaining NaN/Inf raises. Negatives raise
        unless clip_negative=True.
    ranks, n_runs, n_genes : sequence/int
        Defaults: K=4..9, ten random restarts per K, top 50 genes per component.
        Keep the restart with the smallest objective (not a consensus fit).
    min_cells_per_sample : int
        Exclude samples with fewer cells (default 50), recording every sample
        in sample_filter.tsv. The effective minimum is also bounded by max K.
    solver, beta_loss, basis_scale : str
        sklearn NMF controls. Default MU/KL; basis_scale='l1' rescales basis
        columns for cross-program tie breaking and adjusts activities to retain
        W@H. 'none' retains raw loadings as in the R score comparison.
    max_iter, tol, random_state : int/float
        Optimizer controls and reproducible seed. Nonconvergence raises.
    n_jobs, threads_per_job : int
        Parallel samples and BLAS threads. Product cannot exceed 60% of logical
        CPUs (minimum one CPU); workers hold whole individual sample matrices.
    output_dir : path
        Parent of 'sample_nmf'. Existing compatible sample-by-K NPZ files are
        validated and reused; only missing or invalid combinations are fitted.
        Downstream stages overwrite their existing output directories.

    Returns
    -------
    NMFResult
        Programs and gene universe; NPZ files contain gene basis and cell
        activity with IDs, TSV files contain signatures and fit diagnostics.
    """
    ranks = tuple(ranks)
    _positive_int(min_cells_per_sample, "min_cells_per_sample")
    for value, name in [(n_runs, "n_runs"), (n_genes, "n_genes"), (max_iter, "max_iter"), (n_jobs, "n_jobs")]:
        _positive_int(value, name)
    if not ranks or len(set(ranks)) != len(ranks):
        raise ValueError("ranks must be nonempty and unique")
    for k in ranks:
        _positive_int(k, "rank")
    if basis_scale not in {"l1", "none"} or tol <= 0 or random_state < 0:
        raise ValueError("Invalid basis_scale, tol or random_state")
    budget = max(1, int((os.cpu_count() or 1) * 0.6))
    threads = threads_per_job if threads_per_job is not None else max(1, budget // n_jobs)
    _positive_int(threads, "threads_per_job")
    if n_jobs * threads > budget:
        raise ValueError(f"CPU request exceeds {budget} logical CPUs")
    for key in [sample_key] + ([study_key] if study_key is not None else []):
        if key not in adata.obs:
            raise ValueError(f"Missing obs column: {key}")
    genes = tuple(adata.var_names.tolist())
    if not genes or len(set(genes)) != len(genes) or any(not isinstance(g, str) or not g.strip() for g in genes):
        raise ValueError("var_names must be nonempty unique gene identifiers")
    if not adata.obs_names.is_unique:
        raise ValueError("obs_names must be unique")
    if layer is not None and layer not in adata.layers:
        raise ValueError(f"Missing layer: {layer}")
    matrix = adata.X if layer is None else adata.layers[layer]
    if matrix is None:
        raise ValueError("Expression matrix is None")
    matrix = sparse.csr_matrix(matrix, dtype=np.float64) if sparse.issparse(matrix) else np.array(matrix, dtype=np.float64, copy=True)
    all_nan = _all_nan_rows(matrix)
    excluded_cell_ids = tuple(np.asarray(adata.obs_names)[all_nan].astype(str).tolist())
    keep_cells = ~all_nan
    matrix = matrix[keep_cells]
    obs = adata.obs.iloc[np.flatnonzero(keep_cells)]
    obs_names = np.asarray(adata.obs_names)[keep_cells]
    if matrix.shape[0] == 0:
        raise ValueError("No cells remain after removing all-NaN expression rows")
    for key in [sample_key] + ([study_key] if study_key is not None else []):
        if obs[key].isna().any():
            raise ValueError(f"Missing values in obs column after all-NaN cell removal: {key}")
        if any(not isinstance(v, str) or not v.strip() for v in obs[key]):
            raise ValueError(f"obs[{key!r}] must contain nonempty strings after all-NaN cell removal")
    values = cast(np.ndarray, matrix.data if sparse.issparse(matrix) else matrix)
    if not np.isfinite(values).all():
        raise ValueError("Expression contains NaN/Inf outside cells whose entire row is NaN")
    n_negative = int(np.count_nonzero(values < 0))
    if n_negative and not clip_negative:
        raise ValueError("Negative expression: supply nonnegative data or explicitly set clip_negative=True")
    if clip_negative:
        np.maximum(values, 0, out=values)
    sample_ids = sorted(set(obs[sample_key].tolist()))
    if not sample_ids:
        raise ValueError("No cells")
    specs = []
    sample_filter = []
    for index, sid in enumerate(sample_ids):
        mask = np.asarray(obs[sample_key] == sid)
        studies = set(obs.loc[mask, study_key]) if study_key is not None else set()
        if len(studies) > 1:
            raise ValueError(f"Sample {sid} spans multiple studies; use a globally unique sample ID")
        n_cells = int(mask.sum())
        reason = "below_min_cells" if n_cells < min_cells_per_sample else "below_max_rank" if n_cells < max(ranks) else "retained"
        sample_filter.append(dict(sample_id=sid, n_cells=n_cells, retained=reason == "retained", reason=reason))
        if reason == "retained":
            specs.append((index, sid, mask, next(iter(studies)) if studies else None))
    params = dict(sample_key=sample_key, layer=layer, ranks=ranks, n_runs=n_runs,
                  n_genes=n_genes, study_key=study_key, clip_negative=clip_negative,
                  min_cells_per_sample=min_cells_per_sample,
                  n_negative=n_negative, solver=solver, beta_loss=beta_loss,
                  basis_scale=basis_scale, max_iter=max_iter, tol=tol,
                  random_state=random_state, n_jobs=n_jobs, threads_per_job=threads,
                  versions={p: version(p) for p in ["numpy", "scipy", "scikit-learn", "anndata"]})
    if excluded_cell_ids:
        params["n_excluded_all_nan_cells"] = len(excluded_cell_ids)
    path, logger = _stage(
        output_dir, "sample_nmf", params, resume=True,
        resume_ignore_keys=("n_jobs", "threads_per_job"),
    )
    _replace_tsv(
        path / "excluded_all_nan_cells.tsv",
        (
            dict(cell_id=cell_id, sample_id=adata.obs.loc[cell_id, sample_key])
            for cell_id in excluded_cell_ids
        ),
        ["cell_id", "sample_id"],
    )
    _replace_tsv(path / "sample_filter.tsv", sample_filter, ["sample_id", "n_cells", "retained", "reason"])
    if not specs:
        _finish(logger)
        raise ValueError("No samples pass the cell-count filter; see sample_filter.tsv")
    logger.info("%d samples, %d cells, %d genes; checking checkpoints", len(specs), *matrix.shape)
    try:
        previous_diagnostics = _existing_diagnostics(path)
        loaded_programs: list[Program] = []
        loaded_diagnostics: list[dict[str, Any]] = []
        pending_specs = []
        sample_info = []
        completed_count = 0
        for index, sid, mask, study in specs:
            sample_matrix = matrix[mask]
            active = np.asarray(sample_matrix.sum(axis=0)).ravel() > 0
            local_genes = tuple(np.asarray(genes)[active].tolist())
            cell_ids = tuple(obs_names[mask].tolist())
            missing_ranks = []
            for k in ranks:
                fit_path = path / f"S{index + 1:05d}_K{k}.npz"
                if fit_path.is_file():
                    try:
                        loaded, diagnostics = _load_completed_fit(
                            fit_path, genes=local_genes, cell_ids=cell_ids,
                            sample_id=sid, study_id=study, index=index, k=k,
                            n_genes=n_genes, previous_diagnostics=previous_diagnostics,
                        )
                        loaded_programs.extend(loaded)
                        loaded_diagnostics.extend(diagnostics)
                        completed_count += 1
                        continue
                    except (OSError, ValueError, KeyError) as error:
                        logger.warning("Recomputing invalid checkpoint %s: %s", fit_path.name, error)
                missing_ranks.append(k)
            if missing_ranks:
                pending_specs.append((index, sid, mask, study, tuple(missing_ranks)))
            sample_info.append(dict(
                sample_id=sid, study_id=study, n_cells=len(cell_ids),
                n_active_genes=len(local_genes), n_programs=sum(ranks),
            ))
        logger.info(
            "Found %d completed sample-by-K fits; computing %d missing fits across %d samples",
            completed_count, sum(len(spec[-1]) for spec in pending_specs), len(pending_specs),
        )
        results = cast(Any, Parallel)(n_jobs=n_jobs, pre_dispatch=n_jobs)(
            delayed(_fit_sample)(matrix[mask], genes, tuple(obs_names[mask].tolist()),
                                 sid, study, index, missing_ranks, n_runs, n_genes, max_iter,
                                 tol, solver, beta_loss, basis_scale, random_state, threads, path)
            for index, sid, mask, study, missing_ranks in pending_specs)
        new_programs = [p for ps, _, _ in results for p in ps]
        programs = {p.program_id: p for p in loaded_programs + new_programs}
        expected_programs = len(specs) * sum(ranks)
        if len(programs) != expected_programs:
            raise RuntimeError(f"Expected {expected_programs} programs after resume, found {len(programs)}")
        diagnostics = loaded_diagnostics + [row for _, _, rows in results for row in rows]
        background = tuple(sorted({g for p in programs.values() for g in p.genes}))
        _replace_tsv(path / "sample_summary.tsv", sample_info, ["sample_id", "study_id", "n_cells", "n_active_genes", "n_programs"])
        _replace_tsv(path / "fit_diagnostics.tsv", sorted(diagnostics, key=lambda row: (row["sample_id"], row["k"], row["run"])), ["sample_id", "k", "run", "seed", "objective", "n_iter"])
        _replace_tsv(path / "program_genes.tsv", (dict(program_id=p.program_id, sample_id=p.sample_id, k=p.k, component=p.component, gene_rank=i + 1, gene=g, weight=float(p.weights[p.genes.index(g)])) for p in sorted(programs.values(), key=lambda item: item.program_id) for i, g in enumerate(p.top_genes)), ["program_id", "sample_id", "k", "component", "gene_rank", "gene", "weight"])
        _replace_tsv(path / "background_genes.tsv", (dict(gene=g) for g in background), ["gene"])
        logger.info("Retained %d programs", len(programs))
        return NMFResult(programs, background, sample_info, path)
    except Exception:
        logger.exception("Failed; output is incomplete")
        raise
    finally:
        _finish(logger)


def _overlaps(programs: Sequence[Program]) -> np.ndarray:
    # Sparse incidence multiplication avoids Python pairwise set comparisons.
    universe = {g: i for i, g in enumerate(sorted({g for p in programs for g in p.top_genes}))}
    rows, cols = [], []
    for i, p in enumerate(programs):
        rows.extend([i] * len(p.top_genes))
        cols.extend(universe[g] for g in p.top_genes)
    incidence = sparse.csr_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)), shape=(len(programs), len(universe)))
    return (incidence @ incidence.T).toarray()


def _consensus(ids: Sequence[str], programs: Mapping[str, Program], n_genes: int) -> tuple[str, ...]:
    counts = Counter(g for pid in ids for g in programs[pid].top_genes)
    max_scores = {g: 0.0 for g in counts}
    for pid in ids:
        p = programs[pid]
        for gene, score in zip(p.genes, p.weights):
            if gene in max_scores:
                max_scores[gene] = max(max_scores[gene], float(score))
    return tuple(sorted(counts, key=lambda g: (-counts[g], -max_scores[g], g))[:n_genes])


def derive_metaprograms(
    nmf: NMFResult, *, intra_min: int = 35, min_supporting_programs: int = 1,
    require_different_k: bool = True, inter_min: int = 10,
    inter_filter: bool = True, inter_reference: str = "all",
    redundancy_overlap: int = 10, remove_at_threshold: bool = True,
    founder_overlap: int = 10, min_founder_neighbors: int = 6,
    member_overlap: int = 10, output_dir: PathLike = DEFAULT_OUTPUT,
) -> MPResult:
    """Select robust programs, then greedily construct every automatic MP.

    nmf supplies all sample programs and full loadings for score tie breaking.
    intra_min is minimum same-sample overlap; min_supporting_programs counts
    other programs (not distinct K values); require_different_k restricts these
    matches to other ranks. inter_reference='all' compares against all other-
    sample programs (Method); 'robust' uses only intra-robust ones (R code).
    inter_filter enables inter_min. Rank candidates by maximum inter overlap;
    remove_at_threshold=True rejects overlap >= redundancy_overlap, while False
    rejects only > threshold (R code). All overlaps are absolute gene counts.
    A founder needs >= min_founder_neighbors other remaining programs with
    overlap >= founder_overlap (6 implements the R rule >5). Grow by highest
    overlap with the evolving MP while overlap >= member_overlap. Consensus
    ranks by occurrence count, then maximum basis score, then gene name.
    Tied founders/candidates use program ID. No final biological/QC filtering
    is performed here. Unclustered robust programs are explicitly retained in
    the audit. Saves genes, memberships, selection reasons and similarities.
    """
    ps = sorted(nmf.programs.values(), key=lambda p: p.program_id)
    if not ps:
        raise ValueError("No input programs")
    sizes = {len(p.top_genes) for p in ps}
    if len(sizes) != 1:
        raise ValueError("Programs must have equal signature sizes")
    n_genes = sizes.pop()
    for value, name in [(intra_min, "intra_min"), (inter_min, "inter_min"), (redundancy_overlap, "redundancy_overlap"), (founder_overlap, "founder_overlap"), (member_overlap, "member_overlap")]:
        _positive_int(value, name)
        if value > n_genes:
            raise ValueError(f"{name} exceeds signature size")
    _positive_int(min_supporting_programs, "min_supporting_programs")
    _positive_int(min_founder_neighbors, "min_founder_neighbors")
    if inter_reference not in {"all", "robust"}:
        raise ValueError("inter_reference must be all or robust")
    params = dict(intra_min=intra_min, min_supporting_programs=min_supporting_programs,
                  require_different_k=require_different_k, inter_min=inter_min,
                  inter_filter=inter_filter, inter_reference=inter_reference,
                  redundancy_overlap=redundancy_overlap, remove_at_threshold=remove_at_threshold,
                  founder_overlap=founder_overlap, min_founder_neighbors=min_founder_neighbors,
                  member_overlap=member_overlap, n_genes=n_genes, nmf_dir=nmf.output_dir)
    path, logger = _stage(output_dir, "all_mps", params, overwrite=True)
    ov = _overlaps(ps)
    support = [sum(i != j and p.sample_id == q.sample_id and (not require_different_k or p.k != q.k) and ov[i, j] >= intra_min for j, q in enumerate(ps)) for i, p in enumerate(ps)]
    intra_ok = np.asarray(support) >= min_supporting_programs
    inter = [max((int(ov[i, j]) for j, q in enumerate(ps) if p.sample_id != q.sample_id and (inter_reference == "all" or intra_ok[j])), default=0) for i, p in enumerate(ps)]
    selected: list[int] = []
    reasons = ["not_intra_robust" if not intra_ok[i] else "not_inter_robust" for i in range(len(ps))]
    for sample in sorted({p.sample_id for p in ps}):
        candidates = sorted((i for i, p in enumerate(ps) if p.sample_id == sample and intra_ok[i] and (not inter_filter or inter[i] >= inter_min)), key=lambda i: (-inter[i], ps[i].program_id))
        kept: list[int] = []
        for i in candidates:
            redundant = any(ov[i, j] >= redundancy_overlap if remove_at_threshold else ov[i, j] > redundancy_overlap for j in kept)
            reasons[i] = "redundant" if redundant else "selected"
            if not redundant:
                kept.append(i)
        selected.extend(kept)
    remaining = set(selected)
    mp_genes: dict[str, tuple[str, ...]] = {}
    members: dict[str, tuple[str, ...]] = {}
    while remaining:
        neighbors = {i: sum(ov[i, j] >= founder_overlap for j in remaining if i != j) for i in remaining}
        founder = min(remaining, key=lambda i: (-neighbors[i], ps[i].program_id))
        if neighbors[founder] < min_founder_neighbors:
            break
        cluster = [ps[founder].program_id]
        signature = ps[founder].top_genes
        remaining.remove(founder)
        while remaining:
            overlaps = {i: len(set(signature).intersection(ps[i].top_genes)) for i in remaining}
            candidate = min(remaining, key=lambda i: (-overlaps[i], ps[i].program_id))
            if overlaps[candidate] < member_overlap:
                break
            cluster.append(ps[candidate].program_id)
            remaining.remove(candidate)
            signature = _consensus(cluster, nmf.programs, n_genes)
        mp_id = f"MP_{len(mp_genes) + 1:03d}"
        mp_genes[mp_id], members[mp_id] = signature, tuple(cluster)
    robust_ids = tuple(ps[i].program_id for i in selected)
    unassigned = tuple(ps[i].program_id for i in sorted(remaining))
    result = MPResult(mp_genes, members, nmf, robust_ids, unassigned, path)
    _save_mps(result)
    _tsv(path / "program_selection.tsv", (dict(program_id=p.program_id, sample_id=p.sample_id, k=p.k, supporting_programs=support[i], max_inter_overlap=inter[i], reason=reasons[i]) for i, p in enumerate(ps)), ["program_id", "sample_id", "k", "supporting_programs", "max_inter_overlap", "reason"])
    sub = ov[np.ix_(selected, selected)]
    np.savez_compressed(path / "robust_similarity.npz", program_ids=np.asarray(robust_ids, dtype=str), overlap=sub, jaccard=sub / (2 * n_genes - sub))
    logger.info("%d raw, %d robust, %d MPs, %d unassigned", len(ps), len(selected), len(mp_genes), len(unassigned))
    _finish(logger)
    return result


def _save_mps(result: MPResult) -> None:
    _tsv(result.output_dir / "mp_genes.tsv", (dict(mp_id=mid, gene_rank=i + 1, gene=g) for mid, genes in result.genes.items() for i, g in enumerate(genes)), ["mp_id", "gene_rank", "gene"])
    _tsv(result.output_dir / "mp_members.tsv", (dict(mp_id=mid, program_id=pid, sample_id=result.nmf.programs[pid].sample_id, study_id=result.nmf.programs[pid].study_id, k=result.nmf.programs[pid].k) for mid, ids in result.members.items() for pid in ids), ["mp_id", "program_id", "sample_id", "study_id", "k"])
    _json(result.output_dir / "mp_lists.json", dict(genes=result.genes, members=result.members, robust_ids=result.robust_ids, unassigned_ids=result.unassigned_ids))


def filter_metaprograms(
    mps: MPResult, *, min_programs: int = 5, min_samples: int = 2,
    min_studies: Optional[int] = None, min_programs_per_study: int = 1,
    max_single_sample_fraction: float = 1.0,
    max_mito_fraction: Optional[float] = None,
    max_ribo_fraction: Optional[float] = None,
    mito_prefixes: Sequence[str] = ("MT-",),
    ribo_prefixes: Sequence[str] = ("RPL", "RPS"),
    exclude_mp_ids: Sequence[str] = (),
    exclusion_signatures: Optional[Mapping[str, Sequence[str]]] = None,
    exclusion_overlap_fraction: float = 0.2,
    output_dir: PathLike = DEFAULT_OUTPUT,
) -> MPResult:
    """Filter MPs independently of clustering and return a new MPResult.

    min_programs/min_samples count contributing programs/distinct samples.
    min_studies requires study_key during NMF. min_programs_per_study applies
    to EVERY contributing study; use min_studies=2 and this value=2 to reject
    MPs with an under-supported study. No studies are inferred from names.
    max_single_sample_fraction limits the largest sample's program fraction.
    Optional mito/ribo limits compare prefix gene fractions, case-insensitively;
    these are flags/proxies, not pathway enrichment tests. None disables them.
    exclusion_signatures rejects MP overlap / MP size >= the given fraction;
    callers must supply biologically appropriate doublet/cell-type signatures.
    exclude_mp_ids supports explicit reviewed exclusions. All input MPs receive
    a QC row with pass/fail reasons, even if none survive. IDs remain stable.
    """
    for value, name in [(min_programs, "min_programs"), (min_samples, "min_samples"), (min_programs_per_study, "min_programs_per_study")]:
        _positive_int(value, name)
    if min_studies is not None:
        _positive_int(min_studies, "min_studies")
    for value in [max_single_sample_fraction, max_mito_fraction, max_ribo_fraction, exclusion_overlap_fraction]:
        if value is not None and not 0 <= value <= 1:
            raise ValueError("Fraction thresholds must be in [0, 1]")
    if set(exclude_mp_ids) - set(mps.genes):
        raise ValueError("Unknown exclude_mp_ids")
    study_available = all(p.study_id is not None for p in mps.nmf.programs.values())
    if (min_studies is not None or min_programs_per_study > 1) and not study_available:
        raise ValueError("Study filtering requires study_key in run_sample_nmf")
    params = dict(min_programs=min_programs, min_samples=min_samples, min_studies=min_studies,
                  min_programs_per_study=min_programs_per_study,
                  max_single_sample_fraction=max_single_sample_fraction,
                  max_mito_fraction=max_mito_fraction, max_ribo_fraction=max_ribo_fraction,
                  mito_prefixes=mito_prefixes, ribo_prefixes=ribo_prefixes,
                  exclude_mp_ids=exclude_mp_ids, exclusion_signatures=exclusion_signatures,
                  exclusion_overlap_fraction=exclusion_overlap_fraction, source=mps.output_dir)
    path, logger = _stage(output_dir, "filtered_mps", params, overwrite=True)
    qc: list[dict[str, Any]] = []
    keep = []
    for mid, genes in mps.genes.items():
        programs = [mps.nmf.programs[pid] for pid in mps.members[mid]]
        samples = Counter(p.sample_id for p in programs)
        studies = Counter(p.study_id for p in programs) if study_available else Counter()
        mito = sum(g.upper().startswith(tuple(p.upper() for p in mito_prefixes)) for g in genes) / len(genes)
        ribo = sum(g.upper().startswith(tuple(p.upper() for p in ribo_prefixes)) for g in genes) / len(genes)
        largest = max(samples.values()) / len(programs)
        reasons = []
        checks = [(len(programs) < min_programs, "too_few_programs"), (len(samples) < min_samples, "too_few_samples"),
                  (min_studies is not None and len(studies) < min_studies, "too_few_studies"),
                  (bool(studies) and min(studies.values()) < min_programs_per_study, "under_supported_study"),
                  (largest > max_single_sample_fraction, "sample_dominance"),
                  (max_mito_fraction is not None and mito > max_mito_fraction, "mitochondrial"),
                  (max_ribo_fraction is not None and ribo > max_ribo_fraction, "ribosomal"),
                  (mid in exclude_mp_ids, "manual_exclusion")]
        reasons.extend(reason for condition, reason in checks if condition)
        for name, signature in (exclusion_signatures or {}).items():
            if len(set(genes).intersection(signature)) / len(genes) >= exclusion_overlap_fraction:
                reasons.append(f"signature:{name}")
        passed = not reasons
        if passed:
            keep.append(mid)
        qc.append(dict(mp_id=mid, n_programs=len(programs), n_samples=len(samples),
                       n_studies=len(studies) if study_available else None,
                       sample_program_counts=json.dumps(samples, sort_keys=True),
                       study_program_counts=json.dumps(studies, sort_keys=True),
                       largest_sample_fraction=largest, mito_fraction=mito, ribo_fraction=ribo,
                       passed=passed, reasons=";".join(reasons)))
    result = MPResult({k: mps.genes[k] for k in keep}, {k: mps.members[k] for k in keep}, mps.nmf,
                      mps.robust_ids, mps.unassigned_ids, path, qc)
    _save_mps(result)
    _tsv(path / "mp_qc.tsv", qc, ["mp_id", "n_programs", "n_samples", "n_studies", "sample_program_counts", "study_program_counts", "largest_sample_fraction", "mito_fraction", "ribo_fraction", "passed", "reasons"])
    logger.info("Retained %d of %d MPs", len(keep), len(mps.genes))
    _finish(logger)
    return result


def enrich_metaprograms(
    mps: MPResult, *, libraries: Sequence[str] = ("MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023"),
    organism: str = "Human", background: Optional[Sequence[str]] = None,
    min_gs_size: int = 40, max_gs_size: int = 500, alpha: float = 0.05,
    bh_scope: str = "per_mp", gene_sets: Optional[Mapping[str, Sequence[str]]] = None,
    output_dir: PathLike = DEFAULT_OUTPUT,
) -> list[dict[str, Any]]:
    """Download Enrichr libraries via GSEApy and run local hypergeometric/BH.

    mps contains final MP signatures. libraries are explicit Enrichr library
    names; availability is validated, downloaded sets are saved for provenance.
    gene_sets optionally supplies a dictionary instead (offline reproducibility).
    background defaults to the union of nonzero genes eligible for input NMF;
    a custom background must contain every query gene. Symbols must match the
    library identifiers; no automatic capitalization/mapping is performed.
    Sizes are filtered AFTER intersection with background. Hypergeometric
    P(X>=overlap) and BH use gseapy.stats. Zero-overlap terms are included with
    p=1 in the testing family. bh_scope='per_mp' corrects across all selected
    libraries for each MP; 'global' corrects across all MP-term pairs.
    alpha defines significance strictly q<alpha. Returns all tested terms,
    including nonsignificant ones, and writes all/significant TSV tables.
    """
    _positive_int(min_gs_size, "min_gs_size")
    if max_gs_size < min_gs_size or not 0 < alpha < 1 or bh_scope not in {"per_mp", "global"}:
        raise ValueError("Invalid gene-set size, alpha or bh_scope")
    if not mps.genes:
        raise ValueError("No final MPs to enrich")
    import gseapy as gp
    from gseapy.stats import calc_pvalues, multiple_testing_correction

    universe = set(mps.nmf.background if background is None else background)
    if not universe or any(set(gs) - universe for gs in mps.genes.values()):
        raise ValueError("Background is empty or does not contain every MP gene")
    if gene_sets is None and (not libraries or len(set(libraries)) != len(libraries)):
        raise ValueError("libraries must be nonempty and unique")
    path, logger = _stage(output_dir, "enrichment", dict(libraries=libraries, organism=organism,
                         background_size=len(universe), min_gs_size=min_gs_size,
                         max_gs_size=max_gs_size, alpha=alpha, bh_scope=bh_scope,
                         gseapy_version=version("gseapy"), source=mps.output_dir), overwrite=True)
    _tsv(path / "background_genes.tsv", (dict(gene=g) for g in sorted(universe)), ["gene"])
    raw: dict[str, list[str]] = {}
    if gene_sets is not None:
        raw = {name: list(genes) for name, genes in gene_sets.items()}
    else:
        available = set(gp.get_library_name(organism=organism))
        missing = set(libraries) - available
        if missing:
            raise ValueError(f"Unavailable Enrichr libraries: {sorted(missing)}; choose explicit available versions")
        for library in libraries:
            logger.info("Fetching %s", library)
            downloaded = gp.get_library(name=library, organism=organism, min_size=0, max_size=1000000)
            raw.update({f"{library}::{term}": list(genes) for term, genes in downloaded.items()})
    _json(path / "downloaded_gene_sets.json", raw)
    sets = {name: sorted(set(genes) & universe) for name, genes in raw.items()}
    sets = {name: genes for name, genes in sets.items() if min_gs_size <= len(genes) <= max_gs_size}
    if not sets:
        raise ValueError("No gene sets remain after background/size filtering; check identifiers")
    _json(path / "tested_gene_sets.json", sets)
    rows: list[dict[str, Any]] = []
    for mid, genes in mps.genes.items():
        query = set(genes)
        # GSEApy omits zero-hit categories; restore them before BH correction.
        columns = list(cast(Any, calc_pvalues)(query, sets, background=universe))
        pvalues = dict(zip(columns[0], columns[1])) if columns else {}
        for term, members in sorted(sets.items()):
            hits = sorted(query.intersection(members))
            rows.append(dict(mp_id=mid, term=term, n_background=len(universe),
                             n_query=len(query), n_gene_set=len(members), n_overlap=len(hits),
                             gene_ratio=len(hits) / len(query),
                             fold_enrichment=(len(hits) / len(query)) / (len(members) / len(universe)),
                             p_value=float(pvalues.get(term, 1.0)), q_value=1.0,
                             significant=False, overlap_genes=";".join(hits)))
    families = [list(range(len(rows)))] if bh_scope == "global" else [[i for i, row in enumerate(rows) if row["mp_id"] == mid] for mid in mps.genes]
    for indices in families:
        qvalues, _ = multiple_testing_correction([rows[i]["p_value"] for i in indices], alpha=alpha, method="benjamini-hochberg")
        for i, q in zip(indices, qvalues):
            rows[i]["q_value"], rows[i]["significant"] = float(q), bool(q < alpha)
    rows.sort(key=lambda row: (row["mp_id"], row["q_value"], row["term"]))
    columns_out = list(rows[0])
    _tsv(path / "enrichment_all.tsv", rows, columns_out)
    _tsv(path / "enrichment_significant.tsv", (row for row in rows if row["significant"]), columns_out)
    logger.info("Tested %d pairs; %d significant", len(rows), sum(row["significant"] for row in rows))
    _finish(logger)
    return rows


def plot_metaprograms(
    mps: MPResult, *, enrichment: Optional[Sequence[Mapping[str, Any]]] = None,
    top_terms: int = 10, alpha: float = 0.05, figsize: tuple[float, float] = (8, 6),
    dpi: int = 180, program_scope: str = "all", order_mode: str = "mp",
    linkage_method: str = "average", jaccard_percent_limits: tuple[float, float] = (2, 25),
    include_unassigned_programs: bool = True,
    output_dir: PathLike = DEFAULT_OUTPUT,
) -> list[Path]:
    """Save clustered program Jaccard heatmap, MP support and enrichment dots.

    mps can be the unfiltered or filtered result. program_scope='all' includes
    every per-sample, per-K NMF program; 'robust' restricts to robust programs.
    order_mode='global' clusters the complete selected matrix by Jaccard and
    uses final-MP membership only as annotation, without affecting order.
    order_mode='mp' places final MPs consecutively in ascending MP-number order,
    clusters programs within each MP, and appends
    clustered unassigned programs. The same final order is used for heatmap
    rows, columns and both annotation strips. 'global' performs one clustering
    across the selected scope and therefore may split members of a named MP.
    include_unassigned_programs=False restricts the matrix to programs belonging
    to final MPs. Here 'unassigned programs' means orphan robust programs, not
    final MPs lacking a higher-level functional annotation in the paper.
    Jaccard is displayed as percentage with a paper-like pale-to-purple scale;
    jaccard_percent_limits clips only display colors, not clustering distances.
    Dashed boxes and centered axis labels mark final-MP blocks. The order and
    membership are saved in program_jaccard_order.tsv.
    Support compares distinct samples and component counts (not cell counts).
    enrichment is the list returned by enrich_metaprograms; plot up to top_terms
    significant terms per MP using gene ratio, overlap size and -log10(BH q).
    Empty results get explicit placeholder panels rather than plotting errors.
    figsize/dpi control every PNG/PDF pair; returns their paths.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _positive_int(top_terms, "top_terms")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0, 1)")
    if program_scope not in {"all", "robust"}:
        raise ValueError("program_scope must be 'all' or 'robust'")
    if order_mode not in {"mp", "global"}:
        raise ValueError("order_mode must be 'mp' or 'global'")
    if len(jaccard_percent_limits) != 2 or not 0 <= jaccard_percent_limits[0] < jaccard_percent_limits[1] <= 100:
        raise ValueError("jaccard_percent_limits must be increasing values in [0, 100]")
    path, logger = _stage(output_dir, "plots", dict(source=mps.output_dir, top_terms=top_terms,
                          alpha=alpha, figsize=figsize, dpi=dpi,
                          program_scope=program_scope, order_mode=order_mode,
                          linkage_method=linkage_method,
                          jaccard_percent_limits=jaccard_percent_limits,
                          include_unassigned_programs=include_unassigned_programs), overwrite=True)
    outputs: list[Path] = []

    def save(fig: Any, name: str, tight: bool = True) -> None:
        if tight:
            fig.tight_layout()
        for suffix in ["pdf", "png"]:
            dest = path / f"{name}.{suffix}"
            fig.savefig(dest, dpi=dpi, bbox_inches="tight")
            outputs.append(dest)
        plt.close(fig)

    with plt.rc_context({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42}):
        from matplotlib.patches import Rectangle
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        membership = {pid: mid for mid, group in mps.members.items() for pid in group}
        ids = sorted(mps.nmf.programs) if program_scope == "all" else sorted(mps.robust_ids)
        if not include_unassigned_programs:
            ids = [pid for pid in ids if pid in membership]
        fig = plt.figure(figsize=figsize)
        if ids:
            ps = [mps.nmf.programs[pid] for pid in ids]
            ov = _overlaps(ps)
            sim = ov / (2 * len(ps[0].top_genes) - ov)
            distance = np.clip(1.0 - sim, 0.0, 1.0)
            np.fill_diagonal(distance, 0.0)
            def cluster_indices(indices: list[int]) -> list[int]:
                # Cluster one display block without changing membership grouping.
                if len(indices) < 2:
                    return indices
                local = distance[np.ix_(indices, indices)]
                local_order = leaves_list(linkage(squareform(local, checks=True), method=linkage_method))
                return [indices[int(i)] for i in local_order]

            block_sizes: list[int] = []
            block_labels: list[str] = []
            if order_mode == "mp":
                order_list: list[int] = []
                populated_mids = sorted([
                    mid for mid in mps.members
                    if any(membership.get(pid) == mid for pid in ids)
                ], key=_mp_sort_key)
                for mid in populated_mids:
                    block = cluster_indices([i for i, pid in enumerate(ids) if membership.get(pid) == mid])
                    if block:
                        order_list.extend(block)
                        block_sizes.append(len(block))
                        block_labels.append(mid)
                unassigned = cluster_indices([i for i, pid in enumerate(ids) if pid not in membership])
                if include_unassigned_programs and unassigned:
                    order_list.extend(unassigned)
                    block_sizes.append(len(unassigned))
                    block_labels.append("Unassigned")
                order = np.asarray(order_list)
            else:
                order = np.asarray(cluster_indices(list(range(len(ids)))))
            ids = [ids[int(i)] for i in order]
            sim = sim[np.ix_(order, order)]
            ax = fig.add_subplot(111)
            im = ax.imshow(100 * sim, vmin=jaccard_percent_limits[0], vmax=jaccard_percent_limits[1],
                           cmap="magma_r", interpolation="nearest", rasterized=True,
                           origin="upper")
            colorbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            colorbar.set_label("Similarity (Jaccard index, %)")
            if order_mode == "mp":
                ends = np.cumsum(block_sizes)
                starts = np.r_[0, ends[:-1]]
                centers = (starts + ends - 1) / 2
                ax.set_xticks(centers, block_labels, rotation=90)
                ax.set_yticks(centers, block_labels)
                for start, size, label in zip(starts, block_sizes, block_labels):
                    if label != "Unassigned":
                        ax.add_patch(Rectangle((start - 0.5, start - 0.5), size, size,
                                              fill=False, edgecolor="black", linewidth=1.2,
                                              linestyle=(0, (5, 4))))
            else:
                ax.set_xticks([])
                ax.set_yticks([])
            _tsv(path / "program_jaccard_order.tsv",
                 (dict(cluster_order=i + 1, program_id=pid,
                       sample_id=mps.nmf.programs[pid].sample_id,
                       k=mps.nmf.programs[pid].k,
                       component=mps.nmf.programs[pid].component,
                       final_mp=membership.get(pid, "Unassigned"))
                  for i, pid in enumerate(ids)),
                 ["cluster_order", "program_id", "sample_id", "k", "component", "final_mp"])
        else:
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, "No metaprograms", ha="center", va="center")
            ax.set_axis_off()
        scope_title = f"{program_scope} sample NMF programs (global Jaccard clustering)" if order_mode == "global" else f"{program_scope} sample NMF programs grouped by final MP"
        if not include_unassigned_programs:
            scope_title += " assigned to final MPs"
        fig.suptitle(f"Jaccard similarity of {scope_title}")
        save(fig, "program_jaccard", tight=False)
        fig, ax = plt.subplots(figsize=figsize)
        mids = list(mps.members)
        if mids:
            y = np.arange(len(mids))
            ax.barh(y - 0.2, [len(mps.members[mid]) for mid in mids], height=0.4, label="NMF programs", color="#4477AA")
            ax.barh(y + 0.2, [len({mps.nmf.programs[p].sample_id for p in mps.members[mid]}) for mid in mids], height=0.4, label="Samples", color="#66CCEE")
            ax.set_yticks(y, mids)
            ax.set_xlabel("Count")
            ax.legend(frameon=False)
        else:
            ax.text(0.5, 0.5, "No metaprograms", ha="center")
            ax.set_axis_off()
        save(fig, "mp_support")
        if enrichment is not None:
            for mid in mids:
                rows = sorted((row for row in enrichment if row["mp_id"] == mid and row["q_value"] < alpha), key=lambda row: (row["q_value"], row["term"]))[:top_terms]
                fig, ax = plt.subplots(figsize=figsize)
                if rows:
                    import textwrap
                    y = np.arange(len(rows))
                    dots = ax.scatter([r["gene_ratio"] for r in rows], y, s=[25 + 8 * r["n_overlap"] for r in rows], c=[-math.log10(max(r["q_value"], np.finfo(float).tiny)) for r in rows], cmap="viridis")
                    ax.set_yticks(y, [textwrap.fill(str(r["term"]).split("::")[-1], 45) for r in rows])
                    ax.invert_yaxis()
                    ax.set_xlabel("Overlap / MP genes (dot area reflects overlap count)")
                    fig.colorbar(dots, ax=ax, label="−log10(BH q)")
                else:
                    ax.text(0.5, 0.5, "No significant gene sets", ha="center")
                    ax.set_axis_off()
                ax.set_title(mid)
                save(fig, f"{mid}_enrichment")
    _finish(logger)
    return outputs
