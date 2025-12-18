#!/usr/bin/env python3
"""
Compute Structure-Function (SC-FC) Coupling

This script computes the coupling between structural connectivity (SC) and
functional connectivity (FC) matrices, measuring how well structural
connections predict functional connectivity.

Analyses:
1. Whole-brain SC-FC coupling (Spearman correlation across all edges)
2. Network-specific coupling (within each functional network)
3. Within-network vs between-network coupling
4. Regional (parcel-level) coupling


"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import logging
import sys
from scipy.stats import spearmanr, pearsonr
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="""
Compute Structure-Function (SC-FC) Coupling.

This script computes the correlation between structural and functional
connectivity matrices at multiple scales:
  - Whole-brain coupling
  - Network-level coupling
  - Regional (parcel-level) coupling

Example usage:
  python 03_compute_sc_fc_coupling.py \\
      --participants participants_qc_passed.txt \\
      --sc-dir data/sc_matrices \\
      --fc-dir data/fc_matrices \\
      --atlas atlas/combined_414_labels.csv \\
      --output-dir data/coupling_metrics
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Required arguments
    parser.add_argument(
        '--participants', '-p',
        type=str,
        required=True,
        help='Path to participant list file (one ID per line)'
    )

    parser.add_argument(
        '--sc-dir',
        type=str,
        required=True,
        help='Directory containing SC matrices (expects sub-{ID}/connectome.csv structure)'
    )

    parser.add_argument(
        '--fc-dir',
        type=str,
        required=True,
        help='Directory containing FC matrices (expects sub-{ID}/fc_matrix.csv structure)'
    )

    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        required=True,
        help='Output directory for coupling metrics'
    )

    # Optional arguments
    parser.add_argument(
        '--atlas',
        type=str,
        default=None,
        help='Path to atlas file with network assignments (CSV with parcel_id, network columns)'
    )

    parser.add_argument(
        '--sc-filename',
        type=str,
        default='connectome_sift2_fbc_10M.csv',
        help='SC matrix filename pattern (default: connectome_sift2_fbc_10M.csv)'
    )

    parser.add_argument(
        '--fc-filename',
        type=str,
        default='fc_matrix.csv',
        help='FC matrix filename pattern (default: fc_matrix.csv)'
    )

    parser.add_argument(
        '--n-parcels',
        type=int,
        default=414,
        help='Number of parcels (default: 414)'
    )

    parser.add_argument(
        '--n-cortical',
        type=int,
        default=360,
        help='Number of cortical parcels for network analysis (default: 360)'
    )

    parser.add_argument(
        '--nonzero-only',
        action='store_true',
        default=True,
        help='Only use edges with non-zero SC (default: True)'
    )

    parser.add_argument(
        '--include-zero-edges',
        action='store_true',
        help='Include edges with zero SC in coupling computation'
    )

    parser.add_argument(
        '--compute-regional',
        action='store_true',
        help='Compute regional (parcel-level) coupling'
    )

    parser.add_argument(
        '--reports-dir',
        type=str,
        default=None,
        help='Directory for reports (default: output-dir/reports)'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )

    return parser.parse_args()


def compute_wholebrain_coupling(
    sc_matrix: np.ndarray,
    fc_matrix: np.ndarray,
    nonzero_only: bool = True
) -> dict:
    """
    Compute whole-brain SC-FC coupling.

    Parameters
    ----------
    sc_matrix : np.ndarray
        Structural connectivity matrix
    fc_matrix : np.ndarray
        Functional connectivity matrix
    nonzero_only : bool
        Only use edges with non-zero SC

    Returns
    -------
    dict
        Coupling metrics
    """
    # Get upper triangle indices (exclude diagonal)
    upper_tri_idx = np.triu_indices_from(sc_matrix, k=1)

    # Extract edge weights
    sc_edges = sc_matrix[upper_tri_idx]
    fc_edges = fc_matrix[upper_tri_idx]

    result = {
        'n_edges_total': len(sc_edges),
        'n_edges_nonzero_sc': int(np.sum(sc_edges > 0)),
        'pct_nonzero_sc': 100 * np.sum(sc_edges > 0) / len(sc_edges)
    }

    # Compute coupling on all edges
    if len(sc_edges) > 10:
        rho_all, p_all = spearmanr(sc_edges, fc_edges)
        r_all, _ = pearsonr(sc_edges, fc_edges)
        result['coupling_spearman_all'] = rho_all
        result['coupling_pearson_all'] = r_all
        result['pval_spearman_all'] = p_all
    else:
        result['coupling_spearman_all'] = np.nan
        result['coupling_pearson_all'] = np.nan
        result['pval_spearman_all'] = np.nan

    # Compute coupling on non-zero SC edges only
    if nonzero_only:
        nonzero_mask = sc_edges > 0
        sc_nonzero = sc_edges[nonzero_mask]
        fc_nonzero = fc_edges[nonzero_mask]

        if len(sc_nonzero) > 10:
            rho_nz, p_nz = spearmanr(sc_nonzero, fc_nonzero)
            r_nz, _ = pearsonr(sc_nonzero, fc_nonzero)
            result['coupling_spearman_nonzero'] = rho_nz
            result['coupling_pearson_nonzero'] = r_nz
            result['pval_spearman_nonzero'] = p_nz
        else:
            result['coupling_spearman_nonzero'] = np.nan
            result['coupling_pearson_nonzero'] = np.nan
            result['pval_spearman_nonzero'] = np.nan

    return result


def compute_network_coupling(
    sc_matrix: np.ndarray,
    fc_matrix: np.ndarray,
    network_masks: dict,
    nonzero_only: bool = True
) -> list:
    """
    Compute network-level SC-FC coupling.

    Parameters
    ----------
    sc_matrix : np.ndarray
        Structural connectivity matrix
    fc_matrix : np.ndarray
        Functional connectivity matrix
    network_masks : dict
        Dictionary mapping network names to parcel indices
    nonzero_only : bool
        Only use edges with non-zero SC

    Returns
    -------
    list
        List of coupling metrics per network
    """
    results = []
    all_parcels = np.arange(sc_matrix.shape[0])

    for network, network_idx in network_masks.items():
        n_parcels = len(network_idx)

        # Within-network edges
        within_edges = []
        for i in range(len(network_idx)):
            for j in range(i + 1, len(network_idx)):
                within_edges.append((network_idx[i], network_idx[j]))

        # Between-network edges
        other_parcels = np.setdiff1d(all_parcels, network_idx)
        between_edges = []
        for i in network_idx:
            for j in other_parcels:
                if i < j:
                    between_edges.append((i, j))
                else:
                    between_edges.append((j, i))

        result = {
            'network': network,
            'n_parcels': n_parcels,
            'n_within_edges': len(within_edges),
            'n_between_edges': len(between_edges)
        }

        # Within-network coupling
        if len(within_edges) > 10:
            within_sc = np.array([sc_matrix[i, j] for i, j in within_edges])
            within_fc = np.array([fc_matrix[i, j] for i, j in within_edges])

            if nonzero_only:
                mask = within_sc > 0
                within_sc_use = within_sc[mask]
                within_fc_use = within_fc[mask]
            else:
                within_sc_use = within_sc
                within_fc_use = within_fc

            result['n_within_nonzero_sc'] = int(np.sum(within_sc > 0))

            if len(within_sc_use) > 10:
                rho, _ = spearmanr(within_sc_use, within_fc_use)
                result['within_coupling_spearman'] = rho
            else:
                result['within_coupling_spearman'] = np.nan
        else:
            result['within_coupling_spearman'] = np.nan
            result['n_within_nonzero_sc'] = 0

        # Between-network coupling
        if len(between_edges) > 10:
            between_sc = np.array([sc_matrix[i, j] for i, j in between_edges])
            between_fc = np.array([fc_matrix[i, j] for i, j in between_edges])

            if nonzero_only:
                mask = between_sc > 0
                between_sc_use = between_sc[mask]
                between_fc_use = between_fc[mask]
            else:
                between_sc_use = between_sc
                between_fc_use = between_fc

            result['n_between_nonzero_sc'] = int(np.sum(between_sc > 0))

            if len(between_sc_use) > 10:
                rho, _ = spearmanr(between_sc_use, between_fc_use)
                result['between_coupling_spearman'] = rho
            else:
                result['between_coupling_spearman'] = np.nan
        else:
            result['between_coupling_spearman'] = np.nan
            result['n_between_nonzero_sc'] = 0

        results.append(result)

    return results


def compute_regional_coupling(
    sc_matrix: np.ndarray,
    fc_matrix: np.ndarray,
    nonzero_only: bool = True
) -> list:
    """
    Compute regional (parcel-level) SC-FC coupling.

    For each parcel, compute the correlation between its SC and FC
    connectivity profiles (connections to all other parcels).

    Parameters
    ----------
    sc_matrix : np.ndarray
        Structural connectivity matrix
    fc_matrix : np.ndarray
        Functional connectivity matrix
    nonzero_only : bool
        Only use edges with non-zero SC

    Returns
    -------
    list
        List of coupling metrics per parcel
    """
    n_parcels = sc_matrix.shape[0]
    results = []

    for parcel_idx in range(n_parcels):
        # Get connectivity profile (exclude self-connection)
        sc_profile = np.delete(sc_matrix[parcel_idx, :], parcel_idx)
        fc_profile = np.delete(fc_matrix[parcel_idx, :], parcel_idx)

        result = {
            'parcel_idx': parcel_idx,
            'n_connections': len(sc_profile),
            'n_nonzero_sc': int(np.sum(sc_profile > 0))
        }

        if nonzero_only:
            mask = sc_profile > 0
            sc_use = sc_profile[mask]
            fc_use = fc_profile[mask]
        else:
            sc_use = sc_profile
            fc_use = fc_profile

        if len(sc_use) > 10:
            rho, pval = spearmanr(sc_use, fc_use)
            result['coupling_spearman'] = rho
            result['pval_spearman'] = pval
        else:
            result['coupling_spearman'] = np.nan
            result['pval_spearman'] = np.nan

        results.append(result)

    return results


def generate_report(
    wholebrain_df: pd.DataFrame,
    network_df: pd.DataFrame,
    output_file: Path,
    n_participants: int
):
    """Generate markdown summary report."""
    with open(output_file, 'w') as f:
        f.write("# SC-FC Coupling Analysis Report\n\n")
        f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"**Participants:** {n_participants}\n\n")
        f.write("---\n\n")

        # Whole-brain coupling
        f.write("## Whole-Brain SC-FC Coupling\n\n")

        if 'coupling_spearman_nonzero' in wholebrain_df.columns:
            mean_coupling = wholebrain_df['coupling_spearman_nonzero'].mean()
            std_coupling = wholebrain_df['coupling_spearman_nonzero'].std()
            f.write("### Non-zero SC edges only (recommended):\n\n")
            f.write(f"- **Mean coupling (Spearman):** {mean_coupling:.4f} ± {std_coupling:.4f}\n")
            f.write(f"- **Range:** {wholebrain_df['coupling_spearman_nonzero'].min():.4f} to {wholebrain_df['coupling_spearman_nonzero'].max():.4f}\n\n")

        if 'coupling_spearman_all' in wholebrain_df.columns:
            mean_coupling = wholebrain_df['coupling_spearman_all'].mean()
            std_coupling = wholebrain_df['coupling_spearman_all'].std()
            f.write("### All edges:\n\n")
            f.write(f"- **Mean coupling (Spearman):** {mean_coupling:.4f} ± {std_coupling:.4f}\n\n")

        # Edge statistics
        f.write("### Edge Statistics:\n\n")
        f.write(f"- **Total edges:** {wholebrain_df['n_edges_total'].iloc[0]:,}\n")
        f.write(f"- **Non-zero SC edges:** {wholebrain_df['n_edges_nonzero_sc'].mean():.0f} ± {wholebrain_df['n_edges_nonzero_sc'].std():.0f}\n")
        f.write(f"- **Percent non-zero:** {wholebrain_df['pct_nonzero_sc'].mean():.1f}%\n\n")

        # Network-level coupling
        if network_df is not None and len(network_df) > 0:
            f.write("---\n\n## Network-Level Coupling\n\n")

            # Within-network
            network_summary = network_df.groupby('network').agg({
                'within_coupling_spearman': ['mean', 'std'],
                'between_coupling_spearman': ['mean', 'std']
            }).round(4)

            f.write("### Within-Network Coupling:\n\n")
            f.write("| Network | Mean | Std |\n")
            f.write("|---------|------|-----|\n")
            for network in sorted(network_df['network'].unique()):
                try:
                    mean_val = network_summary.loc[network, ('within_coupling_spearman', 'mean')]
                    std_val = network_summary.loc[network, ('within_coupling_spearman', 'std')]
                    f.write(f"| {network} | {mean_val:.4f} | {std_val:.4f} |\n")
                except:
                    pass

            f.write("\n### Between-Network Coupling:\n\n")
            f.write("| Network | Mean | Std |\n")
            f.write("|---------|------|-----|\n")
            for network in sorted(network_df['network'].unique()):
                try:
                    mean_val = network_summary.loc[network, ('between_coupling_spearman', 'mean')]
                    std_val = network_summary.loc[network, ('between_coupling_spearman', 'std')]
                    f.write(f"| {network} | {mean_val:.4f} | {std_val:.4f} |\n")
                except:
                    pass

        f.write(f"\n---\n\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")


def main():
    """Main function."""
    args = parse_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Handle nonzero flag
    nonzero_only = not args.include_zero_edges

    # Convert paths
    participants_file = Path(args.participants)
    sc_dir = Path(args.sc_dir)
    fc_dir = Path(args.fc_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    reports_dir = Path(args.reports_dir) if args.reports_dir else output_dir / 'reports'
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Validate inputs
    if not participants_file.exists():
        logger.error(f"Participants file not found: {participants_file}")
        sys.exit(1)

    if not sc_dir.exists():
        logger.error(f"SC directory not found: {sc_dir}")
        sys.exit(1)

    if not fc_dir.exists():
        logger.error(f"FC directory not found: {fc_dir}")
        sys.exit(1)

    # Load participants
    logger.info(f"Loading participants from {participants_file}")
    participants = pd.read_csv(participants_file, header=None)[0].astype(str).tolist()
    logger.info(f"Found {len(participants)} participants")

    # Load atlas if provided
    network_masks = None
    if args.atlas:
        atlas_file = Path(args.atlas)
        if atlas_file.exists():
            logger.info(f"Loading atlas from {atlas_file}")
            atlas_df = pd.read_csv(atlas_file)

            # Get cortical parcels with network assignments
            cortical_atlas = atlas_df[atlas_df['parcel_id'] <= args.n_cortical].copy()
            networks = cortical_atlas['network'].dropna().unique()

            network_masks = {}
            for network in networks:
                network_parcels = cortical_atlas[cortical_atlas['network'] == network]['parcel_id'].values
                network_masks[network] = network_parcels - 1  # Convert to 0-indexed

            logger.info(f"Found {len(networks)} networks: {', '.join(sorted(networks))}")
        else:
            logger.warning(f"Atlas file not found: {atlas_file}")

    # Process participants
    logger.info("Computing SC-FC coupling...")
    wholebrain_results = []
    network_results = []
    regional_results = []

    for pid in tqdm(participants, desc="Computing coupling", disable=not args.verbose):
        try:
            # Load matrices
            sc_file = sc_dir / f'sub-{pid}' / args.sc_filename
            fc_file = fc_dir / f'sub-{pid}' / args.fc_filename

            if not sc_file.exists() or not fc_file.exists():
                logger.warning(f"Missing files for {pid}")
                continue

            sc_matrix = pd.read_csv(sc_file, header=None).values
            fc_matrix = pd.read_csv(fc_file, header=None).values

            # Whole-brain coupling
            wb_result = compute_wholebrain_coupling(sc_matrix, fc_matrix, nonzero_only)
            wb_result['participant_id'] = pid
            wholebrain_results.append(wb_result)

            # Network-level coupling
            if network_masks:
                net_results = compute_network_coupling(sc_matrix, fc_matrix, network_masks, nonzero_only)
                for nr in net_results:
                    nr['participant_id'] = pid
                network_results.extend(net_results)

            # Regional coupling
            if args.compute_regional:
                reg_results = compute_regional_coupling(sc_matrix, fc_matrix, nonzero_only)
                for rr in reg_results:
                    rr['participant_id'] = pid
                regional_results.extend(reg_results)

        except Exception as e:
            logger.error(f"Error processing {pid}: {e}")

    # Save results
    wholebrain_df = pd.DataFrame(wholebrain_results)
    wholebrain_file = output_dir / 'sc_fc_coupling_wholebrain.csv'
    wholebrain_df.to_csv(wholebrain_file, index=False)
    logger.info(f"Whole-brain coupling saved: {wholebrain_file}")

    network_df = None
    if network_results:
        network_df = pd.DataFrame(network_results)
        network_file = output_dir / 'sc_fc_coupling_network.csv'
        network_df.to_csv(network_file, index=False)
        logger.info(f"Network coupling saved: {network_file}")

    if regional_results:
        regional_df = pd.DataFrame(regional_results)
        regional_file = output_dir / 'sc_fc_coupling_regional.csv'
        regional_df.to_csv(regional_file, index=False)
        logger.info(f"Regional coupling saved: {regional_file}")

    # Generate report
    report_file = reports_dir / 'SC_FC_COUPLING_REPORT.md'
    generate_report(wholebrain_df, network_df, report_file, len(participants))
    logger.info(f"Report saved: {report_file}")

    # Print summary
    print("\n" + "="*60)
    print("SC-FC COUPLING ANALYSIS SUMMARY")
    print("="*60)
    print(f"Participants processed: {len(wholebrain_results)}/{len(participants)}")
    print(f"Non-zero SC edges only: {nonzero_only}")

    if 'coupling_spearman_nonzero' in wholebrain_df.columns:
        mean_coupling = wholebrain_df['coupling_spearman_nonzero'].mean()
        std_coupling = wholebrain_df['coupling_spearman_nonzero'].std()
        print(f"\nWhole-brain coupling (Spearman): {mean_coupling:.4f} ± {std_coupling:.4f}")

    print("="*60)
    print("✅ SC-FC COUPLING ANALYSIS COMPLETE")


if __name__ == '__main__':
    main()
