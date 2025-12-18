#!/usr/bin/env python3
"""
Compute Functional Connectivity Matrices from Time Series

This script:
1. Loads parcellated fMRI time series
2. Computes Pearson correlation between all parcel pairs
3. Applies Fisher z-transformation (optional)
4. Saves FC matrices for all participants


"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import logging
import sys
from tqdm import tqdm

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
Compute Functional Connectivity Matrices from Time Series.

This script computes Pearson correlation matrices from parcellated
fMRI time series data, with optional Fisher z-transformation.

Example usage:
  python 02_compute_fc_matrices.py \\
      --participants participants_qc_passed.txt \\
      --input-dir data/fc_timeseries \\
      --output-dir data/fc_matrices \\
      --fisher-z
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
        '--input-dir', '-i',
        type=str,
        required=True,
        help='Directory containing FC time series (expects sub-{ID}/timeseries.csv structure)'
    )

    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        required=True,
        help='Output directory for FC matrices'
    )

    # Optional arguments
    parser.add_argument(
        '--input-filename',
        type=str,
        default='fMRI.combined_414.csv',
        help='Input time series filename pattern (default: fMRI.combined_414.csv)'
    )

    parser.add_argument(
        '--output-filename',
        type=str,
        default='fc_matrix.csv',
        help='Output matrix filename pattern (default: fc_matrix.csv)'
    )

    parser.add_argument(
        '--n-parcels',
        type=int,
        default=414,
        help='Expected number of parcels (default: 414)'
    )

    parser.add_argument(
        '--n-timepoints',
        type=int,
        default=490,
        help='Expected number of fMRI timepoints (default: 490)'
    )

    parser.add_argument(
        '--fisher-z',
        action='store_true',
        help='Apply Fisher z-transformation to correlations'
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

    parser.add_argument(
        '--n-jobs',
        type=int,
        default=1,
        help='Number of parallel jobs (default: 1)'
    )

    return parser.parse_args()


def compute_fc_matrix(
    timeseries: np.ndarray,
    fisher_z: bool = False
) -> np.ndarray:
    """
    Compute functional connectivity matrix from time series.

    Parameters
    ----------
    timeseries : np.ndarray
        Time series data of shape (n_parcels, n_timepoints)
    fisher_z : bool
        Whether to apply Fisher z-transformation

    Returns
    -------
    np.ndarray
        FC correlation matrix of shape (n_parcels, n_parcels)
    """
    # Compute Pearson correlation
    fc_matrix = np.corrcoef(timeseries)

    # Handle NaN values (from zero-variance parcels)
    fc_matrix = np.nan_to_num(fc_matrix, nan=0.0)

    # Enforce symmetry
    fc_matrix = (fc_matrix + fc_matrix.T) / 2.0

    # Ensure diagonal is 1.0
    np.fill_diagonal(fc_matrix, 1.0)

    # Apply Fisher z-transformation if requested
    if fisher_z:
        # Clip to avoid inf values at r=1 or r=-1
        fc_matrix_clipped = np.clip(fc_matrix, -0.9999, 0.9999)
        fc_matrix = np.arctanh(fc_matrix_clipped)
        # Reset diagonal to 0 (arctanh(1) = inf)
        np.fill_diagonal(fc_matrix, 0.0)

    return fc_matrix


def process_participant(
    participant_id: str,
    input_dir: Path,
    output_dir: Path,
    input_filename: str,
    output_filename: str,
    n_parcels: int,
    n_timepoints: int,
    fisher_z: bool = False
) -> dict:
    """
    Process a single participant's time series data.

    Returns
    -------
    dict
        Processing results and statistics
    """
    result = {
        'participant_id': participant_id,
        'success': False,
        'error': None
    }

    try:
        # Load time series
        ts_file = input_dir / f'sub-{participant_id}' / input_filename
        if not ts_file.exists():
            result['error'] = f"File not found: {ts_file}"
            return result

        ts_df = pd.read_csv(ts_file)

        # Extract numeric data (skip first column which is labels)
        ts_data = ts_df.iloc[:, 1:].values

        # Verify dimensions
        if ts_data.shape != (n_parcels, n_timepoints):
            result['error'] = f"Wrong dimensions: {ts_data.shape}, expected ({n_parcels}, {n_timepoints})"
            return result

        # Compute FC matrix
        fc_matrix = compute_fc_matrix(ts_data, fisher_z=fisher_z)

        # Verify output
        if fc_matrix.shape != (n_parcels, n_parcels):
            result['error'] = f"FC matrix wrong shape: {fc_matrix.shape}"
            return result

        if not np.allclose(fc_matrix, fc_matrix.T, rtol=1e-10):
            result['error'] = "FC matrix not symmetric"
            return result

        # Save FC matrix
        out_dir = output_dir / f'sub-{participant_id}'
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / output_filename

        pd.DataFrame(fc_matrix).to_csv(out_file, index=False, header=False)

        # Collect statistics (upper triangle excluding diagonal)
        upper_tri_idx = np.triu_indices_from(fc_matrix, k=1)
        upper_tri_vals = fc_matrix[upper_tri_idx]

        result.update({
            'success': True,
            'min_corr': float(np.min(upper_tri_vals)),
            'max_corr': float(np.max(upper_tri_vals)),
            'mean_corr': float(np.mean(upper_tri_vals)),
            'median_corr': float(np.median(upper_tri_vals)),
            'std_corr': float(np.std(upper_tri_vals)),
            'n_positive': int(np.sum(upper_tri_vals > 0)),
            'n_negative': int(np.sum(upper_tri_vals < 0)),
            'output_file': str(out_file)
        })

    except Exception as e:
        result['error'] = str(e)

    return result


def generate_report(
    stats_df: pd.DataFrame,
    output_file: Path,
    n_total: int,
    fisher_z: bool
):
    """Generate markdown report."""
    success_df = stats_df[stats_df['success']]
    n_success = len(success_df)

    with open(output_file, 'w') as f:
        f.write("# FC Matrix Computation Report\n\n")
        f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"**Fisher z-transform:** {'Yes' if fisher_z else 'No'}\n")
        f.write(f"**Total Participants:** {n_total}\n")
        f.write(f"**Successful:** {n_success}\n")
        f.write(f"**Failed:** {n_total - n_success}\n\n")

        if n_success > 0:
            f.write("---\n\n## Correlation Statistics\n\n")
            f.write(f"- **Mean correlation:** {success_df['mean_corr'].mean():.4f} ± {success_df['mean_corr'].std():.4f}\n")
            f.write(f"- **Median correlation:** {success_df['median_corr'].mean():.4f} ± {success_df['median_corr'].std():.4f}\n")
            f.write(f"- **Min range:** {success_df['min_corr'].min():.4f} to {success_df['min_corr'].max():.4f}\n")
            f.write(f"- **Max range:** {success_df['max_corr'].min():.4f} to {success_df['max_corr'].max():.4f}\n\n")

            f.write("## Edge Distribution\n\n")
            n_edges = 85491  # (414 * 413) / 2
            pct_pos = 100 * success_df['n_positive'].mean() / n_edges
            pct_neg = 100 * success_df['n_negative'].mean() / n_edges
            f.write(f"- **Positive correlations:** {success_df['n_positive'].mean():.0f} ({pct_pos:.1f}%)\n")
            f.write(f"- **Negative correlations:** {success_df['n_negative'].mean():.0f} ({pct_neg:.1f}%)\n")

        if n_total - n_success > 0:
            f.write("\n---\n\n## Failed Participants\n\n")
            failed_df = stats_df[~stats_df['success']]
            for _, row in failed_df.iterrows():
                f.write(f"- **{row['participant_id']}**: {row['error']}\n")

        f.write(f"\n---\n\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")


def main():
    """Main function."""
    args = parse_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Convert paths
    participants_file = Path(args.participants)
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    reports_dir = Path(args.reports_dir) if args.reports_dir else output_dir / 'reports'
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Validate inputs
    if not participants_file.exists():
        logger.error(f"Participants file not found: {participants_file}")
        sys.exit(1)

    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        sys.exit(1)

    # Load participants
    logger.info(f"Loading participants from {participants_file}")
    participants = pd.read_csv(participants_file, header=None)[0].astype(str).tolist()
    logger.info(f"Found {len(participants)} participants")

    # Process participants
    logger.info("Computing FC matrices...")
    results = []

    for pid in tqdm(participants, desc="Computing FC", disable=not args.verbose):
        result = process_participant(
            participant_id=pid,
            input_dir=input_dir,
            output_dir=output_dir,
            input_filename=args.input_filename,
            output_filename=args.output_filename,
            n_parcels=args.n_parcels,
            n_timepoints=args.n_timepoints,
            fisher_z=args.fisher_z
        )
        results.append(result)

    # Save statistics
    stats_df = pd.DataFrame(results)
    stats_file = reports_dir / 'fc_computation_stats.csv'
    stats_df.to_csv(stats_file, index=False)
    logger.info(f"Statistics saved: {stats_file}")

    # Generate report
    report_file = reports_dir / 'FC_COMPUTATION_REPORT.md'
    generate_report(stats_df, report_file, len(participants), args.fisher_z)
    logger.info(f"Report saved: {report_file}")

    # Print summary
    n_success = stats_df['success'].sum()
    n_failed = len(stats_df) - n_success

    print("\n" + "="*60)
    print("FC MATRIX COMPUTATION SUMMARY")
    print("="*60)
    print(f"Total participants:  {len(participants)}")
    print(f"Successful:          {n_success}")
    print(f"Failed:              {n_failed}")
    print(f"Fisher z-transform:  {'Yes' if args.fisher_z else 'No'}")
    print("="*60)

    if n_success > 0:
        success_df = stats_df[stats_df['success']]
        print(f"\nMean correlation: {success_df['mean_corr'].mean():.4f} ± {success_df['mean_corr'].std():.4f}")

    if n_failed > 0:
        print(f"\n⚠️  {n_failed} participants failed - see report for details")
    else:
        print("\n✅ ALL PARTICIPANTS PROCESSED SUCCESSFULLY")


if __name__ == '__main__':
    main()
