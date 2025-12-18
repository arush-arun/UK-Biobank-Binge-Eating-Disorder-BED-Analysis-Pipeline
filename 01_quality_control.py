#!/usr/bin/env python3
"""
Quality Control for Structural and Functional Connectivity Data

This script performs comprehensive QC checks on connectivity data:
1. SC matrices: dimensions, symmetry, NaN values, value ranges
2. FC time series: dimensions, NaN values, value ranges
3. Generates detailed QC report with pass/fail for each participant

Author: Arush Honnesetty
License: MIT
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import logging
import sys

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
Quality Control for Structural and Functional Connectivity Data.

This script performs comprehensive QC checks on:
  - Structural Connectivity (SC) matrices
  - Functional Connectivity (FC) time series

Example usage:
  python 01_quality_control.py \\
      --participants participants.txt \\
      --sc-dir data/sc_matrices \\
      --fc-dir data/fc_timeseries \\
      --output-dir reports/qc
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
        '--sc-dir', '-s',
        type=str,
        required=True,
        help='Directory containing SC matrices (expects sub-{ID}/connectome.csv structure)'
    )

    parser.add_argument(
        '--fc-dir', '-f',
        type=str,
        required=True,
        help='Directory containing FC time series (expects sub-{ID}/timeseries.csv structure)'
    )

    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        required=True,
        help='Output directory for QC reports'
    )

    # Optional arguments
    parser.add_argument(
        '--sc-filename',
        type=str,
        default='connectome_sift2_fbc_10M.csv',
        help='SC matrix filename pattern (default: connectome_sift2_fbc_10M.csv)'
    )

    parser.add_argument(
        '--fc-filename',
        type=str,
        default='fMRI.combined_414.csv',
        help='FC time series filename pattern (default: fMRI.combined_414.csv)'
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
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )

    return parser.parse_args()


def qc_sc_matrix(sc_file: Path, n_parcels: int) -> dict:
    """
    Perform QC on a structural connectivity matrix.

    Parameters
    ----------
    sc_file : Path
        Path to SC matrix CSV file
    n_parcels : int
        Expected number of parcels

    Returns
    -------
    dict
        QC results dictionary
    """
    qc_result = {
        'file_exists': False,
        'correct_dimensions': False,
        'is_symmetric': False,
        'no_nan': False,
        'no_inf': False,
        'no_negative': False,
        'non_zero_edges': False,
        'diagonal_valid': False,
        'min_value': np.nan,
        'max_value': np.nan,
        'mean_value': np.nan,
        'median_value': np.nan,
        'pass_qc': False
    }

    if not sc_file.exists():
        return qc_result

    qc_result['file_exists'] = True

    try:
        # Load SC matrix
        sc_matrix = pd.read_csv(sc_file, header=None).values

        # Check dimensions
        if sc_matrix.shape == (n_parcels, n_parcels):
            qc_result['correct_dimensions'] = True

        # Check symmetry
        if np.allclose(sc_matrix, sc_matrix.T, rtol=1e-5):
            qc_result['is_symmetric'] = True

        # Check for NaN
        if not np.isnan(sc_matrix).any():
            qc_result['no_nan'] = True

        # Check for Inf
        if not np.isinf(sc_matrix).any():
            qc_result['no_inf'] = True

        # Check for negative values
        if (sc_matrix >= 0).all():
            qc_result['no_negative'] = True

        # Check diagonal
        diagonal_vals = np.diag(sc_matrix)
        if (diagonal_vals >= 0).all():
            qc_result['diagonal_valid'] = True

        # Check for non-zero edges
        if np.count_nonzero(sc_matrix) > 0:
            qc_result['non_zero_edges'] = True

        # Compute statistics
        qc_result['min_value'] = float(np.min(sc_matrix))
        qc_result['max_value'] = float(np.max(sc_matrix))
        qc_result['mean_value'] = float(np.mean(sc_matrix))
        qc_result['median_value'] = float(np.median(sc_matrix))

        # Overall pass/fail
        qc_result['pass_qc'] = all([
            qc_result['file_exists'],
            qc_result['correct_dimensions'],
            qc_result['is_symmetric'],
            qc_result['no_nan'],
            qc_result['no_inf'],
            qc_result['no_negative'],
            qc_result['non_zero_edges'],
            qc_result['diagonal_valid']
        ])

    except Exception as e:
        logger.error(f"Error loading SC matrix {sc_file}: {e}")

    return qc_result


def qc_fc_timeseries(fc_file: Path, n_parcels: int, n_timepoints: int) -> dict:
    """
    Perform QC on functional connectivity time series.

    Parameters
    ----------
    fc_file : Path
        Path to FC time series CSV file
    n_parcels : int
        Expected number of parcels
    n_timepoints : int
        Expected number of timepoints

    Returns
    -------
    dict
        QC results dictionary
    """
    qc_result = {
        'file_exists': False,
        'correct_dimensions': False,
        'correct_parcels': False,
        'correct_timepoints': False,
        'no_nan': False,
        'no_inf': False,
        'reasonable_range': False,
        'min_value': np.nan,
        'max_value': np.nan,
        'mean_value': np.nan,
        'std_value': np.nan,
        'pass_qc': False
    }

    if not fc_file.exists():
        return qc_result

    qc_result['file_exists'] = True

    try:
        # Load FC time series
        fc_df = pd.read_csv(fc_file)

        # Check dimensions
        actual_parcels = len(fc_df)
        actual_timepoints = len(fc_df.columns) - 1  # Subtract label column

        if actual_parcels == n_parcels:
            qc_result['correct_parcels'] = True

        if actual_timepoints == n_timepoints:
            qc_result['correct_timepoints'] = True

        qc_result['correct_dimensions'] = (
            qc_result['correct_parcels'] and qc_result['correct_timepoints']
        )

        # Extract numeric data
        fc_data = fc_df.iloc[:, 1:].values

        # Check for NaN
        if not np.isnan(fc_data).any():
            qc_result['no_nan'] = True

        # Check for Inf
        if not np.isinf(fc_data).any():
            qc_result['no_inf'] = True

        # Compute statistics
        qc_result['min_value'] = float(np.min(fc_data))
        qc_result['max_value'] = float(np.max(fc_data))
        qc_result['mean_value'] = float(np.mean(fc_data))
        qc_result['std_value'] = float(np.std(fc_data))

        # Check reasonable range
        if qc_result['min_value'] >= 0 and qc_result['min_value'] < qc_result['max_value']:
            qc_result['reasonable_range'] = True

        # Overall pass/fail
        qc_result['pass_qc'] = all([
            qc_result['file_exists'],
            qc_result['correct_dimensions'],
            qc_result['no_nan'],
            qc_result['no_inf'],
            qc_result['reasonable_range']
        ])

    except Exception as e:
        logger.error(f"Error loading FC time series {fc_file}: {e}")

    return qc_result


def generate_qc_summary(
    qc_sc_df: pd.DataFrame,
    qc_fc_df: pd.DataFrame,
    output_file: Path,
    n_total: int
):
    """Generate markdown QC summary report."""

    n_sc_pass = qc_sc_df['pass_qc'].sum()
    n_fc_pass = qc_fc_df['pass_qc'].sum()

    # Participants passing both
    sc_pass_ids = set(qc_sc_df[qc_sc_df['pass_qc']]['participant_id'].astype(str))
    fc_pass_ids = set(qc_fc_df[qc_fc_df['pass_qc']]['participant_id'].astype(str))
    both_pass_ids = sorted(sc_pass_ids & fc_pass_ids)
    n_both_pass = len(both_pass_ids)

    with open(output_file, 'w') as f:
        f.write("# Quality Control Summary\n\n")
        f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"**Total Participants:** {n_total}\n\n")
        f.write("---\n\n")

        f.write("## Overall Results\n\n")
        f.write("| Data Type | Pass | Fail | Pass Rate |\n")
        f.write("|-----------|------|------|-----------|\n")
        f.write(f"| SC Matrices | {n_sc_pass}/{n_total} | {n_total-n_sc_pass} | {100*n_sc_pass/n_total:.1f}% |\n")
        f.write(f"| FC Time Series | {n_fc_pass}/{n_total} | {n_total-n_fc_pass} | {100*n_fc_pass/n_total:.1f}% |\n")
        f.write(f"| **Both Pass** | {n_both_pass}/{n_total} | {n_total-n_both_pass} | {100*n_both_pass/n_total:.1f}% |\n\n")

        # SC QC details
        f.write("---\n\n## SC Matrix QC Checks\n\n")
        for check in ['file_exists', 'correct_dimensions', 'is_symmetric', 'no_nan',
                      'no_inf', 'no_negative', 'non_zero_edges', 'diagonal_valid']:
            n_pass = qc_sc_df[check].sum()
            status = "✅" if n_pass == n_total else "⚠️"
            f.write(f"- **{check}**: {status} {n_pass}/{n_total} pass\n")

        # FC QC details
        f.write("\n---\n\n## FC Time Series QC Checks\n\n")
        for check in ['file_exists', 'correct_dimensions', 'no_nan', 'no_inf', 'reasonable_range']:
            n_pass = qc_fc_df[check].sum()
            status = "✅" if n_pass == n_total else "⚠️"
            f.write(f"- **{check}**: {status} {n_pass}/{n_total} pass\n")

        f.write(f"\n---\n\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    return both_pass_ids


def main():
    """Main function."""
    args = parse_arguments()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Convert paths
    participants_file = Path(args.participants)
    sc_dir = Path(args.sc_dir)
    fc_dir = Path(args.fc_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

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

    # QC SC matrices
    logger.info("Running QC on SC matrices...")
    qc_sc_results = []
    for i, pid in enumerate(participants, 1):
        if args.verbose or i % 50 == 0:
            logger.info(f"[{i}/{len(participants)}] QC SC: sub-{pid}")

        sc_file = sc_dir / f'sub-{pid}' / args.sc_filename
        result = qc_sc_matrix(sc_file, args.n_parcels)
        result['participant_id'] = pid
        qc_sc_results.append(result)

    # QC FC time series
    logger.info("Running QC on FC time series...")
    qc_fc_results = []
    for i, pid in enumerate(participants, 1):
        if args.verbose or i % 50 == 0:
            logger.info(f"[{i}/{len(participants)}] QC FC: sub-{pid}")

        fc_file = fc_dir / f'sub-{pid}' / args.fc_filename
        result = qc_fc_timeseries(fc_file, args.n_parcels, args.n_timepoints)
        result['participant_id'] = pid
        qc_fc_results.append(result)

    # Save results
    qc_sc_df = pd.DataFrame(qc_sc_results)
    qc_fc_df = pd.DataFrame(qc_fc_results)

    qc_sc_file = output_dir / 'qc_sc_matrices.csv'
    qc_fc_file = output_dir / 'qc_fc_timeseries.csv'

    qc_sc_df.to_csv(qc_sc_file, index=False)
    qc_fc_df.to_csv(qc_fc_file, index=False)

    logger.info(f"SC QC results saved: {qc_sc_file}")
    logger.info(f"FC QC results saved: {qc_fc_file}")

    # Generate summary
    summary_file = output_dir / 'QC_SUMMARY.md'
    passed_ids = generate_qc_summary(qc_sc_df, qc_fc_df, summary_file, len(participants))
    logger.info(f"QC summary saved: {summary_file}")

    # Save passed participants
    passed_file = output_dir / 'participants_qc_passed.txt'
    with open(passed_file, 'w') as f:
        for pid in passed_ids:
            f.write(f"{pid}\n")
    logger.info(f"Passed participants saved: {passed_file}")

    # Print summary
    n_sc_pass = qc_sc_df['pass_qc'].sum()
    n_fc_pass = qc_fc_df['pass_qc'].sum()

    print("\n" + "="*60)
    print("QC SUMMARY")
    print("="*60)
    print(f"Total participants: {len(participants)}")
    print(f"SC matrices pass:   {n_sc_pass}/{len(participants)} ({100*n_sc_pass/len(participants):.1f}%)")
    print(f"FC time series pass: {n_fc_pass}/{len(participants)} ({100*n_fc_pass/len(participants):.1f}%)")
    print(f"Both pass:          {len(passed_ids)}/{len(participants)} ({100*len(passed_ids)/len(participants):.1f}%)")
    print("="*60)

    if len(passed_ids) == len(participants):
        print("✅ ALL PARTICIPANTS PASSED QC")
    else:
        print(f"⚠️  {len(participants) - len(passed_ids)} participants failed QC")


if __name__ == '__main__':
    main()
