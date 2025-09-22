#!/usr/bin/env python3
"""
UK Biobank BED Data Organization and Analysis - Refactored Version

This module provides a comprehensive, object-oriented approach to processing
UK Biobank eating disorder data with proper class structure and modular design.

Usage:
    python3 data_org_BED_refactored.py

Author: Refactored for better maintainability and modularity
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ast
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


class UKBiobankDataProcessor:
    """
    Main class for processing UK Biobank eating disorder data.
    Handles data loading, filtering, group assignment, and exclusion criteria.
    """

    def __init__(self, config: Dict = None):
        """
        Initialize the data processor with configuration.

        Args:
            config: Dictionary containing file paths and settings
        """
        self.config = config or self._get_default_config()
        self.raw_data = None
        self.processed_data = None
        self.exclusion_stats = {}

    def _get_default_config(self) -> Dict:
        """Get default configuration for file paths and settings."""
        return {
            'input_file': "/home/uqahonne/uq/ukb/data/bed_data/updated_13052025/updated_bed_db_30052025.csv",
            'output_dir': "/home/uqahonne/uq/ukb/scripts_updated",
            'figures_dir': "/home/uqahonne/uq/ukb/scripts_updated/figures",
            'dmri_filtered_file': "filtered_participants_dmri_bed_18092025.csv",
            'bed_filtered_file': "filtered_inhouse_BED_only_dmri_18092025.csv"
        }

    def load_data(self, file_path: str = None) -> pd.DataFrame:
        """
        Load the UK Biobank dataset.

        Args:
            file_path: Path to the CSV file (optional, uses config if not provided)

        Returns:
            Loaded DataFrame
        """
        file_path = file_path or self.config['input_file']

        print(f"📂 Loading data from: {file_path}")
        try:
            self.raw_data = pd.read_csv(file_path)
            print(f"✅ Successfully loaded {len(self.raw_data):,} participants")
            print(f"📊 Dataset shape: {self.raw_data.shape}")
            print(f"📋 Columns: {len(self.raw_data.columns)}")
            return self.raw_data
        except Exception as e:
            print(f"❌ Error loading data: {e}")
            raise

    def filter_dmri_participants(self) -> pd.DataFrame:
        """
        Filter participants to include only those with dMRI data available.

        Returns:
            Filtered DataFrame with dMRI participants only
        """
        if self.raw_data is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        print("\n🔍 FILTERING FOR dMRI AVAILABILITY")
        print("=" * 40)

        # Check column names
        print(f"Available columns: {len(self.raw_data.columns)}")
        second_col = self.raw_data.columns[1]
        print(f"Using column for dMRI filtering: {second_col}")

        # Filter for dMRI availability
        initial_count = len(self.raw_data)
        dmri_filtered = self.raw_data[
            self.raw_data[second_col].notna() &
            (self.raw_data[second_col].astype(str).str.strip() != "")
        ].copy()

        # Save intermediate result
        output_path = os.path.join(self.config['output_dir'], self.config['dmri_filtered_file'])
        dmri_filtered.to_csv(output_path, index=False)

        print(f"Initial participants: {initial_count:,}")
        print(f"With dMRI data: {len(dmri_filtered):,}")
        print(f"Exclusion rate: {(initial_count - len(dmri_filtered))/initial_count*100:.1f}%")
        print(f"💾 Saved to: {output_path}")

        self.processed_data = dmri_filtered
        return dmri_filtered


class ExclusionCriteriaHandler:
    """
    Handles all exclusion criteria for the study including substance use disorders.
    """

    @staticmethod
    def has_substance_use_disorder(row: pd.Series) -> bool:
        """
        Check if participant should be excluded due to substance use disorders.

        Args:
            row: Participant data row

        Returns:
            True if participant should be EXCLUDED
        """
        # Criterion 1: participant.p20457 coded as 1
        if row.get("participant.p20457") == 1:
            return True

        # Criterion 2: participant.p41202 contains ICD10 codes F10-F19
        icd_codes = row.get("participant.p41202")
        if pd.notna(icd_codes):
            # Handle both string and list formats
            if isinstance(icd_codes, str):
                if icd_codes.startswith('['):
                    try:
                        icd_codes = ast.literal_eval(icd_codes)
                    except:
                        return False
                else:
                    icd_codes = [icd_codes]

            if isinstance(icd_codes, list):
                for code in icd_codes:
                    if isinstance(code, str) and code.startswith('F1') and len(code) >= 3:
                        # Check if it's F10, F11, F12, ..., F19
                        try:
                            f_number = int(code[1:3])
                            if 10 <= f_number <= 19:
                                return True
                        except:
                            continue

        return False

    def apply_exclusions(self, df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Apply all exclusion criteria to the dataset.

        Args:
            df: Input DataFrame
            verbose: Whether to print detailed information

        Returns:
            DataFrame after exclusions
        """
        if verbose:
            print("\n🚫 APPLYING EXCLUSION CRITERIA")
            print("=" * 40)

        initial_count = len(df)

        # Apply substance use disorder exclusion
        excluded_mask = df.apply(self.has_substance_use_disorder, axis=1)
        excluded_count = excluded_mask.sum()
        df_filtered = df[~excluded_mask].copy()

        if verbose:
            print(f"Substance use disorder exclusions:")
            print(f"  Initial participants: {initial_count:,}")
            print(f"  Excluded: {excluded_count:,}")
            print(f"  Remaining: {len(df_filtered):,}")
            print(f"  Exclusion rate: {excluded_count/initial_count*100:.1f}%")

        return df_filtered


class GroupAssigner:
    """
    Handles group assignment for eating disorder classifications.
    """

    def __init__(self):
        """Initialize with eating disorder classification criteria."""
        self.criteria = {
            'BED': {'p29000_code': 19, 'description': 'Binge Eating Disorder'},
            'AN': {'p29000_code': 17, 'description': 'Anorexia Nervosa'},
            'BN': {'p29000_code': 18, 'description': 'Bulimia Nervosa'},
            'other_eating_disorder': {'p29000_code': 20, 'description': 'Other Eating Disorder'}
        }

    def _prepare_list_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare list columns for proper evaluation.

        Args:
            df: Input DataFrame

        Returns:
            DataFrame with properly formatted list columns
        """
        df = df.copy()

        # Convert p29000 to list format
        df["participant.p29000"] = df["participant.p29000"].apply(
            lambda x: ast.literal_eval(x) if pd.notna(x) else []
        )

        # Convert p29140 to list format
        df["participant.p29140"] = df["participant.p29140"].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith("[") else x
        )

        return df

    def assign_group(self, row: pd.Series) -> Optional[str]:
        """
        Assign group based on eating disorder criteria.

        Args:
            row: Participant data row

        Returns:
            Group assignment string or None
        """
        p29000_codes = row.get("participant.p29000", [])

        # Check for specific eating disorders
        if 19 in p29000_codes:
            return "BED"
        elif 17 in p29000_codes:
            return "AN"
        elif 18 in p29000_codes:
            return "BN"
        elif 20 in p29000_codes:
            return "other_eating_disorder"

        # Check for filtered criteria (subclinical BED)
        criteria_met = (
            row.get("participant.p29132") in [1, 2] and
            row.get("participant.p29137") == 1 and
            row.get("participant.p29134") == -4 and
            row.get("participant.p29135") == 0 and
            row.get("participant.p29133") == 0 and
            row.get("participant.p29140") == [0] and
            row.get("participant.p29144") == 0
        )

        if criteria_met:
            return "Filtered"

        return None

    def process_groups(self, df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
        """
        Process group assignments for the entire dataset.

        Args:
            df: Input DataFrame
            verbose: Whether to print detailed information

        Returns:
            DataFrame with group assignments
        """
        if verbose:
            print("\n🏷️  ASSIGNING GROUPS")
            print("=" * 40)

        # Prepare list columns
        df = self._prepare_list_columns(df)

        # Assign groups
        df["Group"] = df.apply(self.assign_group, axis=1)

        # Remove 'other_eating_disorder' group
        df = df[df["Group"] != "other_eating_disorder"]

        # Print group statistics
        if verbose:
            group_counts = df['Group'].value_counts()
            print("Group assignments:")
            for group, count in group_counts.items():
                if group in self.criteria:
                    desc = self.criteria[group]['description']
                    print(f"  {group}: {count:,} participants ({desc})")
                else:
                    print(f"  {group}: {count:,} participants")

            print(f"\nTotal participants with group assignments: {len(df):,}")

        return df


class TrajectoryAnalyzer:
    """
    Handles trajectory analysis and visualization for eating disorder progression.
    """

    def __init__(self, figures_dir: str):
        """
        Initialize trajectory analyzer.

        Args:
            figures_dir: Directory to save figures
        """
        self.figures_dir = Path(figures_dir)
        self.figures_dir.mkdir(exist_ok=True)

        # Color scheme for groups
        self.colors = {
            'BED': '#d62728',
            'AN': '#2ca02c',
            'BN': '#1f77b4',
            'Filtered': '#ff7f0e'
        }

    def prepare_trajectory_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Prepare data for trajectory analysis.

        Args:
            df: Input DataFrame

        Returns:
            DataFrame ready for trajectory analysis
        """
        trajectory_data = df[
            df["Group"].notna() &
            df["participant.p29138"].isin(range(0, 100)) &  # Age at onset
            df["participant.p29139"].isin(range(0, 100)) &  # Age at last episode
            df["participant.p21003_i2"].isin(range(0, 100))  # Current age
        ].copy()

        print(f"📊 Trajectory analysis sample: {len(trajectory_data):,} participants")
        return trajectory_data

    def create_boxplot_visualization(self, trajectory_data: pd.DataFrame) -> None:
        """
        Create boxplot visualization of age trajectories.

        Args:
            trajectory_data: Prepared trajectory data
        """
        # Reshape for plotting
        trajectory_long = trajectory_data.melt(
            id_vars="Group",
            value_vars=["participant.p29138", "participant.p29139", "participant.p21003_i2"],
            var_name="Stage",
            value_name="Age"
        )

        # Rename stages
        stage_mapping = {
            "participant.p29138": "Onset",
            "participant.p29139": "Last Binge",
            "participant.p21003_i2": "Current Age at Assessment"
        }
        trajectory_long["Stage"] = trajectory_long["Stage"].map(stage_mapping)

        # Create plot
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=trajectory_long, x="Stage", y="Age", hue="Group")
        plt.title("Trajectory: Onset, Last Binge, and Current Age by Diagnosis Group")
        plt.ylabel("Age (years)")
        plt.tight_layout()

        output_path = self.figures_dir / "age_onset_allgroups.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"💾 Boxplot saved to: {output_path}")
        plt.show()

    def create_individual_trajectories(self, trajectory_data: pd.DataFrame) -> None:
        """
        Create individual trajectory/spaghetti plot.

        Args:
            trajectory_data: Prepared trajectory data
        """
        plt.figure(figsize=(14, 8))

        for group in trajectory_data['Group'].dropna().unique():
            group_data = trajectory_data[trajectory_data['Group'] == group]
            group_color = self.colors.get(group, '#7f7f7f')

            first_line = True
            for _, participant in group_data.iterrows():
                ages = [
                    participant['participant.p29138'],
                    participant['participant.p29139'],
                    participant['participant.p21003_i2']
                ]
                stages = [0, 1, 2]

                if first_line:
                    plt.plot(stages, ages, color=group_color, alpha=0.3,
                            linewidth=1, label=group)
                    first_line = False
                else:
                    plt.plot(stages, ages, color=group_color, alpha=0.3, linewidth=1)

        plt.xticks([0, 1, 2], ['Onset', 'Last Binge', 'Current Age'])
        plt.ylabel('Age (years)')
        plt.xlabel('Disease Stage')
        plt.title('Individual Disease Trajectories by Eating Disorder Group')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        output_path = self.figures_dir / "individual_trajectories.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"💾 Individual trajectories saved to: {output_path}")
        plt.show()

    def create_duration_analysis(self, trajectory_data: pd.DataFrame) -> None:
        """
        Create duration analysis visualization.

        Args:
            trajectory_data: Prepared trajectory data
        """
        # Calculate duration metrics
        duration_data = trajectory_data.copy()
        duration_data['Duration_Active'] = (
            duration_data['participant.p29139'] - duration_data['participant.p29138']
        )
        duration_data['Years_Since_Last'] = (
            duration_data['participant.p21003_i2'] - duration_data['participant.p29139']
        )

        # Remove negative durations (data quality issues)
        duration_data = duration_data[
            (duration_data['Duration_Active'] >= 0) &
            (duration_data['Years_Since_Last'] >= 0)
        ]

        # Create plots
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Duration of active eating disorder
        sns.boxplot(data=duration_data, x='Group', y='Duration_Active', ax=axes[0])
        axes[0].set_title('Duration of Active Eating Disorder')
        axes[0].set_ylabel('Years')
        axes[0].tick_params(axis='x', rotation=45)

        # Years since last episode
        sns.boxplot(data=duration_data, x='Group', y='Years_Since_Last', ax=axes[1])
        axes[1].set_title('Years Since Last Episode')
        axes[1].set_ylabel('Years')
        axes[1].tick_params(axis='x', rotation=45)

        plt.tight_layout()

        output_path = self.figures_dir / "duration_analysis.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"💾 Duration analysis saved to: {output_path}")
        plt.show()

        # Print summary statistics
        self._print_duration_summary(duration_data)

    def _print_duration_summary(self, duration_data: pd.DataFrame) -> None:
        """
        Print summary statistics for duration analysis.

        Args:
            duration_data: Duration analysis data
        """
        print("\n📊 DURATION ANALYSIS SUMMARY")
        print("=" * 40)

        for group in duration_data['Group'].dropna().unique():
            group_data = duration_data[duration_data['Group'] == group]
            print(f"\n{group} (n={len(group_data):,}):")
            print(f"  Duration Active - Mean: {group_data['Duration_Active'].mean():.1f} years, "
                  f"Median: {group_data['Duration_Active'].median():.1f} years")
            print(f"  Years Since Last - Mean: {group_data['Years_Since_Last'].mean():.1f} years, "
                  f"Median: {group_data['Years_Since_Last'].median():.1f} years")


class BEDDataOrganizer:
    """
    Main orchestrator class that coordinates all processing steps.
    """

    def __init__(self, config: Dict = None):
        """
        Initialize the BED data organizer.

        Args:
            config: Configuration dictionary
        """
        # Initialize data processor to get default config if none provided
        self.data_processor = UKBiobankDataProcessor(config)
        self.config = self.data_processor.config  # Use the processor's config (includes defaults)

        self.exclusion_handler = ExclusionCriteriaHandler()
        self.group_assigner = GroupAssigner()
        self.trajectory_analyzer = TrajectoryAnalyzer(self.config['figures_dir'])

    def run_complete_analysis(self, save_intermediate: bool = True) -> pd.DataFrame:
        """
        Run the complete analysis pipeline.

        Args:
            save_intermediate: Whether to save intermediate files

        Returns:
            Final processed DataFrame
        """
        print("🚀 STARTING UK BIOBANK BED DATA ANALYSIS")
        print("=" * 60)

        # Step 1: Load and filter for dMRI
        self.data_processor.load_data()
        dmri_data = self.data_processor.filter_dmri_participants()

        # Step 2: Apply exclusion criteria
        excluded_data = self.exclusion_handler.apply_exclusions(dmri_data)

        # Step 3: Assign groups
        grouped_data = self.group_assigner.process_groups(excluded_data)

        # Step 4: Save filtered BED data
        if save_intermediate:
            filtered_bed = grouped_data[grouped_data["Group"] == "Filtered"]
            output_path = os.path.join(
                self.config['output_dir'],
                self.config['bed_filtered_file']
            )
            filtered_bed.to_csv(output_path, index=False)
            print(f"\n💾 Filtered BED data saved to: {output_path}")
            print(f"   Contains {len(filtered_bed):,} participants")

        # Step 5: Trajectory analysis
        print(f"\n📈 TRAJECTORY ANALYSIS")
        print("=" * 40)
        trajectory_data = self.trajectory_analyzer.prepare_trajectory_data(grouped_data)

        if len(trajectory_data) > 0:
            self.trajectory_analyzer.create_boxplot_visualization(trajectory_data)
            self.trajectory_analyzer.create_individual_trajectories(trajectory_data)
            self.trajectory_analyzer.create_duration_analysis(trajectory_data)
        else:
            print("⚠️  No valid trajectory data found for visualization")

        # Step 6: Final summary
        self._print_final_summary(grouped_data)

        return grouped_data

    def _print_final_summary(self, final_data: pd.DataFrame) -> None:
        """
        Print final analysis summary.

        Args:
            final_data: Final processed data
        """
        print(f"\n🎯 ANALYSIS COMPLETE")
        print("=" * 60)
        print("FOR INTERACTIVE CONTROL SELECTION:")
        print("python3 interactive_control_matching.py")
        print("\nNext steps:")
        print("- Review generated visualizations in figures/ directory")
        print("- Use interactive_control_matching.py for control selection")
        print("- Proceed with statistical analysis")
        print("=" * 60)


def main():
    """
    Main function to run the BED data organization analysis.
    """
    try:
        # Initialize and run analysis
        organizer = BEDDataOrganizer()
        final_data = organizer.run_complete_analysis()

        print(f"\n✅ Analysis completed successfully!")
        print(f"Final dataset contains {len(final_data):,} participants")

    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        raise


if __name__ == "__main__":
    main()