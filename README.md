# UK Biobank Binge Eating Disorder (BED) Analysis Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive pipeline for analyzing Binge Eating Disorder (BED) data from UK Biobank with advanced control matching and neuroimaging data preparation.

## Overview

This project provides a complete workflow for:
- Processing UK Biobank eating disorder data with diffusion MRI filtering
- Implementing sophisticated control matching algorithms
- Applying evidence-based exclusion criteria for substance use and psychiatric comorbidities
- Generating comprehensive visualizations and quality assessments
- Creating publication-ready matched case-control datasets

##  Key Features

### Data Processing
- **Automated dMRI filtering**: Selects participants with available diffusion MRI data
- **Multi-level group assignment**: Clinical BED, AN, BN, and subclinical BED cases
- **Robust exclusion criteria**: Substance use disorders, psychiatric comorbidities
- **Data quality assurance**: Safe parsing of UK Biobank array fields

### Advanced Control Matching
- **5 Matching Algorithms**: Exact, Propensity Score, Nearest Neighbor, Optimal, Flexible
- **Logit-based Propensity Scores**: Austin (2011) standard with 0.2×std caliper
- **Quality Assessment**: Standardized Mean Difference (SMD), balance metrics
- **Interactive Selection**: User-guided method selection with performance comparison

### 📈 Visualization & Analysis
- **Propensity Score Distributions**: Before/after matching comparisons
- **Disease Trajectory Analysis**: Age onset, duration, recovery patterns
- **Quality Dashboards**: Comprehensive matching assessment
- **Publication-ready Figures**: High-resolution plots with professional styling

## Quick Start

### Prerequisites
```bash
pip install pandas numpy scikit-learn scipy matplotlib seaborn
```

### Basic Usage

#### 1. Process Fresh UK Biobank Data
```python
from data_org_BED_refactored import BEDDataOrganizer

# Initialize and run complete analysis
organizer = BEDDataOrganizer()
processed_data = organizer.run_complete_analysis()
```

#### 2. Interactive Control Matching
```python
python3 interactive_control_matching_refactored.py
```

#### 3. Quick Field Lookup
```python
python3 web_field_lookup.py 29000
```

## 📁 Project Structure

```
ukb-bed-analysis/
├── data_org_BED_refactored.py              # Main data processing pipeline
├── interactive_control_matching.py   # Interactive matching system
├── control_selection_pattern_matching.py    # Advanced exclusion criteria
├── web_field_lookup.py                     # UK Biobank field descriptions
├── matching_visualization_dashboard.py      # Visualization system
├── summary_BED_ukbiobank_19092025.md       # Methodology documentation
└── README.md                               # This file
```

## Methodology

### Sample Selection Criteria

#### Inclusion Criteria
- **dMRI Data Available**: `participant.p20218_i2` (Multiband diffusion brain images)
- **Case Definitions**:
  - **Clinical BED**: `participant.p29000` code 19
  - **Clinical AN**: `participant.p29000` code 17
  - **Clinical BN**: `participant.p29000` code 18
  - **Subclinical BED**: Meeting all 7 specific criteria from mental health questionnaire

#### Exclusion Criteria (Hierarchical)

**MINIMAL (Essential)**:
- Substance use disorders (`p20457=1`, `p41202` F10-F19)
- Binge eating disorders in controls (`p20544` codes 13,16, `p41202` F50)
- Eating disorders in controls (`p29000` codes 17-20)

**MODERATE (+Additional)**:
- Severe depression/bipolar (`p20126`)

**STRICT (+Comprehensive)**:
- High mood symptoms (`p1930`, `p1980`)

### Sample Composition
From **502,128** initial participants → **82,342** with dMRI → **81,876** final sample:
- **AN**: 306 participants
- **BN**: 115 participants
- **BED**: 106 participants
- **Subclinical BED**: 125 participants
- **Controls**: 81,224 participants

##  Matching Algorithms

### 1. Propensity Score Matching (Recommended)
- **Logit transformation**: Austin (2011) standard
- **Caliper**: 0.2 × std(logit propensity scores)
- **Convergence**: 5000 max iterations with monitoring
- **Quality metrics**: SMD < 0.25 for good balance

### 2. Exact Matching
- Direct demographic matching with optional tolerance
- Best for: Clear matching criteria, small datasets

### 3. Nearest Neighbor Matching
- Euclidean distance with standardized features
- Best for: Continuous variables, medium datasets

### 4. Optimal Matching
- Hungarian algorithm for globally optimal 1:1 pairs
- Best for: When global optimization is priority

### 5. Flexible Matching
- Multi-round with progressively relaxed criteria
- Best for: Maximizing sample retention

##  Usage Examples

### Process Data and Run Matching
```python
# Complete pipeline
from data_org_BED_refactored import BEDDataOrganizer
from interactive_control_matching_refactored import InteractiveControlMatcher

# Process data
organizer = BEDDataOrganizer()
df = organizer.run_complete_analysis()

# Interactive matching
matcher = InteractiveControlMatcher(df, 'Filtered')
results = matcher.run_interactive_matching()
```

### Programmatic Matching
```python
# Direct propensity score matching
matcher = InteractiveControlMatcher(df, 'Filtered')
matcher.prepare_control_pool(exclusion_level='minimal')

# Get methods manager
from interactive_control_matching_refactored import MatchingMethodsManager
methods = MatchingMethodsManager(matcher.bed_cases, matcher.prepared_pool)

# Run propensity score matching
controls = methods.run_propensity_matching(
    matching_vars=['participant.p21003_i2', 'participant.p31'],
    ratio=2,
    caliper=0.1,
    verbose=True
)
```

### Visualization
```python
# Propensity score distributions
matcher.visualize_propensity_scores(
    cases_with_ps,
    controls_with_ps,
    matched_controls,
    caliper_logit
)
```

## 📊 Output Files

- **Matched datasets**: `matched_[method]_[date].csv`
- **Controls only**: `controls_[method]_[date].csv`
- **Visualizations**: `figures/` directory with high-resolution plots
- **Quality reports**: Console output with detailed metrics

## 🔧 Configuration

### Default Paths
```python
config = {
    'input_file': "/path/to/uk_biobank_data.csv",
    'output_dir': "/path/to/output",
    'figures_dir': "/path/to/figures"
}
```

### Matching Parameters
```python
preferences = {
    'ratio': 2,                    # 2:1 control:case ratio
    'exclusion_level': 'minimal',  # Exclusion stringency
    'matching_vars': ['participant.p21003_i2', 'participant.p31']  # Age + Sex
}
```

## 📚 Documentation

- **Methodology**: See `summary_BED_ukbiobank_19092025.md` for complete methodology
- **Field Descriptions**: Use `web_field_lookup.py` for UK Biobank field details
- **API Reference**: Comprehensive docstrings throughout codebase

## 🛠️ Troubleshooting

### Common Issues

**Array Comparison Errors**:
```python
# Fixed with safe parsing in refactored version
df["participant.p29000"] = df["participant.p29000"].apply(safe_eval)
```

**Convergence Issues**:
```python
# Increased max_iter based on literature
LogisticRegression(max_iter=5000, solver='lbfgs')
```

**Memory Issues**:
```python
# Chunked processing for large datasets
chunk_size = 5000  # Adjustable based on available memory
```

## Quality Metrics

### Matching Quality Assessment
- **Balance Score** (35%): Demographic balance between groups
- **Quality Score** (45%): Matching efficiency and ratio achievement
- **Efficiency Score** (20%): Sample utilization and computation time
- **Composite Score**: Weighted combination for method selection

### Balance Thresholds
- **SMD < 0.1**:  Excellent balance
- **SMD < 0.25**: Good balance
- **SMD ≥ 0.25**: Moderate balance



## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **UK Biobank**: Data access and research infrastructure
- **Austin (2011)**: Propensity score methodology standards
- **Stuart (2010)**: Matching methods framework
- **Research Team**: Methodology development and validation

##  Contact

For questions about methodology or implementation:
- Open an issue on GitHub
- Email: [your-email@institution.edu]

## 🔗 References

1. Austin, P.C. (2011). An introduction to propensity score methods for reducing the effects of confounding in observational studies. *Multivariate Behavioral Research*, 46(3), 399-424.

2. Stuart, E.A. (2010). Matching methods for causal inference: A review and a look forward. *Statistical Science*, 25(1), 1-21.

3. Brookhart, M.A., et al. (2013). Variable selection for propensity score models. *American Journal of Epidemiology*, 163(12), 1149-1156.

---

**Last Updated**: September 19, 2025
**Version**: 1.0.0
**Python**: 3.8+