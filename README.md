# TCGA BRCA Clinical Dashboard

An interactive dashboard for analyzing TCGA breast cancer clinical data using Streamlit.

## Features

- Interactive filtering by age, gender, and ethnicity
- Key patient metrics and statistics
- Visualizations for:
  - Age distribution
  - Gender distribution
  - Receptor status (ER, PR, HER2)
  - Race distribution
  - Kaplan-Meier survival analysis

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/tcga-analytics.git
cd tcga-analytics
```

2. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

1. Place your TCGA clinical data file in the `data/` directory
2. Run the dashboard:
```bash
streamlit run clinical_dashboard.py
```

3. Open your browser and navigate to `http://localhost:8501`

## Data Sources

This dashboard uses TCGA BRCA clinical data from the [cBioPortal](https://www.cbioportal.org/).

## License

MIT License