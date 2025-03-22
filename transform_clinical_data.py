import pandas as pd

# Read the clinical data
print("Reading clinical data...")
clinical_data = pd.read_csv('data/brca_tcga_pub2015_clinical.csv')

# Check for duplicates
print("\nChecking for duplicates...")
duplicates = clinical_data.groupby(['patientId', 'clinicalAttributeId']).size().reset_index(name='count')
duplicates = duplicates[duplicates['count'] > 1]
if len(duplicates) > 0:
    print(f"Found {len(duplicates)} duplicate patient-attribute combinations")
    print("\nExample duplicates:")
    print(duplicates.head())

# Keep first occurrence of each patient-attribute combination
print("\nRemoving duplicates (keeping first occurrence)...")
clinical_data_dedup = clinical_data.drop_duplicates(subset=['patientId', 'clinicalAttributeId'], keep='first')

# Convert from long to wide format using pivot
print("Transforming data to wide format...")
clinical_data_wide = clinical_data_dedup.pivot(
    index='patientId',
    columns='clinicalAttributeId',
    values='value'
)

# Reset index to make patientId a column
clinical_data_wide.reset_index(inplace=True)

# Save the transformed data
print("\nSaving transformed data...")
output_file = 'data/brca_tcga_pub2015_clinical_wide.csv'
clinical_data_wide.to_csv(output_file, index=False)

print(f"\nTransformation complete. Data saved to {output_file}")
print(f"Original shape: {clinical_data.shape}")
print(f"After deduplication: {clinical_data_dedup.shape}")
print(f"Transformed shape: {clinical_data_wide.shape}")

# Print some basic statistics about the transformed data
print("\nDataset statistics:")
print(f"Number of patients: {len(clinical_data_wide)}")
print(f"Number of clinical attributes: {len(clinical_data_wide.columns) - 1}")  # -1 for patientId column 