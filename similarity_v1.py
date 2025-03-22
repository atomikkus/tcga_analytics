import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cosine, jaccard
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

class PatientSimilarity:
    def __init__(self, clinical_weight=0.5, mutation_weight=0.5):
        self.clinical_weight = clinical_weight
        self.mutation_weight = mutation_weight
        self.scaler = StandardScaler()
        
    def prepare_clinical_features(self, clinical_df):
        """Prepare clinical features for similarity calculation."""
        # Select relevant clinical features
        clinical_features = [
            'AGE',
            'SEX',
            'RACE',
            'ETHNICITY',
            'ER_STATUS_BY_IHC',
            'PR_STATUS_BY_IHC',
            'OS_MONTHS',
            'OS_STATUS'
        ]
        
        # Create a copy of the dataframe with selected features
        df = clinical_df[clinical_features].copy()
        
        # Handle missing values
        df['AGE'] = pd.to_numeric(df['AGE'], errors='coerce')
        df['OS_MONTHS'] = pd.to_numeric(df['OS_MONTHS'], errors='coerce')
        
        # Fill missing values
        df['AGE'].fillna(df['AGE'].median(), inplace=True)
        df['OS_MONTHS'].fillna(df['OS_MONTHS'].median(), inplace=True)
        
        # One-hot encode categorical variables
        categorical_features = ['SEX', 'RACE', 'ETHNICITY', 'ER_STATUS_BY_IHC', 'PR_STATUS_BY_IHC', 'OS_STATUS']
        df_encoded = pd.get_dummies(df, columns=categorical_features)
        
        # Scale numerical features
        numerical_features = ['AGE', 'OS_MONTHS']
        df_encoded[numerical_features] = self.scaler.fit_transform(df_encoded[numerical_features])
        
        return df_encoded
    
    def prepare_mutation_features(self, mutations_df):
        """Create a binary mutation profile for each patient."""
        # Create gene mutation profile
        mutation_profile = pd.crosstab(
            mutations_df['patientId'], 
            mutations_df['hugoGeneSymbol']
        ).astype(bool).astype(int)
        
        return mutation_profile
    
    def fit(self, clinical_df, mutations_df):
        """Prepare both clinical and mutation data."""
        self.clinical_features = self.prepare_clinical_features(clinical_df)
        self.mutation_features = self.prepare_mutation_features(mutations_df)
        
        # Ensure we have the same patients in both datasets
        common_patients = list(set(self.clinical_features.index) & set(self.mutation_features.index))
        self.clinical_features = self.clinical_features.loc[common_patients]
        self.mutation_features = self.mutation_features.loc[common_patients]
        
        self.patient_ids = common_patients
    
    def calculate_similarity(self, patient_id_1, patient_id_2):
        """Calculate similarity between two patients."""
        # Clinical similarity (cosine similarity)
        clinical_sim = 1 - cosine(
            self.clinical_features.loc[patient_id_1],
            self.clinical_features.loc[patient_id_2]
        )
        
        # Mutation similarity (Jaccard similarity)
        mutation_sim = 1 - jaccard(
            self.mutation_features.loc[patient_id_1],
            self.mutation_features.loc[patient_id_2]
        )
        
        # Weighted combination
        total_similarity = (
            self.clinical_weight * clinical_sim +
            self.mutation_weight * mutation_sim
        )
        
        return {
            'total_similarity': total_similarity,
            'clinical_similarity': clinical_sim,
            'mutation_similarity': mutation_sim
        }
    
    def find_similar_patients(self, patient_id, n_matches=5):
        """Find the most similar patients to a given patient."""
        if patient_id not in self.patient_ids:
            raise ValueError(f"Patient {patient_id} not found in the dataset")
        
        similarities = []
        for other_id in self.patient_ids:
            if other_id != patient_id:
                sim = self.calculate_similarity(patient_id, other_id)
                similarities.append((other_id, sim))
        
        # Sort by total similarity (descending)
        similarities.sort(key=lambda x: x[1]['total_similarity'], reverse=True)
        return similarities[:n_matches]

def main():
    st.title("Patient Similarity Analysis")
    
    # Load data
    try:
        clinical_df = pd.read_csv('data/brca_tcga_pub2015_clinical_wide.csv', index_col='patientId')
        mutations_df = pd.read_csv('data/brca_tcga_pub2015_mutations.csv')
        
        # Initialize similarity model
        similarity_model = PatientSimilarity(clinical_weight=0.5, mutation_weight=0.5)
        similarity_model.fit(clinical_df, mutations_df)
        
        # Sidebar controls
        st.sidebar.header("Settings")
        clinical_weight = st.sidebar.slider(
            "Clinical Data Weight",
            min_value=0.0,
            max_value=1.0,
            value=0.5,
            step=0.1
        )
        mutation_weight = 1 - clinical_weight
        st.sidebar.write(f"Mutation Data Weight: {mutation_weight:.1f}")
        
        # Patient selection
        selected_patient = st.selectbox(
            "Select a patient to analyze",
            options=similarity_model.patient_ids
        )
        
        n_matches = st.slider(
            "Number of similar patients to find",
            min_value=1,
            max_value=20,
            value=5
        )
        
        if st.button("Find Similar Patients"):
            # Update weights
            similarity_model.clinical_weight = clinical_weight
            similarity_model.mutation_weight = mutation_weight
            
            # Display query patient information
            st.header("Query Patient Information")
            
            # Clinical information for query patient
            st.subheader("Clinical Information")
            query_clinical = clinical_df.loc[selected_patient]
            col1, col2 = st.columns(2)
            
            with col1:
                st.write("Demographics:")
                st.write(f"- Age: {query_clinical['AGE']}")
                st.write(f"- Sex: {query_clinical['SEX']}")
                st.write(f"- Race: {query_clinical['RACE']}")
                st.write(f"- Ethnicity: {query_clinical['ETHNICITY']}")
            
            with col2:
                st.write("Clinical Status:")
                st.write(f"- ER Status: {query_clinical['ER_STATUS_BY_IHC']}")
                st.write(f"- PR Status: {query_clinical['PR_STATUS_BY_IHC']}")
                st.write(f"- Overall Survival: {query_clinical['OS_MONTHS']:.1f} months")
                st.write(f"- Survival Status: {query_clinical['OS_STATUS']}")
            
            # Mutation information for query patient
            st.subheader("Mutation Profile")
            query_mutations = mutations_df[mutations_df['patientId'] == selected_patient]
            
            # Create mutation summary
            mutation_summary = query_mutations.groupby('mutationType').size().reset_index()
            mutation_summary.columns = ['Mutation Type', 'Count']
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.write(f"Total Mutations: {len(query_mutations)}")
                st.write("Mutation Types:")
                st.dataframe(mutation_summary, hide_index=True)
            
            with col2:
                # Create pie chart of mutation types
                fig_mutation_types = px.pie(
                    mutation_summary,
                    values='Count',
                    names='Mutation Type',
                    title='Distribution of Mutation Types'
                )
                fig_mutation_types.update_layout(height=300)
                st.plotly_chart(fig_mutation_types, use_container_width=True)
            
            # Show all mutations in an expandable section
            with st.expander("View All Mutations"):
                st.dataframe(
                    query_mutations[['hugoGeneSymbol', 'mutationType', 'proteinChange']].sort_values('hugoGeneSymbol'),
                    hide_index=True
                )
            
            st.markdown("---")  # Add a divider
            
            # Find similar patients
            similar_patients = similarity_model.find_similar_patients(
                selected_patient,
                n_matches=n_matches
            )
            
            # Display results
            st.header("Similar Patients")
            
            # Create visualization
            similarity_data = {
                'Patient ID': [p[0] for p in similar_patients],
                'Total Similarity': [p[1]['total_similarity'] for p in similar_patients],
                'Clinical Similarity': [p[1]['clinical_similarity'] for p in similar_patients],
                'Mutation Similarity': [p[1]['mutation_similarity'] for p in similar_patients]
            }
            
            df_similarities = pd.DataFrame(similarity_data)
            
            # Plot similarities
            fig = go.Figure()
            
            # Add traces for each similarity type
            fig.add_trace(go.Bar(
                name='Clinical Similarity',
                x=df_similarities['Patient ID'],
                y=df_similarities['Clinical Similarity'],
                marker_color='rgb(55, 83, 109)'
            ))
            
            fig.add_trace(go.Bar(
                name='Mutation Similarity',
                x=df_similarities['Patient ID'],
                y=df_similarities['Mutation Similarity'],
                marker_color='rgb(26, 118, 255)'
            ))
            
            fig.add_trace(go.Scatter(
                name='Total Similarity',
                x=df_similarities['Patient ID'],
                y=df_similarities['Total Similarity'],
                mode='lines+markers',
                line=dict(color='rgb(219, 64, 82)', width=2)
            ))
            
            fig.update_layout(
                title=f'Similarity Scores for Patient {selected_patient}',
                xaxis_title='Patient ID',
                yaxis_title='Similarity Score',
                barmode='group',
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Display detailed information
            st.subheader("Detailed Patient Information")
            
            # Create comparison table
            comparison_data = clinical_df.loc[[selected_patient] + list(df_similarities['Patient ID'])]
            st.write("Clinical Comparison:")
            st.dataframe(comparison_data[['AGE', 'SEX', 'RACE', 'ETHNICITY', 'ER_STATUS_BY_IHC', 'PR_STATUS_BY_IHC', 'OS_MONTHS', 'OS_STATUS']])
            
            # Show mutation overlap
            st.write("Common Mutations:")
            for patient_id in df_similarities['Patient ID']:
                common_mutations = set(
                    mutations_df[mutations_df['patientId'] == selected_patient]['hugoGeneSymbol']
                ) & set(
                    mutations_df[mutations_df['patientId'] == patient_id]['hugoGeneSymbol']
                )
                st.write(f"Patient {patient_id}: {len(common_mutations)} common mutations")
                if common_mutations:
                    st.write(", ".join(sorted(common_mutations)))
    except FileNotFoundError as e:
        st.error("Error: Could not find the data files. Please ensure both clinical and mutation data files are present in the data directory.")
        st.stop()
    except Exception as e:
        st.error(f"Error loading or processing data: {str(e)}")
        st.stop()

if __name__ == "__main__":
    main() 