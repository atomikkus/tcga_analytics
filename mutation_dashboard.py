import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
import io

# Set page config
st.set_page_config(
    page_title="TCGA Mutation Analysis",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Add custom CSS
st.markdown("""
    <style>
    .main {
        padding: 0rem 1rem;
    }
    .stPlotlyChart {
        background-color: #ffffff;
        border-radius: 5px;
        padding: 10px;
        margin-bottom: 20px;
    }
    </style>
""", unsafe_allow_html=True)

# Common plot layout settings
plot_layout = dict(
    height=400,
    margin=dict(l=50, r=30, t=50, b=50),
    paper_bgcolor='black',
    plot_bgcolor='black',
    font=dict(size=12)
)

# Required columns for the mutation data
REQUIRED_COLUMNS = [
    'hugoGeneSymbol',
    'patientId',
    'mutationType',
    'proteinChange',
    'proteinPosStart'
]

# Load mutation data
@st.cache_data
def load_data(file_path=None, uploaded_file=None):
    try:
        if uploaded_file is not None:
            # Read the uploaded file
            df = pd.read_csv(uploaded_file)
        elif file_path is not None:
            # Read from local file path
            df = pd.read_csv(file_path)
        else:
            return None
        
        # Validate required columns
        missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
        if missing_cols:
            st.error(f"Missing required columns: {', '.join(missing_cols)}")
            return None
            
        return df
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return None

# File upload section in sidebar
st.sidebar.header('Data Input')
uploaded_file = st.sidebar.file_uploader(
    "Upload mutation data file (CSV)",
    type=['csv'],
    help="Upload a CSV file containing mutation data. Required columns: " + ", ".join(REQUIRED_COLUMNS)
)

# Load either uploaded file or default file
if uploaded_file is not None:
    mutations_df = load_data(uploaded_file=uploaded_file)
else:
    st.sidebar.write("No file uploaded. Using default TCGA BRCA mutation data.")
    mutations_df = load_data(file_path='data/brca_tcga_pub2015_mutations.csv')

if mutations_df is None:
    st.error("Failed to load mutation data. Please check the data file and try again.")
    st.stop()

# Display data summary
st.sidebar.write("Data Summary:")
st.sidebar.write(f"Total mutations: {len(mutations_df)}")
st.sidebar.write(f"Total genes: {mutations_df['hugoGeneSymbol'].nunique()}")
st.sidebar.write(f"Total patients: {mutations_df['patientId'].nunique()}")

# Get unique genes for selection
unique_genes = sorted(mutations_df['hugoGeneSymbol'].unique())

# Main dashboard
st.title('Mutation Analysis Dashboard')

# Sidebar for gene selection
st.sidebar.header('Gene Selection')

# Search box for genes
gene_search = st.sidebar.text_input(
    "Search genes",
    "",
    help="Type to search for specific genes"
)

# Filter genes based on search
if gene_search:
    filtered_genes = [gene for gene in unique_genes if gene_search.upper() in gene.upper()]
else:
    filtered_genes = unique_genes

# Multi-select for genes
selected_genes = st.sidebar.multiselect(
    'Select genes to analyze',
    options=filtered_genes,
    default=filtered_genes[:5] if len(filtered_genes) > 5 else filtered_genes
)

if not selected_genes:
    st.warning("Please select at least one gene to analyze.")
    st.stop()

# Filter data for selected genes
filtered_df = mutations_df[mutations_df['hugoGeneSymbol'].isin(selected_genes)]

# Main content
st.header('Mutation Analysis Results')

# Top metrics row
col1, col2, col3 = st.columns(3)

with col1:
    total_mutations = len(filtered_df)
    st.metric("Total Mutations", total_mutations)

with col2:
    unique_patients = filtered_df['patientId'].nunique()
    st.metric("Unique Patients", unique_patients)

with col3:
    mutations_per_patient = total_mutations / unique_patients if unique_patients > 0 else 0
    st.metric("Avg. Mutations per Patient", f"{mutations_per_patient:.2f}")

# First row of visualizations
col1, col2 = st.columns(2)

with col1:
    # Gene frequency plot
    gene_counts = filtered_df['hugoGeneSymbol'].value_counts()
    fig_gene_freq = px.bar(
        x=gene_counts.index,
        y=gene_counts.values,
        title="Mutation Frequency by Gene",
        labels={'x': 'Gene', 'y': 'Number of Mutations'}
    )
    fig_gene_freq.update_layout(**plot_layout)
    st.plotly_chart(fig_gene_freq, use_container_width=True)

with col2:
    # Mutation type distribution
    mutation_type_counts = filtered_df['mutationType'].value_counts()
    fig_mutation_types = px.pie(
        values=mutation_type_counts.values,
        names=mutation_type_counts.index,
        title="Distribution of Mutation Types"
    )
    fig_mutation_types.update_layout(
        **plot_layout,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    st.plotly_chart(fig_mutation_types, use_container_width=True)

# Second row - Top protein changes by gene
st.subheader("Top Protein Changes by Gene")

for gene in selected_genes:
    gene_data = filtered_df[filtered_df['hugoGeneSymbol'] == gene]
    
    if len(gene_data) > 0:
        st.write(f"### {gene}")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Top protein changes
            protein_changes = gene_data['proteinChange'].value_counts().head(10)
            fig_protein = px.bar(
                x=protein_changes.index,
                y=protein_changes.values,
                title=f"Top 10 Protein Changes in {gene}",
                labels={'x': 'Protein Change', 'y': 'Frequency'}
            )
            fig_protein.update_layout(
                **plot_layout,
                xaxis_tickangle=45
            )
            st.plotly_chart(fig_protein, use_container_width=True)
        
        with col2:
            # Mutation types for this gene
            mutation_types = gene_data['mutationType'].value_counts()
            fig_variants = px.pie(
                values=mutation_types.values,
                names=mutation_types.index,
                title=f"Mutation Types in {gene}"
            )
            fig_variants.update_layout(
                **plot_layout,
                legend=dict(
                    orientation="h",
                    yanchor="bottom",
                    y=1.02,
                    xanchor="right",
                    x=1
                )
            )
            st.plotly_chart(fig_variants, use_container_width=True)

# Additional analysis - Mutation position distribution
st.subheader("Mutation Position Distribution")

for gene in selected_genes:
    gene_data = filtered_df[filtered_df['hugoGeneSymbol'] == gene]
    
    if len(gene_data) > 0:
        # Convert position to numeric, handling any errors
        gene_data['proteinPosStart'] = pd.to_numeric(gene_data['proteinPosStart'], errors='coerce')
        
        # Create a count of mutations at each position
        position_counts = gene_data['proteinPosStart'].value_counts().reset_index()
        position_counts.columns = ['Position', 'Count']
        
        fig_position = px.scatter(
            position_counts,
            x='Position',
            y='Count',
            title=f"Mutation Position Distribution in {gene}",
            labels={'Position': 'Protein Position', 'Count': 'Number of Mutations'},
            size='Count',  # Size of points based on mutation count
            hover_data={'Position': True, 'Count': True}
        )
        
        fig_position.update_layout(
            **plot_layout,
            xaxis_title="Protein Position",
            yaxis_title="Number of Mutations",
            showlegend=False,
            hovermode='closest'
        )
        
        # Update marker properties
        fig_position.update_traces(
            marker=dict(
                sizeref=2.*max(position_counts['Count'])/(40.**2),  # Adjust point size scaling
                sizemin=4,  # Minimum point size
                opacity=0.7
            )
        )
        
        st.plotly_chart(fig_position, use_container_width=True)
        
        # Add a data table with hotspot positions (positions with multiple mutations)
        hotspots = position_counts[position_counts['Count'] > 1].sort_values('Count', ascending=False).head(10)
        if not hotspots.empty:
            st.write(f"### Top Mutation Hotspots in {gene}")
            st.write("Positions with multiple mutations:")
            st.dataframe(
                hotspots,
                column_config={
                    "Position": st.column_config.NumberColumn("Position", help="Protein position"),
                    "Count": st.column_config.NumberColumn("Mutations", help="Number of mutations at this position")
                },
                hide_index=True
            )

# Download section
st.sidebar.header("Download Data")
if st.sidebar.button("Download Selected Gene Data"):
    csv = filtered_df.to_csv(index=False)
    st.sidebar.download_button(
        label="Click to Download",
        data=csv,
        file_name="selected_mutations.csv",
        mime="text/csv"
    ) 