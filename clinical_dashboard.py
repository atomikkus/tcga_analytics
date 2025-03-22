import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from lifelines import KaplanMeierFitter

# Set page config
st.set_page_config(
    page_title="TCGA BRCA Clinical Dashboard",
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

# Load data
@st.cache_data
def load_data():
    try:
        df = pd.read_csv('data/brca_tcga_pub2015_clinical_wide.csv')
        st.write("Data loaded successfully!")
        st.write(f"Number of rows: {len(df)}")
        st.write("Columns:", df.columns.tolist())
        
        # Convert AGE to numeric
        df['AGE'] = pd.to_numeric(df['AGE'], errors='coerce')
        # Convert survival data
        df['OS_MONTHS'] = pd.to_numeric(df['OS_MONTHS'], errors='coerce')
        return df
    except Exception as e:
        st.error(f"Error loading data: {str(e)}")
        return None

# Load the data
df = load_data()

if df is None:
    st.error("Failed to load data. Please check the data file and try again.")
    st.stop()

# Sidebar filters
st.sidebar.header('Filters')

# Age range filter
age_range = st.sidebar.slider(
    'Age Range',
    float(df['AGE'].min()),
    float(df['AGE'].max()),
    (float(df['AGE'].min()), float(df['AGE'].max()))
)

# Gender filter
genders = ['All'] + sorted(df['SEX'].dropna().unique().tolist())
selected_gender = st.sidebar.selectbox('Gender', genders)

# Ethnicity filter
ethnicities = ['All'] + sorted(df['ETHNICITY'].dropna().unique().tolist())
selected_ethnicity = st.sidebar.selectbox('Ethnicity', ethnicities)

# Apply filters
mask = (df['AGE'].between(age_range[0], age_range[1]))
if selected_gender != 'All':
    mask &= (df['SEX'] == selected_gender)
if selected_ethnicity != 'All':
    mask &= (df['ETHNICITY'] == selected_ethnicity)

filtered_df = df[mask]

# Main dashboard
st.title('TCGA BRCA Clinical Dashboard')

# Add spacing after title
st.markdown("<br>", unsafe_allow_html=True)

# Top metrics
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric("Total Patients", len(filtered_df))
with col2:
    st.metric("Median Age", f"{filtered_df['AGE'].median():.1f}")
with col3:
    female_pct = (filtered_df['SEX'] == 'Female').mean() * 100
    st.metric("Female %", f"{female_pct:.1f}%")
with col4:
    survival_rate = (filtered_df['OS_STATUS'] == '0:LIVING').mean() * 100
    st.metric("Overall Survival Rate", f"{survival_rate:.1f}%")

# Add spacing after metrics
st.markdown("<br>", unsafe_allow_html=True)

# First row of visualizations
col1, col2 = st.columns(2)

with col1:
    st.subheader("Age Distribution")
    fig_age = px.histogram(
        filtered_df,
        x="AGE",
        nbins=30,
        color_discrete_sequence=['#1f77b4']
    )
    fig_age.update_layout(
        **plot_layout,
        bargap=0.1,
        xaxis_title="Age (years)",
        yaxis_title="Count",
        showlegend=False
    )
    st.plotly_chart(fig_age, use_container_width=True)

with col2:
    st.subheader("Gender Distribution")
    gender_counts = filtered_df['SEX'].value_counts()
    fig_gender = px.pie(
        values=gender_counts.values,
        names=gender_counts.index
    )
    fig_gender.update_layout(
        **plot_layout,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    st.plotly_chart(fig_gender, use_container_width=True)

# Add spacing between rows
st.markdown("<br>", unsafe_allow_html=True)

# Second row of visualizations
col1, col2 = st.columns(2)

with col1:
    st.subheader("Receptor Status")
    receptor_data = pd.DataFrame({
        'Receptor': ['ER', 'PR', 'HER2'],
        'Positive': [
            (filtered_df['ER_STATUS_BY_IHC'] == 'Positive').mean() * 100,
            (filtered_df['PR_STATUS_BY_IHC'] == 'Positive').mean() * 100,
            (filtered_df['HER2_FISH_STATUS'] == 'Positive').mean() * 100
        ]
    })
    fig_receptor = px.bar(
        receptor_data,
        x='Receptor',
        y='Positive',
        text=receptor_data['Positive'].round(1).astype(str) + '%'
    )
    fig_receptor.update_layout(
        **plot_layout,
        yaxis_title="Positive %",
        showlegend=False,
        xaxis_title="",
        yaxis=dict(range=[0, 100])
    )
    fig_receptor.update_traces(textposition='outside')
    st.plotly_chart(fig_receptor, use_container_width=True)

with col2:
    st.subheader("Race Distribution")
    race_counts = filtered_df['RACE'].value_counts()
    fig_race = px.pie(
        values=race_counts.values,
        names=race_counts.index
    )
    fig_race.update_layout(
        **plot_layout,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    st.plotly_chart(fig_race, use_container_width=True)

# Add spacing before survival analysis
st.markdown("<br>", unsafe_allow_html=True)

# Kaplan-Meier plot
st.subheader("Kaplan-Meier Survival Analysis")

# Prepare survival data
survival_df = df.dropna(subset=['OS_MONTHS', 'OS_STATUS'])
survival_df['OS_STATUS'] = survival_df['OS_STATUS'].map({'0:LIVING': 0, '1:DECEASED': 1})

kmf = KaplanMeierFitter()
kmf.fit(
    survival_df['OS_MONTHS'],
    survival_df['OS_STATUS'],
    label='Overall'
)

fig_km = go.Figure()
fig_km.add_trace(
    go.Scatter(
        x=kmf.timeline,
        y=kmf.survival_function_.values.flatten(),
        mode='lines',
        name='Overall Survival',
        line=dict(color='red', width=2)
    )
)

# Create a modified layout for the KM plot with greater height
km_plot_layout = plot_layout.copy()
km_plot_layout['height'] = 600  # Override the height for this specific plot

fig_km.update_layout(
    **km_plot_layout,
    xaxis_title="Months",
    yaxis_title="Survival Probability",
    showlegend=True,
    legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1
    ),
    yaxis=dict(range=[0, 1])
)

st.plotly_chart(fig_km, use_container_width=True) 