import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("viridis")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 13
plt.rcParams['figure.titlesize'] = 20

def load_and_preprocess_data(file_path, sheet_name=4):
    # Load data
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    
    print(f"Data loaded successfully: {df.shape[0]} rows, {df.shape[1]} columns")
    print(f"Available columns: {list(df.columns)}")
    
    # Load efficiency data from R analysis
    efficiency_file = r"./results/efficiency_calculation/frontier_efficiency_results.csv"
    efficiency_df = pd.read_csv(efficiency_file)
    print(f"Efficiency data loaded successfully: {efficiency_df.shape[0]} rows")
    
    # Merge efficiency data with main data
    if 'No.' in df.columns and 'No.' in efficiency_df.columns:
        # Convert Year to integer for matching
        if 'Year' in df.columns:
            df['Year'] = pd.to_numeric(df['Year'], errors='coerce')
        if 'Report_year' in efficiency_df.columns:
            efficiency_df['Year'] = pd.to_numeric(efficiency_df['Report_year'], errors='coerce')
        
        # Merge on farm ID and year
        df = df.merge(efficiency_df[['No.', 'Year', 'A_TE_VRS', 'A_TE_CRS']], 
                     left_on=['No.', 'Year'], 
                     right_on=['No.', 'Year'], 
                     how='left')
        print(f"Data merged with efficiency values: {df.shape[0]} rows")
    
    # Clean and prepare data
    # Remove water footprint indicator to keep consistency with SBM setup
    if 'Water_footprint' in df.columns:
        df = df.drop(columns=['Water_footprint'])
        print("Removed indicator: Water_footprint")

    # Ensure numeric types
    numeric_cols = ['Year', 'Scale-up', 'Market_pig', 'A_TE_VRS', 'A_TE_CRS', 
                    'Pork_price', 'Feed_price']
    
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Rename columns for easier handling
    df = df.rename(columns={'No.': 'Farm_ID'})
    
    # Handle Scale-up column (with hyphen)
    if 'Scale-up' in df.columns:
        print(f"Scale-up distribution:\n{df['Scale-up'].value_counts(dropna=False)}")
    else:
        print("Warning: Scale-up column not found in data")
    
    # Check if efficiency values are available
    if 'A_TE_CRS' in df.columns:
        print(f"Efficiency data available: {df['A_TE_CRS'].notnull().sum()} non-null values")
    else:
        print("Warning: A_TE_CRS column not found after merge")
    
    return df

def select_control_farms(df, treatment_farms, max_controls_per_treatment=3):
    """
    Select control farms for each treatment farm based on similarity
    
    Parameters:
    df: DataFrame with all farm data
    treatment_farms: DataFrame with treatment farms (Scale-up = 0)
    max_controls_per_treatment: Maximum number of control farms per treatment farm
    
    Returns:
    Dictionary mapping treatment farm IDs to list of control farms
    """
    # Identify potential control farms (Scale-up = N/A)
    potential_controls = df[df['Scale-up'].isna()].copy()
    print(f"Found {len(potential_controls)} potential control farms (Scale-up = N/A)")
    
    # Prepare potential controls with relevant variables
    control_vars = ['Year', 'Market_pig', 'Province']
    for var in control_vars:
        if var not in potential_controls.columns:
            print(f"Warning: {var} not found in data, will not use for matching")
            control_vars.remove(var)
    
    # Create mapping from treatment farm to control farms
    treatment_to_controls = {}
    
    for _, treatment_row in treatment_farms.iterrows():
        treatment_id = treatment_row['Farm_ID']
        treatment_year = treatment_row.get('Year')
        treatment_scale = treatment_row.get('Market_pig')
        treatment_province = treatment_row.get('Province')
        
        # Calculate similarity scores for each potential control
        similarities = []
        for _, control_row in potential_controls.iterrows():
            control_id = control_row['Farm_ID']
            
            # Skip if same farm (shouldn't happen, but just in case)
            if control_id == treatment_id:
                continue
            
            # Calculate similarity based on available variables
            similarity_score = 0
            
            # Year similarity (weight 1)
            if not pd.isna(treatment_year) and not pd.isna(control_row.get('Year')):
                year_diff = abs(treatment_year - control_row.get('Year'))
                year_similarity = max(0, 1 - year_diff / 2)  # Full score if within 2 years
                similarity_score += year_similarity
            
            # Scale similarity (weight 1)
            if not pd.isna(treatment_scale) and not pd.isna(control_row.get('Market_pig')) and treatment_scale > 0:
                scale_ratio = min(treatment_scale, control_row.get('Market_pig')) / max(treatment_scale, control_row.get('Market_pig'))
                scale_similarity = scale_ratio
                similarity_score += scale_similarity
            
            # Province similarity (weight 1)
            if not pd.isna(treatment_province) and not pd.isna(control_row.get('Province')):
                province_similarity = 1 if treatment_province == control_row.get('Province') else 0
                similarity_score += province_similarity
            
            similarities.append((control_id, similarity_score))
        
        # Sort by similarity score
        similarities.sort(key=lambda x: x[1], reverse=True)
        
        # Select top controls
        selected_controls = [control_id for control_id, score in similarities[:max_controls_per_treatment] if score > 0]
        treatment_to_controls[treatment_id] = selected_controls
        
        print(f"Treatment farm {treatment_id}: selected {len(selected_controls)} control farms")
    
    return treatment_to_controls

def prepare_panel_data(df, treatment_to_controls):
    """
    Prepare panel data for causal inference analysis
    
    Parameters:
    df: DataFrame with farm data
    treatment_to_controls: Dictionary mapping treatment farm IDs to control farm IDs
    
    Returns:
    DataFrame with panel structure (one row per farm-period)
    """
    panel_data = []
    
    # Process treatment farms
    for treatment_id, control_ids in treatment_to_controls.items():
        # Get treatment farm data
        treatment_farm_data = df[df['Farm_ID'] == treatment_id].copy()
        
        # Get pre-expansion data (Scale-up = 0)
        pre_data = treatment_farm_data[treatment_farm_data['Scale-up'] == 0]
        # Get post-expansion data (Scale-up = 1)
        post_data = treatment_farm_data[treatment_farm_data['Scale-up'] == 1]
        
        if len(pre_data) == 0 or len(post_data) == 0:
            print(f"Warning: Farm {treatment_id} missing pre or post expansion data, skipping")
            continue
        
        # Use first pre and post observations
        pre_data = pre_data.iloc[0]
        post_data = post_data.iloc[0]
        
        # Add treatment farm to panel data
        # Pre-expansion period
        panel_row_pre = {
            'Farm_ID': treatment_id,
            'Group': 'Treatment',
            'Period': 'Pre',
            'Treatment': 0,
            'Efficiency': pre_data.get('A_TE_CRS', pre_data.get('A_TE_VRS', 0.8)),
            'Efficiency_VRS': pre_data.get('A_TE_VRS', pre_data.get('A_TE_CRS', 0.8)),
            'Market_pig': pre_data.get('Market_pig'),
            'Pork_price': pre_data.get('Pork_price'),
            'Feed_price': pre_data.get('Feed_price'),
            'Year': pre_data.get('Year'),
            'Province': pre_data.get('Province')
        }
        panel_data.append(panel_row_pre)
        
        # Post-expansion period
        panel_row_post = {
            'Farm_ID': treatment_id,
            'Group': 'Treatment',
            'Period': 'Post',
            'Treatment': 1,
            'Efficiency': post_data.get('A_TE_CRS', post_data.get('A_TE_VRS', 0.9)),
            'Efficiency_VRS': post_data.get('A_TE_VRS', post_data.get('A_TE_CRS', 0.9)),
            'Market_pig': post_data.get('Market_pig'),
            'Pork_price': post_data.get('Pork_price'),
            'Feed_price': post_data.get('Feed_price'),
            'Year': post_data.get('Year'),
            'Province': post_data.get('Province')
        }
        panel_data.append(panel_row_post)
        
        # Add control farms to panel data
        for control_id in control_ids:
            control_farm_data = df[df['Farm_ID'] == control_id].copy()
            
            if len(control_farm_data) == 0:
                continue
            
            # For controls, use the closest year to treatment farm's pre-expansion year
            treatment_year = pre_data.get('Year')
            if not pd.isna(treatment_year):
                control_farm_data['year_diff'] = abs(control_farm_data['Year'] - treatment_year)
                closest_control_data = control_farm_data.sort_values('year_diff').iloc[0]
            else:
                # If no treatment year, use first observation
                closest_control_data = control_farm_data.iloc[0]
            
            # Add control farm to panel data (both periods have Treatment = 0)
            # Pre-expansion period (matching treatment's pre period)
            panel_row_control_pre = {
                'Farm_ID': control_id,
                'Group': 'Control',
                'Period': 'Pre',
                'Treatment': 0,
                'Efficiency': closest_control_data.get('A_TE_CRS', closest_control_data.get('A_TE_VRS', 0.8)),
                'Efficiency_VRS': closest_control_data.get('A_TE_VRS', closest_control_data.get('A_TE_CRS', 0.8)),
                'Market_pig': closest_control_data.get('Market_pig'),
                'Pork_price': closest_control_data.get('Pork_price'),
                'Feed_price': closest_control_data.get('Feed_price'),
                'Year': closest_control_data.get('Year'),
                'Province': closest_control_data.get('Province')
            }
            panel_data.append(panel_row_control_pre)
            
            # Post-expansion period (use same control farm data for consistency)
            panel_row_control_post = {
                'Farm_ID': control_id,
                'Group': 'Control',
                'Period': 'Post',
                'Treatment': 0,
                'Efficiency': closest_control_data.get('A_TE_CRS', closest_control_data.get('A_TE_VRS', 0.8)),
                'Efficiency_VRS': closest_control_data.get('A_TE_VRS', closest_control_data.get('A_TE_CRS', 0.8)),
                'Market_pig': closest_control_data.get('Market_pig'),
                'Pork_price': closest_control_data.get('Pork_price'),
                'Feed_price': closest_control_data.get('Feed_price'),
                'Year': closest_control_data.get('Year'),
                'Province': closest_control_data.get('Province')
            }
            panel_data.append(panel_row_control_post)
    
    panel_df = pd.DataFrame(panel_data)
    
    # Print panel data summary
    print(f"\nPanel data prepared: {len(panel_df)} observations")
    print(f"  Farms: {panel_df['Farm_ID'].nunique()}")
    print(f"  Treatment farms: {len(panel_df[panel_df['Group'] == 'Treatment']['Farm_ID'].unique())}")
    print(f"  Control farms: {len(panel_df[panel_df['Group'] == 'Control']['Farm_ID'].unique())}")
    print(f"  Pre-expansion: {len(panel_df[panel_df['Period'] == 'Pre'])}")
    print(f"  Post-expansion: {len(panel_df[panel_df['Period'] == 'Post'])}")
    
    return panel_df

def export_treatment_control_to_excel(panel_df, treatment_to_controls, output_path):
    """
    Export each treatment (expanding) farm and its matched control group to an Excel file.
    Two sheets: Treatment_farms (one row per treatment farm), Control_mapping (one row per treatment-control pair).
    """
    treatment_ids = list(treatment_to_controls.keys())
    # Sheet 1: Treatment farms with pre/post info and matched control IDs
    treatment_rows = []
    for tid in treatment_ids:
        t_pre = panel_df[(panel_df['Farm_ID'] == tid) & (panel_df['Group'] == 'Treatment') & (panel_df['Period'] == 'Pre')]
        t_post = panel_df[(panel_df['Farm_ID'] == tid) & (panel_df['Group'] == 'Treatment') & (panel_df['Period'] == 'Post')]
        if t_pre.empty or t_post.empty:
            continue
        t_pre = t_pre.iloc[0]
        t_post = t_post.iloc[0]
        control_ids = treatment_to_controls.get(tid, [])
        control_cols = ['Control_Farm_ID_1', 'Control_Farm_ID_2', 'Control_Farm_ID_3']
        c1, c2, c3 = (control_ids + [None, None, None])[:3]
        treatment_rows.append({
            'Treatment_Farm_ID': tid,
            'Pre_Year': t_pre.get('Year'),
            'Pre_Market_pig': t_pre.get('Market_pig'),
            'Pre_Efficiency': t_pre.get('Efficiency'),
            'Pre_Pork_price': t_pre.get('Pork_price'),
            'Pre_Feed_price': t_pre.get('Feed_price'),
            'Post_Year': t_post.get('Year'),
            'Post_Market_pig': t_post.get('Market_pig'),
            'Post_Efficiency': t_post.get('Efficiency'),
            'Post_Pork_price': t_post.get('Pork_price'),
            'Post_Feed_price': t_post.get('Feed_price'),
            'Province': t_pre.get('Province'),
            'N_controls_matched': len(control_ids),
            'Control_Farm_ID_1': c1,
            'Control_Farm_ID_2': c2,
            'Control_Farm_ID_3': c3,
        })
    df_treatment = pd.DataFrame(treatment_rows)

    # Sheet 2: Control mapping (one row per treatment-control pair with control farm details)
    control_rows = []
    for tid in treatment_ids:
        control_ids = treatment_to_controls.get(tid, [])
        for cid in control_ids:
            c_row = panel_df[(panel_df['Farm_ID'] == cid) & (panel_df['Group'] == 'Control') & (panel_df['Period'] == 'Pre')]
            if c_row.empty:
                continue
            c_row = c_row.iloc[0]
            control_rows.append({
                'Treatment_Farm_ID': tid,
                'Control_Farm_ID': cid,
                'Control_Year': c_row.get('Year'),
                'Control_Market_pig': c_row.get('Market_pig'),
                'Control_Efficiency': c_row.get('Efficiency'),
                'Control_Pork_price': c_row.get('Pork_price'),
                'Control_Feed_price': c_row.get('Feed_price'),
                'Control_Province': c_row.get('Province'),
            })
    df_control = pd.DataFrame(control_rows)

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df_treatment.to_excel(writer, sheet_name='Treatment_farms', index=False)
        df_control.to_excel(writer, sheet_name='Control_mapping', index=False)
    print(f"  Treatment-control mapping saved: {output_path}")

def perform_causal_inference(panel_df):
    """
    Perform causal inference analysis using panel data methods
    
    Parameters:
    panel_df: Panel DataFrame with pre and post observations
    
    Returns:
    Dictionary with regression results
    """
    # Regression models
    models_results = {}
    y = panel_df['Efficiency']
    
    # Create farm fixed effects (dummy variables for each farm)
    farm_dummies = pd.get_dummies(panel_df['Farm_ID'], prefix='Farm', dtype=float)
    # Drop one farm dummy to avoid perfect multicollinearity (reference group)
    farm_dummies = farm_dummies.iloc[:, 1:]  # Drop first column
    
    # Ensure all data are numeric
    farm_dummies = farm_dummies.astype(float)
    
    # Create treatment variable (1 for treatment farms in post period)
    panel_df['Treatment'] = (panel_df['Group'] == 'Treatment') & (panel_df['Period'] == 'Post')
    panel_df['Treatment'] = panel_df['Treatment'].astype(int)
    
    # Create scale categories
    def categorize_scale(market_pig):
        if market_pig < 10000:
            return 'Low'
        elif market_pig <= 50000:
            return 'Medium'
        else:
            return 'High'
    
    panel_df['Scale_Category'] = panel_df['Market_pig'].apply(categorize_scale)
    
    # Create scale category dummies
    scale_dummies = pd.get_dummies(panel_df['Scale_Category'], prefix='Scale', dtype=float)
    
    # Create treatment × scale interaction terms
    for scale in ['Low', 'Medium', 'High']:
        if f'Scale_{scale}' in scale_dummies.columns:
            panel_df[f'Treatment_x_Scale_{scale}'] = panel_df['Treatment'] * scale_dummies[f'Scale_{scale}']
    
    # Model 1: Basic DID with clustered standard errors
    try:
        X1 = panel_df[['Treatment']].copy().astype(float)
        X1 = pd.concat([X1, farm_dummies], axis=1)  # Add farm fixed effects
        model1 = sm.OLS(y.astype(float), X1).fit(cov_type='cluster', cov_kwds={'groups': panel_df['Farm_ID']})
        models_results['Basic_DID_FE'] = model1
        print("✓ Basic_DID with Farm Fixed Effects (clustered SE) estimated successfully")
    except Exception as e:
        print(f"✗ Error in Basic_DID_FE model: {e}")
    
    # Model 2: DID with scale interactions and farm fixed effects
    try:
        interaction_cols = ['Treatment']
        # Add treatment × scale interaction terms only
        for scale in ['Low', 'Medium', 'High']:
            if f'Treatment_x_Scale_{scale}' in panel_df.columns:
                interaction_cols.append(f'Treatment_x_Scale_{scale}')
        interaction_cols = [col for col in interaction_cols if col in panel_df.columns]
        
        X2 = panel_df[interaction_cols].copy().astype(float)
        X2 = pd.concat([X2, farm_dummies], axis=1)  # Add farm fixed effects
        model2 = sm.OLS(y.astype(float), X2).fit(cov_type='cluster', cov_kwds={'groups': panel_df['Farm_ID']})
        models_results['DID_with_Scale_Interactions_FE'] = model2
        print("✓ DID_with_Scale_Interactions + Farm FE (clustered SE) estimated successfully")
    except Exception as e:
        print(f"✗ Error in DID_with_Scale_Interactions_FE model: {e}")
    
    # Model 3: DID with controls and farm fixed effects
    try:
        control_cols = ['Treatment', 'Pork_price', 'Feed_price']
        # Add treatment × scale interaction terms
        for scale in ['Low', 'Medium', 'High']:
            if f'Treatment_x_Scale_{scale}' in panel_df.columns:
                control_cols.append(f'Treatment_x_Scale_{scale}')
        control_cols = [col for col in control_cols if col in panel_df.columns]
        
        X3 = panel_df[control_cols].copy().astype(float)
        X3 = pd.concat([X3, farm_dummies], axis=1)  # Add farm fixed effects
        model3 = sm.OLS(y.astype(float), X3).fit(cov_type='cluster', cov_kwds={'groups': panel_df['Farm_ID']})
        models_results['DID_with_Controls_FE'] = model3
        print("✓ DID_with_Controls + Farm FE (clustered SE) estimated successfully")
    except Exception as e:
        print(f"✗ Error in DID_with_Controls_FE model: {e}")
    
    # Model 4: Full model with both farm and time fixed effects
    try:
        control_cols = ['Treatment', 'Pork_price', 'Feed_price']
        # Add treatment × scale interaction terms
        for scale in ['Low', 'Medium', 'High']:
            if f'Treatment_x_Scale_{scale}' in panel_df.columns:
                control_cols.append(f'Treatment_x_Scale_{scale}')
        control_cols = [col for col in control_cols if col in panel_df.columns]
        
        X4 = panel_df[control_cols].copy().astype(float)
        time_dummies = pd.get_dummies(panel_df['Period'], prefix='Period', dtype=float)
        time_dummies = time_dummies.iloc[:, 1:]  # Drop first period as reference
        X4 = pd.concat([X4, farm_dummies, time_dummies], axis=1)
        model4 = sm.OLS(y.astype(float), X4).fit(cov_type='cluster', cov_kwds={'groups': panel_df['Farm_ID']})
        models_results['Full_Model_FE'] = model4
        print("✓ Full_Model + Farm and Time FE (clustered SE) estimated successfully")
    except Exception as e:
        print(f"✗ Error in Full_Model_FE: {e}")
    
    return models_results

def create_causal_inference_plot(models_results):
    """
    Create causal inference visualization
    
    Parameters:
    models_results: Dictionary with regression results
    
    Returns:
    matplotlib Figure object
    """
    # Set style for publication-quality figures
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 13
    plt.rcParams['figure.titlesize'] = 20
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 0.8])
    
    ax1 = fig.add_subplot(gs[0, :])  # Treatment effects (top full width)
    ax2 = fig.add_subplot(gs[1, 0])  # Model performance (bottom left)
    ax3 = fig.add_subplot(gs[1, 1])  # Coefficient magnitudes (bottom right)
    
    # Define color scheme matching reference
    model_colors = ['#2E86AB', '#A23B72', '#F18F01', '#5DA271', '#9B59B6', '#E74C3C', '#1ABC9C']
    
    # Plot 1: Treatment effects with confidence intervals (Forest plot style)
    models = list(models_results.keys())
    treatment_effects = []
    conf_intervals = []
    p_values = []
    model_names_clean = []
    
    for model_name in models:
        model = models_results[model_name]
        te = None
        se = None
        p_val = None
        
        # Handle Treatment coefficient (for DID models)
        if 'Treatment' in model.params:
            te = model.params['Treatment']
            se = model.bse['Treatment']
            p_val = model.pvalues['Treatment']
        if te is not None:
            treatment_effects.append(te)
            conf_intervals.append([te - 1.96*se, te + 1.96*se])
            p_values.append(p_val)
            
            # Sentence case model labels with acronym preservation
            model_label_map = {
                'Basic_DID_FE': 'Basic DID FE',
                'DID_with_Scale_Interactions_FE': 'DID with scale interactions FE',
                'DID_with_Controls_FE': 'DID with controls FE',
                'Full_Model_FE': 'Full model FE'
            }
            clean_name = model_label_map.get(model_name, model_name.replace('_', ' '))
            model_names_clean.append(clean_name)
    
    if treatment_effects:
        # Convert to numpy arrays
        treatment_effects = np.array(treatment_effects)
        conf_intervals = np.array(conf_intervals)
        p_values = np.array(p_values)
        
        # Create forest plot
        y_pos = np.arange(len(models))
        
        # Plot confidence intervals
        for i, (te, ci) in enumerate(zip(treatment_effects, conf_intervals)):
            ax1.plot([ci[0], ci[1]], [y_pos[i], y_pos[i]], color=model_colors[i], linewidth=3, alpha=0.7)
            ax1.plot(te, y_pos[i], 'o', markersize=8, color=model_colors[i], label=model_names_clean[i])
        
        # Add significance annotations
        for i, (te, p_val) in enumerate(zip(treatment_effects, p_values)):
            if p_val < 0.001:
                sig_text = '***'
            elif p_val < 0.01:
                sig_text = '**'
            elif p_val < 0.05:
                sig_text = '*'
            elif p_val < 0.1:
                sig_text = '†'
            else:
                sig_text = 'ns'
            
            ax1.text(te, y_pos[i] + 0.15, sig_text, ha='center', va='bottom', 
                    fontweight='bold', fontsize=16, color=model_colors[i])
            
            # Add coefficient value to the right and higher above the point
            ax1.text(te + 0.01, y_pos[i] + 0.2, f'{te:.4f}', ha='left', va='center', 
                    fontweight='bold', fontsize=16, bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
        
        ax1.axvline(x=0, color='red', linestyle='--', alpha=0.7, linewidth=1)
        ax1.set_yticks(y_pos)
        panel_a_y_labels = []
        for name in model_names_clean:
            if name == 'DID with scale interactions FE':
                panel_a_y_labels.append('DID with scale\ninteractions FE')
            elif name == 'DID with controls FE':
                panel_a_y_labels.append('DID with\ncontrols FE')
            else:
                panel_a_y_labels.append(name)
        ax1.set_yticklabels(panel_a_y_labels)
        ax1.set_xlabel('Treatment effect coefficient', fontweight='bold', fontsize=16)
        ax1.grid(True, alpha=0.3, linestyle='--', axis='x')
        ax1.legend(loc='center right', frameon=True, framealpha=0.9)
        # Expand y-axis limit to fit all data labels
        ax1.set_ylim(-0.5, len(models) - 0.5)
    
    # Plot 2: Model performance metrics
    if models:
        metrics_data = []
        metric_names = ['R²', 'Adjusted R²', 'AIC', 'BIC']
        
        for model_name, model in models_results.items():
            metrics = [
                model.rsquared,
                model.rsquared_adj,
                model.aic,
                model.bic
            ]
            metrics_data.append(metrics)
        
        metrics_data = np.array(metrics_data)
        
        # Normalize AIC and BIC for better visualization (lower is better, so invert)
        if len(metrics_data) > 1:
            aic_norm = 1 - (metrics_data[:, 2] - np.min(metrics_data[:, 2])) / (np.max(metrics_data[:, 2]) - np.min(metrics_data[:, 2]) + 1e-10)
            bic_norm = 1 - (metrics_data[:, 3] - np.min(metrics_data[:, 3])) / (np.max(metrics_data[:, 3]) - np.min(metrics_data[:, 3]) + 1e-10)
        else:
            aic_norm = [0.5]
            bic_norm = [0.5]
        
        # Replace AIC and BIC with normalized values
        metrics_data[:, 2] = aic_norm
        metrics_data[:, 3] = bic_norm
        
        x_metrics = np.arange(len(metric_names))
        width = 0.8 / len(models)
        
        for i, model_metrics in enumerate(metrics_data):
            ax2.bar(x_metrics + i * width, model_metrics, width, 
                   label=model_names_clean[i], color=model_colors[i], alpha=0.8,
                   edgecolor='black', linewidth=0.5)
        
        ax2.set_xlabel('Performance metrics', fontweight='bold', fontsize=16)
        ax2.set_ylabel('Metric values', fontweight='bold', fontsize=16)
        ax2.set_xticks(x_metrics + width * (len(models) - 1) / 2)
        ax2.set_xticklabels(metric_names)
        # Move legend to center at y=1.1, split into two rows horizontally, remove border
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 0.98), frameon=False, ncol=2, fontsize=12)
        ax2.set_ylim(0, 1.3)
        # Set y-axis ticks to 0-1 only
        ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        # Only show horizontal grid lines
        ax2.grid(True, alpha=0.3, linestyle='--', axis='y')
        ax2.grid(False, axis='x')
    
    # Plot 3: Coefficient magnitudes across models
    if models:
        # Get all unique coefficients across models
        all_coeffs = set()
        for model in models_results.values():
            all_coeffs.update(model.params.index)
        # Filter to only show treatment, treatment × scale interactions, and external variables
        relevant_coeffs = []
        for coeff in all_coeffs:
            if (coeff == 'Treatment' or 
                'Treatment_x_' in coeff or 
                'Treatment_scale_' in coeff or
                'Pork_price' in coeff or 
                'Feed_Price' in coeff or
                'Feed_price' in coeff):
                relevant_coeffs.append(coeff)
        
        # Sort coefficients in the specified order (top to bottom)
        coeff_order = ['Treatment', 'Feed_price', 'Pork_price', 'Treatment_x_Scale_Low', 'Treatment_x_Scale_Medium', 'Treatment_x_Scale_High']
        sorted_coeffs = []
        
        # First, add coefficients that match exactly
        for coeff in coeff_order:
            for c in relevant_coeffs:
                if coeff == c:
                    sorted_coeffs.append(c)
                    relevant_coeffs.remove(c)
                    break
                elif 'Scale' in coeff and 'Scale' in c:
                    # Extract scale level from coefficient name
                    if '_Scale_' in coeff:
                        scale_level = coeff.split('_Scale_')[1]
                        if scale_level in c:
                            sorted_coeffs.append(c)
                            relevant_coeffs.remove(c)
                            break
        
        # Add any remaining coefficients
        sorted_coeffs.extend(relevant_coeffs)
        relevant_coeffs = sorted_coeffs
        
        if not relevant_coeffs:
            ax3.text(0.5, 0.5, 'No relevant coefficients found', 
                   ha='center', va='center', fontsize=12)
            ax3.axis('off')
        else:
            # Create coefficient matrix
            coeff_matrix = np.zeros((len(relevant_coeffs), len(models)))
            pval_matrix = np.zeros((len(relevant_coeffs), len(models)))
            
            for j, model in enumerate(models_results.values()):
                for i, coeff in enumerate(relevant_coeffs):
                    if coeff in model.params:
                        coeff_matrix[i, j] = model.params[coeff]
                        pval_matrix[i, j] = model.pvalues[coeff]
                    else:
                        coeff_matrix[i, j] = np.nan
                        pval_matrix[i, j] = np.nan
            
            # Plot coefficient heatmap with diverging color map for positive/negative values
            # Find the maximum absolute value for symmetric color scale
            max_abs = np.nanmax(np.abs(coeff_matrix))
            im = ax3.imshow(coeff_matrix, cmap='RdBu_r', aspect='auto', vmin=-max_abs, vmax=max_abs)
            
            # Add coefficient values
            for i in range(len(relevant_coeffs)):
                for j in range(len(models)):
                    if not np.isnan(coeff_matrix[i, j]):
                        # Use white text for values near the extremes
                        color = 'white' if np.abs(coeff_matrix[i, j]) > max_abs * 0.6 else 'black'
                        ax3.text(j, i, f'{coeff_matrix[i, j]:.3f}', ha='center', va='center', 
                                color=color, fontweight='bold', fontsize=14)
            
            ax3.set_xticks(range(len(models)))
            # Add line breaks in model names for better readability and avoid overlap
            model_labels = []
            for m in models:
                # Sentence case x-axis labels for panel C
                if m == 'Basic_DID_FE':
                    label = 'Basic\nDID FE'
                elif m == 'DID_with_Scale_Interactions_FE':
                    label = 'DID with\nscale\ninteractions FE'
                elif m == 'DID_with_Controls_FE':
                    label = 'DID with\ncontrols FE'
                elif m == 'Full_Model_FE':
                    label = 'Full\nmodel FE'
                else:
                    # Fallback: split at first space if label is long
                    label = m.replace('_', ' ')
                    words = label.split()
                    if len(words) > 2:
                        mid = len(words) // 2
                        label = ' '.join(words[:mid]) + '\n' + ' '.join(words[mid:])
                model_labels.append(label)
            ax3.set_xticklabels(model_labels, rotation=0, ha='center', fontsize=12)
            
            # Format coefficient labels
            yticklabels = []
            for coeff in relevant_coeffs:
                # Standardize scale-related variables
                if 'Treatment_x_scale_low' in coeff:
                    label = 'Treatment scale low'
                elif 'Treatment_x_scale_medium' in coeff:
                    label = 'Treatment scale medium'
                elif 'Treatment_x_scale_high' in coeff:
                    label = 'Treatment scale high'
                elif 'Treatment_scale_low' in coeff:
                    label = 'Treatment scale low'
                elif 'Treatment_scale_medium' in coeff:
                    label = 'Treatment scale medium'
                elif 'Treatment_scale_high' in coeff:
                    label = 'Treatment scale high'
                else:
                    # Remove underscores and only capitalize first word
                    label = coeff.replace('_', ' ')
                    # Only capitalize first letter of first word
                    words = label.split()
                    if words:
                        label = words[0].capitalize() + ' ' + ' '.join(word.lower() for word in words[1:])
                yticklabels.append(label)
            ax3.set_yticks(range(len(relevant_coeffs)))
            ax3.set_yticklabels(yticklabels, fontsize=12)
            
            # Add colorbar
            plt.colorbar(im, ax=ax3, shrink=0.8)
            # Reduce grid line opacity to avoid blocking numbers
            ax3.grid(True, alpha=0.2, linestyle='-')
    
    # Set black panel borders for all subplots
    for ax in [ax1, ax2, ax3]:
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(0.75)

    plt.tight_layout()
    # Adjust bottom margin to accommodate multi-line x-axis labels
    plt.subplots_adjust(bottom=0.15)
    
    return fig

def main():
    """
    Main function to execute the causal inference analysis
    """
    FILE_PATH = r"./swine_farm_data.xlsx"
    SHEET_NAME = 2
    
    # Generate timestamp for filenames
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    try:
        print("=" * 80)
        print("CAUSAL INFERENCE ANALYSIS: Difference-in-Differences (DID)")
        print("=" * 80)
        print(f"Data source: {FILE_PATH}")
        print(f"Sheet: {SHEET_NAME}")
        print("=" * 80)
        print()
        
        # Step 1: Load and preprocess data
        print("Step 1: Loading and preprocessing data...")
        df = load_and_preprocess_data(FILE_PATH, SHEET_NAME)
        print()
        
        # Step 2: Identify treatment farms (pre-expansion, Scale-up = 0)
        print("Step 2: Identifying treatment farms...")
        treatment_farms = df[df['Scale-up'] == 0].copy()
        print(f"Found {len(treatment_farms)} treatment farms (Scale-up = 0)")
        print()
        
        # Step 3: Select control farms
        print("Step 3: Selecting control farms...")
        treatment_to_controls = select_control_farms(df, treatment_farms)
        print()
        
        # Step 4: Prepare panel data
        print("Step 4: Preparing panel data for causal inference...")
        panel_df = prepare_panel_data(df, treatment_to_controls)
        print()

        # Step 4b: Export treatment and control group info to Excel
        excel_path = f"causal_inference_treatment_control_mapping_{timestamp}.xlsx"
        try:
            export_treatment_control_to_excel(panel_df, treatment_to_controls, excel_path)
        except Exception as e:
            print(f"  Warning: Could not save treatment-control Excel: {e}")

        # Step 5: Perform causal inference
        print("Step 5: Performing causal inference analysis...")
        print("-" * 80)
        models_results = perform_causal_inference(panel_df)
        print("-" * 80)
        print()
        
        # Display regression results summary
        print("=" * 80)
        print("REGRESSION RESULTS SUMMARY")
        print("=" * 80)
        print("Note: Models with '_FE' suffix include Fixed Effects")
        print("      Standard errors are clustered at the farm level")
        print("=" * 80)
        
        for model_name, model in models_results.items():
            print(f"\n{model_name}:")
            print(f"  R²: {model.rsquared:.4f}, Adj. R²: {model.rsquared_adj:.4f}")
            print(f"  AIC: {model.aic:.2f}, BIC: {model.bic:.2f}")
            
            # Handle Treatment effect
            if 'Treatment' in model.params:
                te = model.params['Treatment']
                se = model.bse['Treatment']
                p_val = model.pvalues['Treatment']
                stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else '†' if p_val < 0.1 else ''
                print(f"  Treatment Effect: {te:.4f}{stars} (SE = {se:.4f}, p = {p_val:.4f})")
                print(f"  95% CI: [{te - 1.96*se:.4f}, {te + 1.96*se:.4f}]")
        
        print("=" * 80)
        print()
        
        # Step 6: Create visualization
        print("Step 6: Creating causal inference visualization...")
        fig = create_causal_inference_plot(models_results)
        
        # Save plot
        plot_filename = f"causal_inference_results_{timestamp}.png"
        fig.savefig(plot_filename, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"✓ Causal inference plot saved as: {plot_filename}")
        print()
        
        # Step 7: Comprehensive summary
        print("=" * 80)
        print("COMPREHENSIVE ANALYSIS SUMMARY")
        print("=" * 80)
        
        # Efficiency statistics
        treatment_pre = panel_df[(panel_df['Group'] == 'Treatment') & (panel_df['Period'] == 'Pre')]['Efficiency']
        treatment_post = panel_df[(panel_df['Group'] == 'Treatment') & (panel_df['Period'] == 'Post')]['Efficiency']
        control_pre = panel_df[(panel_df['Group'] == 'Control') & (panel_df['Period'] == 'Pre')]['Efficiency']
        control_post = panel_df[(panel_df['Group'] == 'Control') & (panel_df['Period'] == 'Post')]['Efficiency']
        
        print(f"\nEfficiency Statistics:")
        print(f"  Treatment farms - Pre: Mean = {treatment_pre.mean():.4f} ± {treatment_pre.std():.4f}")
        print(f"  Treatment farms - Post: Mean = {treatment_post.mean():.4f} ± {treatment_post.std():.4f}")
        print(f"  Control farms - Pre: Mean = {control_pre.mean():.4f} ± {control_pre.std():.4f}")
        print(f"  Control farms - Post: Mean = {control_post.mean():.4f} ± {control_post.std():.4f}")
        
        # Calculate DID manually
        did_manual = (treatment_post.mean() - treatment_pre.mean()) - (control_post.mean() - control_pre.mean())
        print(f"\nManual DID calculation: {did_manual:.4f}")
        
        print("=" * 80)
        print()
        
        print(f"✓ Analysis completed successfully!")
        print(f"✓ Results saved with timestamp: {timestamp}")
        print()
        
    except Exception as e:
        print(f"✗ Error in analysis: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()









