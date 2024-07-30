import pandas as pd
from pandas.api.types import CategoricalDtype
from plotnine import ggplot, aes, geom_bar, theme, element_text, labs, scale_fill_manual, geom_vline, geom_segment
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, f_oneway



def get_raw_data(sys_argv):
    files_path = sys_argv[1]
    samples_loaded = pd.read_csv(f"{files_path}/samples_loaded.txt", delimiter='\t')    
    superkingdom = pd.read_csv(f"{files_path}/superkingdom2descendents.txt", delimiter='\t')
    species_abundance = pd.read_csv(f"{files_path}/species_abundance.txt", delimiter='\t')
    project = pd.read_csv(f"{files_path}/sample_to_run_info.txt", delimiter='\t', low_memory=False)
    families_table = pd.read_csv(f"{files_path}/families.csv")
    return(samples_loaded, superkingdom, species_abundance, project,families_table)



def filter_project(project):
    selected_project_data = project[project['project_id'] == "PRJEB8094"]
    return selected_project_data



def add_accession_id(samples_loaded, selected_project_data):
    uid_table = pd.DataFrame({
        'loaded_uid': samples_loaded['uid'],
        'run_id': samples_loaded['accession_id']
    })
    
    selected_project_w_accession = pd.merge(selected_project_data, uid_table, how='left')
    selected_project_w_accession = selected_project_w_accession[~selected_project_w_accession['loaded_uid'].isin([25687, 2615, 4124])]
    return selected_project_w_accession



def add_kingdom(species_abundance, superkingdom):
    species_abundance_kingdom = pd.merge(species_abundance, superkingdom, how='left')
    species_abundance_kingdom = species_abundance_kingdom[species_abundance_kingdom['taxon_rank_level'] == 'species']
    return species_abundance_kingdom



def filter_species_abundance(selected_project_w_accession, species_abundance_kingdom):
    temp_uid = pd.DataFrame({'loaded_uid': selected_project_w_accession['loaded_uid']})
    filtered_species_abundance = pd.merge(species_abundance_kingdom, temp_uid, how='right', on='loaded_uid')
    #correct NA values in 'scientific_name' to "other"
    filtered_species_abundance['scientific_name'] = filtered_species_abundance['scientific_name'].fillna('Other')
    desired_columns = ['loaded_uid', 'scientific_name', 'relative_abundance']
    filtered_species_abundance = filtered_species_abundance[desired_columns]
    return filtered_species_abundance



def match_sample_name(selected_project_w_accession,filtered_species_abundance):
    temp_samplename_uid = pd.DataFrame(selected_project_w_accession['sample_name'])
    temp_samplename_uid.columns = ['sample_name']
    temp_samplename_uid['loaded_uid'] = selected_project_w_accession['loaded_uid']
    species_abundance_w_sample_name = pd.merge(filtered_species_abundance, temp_samplename_uid, how='left', on='loaded_uid')
    return species_abundance_w_sample_name



def add_family(species_abundance_w_sample_name,families_table):
    families_table.columns = ["x", "scientific_name", "family"]
    families_table = families_table[["scientific_name", "family"]]
    species_abundance_family = pd.merge(species_abundance_w_sample_name, families_table, how='left', on='scientific_name')
    species_abundance_family['family'] = species_abundance_family['family'].str.replace('.*: ', '', regex=True)
    species_abundance_family['family'] = species_abundance_family['family'] + '_' + species_abundance_family['sample_name']
    return species_abundance_family



def cal_relative_abundance(species_abundance_family,selected_project_w_accession):
    species_abundance_sum = species_abundance_family.groupby('family', as_index=False)['relative_abundance'].sum()
    species_abundance_sum.columns = ['family_sample', 'x']
    species_abundance_sum['samplename'] = species_abundance_sum['family_sample']
    species_abundance_sum['samplename'] = species_abundance_sum['samplename'].str.replace('.*_', '', regex=True)
    repetitions = selected_project_w_accession['sample_name'].value_counts().reset_index()
    repetitions.columns = ['samplename', 'Freq']
    species_abundance_sum = pd.merge(species_abundance_sum, repetitions, on='samplename', how='left')
    # Calculating the relative abundance:
    species_abundance_sum['relative_abundance'] = species_abundance_sum['x'] / species_abundance_sum['Freq']
    species_abundance_sum.columns = ["family", "x", "sample_name", "Freq", "relative_abundance"]
    species_abundance_sum['family'] = species_abundance_sum['family'].str.replace('_.*', '', regex=True)
    return species_abundance_sum



def relative_abundance_plot(species_relative_abundance, sys_argv):
    path_to_save = os.path.join(sys_argv[2], "relative_abundance_plot")
    # Define the order of samples (Experiment (time point 0,7,9) and Control (time point 0,7,9))
    Listofsamples = [
        "P1E0", "P2E0", "P3E0", "P4E0", "P5E0", "P9E0", "P10E0", "P11E0", "P12E0", "P13E0", 
        "P14E0", "P15E0", "P17E0", "P18E0", "P19E0", "P20E0", "P21E0", "P22E0", "P1E7", "P2E7", 
        "P3E7", "P4E7", "P5E7", "P9E7", "P10E7", "P11E7", "P12E7", "P13E7", "P14E7", "P15E7", 
        "P17E7", "P18E7", "P19E7", "P20E7", "P21E7", "P22E7", "P1E90", "P2E90", "P3E90", "P4E90", 
        "P5E90", "P9E90", "P10E90", "P11E90", "P12E90", "P13E90", "P14E90", "P15E90", "P17E90", 
        "P18E90", "P19E90", "P20E90", "P21E90", "P22E90", "P6C0", "P7C0", "P8C0", "P23C0", "P25C0", 
        "P38C0", "P6C7", "P7C7", "P8C7", "P23C7", "P25C7", "P38C7", "P6C90", "P7C90", "P8C90", 
        "P23C90", "P25C90", "P38C90"
    ]
    cat_type = CategoricalDtype(categories=Listofsamples, ordered=True)
    species_relative_abundance['sample_name'] = species_relative_abundance['sample_name'].astype(cat_type)
    # Define the plot colors
    hex_colors = [
    "#011F4B", "#022655", "#03315F", "#033B69", "#044573", "#05507D", "#065A87", "#076491",
    "#087E9B", "#0988A5", "#0A92AF", "#0BAEBF", "#1EC3CB", "#3FD0D7", "#5FE1E3", "#7FF0EF",
    "#8FF3F3", "#9FF6F7", "#AFFAFB", "#C4FBF9", "#D9FCF7", "#EEFDF5", "#F3F3E3", "#F8E9D1",
    "#FDE1BF", "#FED9AD", "#FFD198", "#FFC986", "#FFC174", "#FFB962", "#FFB050", "#FFA83D",
    "#FFA02B", "#FF9819"]
    
    vertical_lines = {
        "P22E0": "P1E7",
        "P22E7": "P1E90",
        "P22E90": "P6C0"
    }
    vertical_lines_pos = {k: Listofsamples.index(k) + 1.45 for k in vertical_lines.keys()}
    max_y = species_relative_abundance['relative_abundance'].max() * 1.165
    
    # Create the plot
    relative_abundance_plot = (
        ggplot(species_relative_abundance, aes(x='sample_name', y='relative_abundance', fill='family')) +
        geom_bar(stat='identity', position='stack') +
        theme(axis_text_x=element_text(rotation=90, hjust=1)) +
        labs(x='Sample Name', y='Relative Abundance (%)', title='Relative Abundance of Bacteria Types by Sample') +
        scale_fill_manual(values=hex_colors) +
        geom_segment(aes(x=vertical_lines_pos["P22E0"], y=0, xend=vertical_lines_pos["P22E0"], yend=max_y), color='black', size=1) +
        geom_segment(aes(x=vertical_lines_pos["P22E7"], y=0, xend=vertical_lines_pos["P22E7"], yend=max_y), color='black', size=1) +
        geom_segment(aes(x=vertical_lines_pos["P22E90"], y=0, xend=vertical_lines_pos["P22E90"], yend=max_y), color='black', size=1)
    )
    
    relative_abundance_plot.save(path_to_save, width=15, height=10, dpi=300)
    # relative_abundance_plot.show()


def add_patient_info(species_relative_abundance):
    species_abundance_avg_tagged = species_relative_abundance.copy()
    # Extract patient information
    species_abundance_avg_tagged['patient'] = species_abundance_avg_tagged['sample_name'].str.replace(r'E.*', '', regex=True)
    species_abundance_avg_tagged['patient'] = species_abundance_avg_tagged['patient'].str.replace(r'C.*', '', regex=True)
    # Extract timepoint information
    species_abundance_avg_tagged['timepoint'] = species_abundance_avg_tagged['sample_name'].str.replace(r'.*E', '', regex=True)
    species_abundance_avg_tagged['timepoint'] = species_abundance_avg_tagged['timepoint'].str.replace(r'.*C', '', regex=True)
    # Read patient info CSV
    patient_info = pd.read_csv("patient_info.csv")
    # Merge the DataFrames
    species_abundance_avg_tagged = pd.merge(species_abundance_avg_tagged, patient_info, how='left')
    return species_abundance_avg_tagged



def get_control_data(species_abundance_avg_tagged):
    control_data = species_abundance_avg_tagged[species_abundance_avg_tagged['Category'] == "Control"]
    return control_data



def get_exposed_data(species_abundance_avg_tagged):
    exposed_data = species_abundance_avg_tagged[species_abundance_avg_tagged['Category'] == "Exposed"]
    return exposed_data



def plot_families_aboundance(data_x, group, sys_argv):
    plt.clf()
    plt.close('all')
    plotgroup = group
    path_to_save = os.path.join(sys_argv[2], f"families_{plotgroup}_plot.png")
    timepoint_order = ['0', '7', '90']
    plt.rcParams.update({'font.size': 20})
    # Define significance levels
    def significance_label(p_value):
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        else:
            return None  # Not significant
    # Perform statistical tests
    families = data_x['family'].unique()
    significance_data = {}
    significant_families = set()  #keep track of significant families
    for family in families:
        family_data = data_x[data_x['family'] == family]
        timepoints = family_data['timepoint'].unique()
        comparisons = [(a, b) for idx, a in enumerate(timepoints) for b in timepoints[idx + 1:]]
        for tp1, tp2 in comparisons:
            group1 = family_data[family_data['timepoint'] == tp1]['relative_abundance']
            group2 = family_data[family_data['timepoint'] == tp2]['relative_abundance']
            stat, p_value = ttest_ind(group1, group2)
            label = significance_label(p_value)
            if label:
                if family not in significance_data:
                    significance_data[family] = []
                significance_data[family].append((tp1, tp2, label))
                significant_families.add(family)  # Add family to significant families set
    sns.set(style="whitegrid")
    g = sns.catplot(
        data=data_x,
        x='timepoint',
        y='relative_abundance',
        kind='box',
        col='family',
        col_wrap=8,
        sharey=False,
        palette='pastel',
        height=8, 
        aspect=0.75,  
        order=timepoint_order,
        flierprops=dict(marker='o', color='black', markersize=6, markerfacecolor='black', markeredgecolor='black')  
    )
    g.set_axis_labels("T.P", "Relative Abundance")
    g.set_titles(template="{col_name}", row_template="{row_name}", size=26, y=1.1)
    g.fig.subplots_adjust(top=0.85, hspace=0.6)  
    for family, p_values in significance_data.items():
        ax = g.axes_dict[family]
        y = data_x[data_x['family'] == family]['relative_abundance'].max()
        h = 0.02 
        increment = 0.03  
        base_y_offset = y + 0.05 
        for i, (tp1, tp2, label) in enumerate(p_values):
            y_offset = base_y_offset + i * increment
            x1, x2 = timepoint_order.index(tp1), timepoint_order.index(tp2)
            ax.plot([x1, x1, x2, x2], [y_offset, y_offset + h, y_offset + h, y_offset], lw=1.5, c='k')
            ax.text((x1 + x2) * .5, y_offset + h, label, ha='center', va='bottom', color='k', fontsize=26)
    g.savefig(path_to_save, bbox_inches="tight")
    plt.close()
    return significant_families



def plot_combined_families(data, family, sys_argv):
    plt.clf()
    plt.close('all')
    path_to_save = os.path.join(sys_argv[2], f"combined_family_{family}_plot.png")
    plt.rcParams.update({'font.size': 20})
    timepoint_order = ['0', '7', '90']
    # Define significance levels
    def significance_label(p_value):
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        else:
            return None  # Not significant
    # Perform statistical tests
    significance_data = {}
    family_data = data[data['family'] == family]
    timepoints = family_data['timepoint'].unique()
    comparisons = [(a, b) for idx, a in enumerate(timepoints) for b in timepoints[idx + 1:]]
    for tp1, tp2 in comparisons:
        group1 = family_data[family_data['timepoint'] == tp1]['relative_abundance']
        group2 = family_data[family_data['timepoint'] == tp2]['relative_abundance']
        stat, p_value = ttest_ind(group1, group2)
        label = significance_label(p_value)
        if label:
            if family not in significance_data:
                significance_data[family] = []
            significance_data[family].append((tp1, tp2, label))
    
    sns.set(style="whitegrid")
    g = sns.catplot(
        data=data,
        x='timepoint',
        y='relative_abundance',
        hue='Category',  # Differentiate between exposed and control
        kind='box',
        palette='pastel',
        height=8, 
        aspect=0.75,  
        order=timepoint_order,
        flierprops=dict(marker='o', color='black', markersize=6, markerfacecolor='black', markeredgecolor='black')  
    )
    g.set_axis_labels("Time Point", "Relative Abundance")
    g.fig.suptitle(f"Family: {family}", fontsize=26, y=1.05)  
    
    for family, p_values in significance_data.items():
        ax = g.ax
        y = data[data['family'] == family]['relative_abundance'].max()
        h = 0.02 
        increment = 0.03  
        base_y_offset = y + 0.05 
        for i, (tp1, tp2, label) in enumerate(p_values):
            y_offset = base_y_offset + i * increment
            x1, x2 = timepoint_order.index(tp1), timepoint_order.index(tp2)
            ax.plot([x1, x1, x2, x2], [y_offset, y_offset + h, y_offset + h, y_offset], lw=1.5, c='k')
            ax.text((x1 + x2) * .5, y_offset + h, label, ha='center', va='bottom', color='k', fontsize=26)
    
    g.savefig(path_to_save, bbox_inches="tight")
    plt.close()

