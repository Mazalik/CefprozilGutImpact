import seaborn as sns
import matplotlib.pyplot as plt
from plotnine import ggplot, aes, geom_bar, theme, element_text, labs, scale_fill_gradientn
import pandas as pd
import functions_analyze_cefprozil_effect as fn
import sys
from pandas.api.types import CategoricalDtype
from scipy.stats import ttest_ind, f_oneway
import os


def main():

    #Getting the raw data and the project data
    samples_loaded, superkingdom, species_abundance, project,families_table = fn.get_raw_data(sys.argv)

    #Keeping only selected project (in our case: "PRJEB8094")
    selected_project_data = fn.filter_project(project)

    #Adding accession ID from 'samples_loaded' to the selected project data
    selected_project_w_accession = fn.add_accession_id(samples_loaded, selected_project_data)

    #Adding super kingdom information to 'species_abundance' and filtering for rank level of species only
    species_abundance_kingdom = fn.add_kingdom(species_abundance, superkingdom)

    #Keeping only selected uids and desired information in 'species_abundance_kingdom'
    filtered_species_abundance = fn.filter_species_abundance(selected_project_w_accession, species_abundance_kingdom)

    #Match and addition of 'sample_name' to filtered_species_abundance from selected_project_w_accession
    species_abundance_w_sample_name = fn.match_sample_name(selected_project_w_accession,filtered_species_abundance)

    #Adding matching families to species_abundance_w_sample_name
    species_abundance_family = fn.add_family(species_abundance_w_sample_name,families_table)

    #Calculating relative aboundance by families
    species_relative_abundance = fn.cal_relative_abundance(species_abundance_family, selected_project_w_accession)
    fn.relative_abundance_plot(species_relative_abundance,sys.argv)

    #2nd analysis - difference between time points (for each family)
    species_abundance_avg_tagged = fn.add_patient_info(species_relative_abundance, sys.argv)
    control_data = fn.get_control_data(species_abundance_avg_tagged)
    exposed_data = fn.get_exposed_data(species_abundance_avg_tagged)
    significant_families_exposed = fn.plot_families_aboundance(exposed_data, "exposed", sys.argv)
    
    #Filter for significant families in the exposed group and comparison to control group
    significant_exposed_data = exposed_data[exposed_data['family'].isin(significant_families_exposed)]
    significant_control_data = control_data[control_data['family'].isin(significant_families_exposed)]
    for family in significant_families_exposed:
        family_exposed_data = significant_exposed_data[significant_exposed_data['family'] == family]
        family_control_data = significant_control_data[significant_control_data['family'] == family]
        combined_data = pd.concat([family_exposed_data, family_control_data])
        fn.plot_combined_families(combined_data, family, sys.argv)

if __name__ == "__main__":
    main()