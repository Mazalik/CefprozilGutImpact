import functions_analyze_cefprozil_effect as fn
import pandas as pd
import os

def test_get_raw_data():
    # Test loading data with correct paths
    sys_argv = ["analyze_cefprozil_effect.py", "C:\\work\\python\\project"] ###update!###
    samples_loaded, superkingdom, species_abundance, project, families_table = fn.get_raw_data(sys_argv)
    assert isinstance(samples_loaded, pd.DataFrame)
    assert isinstance(superkingdom, pd.DataFrame)
    assert isinstance(species_abundance, pd.DataFrame)
    assert isinstance(project, pd.DataFrame)
    assert isinstance(families_table, pd.DataFrame)

def test_filter_project():
    # Test filtering project data
    project = pd.DataFrame({'project_id': ["PRJEB8094", "OTHER_PROJECT"]})
    filtered_project = fn.filter_project(project)
    assert filtered_project['project_id'].unique() == ["PRJEB8094"]

def test_add_accession_id():
    # Test adding accession IDs
    samples_loaded = pd.DataFrame({'uid': [1, 2, 3], 'accession_id': ["A1", "A2", "A3"]})
    selected_project_data = pd.DataFrame({'loaded_uid': [1, 2], 'project_id': ["PRJEB8094", "PRJEB8094"]})
    result = fn.add_accession_id(samples_loaded, selected_project_data)
    assert 'run_id' in result.columns
    assert result['run_id'].tolist() == ["A1", "A2"]

def test_add_kingdom():
    # Test adding kingdom information
    species_abundance = pd.DataFrame({'taxon_rank_level': ['species', 'genus'], 'uid': [1, 2]})
    superkingdom = pd.DataFrame({'uid': [1, 2], 'kingdom': ['Bacteria', 'Archaea']})
    result = fn.add_kingdom(species_abundance, superkingdom)
    assert 'kingdom' in result.columns
    assert result[result['taxon_rank_level'] == 'species']['kingdom'].tolist() == ['Bacteria']

def test_filter_species_abundance():
    # Test filtering species abundance
    selected_project_w_accession = pd.DataFrame({'loaded_uid': [1, 2]})
    species_abundance_kingdom = pd.DataFrame({'loaded_uid': [1, 2, 3], 'scientific_name': ['A', 'B', 'C'], 'relative_abundance': [0.1, 0.2, 0.3]})
    result = fn.filter_species_abundance(selected_project_w_accession, species_abundance_kingdom)
    assert result['scientific_name'].tolist() == ['A', 'B']
    assert 'Other' not in result['scientific_name']

def test_match_sample_name():
    # Test matching sample names
    selected_project_w_accession = pd.DataFrame({'loaded_uid': [1, 2], 'sample_name': ['Sample1', 'Sample2']})
    filtered_species_abundance = pd.DataFrame({'loaded_uid': [1, 2], 'scientific_name': ['A', 'B'], 'relative_abundance': [0.1, 0.2]})
    result = fn.match_sample_name(selected_project_w_accession, filtered_species_abundance)
    assert 'sample_name' in result.columns
    assert result['sample_name'].tolist() == ['Sample1', 'Sample2']

def test_output_files_creation():
    # Define the output folder and the expected files
    output_folder = "output"
    expected_files = ["families_exposed_plot.png", "relative_abundance_plot.png"]
    
    # Check if the output folder exists
    assert os.path.exists(output_folder), f"Output folder '{output_folder}' does not exist."
    
    for file in expected_files:
        file_path = os.path.join(output_folder, file)
        assert os.path.exists(file_path), f"File not found: {file}"


def test_no_null_values():
    samples_loaded, superkingdom, species_abundance, project, families_table = fn.get_raw_data(["analyze_cefprozil_effect.py", "C:\\work\\python\\project"]) ####לעדכן####
    selected_project_data = fn.filter_project(project)
    selected_project_w_accession = fn.add_accession_id(samples_loaded, selected_project_data)
    species_abundance_kingdom = fn.add_kingdom(species_abundance, superkingdom)
    filtered_species_abundance = fn.filter_species_abundance(selected_project_w_accession, species_abundance_kingdom)
    species_abundance_w_sample_name = fn.match_sample_name(selected_project_w_accession, filtered_species_abundance)
    species_abundance_family = fn.add_family(species_abundance_w_sample_name, families_table)
    assert species_abundance_family.isnull().sum().sum() == 0, "Null values found in the processed data"


def test_data_types():
    samples_loaded, superkingdom, species_abundance, project, families_table = fn.get_raw_data(["analyze_cefprozil_effect.py", "C:\\work\\python\\project"])####לעדכן####
    selected_project_data = fn.filter_project(project)
    selected_project_w_accession = fn.add_accession_id(samples_loaded, selected_project_data)
    species_abundance_kingdom = fn.add_kingdom(species_abundance, superkingdom)
    filtered_species_abundance = fn.filter_species_abundance(selected_project_w_accession, species_abundance_kingdom)
    species_abundance_w_sample_name = fn.match_sample_name(selected_project_w_accession, filtered_species_abundance)
    species_abundance_family = fn.add_family(species_abundance_w_sample_name, families_table)
    assert species_abundance_family['relative_abundance'].dtype == float, "Relative abundance column is not of type float"
    assert species_abundance_family['family'].dtype == object, "Family column is not of type object"



