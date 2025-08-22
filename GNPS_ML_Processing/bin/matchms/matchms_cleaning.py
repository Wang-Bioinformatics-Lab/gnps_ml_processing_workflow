import os
from matchms.filtering.default_pipelines import LIBRARY_CLEANING
from matchms.Pipeline import Pipeline, create_workflow
import matchms.filtering as msfilters
import argparse
import shutil

from matchms.yaml_file_functions import load_workflow_from_yaml_file

def main():
    parser = argparse.ArgumentParser(description='Run matchms library cleaning')
    parser.add_argument('--input_mgf_path', type=str, help='Input mgf file')
    parser.add_argument('--cached_compound_name_annotation_path', type=str, help='Path to cache compound name to annotation mapping')
    parser.add_argument('--output_path', type=str, help='Output folder', default="./results_library_cleaning")
    args = parser.parse_args()
    
    results_folder = args.output_path
    os.makedirs(results_folder, exist_ok=True)
    
    # Original (Keeping here to hopefully add compound class back in later)
    # compound_class_csv_file = os.path.join(results_folder, "compound_name_annotation.csv")

    # workflow = create_workflow( yaml_file_name=os.path.join(results_folder, "library_cleaning.yaml"),
    #                             query_filters=LIBRARY_CLEANING + [("derive_annotation_from_compound_name",
    #                                                              {"annotated_compound_names_file": compound_class_csv_file, "mass_tolerance": 0.1})])
    # compound_class_csv_file = os.path.join(results_folder, "compound_name_annotation.csv")

    # Copy the cached compound name annotation file to the results folder
    new_cached_compound_name_annotation_path = "compound_name_annotation.csv"
    shutil.copy(args.cached_compound_name_annotation_path, new_cached_compound_name_annotation_path)

    
    workflow = create_workflow(
        yaml_file_name=os.path.join(results_folder, "matchms_pipeline_settings.yaml"),
        query_filters=[
            msfilters.require_correct_ms_level,
            msfilters.make_charge_int,
            msfilters.add_compound_name,
            msfilters.derive_adduct_from_name,
            msfilters.derive_formula_from_name,
            msfilters.clean_compound_name,
            msfilters.interpret_pepmass,
            msfilters.add_precursor_mz,
            msfilters.add_retention_index,
            msfilters.add_retention_time,
            msfilters.derive_ionmode,
            msfilters.correct_charge,
            msfilters.require_precursor_mz,
            msfilters.harmonize_undefined_inchikey,
            msfilters.harmonize_undefined_inchi,
            msfilters.harmonize_undefined_smiles,
            msfilters.repair_inchi_inchikey_smiles,
            msfilters.clean_adduct,
            msfilters.add_parent_mass,
            # (msfilters.derive_annotation_from_compound_name, {'annotated_compound_names_file': os.path.join(results_folder, new_cached_compound_name_annotation_path)}),
            msfilters.derive_smiles_from_inchi,
            msfilters.derive_inchi_from_smiles,
            msfilters.derive_inchikey_from_inchi,
            (msfilters.repair_smiles_of_salts, {"mass_tolerance": 0.1}),
            (msfilters.repair_parent_mass_is_molar_mass, {"mass_tolerance": 0.1}),
            (msfilters.repair_adduct_and_parent_mass_based_on_smiles, {"mass_tolerance": 0.1}),
            msfilters.repair_not_matching_annotation,
            msfilters.require_valid_annotation,
            (msfilters.require_correct_ionmode, {"ion_mode_to_keep": "both"}),
            (msfilters.require_parent_mass_match_smiles, {"mass_tolerance": 0.1}),
            msfilters.require_matching_adduct_precursor_mz_parent_mass,
            msfilters.require_matching_adduct_and_ionmode,
            msfilters.normalize_intensities,
    ])

    pipeline = Pipeline(workflow,
                        logging_file=os.path.join(results_folder, "library_cleaning_log.log"),
                        logging_level="WARNING")
    

    temp_output = os.path.join(results_folder, "cleaned_spectra.mgf")
    if os.path.exists(temp_output):
        os.remove(temp_output)
        

    pipeline.run(args.input_mgf_path, cleaned_query_file=temp_output)
    
    # Move the temporary output to the final position
    # Check if exists
    if os.path.exists(temp_output):
        print(f"Library cleaning completed. Cleaned spectra saved to {temp_output}")
        shutil.move(temp_output, os.path.join(results_folder, "cleaned_spectra.mgf"))
    else:
        print("No cleaned spectra were generated.")
    
if __name__ == "__main__":
    main()