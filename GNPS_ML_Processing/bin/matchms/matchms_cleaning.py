import os
from matchms.filtering.default_pipelines import LIBRARY_CLEANING
from matchms.Pipeline import Pipeline, create_workflow
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

    # Current
    if not os.path.exists(os.path.join(results_folder, "library_cleaning.yaml")):
        workflow = create_workflow( yaml_file_name=os.path.join(results_folder, "library_cleaning.yaml"),
                                    query_filters=LIBRARY_CLEANING + [("derive_annotation_from_compound_name",
                                                                    {"annotated_compound_names_file":new_cached_compound_name_annotation_path, "mass_tolerance": 0.1})])
    else:
        workflow = load_workflow_from_yaml_file(os.path.join(results_folder, "library_cleaning.yaml"))

    pipeline = Pipeline(workflow,
                        logging_file=os.path.join(results_folder, "library_cleaning_log.log"),
                        logging_level="WARNING")

    temp_output = os.path.join(results_folder, "cleaned_spectra.mgf")
    if os.path.exists(temp_output):
        os.remove(temp_output)

    pipeline.run(args.input_mgf_path, cleaned_query_file=temp_output)
    
    # Move the temporary output to the final position
    shutil.move(temp_output, os.path.join(results_folder, "cleaned_spectra.mgf"))
    
if __name__ == "__main__":
    main()