import json
import pandas as pd
import numpy as np
import os
import sys
import warnings
import argparse
from pathlib import Path
import logging

def integrate_Classyfire(df:pd.DataFrame, cache_path:Path) -> pd.DataFrame:
    """This function integrates Classyfire information into the input DataFrame.
    For each InChIKey, it retrieves the 'kingdom', 'superclass', 'class', 'subclass',
    and 'direct_parent' information from the Classyfire cache and merges it into the DataFrame.

    Parameters:
    df (pd.DataFrame): The input DataFrame containing the data to be enriched.
    cache_path (Path): The path to the Classyfire cache file.

    Returns:
    pd.DataFrame: The enriched DataFrame with Classyfire information.
    """
    logging.info("Integrating Classyfire information...")
    DESIRED_COLUMNS = ['kingdom', 'superclass', 'class', 'subclass', 'direct_parent']  # To Integrate
    REQUIRED_COLUMNS = ['InChIKey_14']  # To Merge

    # Integrate Classyfire information into the DataFrame
    classyfire_dirs = cache_path.glob("**/*.json")
    # File name is the InChIKey[:14]
    data = []
    # Read the JSON files and create a DataFrame
    for file in classyfire_dirs:
        with open(file, 'r') as f:
            d = json.load(f)
            o = {}
            o['InChIKey_14'] = file.stem[:14]
            _kingdom = d.get('kingdom', None)
            _superclass = d.get('superclass', None)
            _class = d.get('class', None)
            _subclass = d.get('subclass', None)
            _direct_parent = d.get('direct_parent', None)

            o['kingdom'] = None
            if _kingdom is not None:
                o['kingdom'] = _kingdom.get('name', '')
                o['kingdom'] = o['kingdom'].replace('"', '')
                o['kingdom'] = o['kingdom'].replace("'", '')
                o['kingdom'] = o['kingdom'].replace(",", '')

            o['superclass'] = None
            if _superclass is not None:
                o['superclass'] = _superclass.get('name', '')
                o['superclass'] = o['superclass'].replace('"', '')
                o['superclass'] = o['superclass'].replace("'", '')
                o['kingdom'] = o['kingdom'].replace(",", '')

            o['class'] = None
            if _class is not None:
                o['class'] = _class.get('name', '')
                o['class'] = o['class'].replace('"', '')
                o['class'] = o['class'].replace("'", '')
                o['kingdom'] = o['kingdom'].replace(",", '')

            o['subclass'] = None
            if _subclass is not None:
                o['subclass'] = _subclass.get('name', '')
                o['subclass'] = o['subclass'].replace('"', '')
                o['subclass'] = o['subclass'].replace("'", '')
                o['kingdom'] = o['kingdom'].replace(",", '')

            o['direct_parent'] = None
            if _direct_parent is not None:
                o['direct_parent'] = _direct_parent.get('name', '')
                o['direct_parent'] = o['direct_parent'].replace('"', '')
                o['direct_parent'] = o['direct_parent'].replace("'", '')
                o['kingdom'] = o['kingdom'].replace(",", '')

            data.append(o)

    # Create a DataFrame from the list of dictionaries
    classyfire_df = pd.DataFrame(data)
    classyfire_df = classyfire_df[REQUIRED_COLUMNS + DESIRED_COLUMNS]
    # Rename desired cols with suffix
    classyfire_df.rename(columns={col: f"classyfire_{col}" for col in DESIRED_COLUMNS}, inplace=True)

    # InChIKey_smiles is the full InChIKey in df
    df['InChIKey_14'] = df['InChIKey_smiles'].str[:14]
    logging.info(f"Found {len(df['InChIKey_14'].unique())} unique InChIKey_14 in the input DataFrame.")
    logging.info(f"Found {len(classyfire_df['InChIKey_14'].unique())} unique InChIKey_14 in the Classyfire cache.")
    # Merge the DataFrames on the InChIKey_14 column
    logging.info(f"Started with {len(df):,} rows in the input DataFrame.")
    df = df.merge(classyfire_df, on='InChIKey_14', how='left', suffixes=('','_classyfire'))
    logging.info(f"After merge, the DataFrame has {len(df):,} rows.") # Should not change

    # For each DESIRED_COLUMNS, count number of nan and non nan values
    for col in DESIRED_COLUMNS:
        full_col_name = f"classyfire_{col}"
        logging.info(f"Column {full_col_name} has {df[full_col_name].isna().sum():,} NaN values and {df[full_col_name].notna().sum():,} non-NaN values.")

    # Drop the InChIKey_14 column
    df.drop(columns=['InChIKey_14'], inplace=True)

    return df

def integrate_ChemInfoService(df:pd.DataFrame, cache_path:Path) -> pd.DataFrame:
    """
    This function integrates ChemInfoService information into the input DataFrame.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame containing the data to be enriched.
    cache_path (Path): The path to the ChemInfoService cache file.

    Returns:
    pd.DataFrame: The enriched DataFrame with ChemInfoService information.
    """
    logging.info("Integrating ChemInfoService information...")
    DESIRED_COLUMNS = ['nplikeness']    # To Integrate
    REQUIRED_COLUMNS = ['InChIKey_14']  # To Merge

    # Integrate nplikeness column into the DataFrame

    np_classifier_dirs = cache_path.glob("**/*.json")
    # File name is the InChIKey[:14]
    data = []
    # Read the JSON files and create a DataFrame
    for file in np_classifier_dirs:
        with open(file, 'r') as f:
            d = json.load(f)
            d['InChIKey_14'] = file.stem[:14]
            data.append(d)

    # Create a DataFrame from the list of dictionaries
    np_classifier_df = pd.DataFrame(data)
    np_classifier_df = np_classifier_df[REQUIRED_COLUMNS + DESIRED_COLUMNS]
    # Rename desired cols with suffix
    np_classifier_df.rename(columns={col: f"np_classifier_{col}" for col in DESIRED_COLUMNS}, inplace=True)

    # InChIKey_smiles is the full InChIKey in df
    
    df['InChIKey_14'] = df['InChIKey_smiles'].str[:14]
    logging.info(f"Found {len(df['InChIKey_14'].unique())} unique InChIKey_14 in the input DataFrame.")
    logging.info(f"Found {len(np_classifier_df['InChIKey_14'].unique())} unique InChIKey_14 in the ChemInfoService cache.")
    # Merge the DataFrames on the InChIKey_14 column
    logging.info(f"Started with {len(df):,} rows in the input DataFrame.")
    df = df.merge(np_classifier_df, on='InChIKey_14', how='left', suffixes=('','_np_classifier'))
    logging.info(f"After merge, the DataFrame has {len(df):,} rows.") # Should not change

    # For each DESIRED_COLUMNS, count number of nan and non nan values
    for col in DESIRED_COLUMNS:
        full_col_name = f"np_classifier_{col}"
        logging.info(f"Column {full_col_name} has {df[full_col_name].isna().sum():,} NaN values and {df[full_col_name].notna().sum():,} non-NaN values.")
    
    # Drop the InChIKey_14 column
    df.drop(columns=['InChIKey_14'], inplace=True)

    return df

def integrate_Npclassifier(df:pd.DataFrame, cache_path:Path) -> pd.DataFrame:
    """This function is no-op until NP-classifier information is desired.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame containing the data to be enriched.
    cache_path (Path): The path to the NP-classifier cache file.
    
    Returns:
    pd.DataFrame: The enriched DataFrame with NP-classifier information.
    """
    # Currently, this function is a no-op
    # If NP-classifier information is desired, implement the integration logic here
    logging.info("Skipping NP-classifier integration.")
    return df

def enrich_with_api_info(df:pd.DataFrame, api_cache_path:Path) -> pd.DataFrame:
    """
    This function enriches the input DataFrame with API information by merging certain columns from 
    the API cache using the first block of the InChiKey as the key.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame containing the data to be enriched.
    api_cache_path (Path): The path to the API cache file.

    Returns:
    pd.DataFrame: The enriched DataFrame with API information.
    """
    # If "Classyfire" is a subdirectory of the API cache path, integrate Classyfire information
    if (api_cache_path / "Classyfire").exists():
        # Integrate Classyfire information
        df = integrate_Classyfire(df, api_cache_path / "Classyfire")
    else:
        logging.warning("Classyfire directory not found in the API cache path.")

    # If "ChemInfoService" is a subdirectory of the API cache path, integrate ChemInfoService information
    if (api_cache_path / "ChemInfoService").exists():
        df = integrate_ChemInfoService(df, api_cache_path / "ChemInfoService")
    else:
        logging.warning("ChemInfoService directory not found in the API cache path.")
    
    # If "Npclassifier" is a subdirectory of the API cache path, integrate Npclassifier information
    if (api_cache_path / "Npclassifier").exists():
        df = integrate_Npclassifier(df, api_cache_path / "Npclassifier")
    else:
        logging.warning("Npclassifier directory not found in the API cache path.")

    return df

def main():
    parser = argparse.ArgumentParser(description="Integrate API information into a CSV file.")
    parser.add_argument("--input_csv_path", type=str, required=True, help="Path to the input CSV file.")
    parser.add_argument("--output_csv_path", type=str, required=True, help="Path to the output CSV file.")
    parser.add_argument("--api_cache_path", type=str, required=True, help="Path to the API information file.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Dump all arguments to logs
    logging.info("Arguments: %s", args)
    for arg, value in vars(args).items():
        logging.info(f"{arg}: {value}")


    input_csv_path = Path(args.input_csv_path)
    output_csv_path = Path(args.output_csv_path)
    api_cache_path = Path(args.api_cache_path)

    if not input_csv_path.exists():
        raise FileNotFoundError(f"Input CSV file not found: {input_csv_path}")
    if not api_cache_path.exists():
        raise FileNotFoundError(f"API cache file not found: {api_cache_path}")
    if not output_csv_path.parent.exists():
        output_csv_path.parent.mkdir(parents=True, exist_ok=True)

    input_df = pd.read_csv(input_csv_path)

    enriched_df = enrich_with_api_info(input_df, api_cache_path)

    enriched_df.to_csv(output_csv_path, index=False)

if __name__ == "__main__":
    main()