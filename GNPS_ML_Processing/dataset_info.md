# GNPS2 ML Data Information

For useful terminology, see this [link](https://www.cs.ucr.edu/~mingxunw/terminology/)
## Data Dict
* 'scan': A numerical index of the entry in the dataset (zero-indexed).
* 'spectrum_id': The ccms spectrum id
* 'collision_energy': Collision energy reported in eV, if applicable
* 'retention_time': Retention time stored in seconds
* 'Adduct':  The molecule or ion that forms as a result of a chemical reaction the substance being analyzed and another molecule or ion during the ionization process
* 'Compound_Source': Contains information about the source of the compound
* 'Precursor_MZ': The mass to charge ratio of the precursor ion
* 'ExactMass': The exact mass of the precusor ion
* 'Charge': The charge of the precursor ion
* 'Ion_Mode': Whether the mass spectrometer was in positive or negative ion mode
* 'Smiles': The spectrum annotation in SMILES format
* 'INCHI': The spectrum annotation in INCHI format
* 'InChIKey_smiles': 
* 'InChIKey_inchi':
* 'msModel': The model of the mass spec instrument (parsed from the mzML/mzXML file, if available)
* 'msManufacturer': The manufacturer of the mass spec instrument (parsed from the mzML/mzXML file, if available)
* 'msDetector': The mass detector used in the instrument (parsed from the mzML/mzXML file, if available)
* 'msMassAnalyzer': The mass analyzer used in the instrument (parsed from the mzML/mzXML file, if available)
* 'msIonisation': The ionization method used (parsed from the mzML/mzXML file, if available)
* 'msDissociationMethod': The dissociation method used (parsed from the mzML/mzXML file, if available)
* 'GNPS_library_membership': The GNPS library the spectrum is in see this [link](https://gnps.ucsd.edu/ProteoSAFe/libraries.jsp) for details
* 'GNPS_Inst': The user specified information about the instrument. Note that this column is not cleaned
* 'ppmBetweenExpAndThMass': The PPM difference between the SMILES annotation and the experimental value
* 'explainable_intensity': The percentage of intensity in the spectra that can be explained by combinatorially enumerating the annotation's chemical formula (not currently in use). 

## Additional Information
Collision energies and retention times we're parsed out of the original mzML and mzXML files. A large number of spectra are backed by MGF files which do not contain this information.

In rare cases some spectra had INCHI values but no SMILES values. In such cases SMILES was imputed using the INCHI value.

SMILES tautamers have been unified to one single tautamer (not currently applied). In addition, unecessary charges are removed, salts have been removed (keeping the larger of the molecules), and stereochemisty is removed.

Adduct values have all been unified in the following format: [M<+/-><Adduct>...<+/-><Adduct>]<Charge><+/->

The GNPS instrument field which contains the user specified information has been propogated to the msModel, msManufacturer, and other ms* columns wherever possible. 

Many columns such as adduct, ionization method, mass analyzer, ion mode, and more have been manually cleaned of spurious entries.
