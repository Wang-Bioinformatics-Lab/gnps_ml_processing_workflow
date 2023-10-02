# GNPS2 ML Data Information

For useful terminology, see this [link](https://www.cs.ucr.edu/~mingxunw/terminology/)
## Data Dict
* 'scan': A numerical index of the entry in the dataset (zero-indexed).
* 'spectrum_id': The ccms spectrum id
* 'collision_energy': Collision energy reported in eV, if applicable
* 'Adduct':  The molecule or ion that forms as a result of a chemical reaction the substance being analyzed and another molecule or ion during the ionization process
* 'Compound_Source': Contains information about the source of the compound
* 'Compund_Name'
* 'Precursor_MZ': The mass to charge ratio of the precursor ion
* 'ExactMass': The exact mass of the precusor ion
* 'Charge': The charge of the precursor ion
* 'Ion_Mode': Whether the mass spectrometer was in positive or negative ion mode
* 'Smiles': The spectrum annotation in SMILES format
* 'INCHI': The spectrum annotation in INCHI format
* 'InChIKey_smiles': The InChIKey generated from the SMILES string
* 'msManufacturer': The manufacturer of the mass spec instrument (parsed from the mzML/mzXML file, if available)
* 'msMassAnalyzer': The mass analyzer used in the instrument (parsed from the mzML/mzXML file, if available)
* 'msIonisation': The ionization method used (parsed from the mzML/mzXML file, if available)
* 'msDissociationMethod': The dissociation method used (parsed from the mzML/mzXML file, if available)
* 'GNPS_library_membership': The GNPS library the spectrum is in see this [link](https://gnps.ucsd.edu/ProteoSAFe/libraries.jsp) for details
* 'ppmBetweenExpAndThMass': The PPM difference between the SMILES annotation and the experimental value

## Additional Information
Collision energies  were parsed out of the original mzML and mzXML files. A large number of spectra are backed by MGF files which do not contain this information.

When SMILES values are unlabeled we impute using the INCHI value, if provided by the user. In any case where INCHI and SMILES conflict, we prefer the INCHI value.

In addition, unecessary charges, salts, and stereochemisty have been removed from annotation.

Adduct values have all been unified in the following format: [M\<+/-\>\<Adduct\>...\<+/-\>\<Adduct\>]\<Charge\><+/-> and common adduct abreviations such as fluoroacetic acid (commonly abrreviated as FA) have been expanded into their chemical formulas.

For inherently charge molecules such as `O=CCCCCCCCc1cccc[n+]1CCCCCC=O`, we validate that the adduct is correct. 291 structures we're changed from an M+H adduct to [M]1+. Any rows with annotations that could be neutralized with protonation were dropped, meaning all annotations are either neutral, or inherently charged.

Vendor and instrument information has been ascribed to each entry with presedence going to vendor information in mzML and mzXML files. If not available, we use the instrument model name to impute it. The GNPS user instrument abbotations have also been propogated to the msMassAnalyzer, msManufacturer, and other ms* columns wherever possible. 

Many columns such as adduct, ionization method, mass analyzer, ion mode, and more have been manually cleaned of spurious entries. For example, in the Ion_Mode column, collision energies have been occasionally listed. These have been removed and the only entries are 'positive' and 'negative'.

### Coming Soon
* SMILES tautamers will be unified to one single tautamer.
* Additional rows
* Additional collision energies
