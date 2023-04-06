import argparse
from bs4 import BeautifulSoup as bs
from lxml import etree
import os
import numpy as np
import obonet
import networkx
import re
import json
import csv

'''
Standard command line usage: nohup python3 GNPS_MGF_maker_lxml.py ALL_GNPS -s -p 10
'''

# Get ids for parsing mzML files
graph = obonet.read_obo('https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo')
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
model_ids = networkx.ancestors(graph, 'MS:1000031')
instrument_vendor_ids = ['MS:1001269'] # Seems that this is not always used
ionization_methods_ids = networkx.ancestors(graph, 'MS:1000008')
detector_ids = networkx.ancestors(graph, 'MS:1000026')
dissociation_ids = networkx.ancestors(graph, 'MS:1000044')
analyzer_ids = networkx.ancestors(graph,'MS:1000443')
retention_id = 'MS:1000894'
collision_energy_id = 'MS:1000045'

def helper(process_num, scan_start, all_spectra_list):
    scan = scan_start
    file_not_found_count = 0
    UnicodeDecodeError_count = 0
    mgf_file_count = 0
    t = all_spectra_list

    parser = etree.XMLParser(huge_tree=True)

    output_mgf = open('./temp/temp_{}.mgf'.format(process_num), "w")
    output_csv = open('./temp/temp_{}.csv'.format(process_num), "w")
    w = csv.DictWriter(output_csv, ['scan', 'spectrum_id','collision_energy','retention_time','Adduct','Compound_Source','Precursor_MZ','ExactMass','Charge','Ion_Mode','Smiles','INCHI','InChIKey_smiles','InChIKey_inchi','msModel','msManufacturer','msDetector','msMassAnalyzer','msIonisation','msDissociationMethod','GNPS_library_membership','GNPS_Inst'])
    if process_num == 0: w.writeheader()    # This way when we merge the csvs we only get one header
    
    for spectrum in t:
        spectrum["new_scan"] = scan
        summary_dict = {}
        
        output_mgf.write("BEGIN IONS\n")
        output_mgf.write("PEPMASS={}\n".format(spectrum["Precursor_MZ"]))
        output_mgf.write("CHARGE={}\n".format(spectrum["Charge"]))
        output_mgf.write("MSLEVEL={}\n".format(2))
        output_mgf.write("TITLE="+str(spectrum.get("spectrum_id")+"\n"))
        output_mgf.write("SCANS={}\n".format(scan))

        peaks = json.loads(spectrum["peaks_json"])
        for peak in peaks:
            output_mgf.write("{} {}\n".format(peak[0], peak[1]))

        output_mgf.write("END IONS\n")
        
        # Get Collision Energy From Original Upload File
        task = spectrum['task']
        source_file = spectrum["source_file"].replace(';','')
        # print(spectrum)
        extension = str(source_file.split('.')[-1].strip())
        source_file = source_file.split('/')[-1]
        source_scan = spectrum.get("scan")
        if source_scan is None or source_scan == -1:
            raise IOError("Expected source scan to be specified but it was not")
        try:
            if extension == 'mzXML':
                # f, _ = request.urlretrieve('ftp://ccms-ftp.ucsd.edu/GNPS_Library_Provenance/{}/{}'.format(task,source_file))
                f = "/home/user/LabData/GNPS_Library_Provenance/{}/{}".format(task,source_file)
                # data_dict = list(mzxml.read(f, read_schema=True))[0]
                
                '''
                data = etree.fromstring(open(f).read())
                spectrum = data.find('.//{*}scan[@num="' + source_scan '"])'
                # mzXML files may or may not have an msInstrumentID independently of whether or not they have a collision energy
                '''
                
                bs_data = bs(open(f).read(), 'xml')
                run = bs_data.find('msRun')
                if run is not None and int(run['scanCount']) > 1:
                    bs_scan = bs_data.find('scan', {"num":source_scan})#data_dict.get("collisionEnergy")
                else:
                    bs_scan = bs_data.find('scan')#data_dict.get("collisionEnergy")
                
                if bs_scan is None:
                    print("WARNING: Unable to find scan {} within {}".format(source_scan,source_file))
                else:
                    ce = bs_scan.get('collisionEnergy')
                    rt = bs_scan.get('retentionTime')
                    if ce is not None: summary_dict["collision_energy"] = ce
                    if rt is not None: summary_dict["retention_time"] = rt
                    msInstrumentID = None
                    msInstrumentID = bs_scan.get('msInstrumentID')   
                    precursorMZ =  bs_scan.find('precursorMz')
                    if precursorMZ is not None:
                        msDissociationMethod =  precursorMZ.get('activationMethod')
                    
                    if msInstrumentID is None:
                        instrument_info = bs_data.find('msInstrument')
                    else:
                        instrument_info = bs_data.find('msInstrument', {'msInstrumentID':None})
                    if instrument_info is not None:
                        msModel        = instrument_info.find('msModel')
                        msManufacturer = instrument_info.find('msManufacturer')
                        msIonisation   = instrument_info.find('msIonisation')
                        msMassAnalyzer = instrument_info.find('msMassAnalyzer')
                        msDetector     = instrument_info.find('msDetector')
                        

                        if msModel is not None: summary_dict["msModel"] = msModel.get('value')
                        if msManufacturer is not None: summary_dict["msManufacturer"] = msManufacturer.get('value')
                        if msIonisation is not None: summary_dict["msIonisation"] = msIonisation.get('value')
                        if msMassAnalyzer is not None: summary_dict["msMassAnalyzer"] = [msMassAnalyzer.get('value')]
                        if msDetector is not None: summary_dict["msDetector"] =msDetector.get('value')
                        if msDissociationMethod is not None: summary_dict['msDissociationMethod'] = msDissociationMethod
                        
            elif extension == 'mzML':
                f = "/home/user/LabData/GNPS_Library_Provenance/{}/{}".format(task,source_file)
                data = etree.parse(f, parser).getroot()
                instrument_conigurations = data.find('.//{*}instrumentConfigurationList')
                num_instrument_configurations = data.find('.//{*}instrumentConfigurationList').attrib['count']
                if num_instrument_configurations == 1: # We'll just shortcut this whole mess if there is only one configuration
                    config_for_specturm = instrument_conigurations.find('.//{*}instrumentConfiguration')
                else:
                    run_list = data.findall('.//{*}run')
                    run = None
                    xml_spectrum=None
                    for r in run_list:
                        spectra_list = data.findall('.//{*}spectrum') # Since etree doesn't support regex, we will have to narrow down the scan numbers regularly
                        xml_spectrum=None
                        for spec in spectra_list:
                            if re.search('scan={}$'.format(source_scan), spec.attrib['id']) is not None:
                                xml_spectrum=spec
                                break
                        if xml_spectrum is not None:
                            run = r
                            break
                    if r is None:
                        raise IOError("Run Not Found.")
                    # A default instrument config is attributed to each run, however this can be superseeded for a particular scan    
                    instrument_config = None
                    
                    scanList = xml_spectrum.find('.//{*}scanList')
                    if scanList.attrib.get('count')!= '1': raise NotImplementedError("Parising has not been implemented for more than one scan.")
                    s = scanList.find('.//{*}scan')
                    instrument_config = s.attrib.get("instrumentConfigurationRef")
                    
                    if instrument_config is None:
                        instrument_config = run.attrib['defaultInstrumentConfigurationRef']
                        assert instrument_config is not None
                    config_for_specturm = instrument_conigurations.find('.//{*}instrumentConfiguration[@id="' + instrument_config + '"]')
                # bs_spectrum = bs_data.find('spectrum', {"id": re.compile('scan={}$'.format(source_scan))})
                
                if config_for_specturm is not None:
                    accessions = [x.attrib['accession'] for x in config_for_specturm.findall('.//{*}cvParam')]
                    
                    # The following are stored in the instrument config
                    for id in model_ids:
                        if id in accessions:
                            summary_dict["msModel"] = id_to_name[id]
                            break
                    for id in instrument_vendor_ids:
                        if id in accessions:
                            summary_dict["msManufacturer"] = id_to_name[id]
                            break
                    for id in ionization_methods_ids:
                        if id in accessions:
                            summary_dict["msIonisation"] = id_to_name[id]
                            break
                    for id in detector_ids:
                        if id in accessions:
                            summary_dict["msDetector"] = id_to_name[id]
                            break
                    summary_dict['msMassAnalyzer'] = [] # we can have multiple mass analyzers for hybrid instruments
                    for id in analyzer_ids:
                        if id in accessions:
                            summary_dict['msMassAnalyzer'].append(id_to_name[id])
                    if len(summary_dict['msMassAnalyzer']) == 0:
                        del summary_dict['msMassAnalyzer']

                # The following values are stored in the spectrum config
                accessions = {x.attrib['accession']: x.attrib.get('value') for x in xml_spectrum.findall('.//{*}cvParam')}
                ce = accessions.get(collision_energy_id)
                if ce is not None: summary_dict["collision_energy"] = ce
                rt = accessions.get(retention_id)
                if rt is not None: summary_dict['retention_time'] = rt
                
                for id in dissociation_ids: # Dissociation is specified with the precursor, not the instrument config
                    if accessions.get(id) is not None:
                        summary_dict['msDissociationMethod'] = id_to_name[id]
                        break
            elif extension.lower() == 'mgf':
                # Right now we are assuming that mgf files do not contain collision energies or retention time
                mgf_file_count +=1
                # if process_num == 0:
                #     t.set_description("(Files not Found: {})(Unicode Decode Error: {})(Unparsed MGFs: {})".format(file_not_found_count, UnicodeDecodeError_count, mgf_file_count))
            
            else:
                print("Warning: Unkown file extension: [{}] ".format(extension))
                
        except KeyError as k:
            # print("Failed to retrieve collisionEnergy or retentionTime from https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID={}".format(spectrum["spectrum_id"]))
            # print("Path to File:", "/home/user/LabData/GNPS_Library_Provenance/{}/{}".format(task,source_file))
            # print(spectrum)
            # raise k 
            pass
        except TypeError as e:
            # Usually this happens when the filetype is not what we expected. There are some MGF files saved as mzXML       
            # print("Warning: An error occured while retrieving collisionEnergy or retentionTime from https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID={}. These field will be left blank.".format(spectrum["spectrum_id"]))
            # print(spectrum["scan"])
            print("TypeError: Path to File:", "/home/user/LabData/GNPS_Library_Provenance/{}/{}".format(task,source_file))
            print(e)

        except FileNotFoundError as e:
            file_not_found_count +=1
            # if process_num == 0:
            #     t.set_description("(Files not Found: {})(Unicode Decode Error: {})(Unparsed MGFs: {})".format(file_not_found_count, UnicodeDecodeError_count, mgf_file_count))

            
        except UnicodeDecodeError as e:
            UnicodeDecodeError_count += 1
            # if process_num == 0:
            #     t.set_description("(Files not Found: {})(Unicode Decode Error: {})(Unparsed MGFs: {})".format(file_not_found_count, UnicodeDecodeError_count, mgf_file_count))
            
        except Exception as e:
            print("Warning: An error occured while retrieving collisionEnergy or retentionTime from https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID={}. These field will be left blank.".format(spectrum["spectrum_id"]))
            print("Path to File:", "/home/user/LabData/GNPS_Library_Provenance/{}/{}".format(task,source_file))
            print(e) 

        if summary_dict.get('msIonisation') is None: summary_dict['msIonisation'] = spectrum.get("Ion_Source")
        summary_dict["spectrum_id"] = spectrum.get("spectrum_id")
        summary_dict["scan"] = scan
        summary_dict['Adduct'] = spectrum.get('Adduct')
        summary_dict['Compound_Source'] = spectrum.get('Compound_Source')
        summary_dict["Precursor_MZ"] = spectrum.get("Precursor_MZ")
        summary_dict["ExactMass"] = spectrum.get("ExactMass")
        summary_dict["Charge"] = spectrum.get("Charge")
        summary_dict["Ion_Mode"] = spectrum.get("Ion_Mode")
        summary_dict["Smiles"] = spectrum.get("Smiles")
        summary_dict["INCHI"] = spectrum.get("INCHI")
        summary_dict["InChIKey_smiles"] = spectrum.get("InChIKey_smiles")
        summary_dict["InChIKey_inchi"] = spectrum.get("InChIKey_inchi")
        summary_dict["GNPS_Inst"] = spectrum.get("Instrument")
        summary_dict["GNPS_library_membership"] = spectrum.get("library_membership")

        w.writerow(summary_dict)
        
        scan += 1
    output_mgf.close()
    output_csv.close()
    return file_not_found_count, UnicodeDecodeError_count, mgf_file_count

def main():
    parser = argparse.ArgumentParser(description='Process some GNPS2 Entries.')
    parser.add_argument('-f', '--file', help="path to params file", required=True)
    args = parser.parse_args()

    if not os.path.isdir('./temp'): os.makedirs('./temp')
    spectra = np.load(args.file, allow_pickle=True)
    # A regex that captures the values from 'params/params_{}_{}.npy'
    regex = r"params\/params_([^_]+)_([^\.]+)\.npy"
    match = re.search(regex, args.file)
    p_idx = int(match.group(1))
    start_scan = int(match.group(2))
    
    helper(p_idx, start_scan, spectra)
    
if __name__ == '__main__':
    main()