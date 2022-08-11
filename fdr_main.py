import os
import shutil
import sys
from sipros_peptides_assembling_mod2 import peptides_assembling, peptides_assembling_lite
import multiprocessing
import sipros_post_module
from functools import partial
from pre_fdr_v3_p import percolator_to_comet, divide_species_auto_p, divide_species_dic_p
from pre_fdr_v3_c import combine_txt_files, clean_comet, divide_species_auto_c, divide_species_dic_c
import re
import glob

# import time
# import matplotlib.pyplot as plt


# IO
path_to_temp_files = './temp'
path_to_species_pin = './temp/species_pin'

## Glboal variables
pep_file_ext = '.pep.txt'
psm_file_ext = '.psm.txt'
decoy_prefix = 'Rev_'
divide = sipros_post_module.divide
FDR_parameter = 1.0
output_result = "./result"


# Clean temp folder
def initialize():
    print("Processing...")
    if os.path.exists(output_result) is True:
        shutil.rmtree(output_result)
    if os.path.exists(output_result) is False:
        os.makedirs(output_result)
    if os.path.exists(path_to_temp_files) is True:
        shutil.rmtree(path_to_temp_files)
    if os.path.exists("./output") is True:
        shutil.rmtree("./output")
    if os.path.exists(path_to_temp_files) is False:
        os.makedirs(path_to_temp_files)
    if os.path.exists(path_to_species_pin) is False:
        os.makedirs(path_to_species_pin)
    if os.path.exists("./output") is False:
        os.makedirs("./output")
    for infile in glob.glob(os.path.join("./", '*.psm.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pep.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pro.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pro2psm.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pro2pep.txt')):
        os.remove(infile)
    return


def clean():
    print("Processing...")
    if os.path.exists("./output") is True:
        shutil.rmtree("./output")
    if os.path.exists(path_to_temp_files) is True:
        shutil.rmtree(path_to_temp_files)
    if os.path.exists('./tmp') is True:
        shutil.rmtree('./tmp')
    for infile in glob.glob(os.path.join("./", '*.psm.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pep.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pro.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pro2psm.txt')):
        os.remove(infile)
    for infile in glob.glob(os.path.join("./", '*.pro2pep.txt')):
        os.remove(infile)
    return


class PSM:
    def __init__(self, filename, file, scan, ParentCharge, rank, MeasuredParentMass, CalculatedParentMass, Massdiff,
                 rescore, PTM_score, IdentifiedPeptide, PSM_Label,
                 Proteins, Proteinname, ProteinCount):
        self.filename = filename
        self.file = file
        self.scan = scan
        self.ParentCharge = ParentCharge
        self.rank = rank
        self.MeasuredParentMass = MeasuredParentMass
        self.CalculatedParentMass = CalculatedParentMass
        self.Massdiff = Massdiff
        self.MassErrorPPM = 'NA'
        self.ScanType = 'NA'
        self.SearchName = 'NA'
        self.ScoringFunction = 'softmax'
        self.rescore = rescore
        self.DeltaZ = 'NA'
        self.DeltaP = 'NA'
        self.PTM_score = PTM_score
        self.IdentifiedPeptide = IdentifiedPeptide
        self.OriginalPeptide = 'NA'
        self.PSM_Label = PSM_Label
        self.Proteins = Proteins
        self.Proteinname = Proteinname
        self.ProteinCount = ProteinCount


class Peptide:
    def __init__(self):
        self.IdentifiedPeptide = ''
        self.ParentCharge = ''
        self.OriginalPeptide = ''
        self.ProteinNames = []
        self.ProteinCount = 0
        self.SpectralCount = 0
        self.BestScore = 0.0
        self.PSMs = []

    def add(self, oPsm):
        self.SpectralCount += 1
        if self.BestScore < oPsm.rescore:
            self.BestScore = oPsm.rescore
        self.PSMs.append('{0}_{1}_{2}_{3}'.format(oPsm.file, oPsm.scan, oPsm.ParentCharge, oPsm.rank))
        self.ScanType = 'NA'
        self.SearchName = 'NA'
        if oPsm.PSM_Label == True:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'

    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide
        self.ProteinNames = oPsm.Proteinname
        self.ProteinCount = oPsm.ProteinCount
        self.SpectralCount = 1
        self.BestScore = oPsm.rescore
        self.PSMs.append('{0}_{1}_{2}_{3}'.format(oPsm.file, oPsm.scan, oPsm.ParentCharge, oPsm.rank))
        self.ScanType = 'NA'
        self.SearchName = 'NA'
        if oPsm.PSM_Label == True:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'


# FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(TP)
    # FDR_accept = True
    if FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True
    return (FDR_accept, float(FDR_value))


def show_Fdr(psm_list, fdr_float, charge_left_given=-1, charge_right_given=-1):
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda psm: (psm.rescore, psm.Massdiff, psm.PTM_score), reverse=True)
    decoy = 0
    target = 0
    best_nums = [0, 0]

    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.PSM_Label:
            target += 1
        elif not oPsm.PSM_Label:
            decoy += 1
        else:
            sys.stderr.write('error 768.\n')
        (FDR_accept, FDR_value) = FDR_calculator(decoy, target)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (best_nums[0] + best_nums[1]) < (decoy + target):
                best_nums = [decoy, target]
                cutoff_probability = oPsm.rescore

    for oPsm in list_sorted:
        if charge_left_given != -1 and (
                oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.rescore >= cutoff_probability:
            psm_filtered_list.append(oPsm)

    return psm_filtered_list


# peptide level filtering
def show_Fdr_Pep(psm_list, fdr_float):
    list_sorted = sorted(psm_list, key=lambda psm: (psm.rescore, psm.Massdiff, psm.PTM_score), reverse=True)

    peptide_set = set()
    decoy = 0
    target = 0
    best_nums = [0, 0]

    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.PSM_Label:
                target += 1
                peptide_set.add(pep_str)
            elif not oPsm.PSM_Label:
                decoy += 1
                peptide_set.add(pep_str)
            else:
                sys.stderr.write('error 768.\n')

        (FDR_accept, FDR_value) = FDR_calculator(decoy, target)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (best_nums[0] + best_nums[1]) < (decoy + target):
                best_nums = [decoy, target]
                cutoff_probability = oPsm.rescore

    peptide = dict()
    for oPsm in list_sorted:
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if oPsm.rescore >= cutoff_probability:
            if pep_str in peptide:
                peptide[pep_str].add(oPsm)
            else:
                oPeptide = Peptide()
                oPeptide.set(oPsm)
                peptide[pep_str] = oPeptide
    return peptide


# remove redundant psm, only one unique spectrum kept
def re_rank(psm_list, consider_charge_bool=False):
    psm_new_list = []
    psm_dict = {}
    if consider_charge_bool:
        for oPsm in psm_list:
            sId = '{0}_{1}_{2}'.format(str(oPsm.file), str(oPsm.scan), str(oPsm.ParentCharge))
            if sId in psm_dict:
                if oPsm.rescore > psm_dict[sId].rescore:
                    psm_dict[sId] = oPsm
                elif oPsm.rescore == psm_dict[sId].rescore:
                    if abs(oPsm.Massdiff) < abs(psm_dict[sId].Massdiff):
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.Massdiff) == abs(psm_dict[sId].Massdiff):
                        # calculate PTM scores
                        if oPsm.PTM_score < psm_dict[sId].PTM_score:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTM_score == psm_dict[sId].PTM_score:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
                            elif oPsm.IdentifiedPeptide.upper() == psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId].add_protein(oPsm.protein_list)

            else:
                psm_dict[sId] = oPsm
    else:
        for oPsm in psm_list:
            sId = '{0}_{1}'.format(str(oPsm.file), str(oPsm.scan))
            if sId in psm_dict:
                if oPsm.rescore > psm_dict[sId].rescore:
                    psm_dict[sId] = oPsm
                elif oPsm.rescore == psm_dict[sId].rescore:
                    if abs(oPsm.Massdiff) < abs(psm_dict[sId].Massdiff):
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.Massdiff) == abs(psm_dict[sId].Massdiff):
                        # calculate PTM scores
                        if oPsm.PTM_score < psm_dict[sId].PTM_score:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTM_score == psm_dict[sId].PTM_score:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
            else:
                psm_dict[sId] = oPsm
    for _key, value in psm_dict.items():
        psm_new_list.append(value)
    return psm_new_list


def read_pin_header(path_to_pin_input):
    with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
        header = pin_input.readline()
        pin_input.close()
    return header


def read_txt(input_txt, method):
    PSMs_ouput = []
    SpecId_rule = r"[a-zA-Z0-9|/.]+"
    SpecId_loc = 0
    scan_loc = 1
    rank_loc = 2
    charge_loc = 3
    MeasuredParentMass_loc = 4
    CalculatedParentMass_loc = 5
    Proteins_loc = 16
    IdentifyPeptide_loc = 13
    if method == "CP":
        score_loc = 21  # P-score
        score_measure = 1
    if method == "C":
        score_loc = 6  # Xcorr
        score_measure = -1
    with open(input_txt, 'r', encoding='utf-8') as txt_input:
        next(txt_input)
        for line in txt_input:
            data = line.split("\t")
            ProteinNames = data[Proteins_loc]
            Proteins = ProteinNames.split(",")
            SpecId = data[SpecId_loc]
            SpecId_info = re.findall(SpecId_rule, SpecId)
            charge = data[charge_loc]
            rank = data[rank_loc]
            file_id = SpecId_info[0]
            for index in range(1, len(SpecId_info) - 3):
                file_id = file_id + '_' + SpecId_info[index]
            filename = file_id
            scan = data[scan_loc]
            score = data[score_loc]
            ProteinCount = 0
            IdentifyPeptide = data[IdentifyPeptide_loc]
            PTM_score = IdentifyPeptide.count('[15.9949]')
            MeasuredParentMass = data[MeasuredParentMass_loc]
            CalculatedParentMass = data[CalculatedParentMass_loc]
            MassDiff = abs(float(CalculatedParentMass) - float(MeasuredParentMass))
            for p in Proteins:
                current = p
                if current.find(decoy_prefix) == -1:
                    PSM_Label = True
                    ProteinCount += 1
                else:
                    PSM_Label = False
            PSMs_ouput.append(
                PSM(filename, file_id, int(scan), int(charge), int(rank), MeasuredParentMass, CalculatedParentMass,
                    MassDiff, float(score) * score_measure,
                    PTM_score, IdentifyPeptide, PSM_Label, Proteins, ProteinNames, ProteinCount))
    return PSMs_ouput


def fdr_control_original(rank_list, input_pin, psm_fdr, pep_fdr, output_dir):
    psm_num_sp = 0
    pep_num_sp = 0
    decoy_pep = 0
    pep_counter = 0
    prefix = input_pin.split('/')[-1]
    print(prefix)
    filter_list = show_Fdr(rank_list, psm_fdr)
    with open(output_dir + '/' + prefix + '.psm.txt', 'w') as f:
        psm_out_list = ['Filename',  # 0
                        'ScanNumber',  # 1
                        'ParentCharge',  # 2
                        'MeasuredParentMass',  # 3
                        'CalculatedParentMass',  # 4
                        'MassErrorDa',  # 5 CalculatedParentMass - MeasuredParentMass
                        'MassErrorPPM',  # 6 MassErrorDa / CalculatedParentMass
                        'ScanType',  # 7
                        'SearchName',  # 8
                        'ScoringFunction',  # 9
                        'Score',  # 10
                        'DeltaZ',  # 11 the difference score between the rank 1 and 2
                        'DeltaP',  # 12
                        'IdentifiedPeptide',  # 13
                        'OriginalPeptide',  # 14
                        'ProteinNames',  # 15
                        'ProteinCount',  # 16
                        'TargetMatch']  # 17
        f.write('\t'.join(psm_out_list) + '\n')
        for psm in filter_list:
            TargetMatch = 'F'
            if psm.PSM_Label:
                TargetMatch = 'T'
                psm_num_sp += 1
            f.write(str(psm.filename) + '\t' + str(psm.scan) + '\t' + str(psm.ParentCharge) + '\t' + str(
                psm.MeasuredParentMass) + '\t' + str(psm.CalculatedParentMass) + '\t' + str(psm.Massdiff) + '\t' + str(
                psm.MassErrorPPM) + '\t' + str(psm.ScanType) + '\t' + str(psm.SearchName) + '\t' + str(
                psm.ScoringFunction) + '\t' + str(-psm.rescore) + '\t' + str(psm.DeltaZ) + '\t' + str(
                psm.DeltaP) + '\t' + str(psm.IdentifiedPeptide) + '\t' + str(psm.OriginalPeptide) + '\t' + str(
                psm.Proteinname) + '\t' + str(psm.ProteinCount) + '\t' + TargetMatch + '\n')
    print('psm: ', psm_num_sp)
    filter_pep_list = show_Fdr_Pep(rank_list, pep_fdr)
    with open(output_dir + '/' + prefix + '.pep.txt', 'w') as f:
        pep_out_list = ['IdentifiedPeptide',  # 0
                        'ParentCharge',  # 1
                        'OriginalPeptide',  # 2
                        'ProteinNames',  # 3
                        'ProteinCount',  # 4
                        'TargetMatch',  # 5
                        'SpectralCount',  # 6 number of PSMs matched to this peptide
                        'BestScore',  # 7 the highest score of those PSMs
                        'PSMs',
                        # 8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                        'ScanType',  # 9
                        'SearchName']  # 10
        f.write('\t'.join(pep_out_list) + '\n')
        for key, pep in filter_pep_list.items():
            if pep.TargetMatch == 'T':
                pep_num_sp += 1
                pep_counter += 1
            else:
                decoy_pep += 1
            f.write(pep.IdentifiedPeptide + '\t' + str(
                pep.ParentCharge) + '\t' + pep.OriginalPeptide + '\t' + '{' + pep.ProteinNames + '}' + '\t' + str(
                pep.ProteinCount) + '\t' + pep.TargetMatch + '\t' + str(pep.SpectralCount) + '\t' + str(
                pep.BestScore) + '\t' + ','.join(pep.PSMs) + '\t' + pep.ScanType + '\t' + pep.SearchName + '\n')
    print('pep:', pep_counter)
    return


def create_new_tab_file(path_to_tab_input, path_to_tab_output):
    header = read_pin_header(path_to_tab_input)
    with open(path_to_tab_output, "wt") as tab_output:
        tab_output.write(header)
        tab_output.close()
    return


def fdr_control_by_species(input_pins_dir, psm_fdr, pep_fdr, output_dir, method):
    target_list = []
    decoy_list = []
    IdentifiedPeptide_id_list = []
    IdentifiedPeptideid_dic = {}
    IdentifiedPeptideid_score_dic = {}
    psm_counter = 0
    pep_counter = 0
    input_dir = os.listdir(input_pins_dir)
    with open(output_dir + '/output_species.psm.txt', 'w') as f:
        psm_out_list = ['Filename',  # 0
                        'ScanNumber',  # 1
                        'ParentCharge',  # 2
                        'MeasuredParentMass',  # 3
                        'CalculatedParentMass',  # 4
                        'MassErrorDa',  # 5 CalculatedParentMass - MeasuredParentMass
                        'MassErrorPPM',  # 6 MassErrorDa / CalculatedParentMass
                        'ScanType',  # 7
                        'SearchName',  # 8
                        'ScoringFunction',  # 9
                        'Score',  # 10
                        'DeltaZ',  # 11 the difference score between the rank 1 and 2
                        'DeltaP',  # 12
                        'IdentifiedPeptide',  # 13
                        'OriginalPeptide',  # 14
                        'ProteinNames',  # 15
                        'ProteinCount',  # 16
                        'TargetMatch']  # 17
        f.write('\t'.join(psm_out_list) + '\n')
    with open(output_dir + '/output_species.pep.txt', 'w') as f:
        pep_out_list = ['IdentifiedPeptide',  # 0
                        'ParentCharge',  # 1
                        'OriginalPeptide',  # 2
                        'ProteinNames',  # 3
                        'ProteinCount',  # 4
                        'TargetMatch',  # 5
                        'SpectralCount',  # 6 number of PSMs matched to this peptide
                        'BestScore',  # 7 the highest score of those PSMs
                        'PSMs',
                        # 8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                        'ScanType',  # 9
                        'SearchName']  # 10
        f.write('\t'.join(pep_out_list) + '\n')
    for specie in input_dir:
        psm_num_sp = 0
        pep_num_sp = 0
        PSMs = read_txt(input_pins_dir + '/' + specie, method)
        print(specie)
        psm_list = sorted(PSMs, key=lambda psm: (psm.rescore, psm.Massdiff, psm.PTM_score), reverse=True)
        rank_list = re_rank(PSMs)
        filter_list = show_Fdr(rank_list, psm_fdr)
        with open(output_dir + '/output_species.psm.txt', "a") as f1:
            for psm in filter_list:
                specie_id = psm.file + "_" + str(psm.scan) + "_" + str(psm.ParentCharge) + "_" + str(psm.rank)
                TargetMatch = 'F'
                if psm.PSM_Label == True:
                    TargetMatch = 'T'
                    psm_num_sp += 1
                    if specie_id not in target_list:
                        psm_counter += 1
                        target_list.append(specie_id)
                        f1.write(str(psm.filename) + '\t' + str(psm.scan) + '\t' + str(psm.ParentCharge) + '\t' + str(
                            psm.MeasuredParentMass) + '\t' + str(psm.CalculatedParentMass) + '\t' + str(
                            psm.Massdiff) + '\t' + str(
                            psm.MassErrorPPM) + '\t' + str(psm.ScanType) + '\t' + str(psm.SearchName) + '\t' + str(
                            psm.ScoringFunction) + '\t' + str(-psm.rescore) + '\t' + str(psm.DeltaZ) + '\t' + str(
                            psm.DeltaP) + '\t' + str(psm.IdentifiedPeptide) + '\t' + str(
                            psm.OriginalPeptide) + '\t' + str(
                            psm.Proteinname) + '\t' + str(psm.ProteinCount) + '\t' + TargetMatch + '\n')
                else:
                    if specie_id not in decoy_list:
                        decoy_list.append(specie_id)
                        f1.write(str(psm.filename) + '\t' + str(psm.scan) + '\t' + str(psm.ParentCharge) + '\t' + str(
                            psm.MeasuredParentMass) + '\t' + str(psm.CalculatedParentMass) + '\t' + str(
                            psm.Massdiff) + '\t' + str(
                            psm.MassErrorPPM) + '\t' + str(psm.ScanType) + '\t' + str(psm.SearchName) + '\t' + str(
                            psm.ScoringFunction) + '\t' + str(-psm.rescore) + '\t' + str(psm.DeltaZ) + '\t' + str(
                            psm.DeltaP) + '\t' + str(psm.IdentifiedPeptide) + '\t' + str(
                            psm.OriginalPeptide) + '\t' + str(
                            psm.Proteinname) + '\t' + str(psm.ProteinCount) + '\t' + TargetMatch + '\n')
        print('Local psm:', psm_num_sp)
        filter_pep_list = show_Fdr_Pep(rank_list, pep_fdr)
        for pep in filter_pep_list.values():
            IdentifiedPep = pep.IdentifiedPeptide
            IdentifiedPep_id = IdentifiedPep + "_" + str(pep.ParentCharge)
            if IdentifiedPep_id not in IdentifiedPeptide_id_list:
                IdentifiedPeptide_id_list.append(IdentifiedPep_id)
                IdentifiedPeptideid_dic[IdentifiedPep_id] = pep
                IdentifiedPeptideid_score_dic[IdentifiedPep_id] = pep.BestScore
                if pep.TargetMatch == 'T':
                    pep_num_sp += 1
            else:
                in_list_pep = IdentifiedPeptideid_dic[IdentifiedPep_id]
                line = pep.ProteinNames + ',' + in_list_pep.ProteinNames
                pep_data = line.split(",")
                pep_data = list(set(pep_data))
                in_list_pep.ProteinCount = len(pep_data)
                in_list_pep.ProteinNames = ','.join(pep_data)
                if pep.TargetMatch == 'T' and in_list_pep.TargetMatch == 'F':
                    in_list_pep.TargetMatch == 'T'
                    pep_num_sp += 1
                if pep.BestScore > in_list_pep.BestScore:
                    in_list_pep.BestScore = pep.BestScore
                in_list_pep.PSMs.extend(pep.PSMs)
                in_list_pep.PSMs = list(set(in_list_pep.PSMs))
                in_list_pep.SpectralCount = len(in_list_pep.PSMs)
                IdentifiedPeptideid_score_dic[IdentifiedPep_id] = in_list_pep.BestScore
        print('Local pep:', pep_num_sp)
    with open(output_dir + '/output_species.pep.txt', 'a') as f2:
        IdentifiedPeptideid_score_dic_sorted = sorted(IdentifiedPeptideid_score_dic.items(), key=lambda item: item[1],
                                                      reverse=1)
        for id in IdentifiedPeptideid_score_dic_sorted:
            pep = IdentifiedPeptideid_dic[id[0]]
            if pep.TargetMatch == 'T':
                pep_counter += 1
            f2.write(pep.IdentifiedPeptide + '\t' + str(
                pep.ParentCharge) + '\t' + pep.OriginalPeptide + '\t' + '{' + pep.ProteinNames + '}' + '\t' + str(
                pep.ProteinCount) + '\t' + pep.TargetMatch + '\t' + str(pep.SpectralCount) + '\t' + str(
                pep.BestScore) + '\t' + ','.join(pep.PSMs) + '\t' + pep.ScanType + '\t' + pep.SearchName + '\n')
    print(input_pins_dir)
    print('Psm: ', psm_counter)
    print('Pep: ', pep_counter)
    return psm_counter, pep_counter


# original fdr control method with binary search
def binary_protein_fdr_search_original(input_pin, psm_fdr, target_protein_fdr, max_run, output_dir, is_report, method):
    pep_fdr_list = []
    run_id_list = []
    protein_fdr_list = []
    res_dic = {}
    run_id = 1
    run_id_list.append(run_id)
    print("run_id:", run_id)
    current_pep_fdr = target_protein_fdr
    print("current_pep_fdr:", current_pep_fdr)
    pep_fdr_list.append(current_pep_fdr)
    PSMs = read_txt(input_pin, method)
    psm_list = sorted(PSMs, key=lambda psm: (psm.rescore, psm.Massdiff, psm.PTM_score), reverse=True)
    rank_list = re_rank(PSMs)
    fdr_control_original(rank_list, input_pin, psm_fdr, current_pep_fdr, output_dir)
    current_protein_fdr, target_proteins_after_filtering = peptides_assembling_lite(output_dir)
    error = current_protein_fdr - target_protein_fdr
    res_dic[current_pep_fdr] = abs(error)
    protein_fdr_list.append(current_protein_fdr)
    past_protein_fdr = current_protein_fdr
    print("Testing the boundary.......")
    if error < 0:
        left_boundary = current_pep_fdr
        while error < 0:
            current_pep_fdr += 0.003
            fdr_control_original(rank_list, input_pin, psm_fdr, current_pep_fdr, output_dir)
            current_protein_fdr, target_proteins_after_filtering = peptides_assembling_lite(output_dir)
            error = current_protein_fdr - target_protein_fdr
        right_boundary = current_pep_fdr
    else:
        right_boundary = current_pep_fdr
        while error > 0:
            current_pep_fdr -= 0.003
            fdr_control_original(rank_list, input_pin, psm_fdr, current_pep_fdr, output_dir)
            current_protein_fdr, target_proteins_after_filtering = peptides_assembling_lite(output_dir)
            error = current_protein_fdr - target_protein_fdr
        left_boundary = current_pep_fdr
    print("Testing Done!")
    print("left_boundary:", left_boundary)
    print("right_boundary:", right_boundary)
    # past_protein_fdr != current_protein_fdr
    # past_protein_fdr not in list(res_dic.values())
    # tolerance = 1 / (target_proteins_after_filtering + 1)
    # if tolerance < 0.0005:
    tolerance = 0.0001
    while abs(error) > tolerance:
        # print(list(res_dic.values()))
        # print(abs(error))
        run_id += 1
        run_id_list.append(run_id)
        print("run_id:", run_id)
        current_pep_fdr = (right_boundary + left_boundary) / 2
        pep_fdr_list.append(current_pep_fdr)
        fdr_control_original(rank_list, input_pin, psm_fdr, current_pep_fdr, output_dir)
        current_protein_fdr, target_proteins_after_filtering = peptides_assembling_lite(output_dir)
        error = current_protein_fdr - target_protein_fdr
        res_dic[current_pep_fdr] = abs(error)
        protein_fdr_list.append(current_protein_fdr)
        past_protein_fdr = current_protein_fdr
        if error < 0:
            left_boundary = current_pep_fdr
        else:
            right_boundary = current_pep_fdr
        if run_id > max_run:
            break
    best_pep_fdr = min(res_dic, key=res_dic.get)
    current_pep_fdr = best_pep_fdr
    fdr_control_original(rank_list, input_pin, psm_fdr, current_pep_fdr, output_dir)
    if is_report:
        current_protein_fdr, target_proteins_after_filtering = peptides_assembling(output_dir)
    else:
        current_protein_fdr, target_proteins_after_filtering = peptides_assembling_lite(output_dir)
    print("best_pep_fdr: ", best_pep_fdr)
    # visualization
    '''
    plt.figure(1)
    plt.scatter(run_id_list, protein_fdr_list)
    plt.xlabel('Run ID')
    plt.ylabel('protein_fdr')
    plt.axhline(y=0.01, color='r', linestyle='-')
    plt.title(input_pin)
    plt.show()
    plt.figure(2)
    plt.xlabel('pep_fdr')
    plt.ylabel('protein_fdr')
    plt.scatter(pep_fdr_list, protein_fdr_list)
    plt.axhline(y=0.01, color='r', linestyle='-')
    plt.title(input_pin)
    plt.show()
    '''
    return current_protein_fdr


# Combine the pep.txt files from different species
def combine_pep_txt(input_dir):
    IdentifiedPeptide_id_list = []
    IdentifiedPeptideid_dic = {}
    IdentifiedPeptideid_score_dic = {}
    pep_counter = 0
    with open('./output_species.pep.txt', 'w') as f:
        pep_out_list = ['IdentifiedPeptide',  # 0
                        'ParentCharge',  # 1
                        'OriginalPeptide',  # 2
                        'ProteinNames',  # 3
                        'ProteinCount',  # 4
                        'TargetMatch',  # 5
                        'SpectralCount',  # 6 number of PSMs matched to this peptide
                        'BestScore',  # 7 the highest score of those PSMs
                        'PSMs',
                        # 8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                        'ScanType',  # 9
                        'SearchName']  # 10
        f.write('\t'.join(pep_out_list) + '\n')
    for sp_pep_txt in os.listdir(input_dir):
        pep_num_sp = 0
        current_file = input_dir + "/" + sp_pep_txt
        with open(current_file, 'r') as txt_input:
            next(txt_input)
            for line in txt_input:
                data = line.replace("{", "").replace("}", "").split("\t")
                IdentifiedPep = data[0]
                IdentifiedPep_id = IdentifiedPep + "_" + str(data[1])
                if IdentifiedPep_id not in IdentifiedPeptide_id_list:
                    IdentifiedPeptide_id_list.append(IdentifiedPep_id)
                    IdentifiedPeptideid_dic[IdentifiedPep_id] = data
                    IdentifiedPeptideid_score_dic[IdentifiedPep_id] = data[7]
                    if data[5] == 'T':
                        pep_num_sp += 1
                else:
                    in_list_pep = IdentifiedPeptideid_dic[IdentifiedPep_id]
                    pr = data[3] + ',' + in_list_pep[3]
                    pep_data = pr.split(",")
                    pep_data = list(set(pep_data))
                    in_list_pep[4] = len(pep_data)
                    in_list_pep[3] = ','.join(pep_data)
                    if data[5] == 'T' and in_list_pep[5] == 'F':
                        in_list_pep[5] == 'T'
                        pep_num_sp += 1
                    if data[7] > in_list_pep[7]:
                        in_list_pep[7] = data[7]
                    PSMs = data[8]
                    PSMs_in_list = in_list_pep[8]
                    PSMs_in_list = PSMs + "," + PSMs_in_list
                    PSMs_in_list_list = list(set(PSMs_in_list.split(",")))
                    in_list_pep[6] = len(PSMs_in_list_list)
                    in_list_pep[8] = str(','.join(PSMs_in_list_list))
                    IdentifiedPeptideid_score_dic[IdentifiedPep_id] = in_list_pep[7]
            print('Local pep:', pep_num_sp)
    with open('./output_species.pep.txt', 'a') as f2:
        IdentifiedPeptideid_score_dic_sorted = sorted(IdentifiedPeptideid_score_dic.items(),
                                                      key=lambda item: item[1],
                                                      reverse=1)
        for id in IdentifiedPeptideid_score_dic_sorted:
            pep = IdentifiedPeptideid_dic[id[0]]
            if pep[5] == 'T':
                pep_counter += 1
            f2.write(pep[0] + '\t' + str(
                pep[1]) + '\t' + pep[2] + '\t' + '{' + pep[3] + '}' + '\t' + str(
                pep[4]) + '\t' + pep[5] + '\t' + str(pep[6]) + '\t' + str(
                pep[7]) + '\t' + pep[8] + '\t' + pep[9] + '\t' + pep[0] + '\n')
    return pep_counter


# Run original fdr control on different species
def select_best_pep_in_sp(method, input_pins_dir, psm_fdr, protein_fdr_sp, max_run, sp):
    current_sp_input_file = input_pins_dir + "/" + sp
    current_sp_output_dir = "./tmp/" + sp
    if os.path.exists(current_sp_output_dir) is False:
        os.makedirs(current_sp_output_dir)
    binary_protein_fdr_search_original(current_sp_input_file, psm_fdr, protein_fdr_sp, max_run, current_sp_output_dir,
                                       False, method)
    pep_txt_file = current_sp_output_dir + "/" + sp + ".pep.txt"
    pep_txt_output = "./output/" + sp + ".pep.txt"
    print(pep_txt_output)
    shutil.copyfile(pep_txt_file, pep_txt_output)


# Species-based fdr control method
def binary_protein_fdr_search_species(input_pins_dir, psm_fdr, pep_fdr_sp, max_run, output_dir, is_report, method):
    '''
    file_dic = {}
    for sp in input_list:
        current_sp_input_file = input_pins_dir + "/" + sp
        file_size = os.path.getsize(current_sp_input_file)
        file_dic[sp] = file_size
    input_list = [i[0] for i in sorted(file_dic.items(), key=lambda kv: (kv[1], kv[0], reverse=True)]
    '''
    input_list = os.listdir(input_pins_dir)
    # print(input_list)
    # print(type(input_list))
    input_list = reversed(input_list)
    print(input_list)
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    pool.map(partial(select_best_pep_in_sp, method, input_pins_dir, psm_fdr, pep_fdr_sp, max_run), input_list)
    pool.close()
    pool.join()
    combine_pep_txt("./output")
    if is_report:
        current_protein_fdr, target_proteins_after_filtering = peptides_assembling(output_dir)
    else:
        current_protein_fdr, target_proteins_after_filtering = peptides_assembling_lite(output_dir)
    print("target_proteins_after_filtering:", target_proteins_after_filtering)
    '''
    input_list = os.listdir(input_pins_dir)
    for sp in input_list:
        select_best_pep_in_sp(input_pins_dir, psm_fdr, pep_fdr_sp, max_run, sp)
    combine_pep_txt("./output")
    current_protein_fdr, target_proteins_after_filtering = peptides_assembling_lite(output_dir)
    print("target_proteins_after_filtering:", target_proteins_after_filtering)
    '''
    return current_protein_fdr, target_proteins_after_filtering


#  binary search for Species-based fdr control method
def binary_protein_fdr_search_species_search(input_pins_dir, psm_fdr, pep_fdr, target_protein_fdr, max_run, output_dir,
                                             method):
    pep_fdr_list = []
    run_id_list = []
    protein_fdr_list = []
    res_dic = {}
    file_dic = {}
    input_list = os.listdir(input_pins_dir)
    for sp in input_list:
        current_sp_input_file = input_pins_dir + "/" + sp
        file_size = os.path.getsize(current_sp_input_file)
        file_dic[sp] = file_size

    psm_counter, pep_counter = fdr_control_by_species(input_pins_dir, psm_fdr, pep_fdr, output_dir, method)
    shutil.copyfile("./output_species.psm.txt", output_result + "/output_species.psm.txt")
    shutil.copyfile("./output_species.pep.txt", output_result + "/output_species.pep.txt")
    print('Psm: ', psm_counter)
    print('Pep: ', pep_counter)
    run_id = 1
    run_id_list.append(run_id)
    print("run_id:", run_id)
    current_protein_fdr_sp = target_protein_fdr
    print("current_protein_fdr_sp:", current_protein_fdr_sp)
    pep_fdr_list.append(current_protein_fdr_sp)
    current_protein_fdr, target_proteins_after_filtering = binary_protein_fdr_search_species(input_pins_dir, psm_fdr,
                                                                                             current_protein_fdr_sp,
                                                                                             max_run,
                                                                                             output_dir, False, method)
    # tolerance = 1 / (target_proteins_after_filtering + 1)
    # if tolerance < 0.0005:
    # tolerance = 0.0005
    error = current_protein_fdr - target_protein_fdr
    res_dic[current_protein_fdr_sp] = abs(error)
    protein_fdr_list.append(current_protein_fdr)
    past_protein_fdr = current_protein_fdr
    print("Testing the boundary.......")
    if error < 0:
        left_boundary = current_protein_fdr_sp
        while error < 0:
            current_protein_fdr_sp += 0.003
            current_protein_fdr, target_proteins_after_filtering = binary_protein_fdr_search_species(input_pins_dir,
                                                                                                     psm_fdr,
                                                                                                     current_protein_fdr_sp,
                                                                                                     max_run,
                                                                                                     output_dir, False,
                                                                                                     method)
            error = current_protein_fdr - target_protein_fdr
        right_boundary = current_protein_fdr_sp
    else:
        right_boundary = current_protein_fdr_sp
        while error > 0:
            current_protein_fdr_sp -= 0.003
            current_protein_fdr, target_proteins_after_filtering = binary_protein_fdr_search_species(input_pins_dir,
                                                                                                     psm_fdr,
                                                                                                     current_protein_fdr_sp,
                                                                                                     max_run,
                                                                                                     output_dir, False,
                                                                                                     method)
            error = current_protein_fdr - target_protein_fdr
        left_boundary = current_protein_fdr_sp
    print("Testing Done!")
    print("left_boundary:", left_boundary)
    print("right_boundary:", right_boundary)
    # past_protein_fdr != current_protein_fdr
    # while past_protein_fdr not in list(res_dic.values()) or abs(error) > tolerance:
    tolerance = 0.0001
    while abs(error) > tolerance:
        print(past_protein_fdr not in list(res_dic.values()))
        # print(abs(error))
        run_id += 1
        run_id_list.append(run_id)
        print("run_id:", run_id)
        current_protein_fdr_sp = (right_boundary + left_boundary) / 2
        pep_fdr_list.append(current_protein_fdr_sp)
        current_protein_fdr, target_proteins_after_filtering = binary_protein_fdr_search_species(input_pins_dir,
                                                                                                 psm_fdr,
                                                                                                 current_protein_fdr_sp,
                                                                                                 max_run, output_dir,
                                                                                                 False, method)
        error = current_protein_fdr - target_protein_fdr
        res_dic[current_protein_fdr_sp] = abs(error)
        protein_fdr_list.append(current_protein_fdr)
        past_protein_fdr = current_protein_fdr
        if error < 0:
            left_boundary = current_protein_fdr_sp
        else:
            right_boundary = current_protein_fdr_sp
        if run_id > max_run:
            break
    best_pep_fdr = min(res_dic, key=res_dic.get)
    current_protein_fdr_sp = best_pep_fdr
    current_protein_fdr, target_proteins_after_filtering = binary_protein_fdr_search_species(input_pins_dir, psm_fdr,
                                                                                             current_protein_fdr_sp,
                                                                                             max_run,
                                                                                             output_dir, True, method)
    # psm_counter, pep_counter = fdr_control_by_species(input_pins_dir, psm_fdr, pep_fdr, output_dir, method)
    shutil.copyfile("./output_species.pro.txt", output_result + "/output_species.pro.txt")
    shutil.copyfile("./output_species.pro2pep.txt", output_result + "/output_species.pro2pep.txt")
    shutil.copyfile("./output_species.pro2psm.txt", output_result + "/output_species.pro2psm.txt")
    print('Psm: ', psm_counter)
    print('Pep: ', pep_counter)
    print("best_pep_fdr: ", best_pep_fdr)
    print("current_protein_fdr: ", current_protein_fdr)
    print("target_proteins_after_filtering: ", target_proteins_after_filtering)
    '''
    plt.figure(1)
    plt.scatter(run_id_list, protein_fdr_list)
    plt.xlabel('Run ID')
    plt.ylabel('protein_fdr')
    plt.axhline(y=0.01, color='r', linestyle='-')
    plt.title(input_pins_dir)
    plt.show()
    plt.figure(2)
    plt.xlabel('pep_fdr')
    plt.ylabel('protein_fdr')
    plt.scatter(pep_fdr_list, protein_fdr_list)
    plt.axhline(y=0.01, color='r', linestyle='-')
    plt.title(input_pins_dir)
    plt.show()
    '''
    return current_protein_fdr


def main(argv):
    print('Configuration File:', argv[1])
    config_dic = {}
    with open(argv[1], 'r') as f:
        for line in f:
            if "#" not in line:
                data = line.split("=")
                config_dic[data[0]] = data[-1]
    # fdr control parameters
    input_txt_dir = config_dic["comet_txt_dir"].replace(" ", "").replace("\n", "")
    global decoy_prefix
    decoy_prefix = config_dic["decoy_prefix"].replace(" ", "").replace("\n", "")
    global output_result
    output_result = config_dic["output_dir"].replace(" ", "").replace("\n", "")
    psm_fdr = float(config_dic["psm_fdr"].replace(" ", "").replace("\n", ""))
    pep_fdr = float(config_dic["pep_fdr"].replace(" ", "").replace("\n", ""))
    target_protein_fdr = float(config_dic["protein_fdr"].replace(" ", "").replace("\n", ""))
    method = config_dic["method"].replace(" ", "").replace("\n", "")
    print("Method: ", method)
    path_to_species_dic = config_dic["sp_dic"].replace(" ", "").replace("\n", "")
    target_file = config_dic["percolator_target"].replace(" ", "").replace("\n", "")
    deocy_file = config_dic["percolator_decoy"].replace(" ", "").replace("\n", "")
    # other parameters
    max_run = 10
    combined_txt_dir = path_to_temp_files + "/"
    postfix = ".txt"
    output_dir_or = './'
    input_pins_dir = './temp/species_pin'
    path_to_converted_output = "./temp/"
    # path_to_cleaned_comet_output = "./temp/cleaned.txt"
    # input_pin_dir = "./hgut_pin/"
    # combined_txt = "./combined.txt"
    # input_pin = 'converted.txt'
    # Main
    initialize()
    # binary_protein_fdr_search_original(input_pin, psm_fdr, target_protein_fdr, max_run, output_dir_or, True, method)
    if method == "C":
        txt_input = combine_txt_files(input_txt_dir, postfix, combined_txt_dir)
        path_to_cleaned_comet_output = clean_comet(txt_input)
        # divide_species_auto_c(path_to_cleaned_comet_output, path_to_species_pin) # U1 Mock
        divide_species_dic_c(path_to_cleaned_comet_output, path_to_species_dic,
                             path_to_species_pin)  # Marine, soil, hgut
        binary_protein_fdr_search_species_search(input_pins_dir, psm_fdr, pep_fdr, target_protein_fdr, max_run,
                                                 output_dir_or, method)
    if method == "CP":
        txt_input = combine_txt_files(input_txt_dir, postfix, combined_txt_dir)
        # combine_pin_files(input_pin_dir, ".pin", combined_txt_dir)
        percolator_to_comet(target_file, deocy_file, txt_input , postfix, path_to_converted_output)
        # divide_species_auto_p("./converted.txt", path_to_species_pin) # U1 Mock
        divide_species_dic_p("./temp/converted.txt", path_to_species_dic, path_to_species_pin) # Marine, soil, hgut
        binary_protein_fdr_search_species_search(input_pins_dir, psm_fdr, psm_fdr, target_protein_fdr, max_run,
                                                 output_dir_or, method)
    return


if __name__ == "__main__":
    main(sys.argv)
    clean()
    print("Done!")
