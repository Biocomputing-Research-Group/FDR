import os
import re
import shutil

# File path to your PIN file
path_to_input = 'U1.pin'
target_fdr = 0.01
decoy_label = "Rev_"

# IO
path_to_temp_files = './temp'
path_to_species_pin = './temp/species_pin'


def initialize():
    print("Processing...")
    if os.path.exists(path_to_temp_files) is True:
        shutil.rmtree(path_to_temp_files)
    if os.path.exists(path_to_temp_files) is False:
        os.makedirs(path_to_temp_files)
    if os.path.exists(path_to_species_pin) is False:
        os.makedirs(path_to_species_pin)
    return


def create_new_txt_file_spid(path_to_pin_input, path_to_pin_output):
    with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
        next(pin_input)
        header = "SpecId\t" + pin_input.readline()
    with open(path_to_pin_output, "w") as pin_output:
        pin_output.write(header)
        pin_output.close()
    return


def read_txt_header(path_to_pin_input):
    with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
        header = pin_input.readline()
        pin_input.close()
    return header


def create_new_txt_file(path_to_pin_input, path_to_pin_output):
    header = read_txt_header(path_to_pin_input).replace("\n", "") + "\tp_score\t" + "posterior_error_prob\n"
    with open(path_to_pin_output, "w") as pin_output:
        pin_output.write(header)
        pin_output.close()
    return


def create_new_pin_file(path_to_pin_input, path_to_pin_output):
    with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
        header = pin_input.readline()
    with open(path_to_pin_output, "w") as pin_output:
        pin_output.write(header)
        pin_output.close()
    return


def combine_txt_files(input_pins_dir, postfix, output_dir):
    input_dir = os.listdir(input_pins_dir)
    for file in input_dir:
        if file.find(postfix) != -1:
            create_new_txt_file_spid(input_pins_dir + file, output_dir + 'combined' + postfix)
    with open(output_dir + 'combined' + postfix, 'a', encoding='utf-8') as f1:
        for file in input_dir:
            if file.find(postfix) != -1:
                with open(input_pins_dir + file, 'r', encoding='utf-8') as f2:
                    previous_scan = -1
                    previous_charge = -1
                    rank = 0
                    title = f2.readline().split("\t")
                    file_name = title[1]
                    next(f2)
                    for line in f2:
                        current = line.split("\t")
                        current_scan = current[0]
                        current_charge = current[2]
                        if current_scan == previous_scan and current_charge == previous_charge:
                            rank += 1
                        else:
                            rank = 1
                        previous_scan = current_scan
                        previous_charge = current_charge
                        SpecId = str(file_name) + "_" + str(current_scan) + "_" + str(current_charge) + "_" + str(rank)
                        current[1] = str(rank)
                        f1.write(SpecId + '\t' + '\t'.join(current))
    print(input_pins_dir)
    return output_dir + 'combined' + postfix


def combine_pin_files(input_pins_dir, postfix, output_dir):
    input_dir = os.listdir(input_pins_dir)
    for file in input_dir:
        if file.find(postfix) != -1:
            create_new_pin_file(input_pins_dir + file, output_dir + 'combined' + postfix)
    with open(output_dir + 'combined' + postfix, 'a', encoding='utf-8') as f1:
        for file in input_dir:
            if file.find(postfix) != -1:
                with open(input_pins_dir + file, 'r', encoding='utf-8') as f2:
                    next(f2)
                    for line in f2:
                        f1.write(line)
    print(input_pins_dir)
    return output_dir + 'combined' + postfix


def clean_comet(path_to_pin_input):
    print("Selecting best Xcorr in each scans...")
    path_to_cleaned_comet_output_tmp = "./temp/cleaned.txt"
    create_new_txt_file(path_to_pin_input, path_to_cleaned_comet_output_tmp)
    with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
        next(pin_input)
        previous_scan = -1
        previous_charge = -1
        with open(path_to_cleaned_comet_output_tmp, "a") as cleaned_comet_output:
            for line in pin_input:
                current = line.split("\t")
                current_scan = current[1]
                current_charge = current[3]
                if current_scan == previous_scan and current_charge == previous_charge:
                    previous_scan = current_scan
                    previous_charge = current_charge
                else:
                    previous_scan = current_scan
                    previous_charge = current_charge
                    cleaned_comet_output.write(line)
    print("Selecting best Xcorr in each scans: Done!")
    return path_to_cleaned_comet_output_tmp


def divide_species_auto_c(path_to_pin_input, path_to_species_output):
    # Extra protein_id
    print("Dividing Comet result by species...")
    target_rule = r"[a-zA-Z0-9|]+"
    decoy_rule = decoy_label + '([a-zA-Z0-9|]+)'
    header = read_txt_header(path_to_pin_input)
    labels = header.split("\t")
    Proteins_loc = labels.index("protein")
    with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
        next(pin_input)
        for line in pin_input:
            data = line.split("\t")
            species_list = []
            proteins = data[Proteins_loc].split(",")
            for p in proteins:
                current = p
                if current.find(decoy_label) == 0:
                    species = re.findall(decoy_rule, current)
                else:
                    species = re.findall(target_rule, current)
                species_list.append(species[0])
            species_list = list(set(species_list))
            for specie in species_list:
                specie_file = path_to_species_output + "/" + specie + ".txt"
                if os.path.exists(specie_file) is False:
                    with open(specie_file, "a") as specie_output:
                        specie_output.write(header)
                with open(specie_file, "a") as specie_output:
                    specie_output.write(line)
    print("Dividing Comet result by species: Done!")
    return path_to_species_output


def divide_species_dic_c(path_to_pin_input, path_to_species_dic, path_to_species_output):
    species_dic = {}
    decoy_rule = decoy_label + '([a-zA-Z0-9|_.]+)'
    print("Reading Species Info...")
    with open(path_to_species_dic, 'r', encoding='utf-8') as sp_dic:
        for line in sp_dic:
            # print(line)
            data = line.split("\t")
            sp_id = data[0]
            sp_name = data[1].replace("\n", "")
            sp_name_list = sp_name.split(",")
            species_dic[sp_id] = sp_name_list
    print("Dividing Comet result by species...")
    header = read_txt_header(path_to_pin_input)
    labels = header.split("\t")
    Proteins_loc = labels.index("protein")
    with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
        next(pin_input)
        for line in pin_input:
            data = line.split("\t")
            species_list = []
            proteins = data[Proteins_loc].split(",")
            for sp in proteins:
                current = sp
                if current.find(decoy_label) == 0:
                    pid = re.findall(decoy_rule, current)
                    for sp in species_dic[pid[0]]:
                        species_list.append(sp)
                else:
                    for sp in species_dic[current]:
                        species_list.append(sp)
            species_list = list(set(species_list))
            for specie in species_list:
                specie_file = path_to_species_output + "/" + specie + ".pin"
                if os.path.exists(specie_file) is False:
                    with open(specie_file, "a") as specie_output:
                        specie_output.write(header)
                with open(specie_file, "a") as specie_output:
                    specie_output.write(line)
    print("Dividing Comet result by species: Done!")
    return path_to_species_output


def percolator_to_comet(target_file, deocy_file, combined_txt, postfix, path_to_converted_output):
    p_PSMId_score_dic = {}
    percolator_input = [target_file, deocy_file]
    percolator_header = read_txt_header(target_file)
    percolator_labels = percolator_header.split("\t")
    p_score_loc = percolator_labels.index("score")
    pep_loc = percolator_labels.index("posterior_error_prob")
    PSMId_loc = percolator_labels.index("PSMId")
    for path_to_pin_input in percolator_input:
        with open(path_to_pin_input, 'r', encoding='utf-8') as pin_input:
            next(pin_input)
            for line in pin_input:
                data = line.split("\t")
                PSMId = data[PSMId_loc]
                p_score = data[p_score_loc]
                pep = data[pep_loc]
                p_PSMId_score_dic[PSMId] = str(p_score) + "\t" + str(pep)
    c_PSMId_record_dic = {}
    with open(combined_txt, 'r', encoding='utf-8') as f2:
        next(f2)
        for line in f2:
            current = line.split("\t")
            PSMId = current[0]
            c_PSMId_record_dic[PSMId] = line
    create_new_txt_file(combined_txt, path_to_converted_output + 'converted' + postfix)
    with open(path_to_converted_output + 'converted' + postfix, 'a', encoding='utf-8') as f3:
        for PSMId in p_PSMId_score_dic.keys():
            record = c_PSMId_record_dic[PSMId]
            data = record.replace("\n", "") + "\t" + p_PSMId_score_dic[PSMId] + "\n"
            f3.write(data)

'''
if __name__ == '__main__':
    initialize()
    # input_pin_dir = "./m2_pin/"
    input_txt_dir = "./u1_txt/"
    combined_txt_dir = path_to_temp_files + "/"
    # target_file = "./target_m1.tsv"
    # deocy_file = "./decoy_m1.tsv"
    # path_to_species_dic = "./humanGut_ProToOTU.tsv"
    combined_txt = "./combined.txt"
    postfix = ".txt"
    path_to_converted_output = "./"
    txt_input = combine_txt_files(input_txt_dir, postfix, combined_txt_dir)
    path_to_cleaned_comet_output = clean_comet(txt_input)
    # divide_species_dic(path_to_cleaned_comet_output, path_to_species_dic, path_to_species_pin)
    divide_species_auto_c(path_to_cleaned_comet_output, path_to_species_pin)
    # combine_pin_files(input_pin_dir, ".pin", combined_txt_dir)
    # divide_species(path_to_sorted_comet_output, path_to_species_dic, path_to_species_pin)
    # percolator_to_comet(target_file, deocy_file, combined_txt, postfix, path_to_converted_output)
    # divide_species_auto("./converted.txt", path_to_species_pin)
    # divide_species_dic("./converted.txt", path_to_species_dic, path_to_species_pin)
    print("Done!")
'''