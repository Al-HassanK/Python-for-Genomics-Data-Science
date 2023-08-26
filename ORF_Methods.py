from MultiFastaAnalyzer import *


### Pre-conditions: takes a DNA sequence and the start position of the wanted reading frame...
### Post-conditions: returns a reading frame in a list...
def Generate_Reading_Frame(sequence, start_position):
    reading_frame = []
    if start_position != 0:
        reading_frame = [sequence[:start_position]]
    
    # Generate the triplets
    for i in range(start_position, len(sequence), 3):
        reading_frame.append(sequence[i:i+3])
    
    return reading_frame
    
### Pre-conditions: takes a DNA sequence...
### Post-conditions: a dictionary of the three reading frames...
def Get_Reading_Frames(sequence):
    reading_frames = {}

    # Fill the dictionary with the three reading frames...
    for i in range(3):
        reading_frames[i+1] = Generate_Reading_Frame(sequence, i)
    
    return reading_frames

### Pre-conditions: takes a reading frame...
### Post-conditions: returns true if it has stop and start codons false other wise...
def Check_Boundary_Codons(reading_frame):
    # Check whether the sequence has start codon...
    first_condition = "ATG" in reading_frame or "atg" in reading_frame
    # Check whether the sequence has stop codon...
    second_condition = False
    for stop_codon in ["TAA", "TAG", "TGA"]:
        if stop_codon in reading_frame or stop_codon.lower() in reading_frame:
            second_condition = True
            break
    
    # an ORF is only determined if it has a start and a stop codon...
    if first_condition and second_condition:
        return True
    else:
        return False

### Pre-conditions: takes a reading frame that contains start and stop condons...
### Post-conditions: Returns the positions of the start condons in the reading frame...
def Get_Start_Codon_Position(reading_frame):
    positions = []
    for pos, codon in enumerate(reading_frame):
        if codon == "ATG" or codon == "atg":
            positions.append(pos)
    
    return positions

### Pre-conditions: takes a reading frame that contains start and stop condons...
### Post-conditions: Returns the position of all possible stop condons in the reading frame...
def Get_Stop_Codon_Position(reading_frame):
    positions = []
    stop_codons = ["TAA", "TAG", "TGA", "taa", "tag", "tga"]
    
    for pos, codon in enumerate(reading_frame):
        if codon in stop_codons:
            positions.append(pos)

    
    return positions

### Pre-conditions: takes a sequence and the reading frame number you want to start search for ORF...
### Post-conditions: Returns a list of all ORFs in that sequence...
def Get_ORF_For_Reading_Frame(sequnece, reading_frame_number):
    reading_frame = Generate_Reading_Frame(sequnece, reading_frame_number - 1)

    orfs = []
    if Check_Boundary_Codons(reading_frame):
        start_codon_positions = Get_Start_Codon_Position(reading_frame)
        for start_codon_position in start_codon_positions:
            if Check_Boundary_Codons(reading_frame[start_codon_position:]):
                stop_codon_position = Get_Stop_Codon_Position(reading_frame[start_codon_position:])[0] + start_codon_position
                orfs.append(''.join(reading_frame[start_codon_position:stop_codon_position+1]))

    return orfs        


### Pre-conditions: takes a path to a Fasta file with multiple records and the reading frame number...
### Post-conditions: returns a dictionary contains each sequence identifier and the all possible ORF for each sequence
# in the specified reading frame...
def Get_ORF_Dictionary(file_path, reading_frame_number):
    records = File_Parser(file_path)
    orf_dict = {}

    for identifier, sequence in records.items():
        orf_dict[identifier] = Get_ORF_For_Reading_Frame(sequence, reading_frame_number)
    
    return orf_dict

### Pre-conditions: takes a path to a Fasta file with multiple records and the reading frame number...
### Post-conditions: returns a list containing pairs of the longest orfs and their identifier and the last element is the 
# longest orf length...
def Find_The_Longest_ORF(file_path, reading_frame_number):
    orf_dict = Get_ORF_Dictionary(file_path, reading_frame_number)
    all_orf_lengths = []
    for _, orf_list in orf_dict.items():
        orf_lengths = [len(orf) for orf in orf_list]
        all_orf_lengths.extend(orf_lengths)
    
    longest_orf_length = max(all_orf_lengths)

    longest_orfs = []
    for identifier, orf_list in orf_dict.items():
        orf_lengths = [len(orf) for orf in orf_list]
        if longest_orf_length in orf_lengths:
            orf_index = orf_lengths.index(longest_orf_length)
            orf = orf_list[orf_index]
            longest_orfs.append((identifier, orf, ))

    longest_orfs.append(longest_orf_length)

    return longest_orfs

### Pre-conditions: takes a path to a Fasta file with multiple records and the reading frame number...
### Post-conditions: returns a list containing pairs of the longest orfs and their identifier and the last element is the 
# longest orf length...
def Find_The_Shortest_ORF(file_path, reading_frame_number):
    orf_dict = Get_ORF_Dictionary(file_path, reading_frame_number)
    all_orf_lengths = []
    for _, orf_list in orf_dict.items():
        orf_lengths = [len(orf) for orf in orf_list]
        all_orf_lengths.extend(orf_lengths)
    
    shortest_orf_length = min(all_orf_lengths)

    shortest_orfs = []
    for identifier, orf_list in orf_dict.items():
        orf_lengths = [len(orf) for orf in orf_list]
        if shortest_orf_length in orf_lengths:
            orf_index = orf_lengths.index(shortest_orf_length)
            orf = orf_list[orf_index]
            shortest_orfs.append((identifier, orf))

    shortest_orfs.append(shortest_orf_length)

    return shortest_orfs

### Pre-conditions: Takes a DNA sequence...
### Post-conditions: Returns a dictionary contains a reading frame and a tuple contains the length of the longest ORF, its start
# position, and the ORF itself...  
def Find_Longest_ORF_In_Sequence(sequence):
    max_orf = {}

    orfs_dict = {f"{i}":Get_ORF_For_Reading_Frame(sequence, i) for i in range(1, 4)}

    for key, orfs in orfs_dict.items():
        max_length = max([len(orf) for orf in orfs])
        for orf in orfs:
            if len(orf) == max_length:
                orf_start_position = sequence.find(orf)
                max_orf[key] = (max_length, orf_start_position, orf)

    return max_orf

def Find_Longest_ORF_In_File(file_path):
    max_orfs = {}

    records = File_Parser(file_path)

    for identifier, sequence in records.items():
        max_orfs_record = Find_Longest_ORF_In_Sequence(sequence)
        max_orfs[identifier] = max([value[0] for _, value in max_orfs_record.items()])

    max_value = max(max_orfs.values())
    max_ids = []

    for identifier, value in max_orfs.items():
        if value == max_value:
            max_ids.append(identifier)
    
    return (max_ids, max_value)

