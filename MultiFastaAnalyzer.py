
### Pre-conditions: takes the header line of a fasta record in the gi type...
### Post-conditions: returns the idetifier of the record...
def Get_Record_Id(line_header):
    end_index = line_header.find(',')

    return line_header[1:end_index]

### Pre-conditions: takes a path to a Fasta file with multiple records...
### Post-conditions: returns a dictionary that has sequence identifier as a key and the value is the sequence itself...
def File_Parser(file_path):
    records = {}
    record_id = ""

    with open(file_path) as infile:
        for line in infile:
            if line[0] == '>':
                record_id = Get_Record_Id(line.rstrip())
                records[record_id] = ""
            
            else:
                records[record_id] += line.rstrip()

    return records


### Pre-conditions: takes the records dectionary returned from the File_Parser function and a file path it will write results to...
### Post-conditions: writes the following: how many records in the file, the lengths of the sequences, 
# and the longest and shortest sequences else a sorted dictionary based on sequence length will be returned...
def Analyze_Sequences_Lengths(records_dict, output_file_path=None):
    number_of_records = len(records_dict)
    lengths_dict = {}

    # Get the lengths of the sequences...
    for key, value in records_dict.items():
        lengths_dict[key] = len(value)
    
    # Sort the dictionary based on the sequence length...
    sorted_lengths_dict = dict(sorted(lengths_dict.items(), key=lambda x:x[1]))

    # Write the data to the file if a file path is assigned else returns the sorted dictionary...
    if output_file_path:
        with open(output_file_path, 'a+') as outfile:
            outfile.write(f"There are {number_of_records} records in the file\n")
            outfile.write(f'{"*"*20}The records sroted based on their sequence length{"*"*20}\n')
            for key, value in sorted_lengths_dict.items():
                outfile.write(f"{key}: {value}\n")

    else:
        return sorted_lengths_dict

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

def Get_Start_Codon_Position(reading_frame):
    try:
        pos = reading_frame.index("ATG")
    except ValueError:
        pos = reading_frame.index("atg")
    
    return pos

def Get_Stop_Codon_Position(reading_frame):
    positions = []
    for stop_codon in ["TAA", "TAG", "TGA"]:
        if stop_codon in reading_frame or stop_codon.lower() in reading_frame:
            try:
                pos = reading_frame.index(stop_codon)
                positions.append(pos)
            except ValueError:
                pos = reading_frame.index(stop_codon.lower())
                positions.append(pos)
    
    return positions

### Pre-conditions: takes a DNA sequence...
### Post-conditions: returns its ORF, returns None if there is no ORF...
def Get_ORF(sequence):
    reading_frames = Get_Reading_Frames(sequence)
    orfs = []
    for _, value in reading_frames.items():
        if Check_Boundary_Codons(value):
            start_codon_pos = Get_Start_Codon_Position(value)
            stop_codon_positions = Get_Stop_Codon_Position(value)
            for pos in stop_codon_positions:
                orfs.append(''.join(value[start_codon_pos:pos+1]))

    return orfs

### Pre-conditions: takes the records dectionary returned from the File_Parser function...
### Post-conditions: returns a dictionary contains the same keys but the values are the all possible ORF for each sequence...
def Get_ORF_Dictionary(fasta_records):
    orf_dict = {}

    for key, value in fasta_records.items():
        orf_dict[key] = Get_ORF(value)
    
    return orf_dict

### Pre-conditions: takes the list containing group of orfs of a sequence...
### Post-conditions: returns a list containing pairs of the longest sequences and their length...
def Find_The_Longest_ORF(orf_list):
    max_len = sorted([len(orf) for orf in orf_list])[-1] # the maximum length between all the orfs...

    # Get the orfs equal to the maximum length 
    longest_orfs = []
    for orf in orf_list:
        if len(orf) == max_len:
            longest_orfs.append((orf, len(orf)))
            max_len = len(orf)
    
    return longest_orfs

### Pre-conditions: takes the orf dectionary returned from the Get_ORF_Dictionary function...
### Post-conditions: returns the longest orf of all the sequences...
def Get_Longest_ORF(orf_dict):
    orf_length_dict = {}

    for key, value in orf_dict.items():
        orf_length_dict[key] = (Find_The_Longest_ORF(value), Find_The_Longest_ORF(value)[0][1])

    return orf_length_dict



### Pre-conidtions: takes the length of the wanted repeat, the sequence, and optional a threshold that is the minmum
# number of occurences each repeat has...
### Post-conditions: returns a dictionary of all repeats and their number of occurences...
def Find_Repeats(repeat_len, sequence, threshold=1):
    repeat_counts = {}
    for i in range(0, len(sequence) - (repeat_len - 1)):
        subsequence = sequence[i:i+repeat_len]
        if sequence.count(subsequence) >= threshold:
            repeat_counts[subsequence] = sequence.count(subsequence)

    return repeat_counts

