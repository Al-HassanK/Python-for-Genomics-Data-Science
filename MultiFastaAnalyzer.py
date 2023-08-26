
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

### Pre-conditions: takes a path to a Fasta file with multiple records...
### Post-conditions: Returns the number of the records in that file...
def Get_Number_Of_Records(file_path):
    return len(File_Parser(file_path))

### Pre-conditions: takes a path to a Fasta file with multiple records...
### Post-conditions: Returns a dictionary contains each sequence identifier and its length...
def Get_Lengths_Of_Sequences(file_path):
    records = File_Parser(file_path)
    sequence_length_dict = {identifier:len(sequence) for identifier, sequence in records.items()}
    
    return sequence_length_dict

### Pre-conditions: takes a path to a Fasta file with multiple records...
### Post-conditions: Returns a tuple consists of two elements the first is a list contains all identifiers with the longest sequence,
# and the second element is the length of the longest sequence...
def Find_Longest_Sequences(file_path):
    sequence_length_dict = Get_Lengths_Of_Sequences(file_path)
    records = File_Parser(file_path)

    longest_sequence_length = max(sequence_length_dict.values())
    longest_sequences_identifiers = []

    for identifer, sequence in records.items():
        if len(sequence) == longest_sequence_length:
            longest_sequences_identifiers.append(identifer)

    return (longest_sequences_identifiers, longest_sequence_length)


### Pre-conditions: takes a path to a Fasta file with multiple records...
### Post-conditions: Returns a tuple consists of two elements the first is a list contains all identifiers with the shortest sequences,
# and the second element is the length of the shortest sequence...
def Find_Shortest_Sequences(file_path):
    sequence_length_dict = Get_Lengths_Of_Sequences(file_path)
    records = File_Parser(file_path)

    shortest_sequence_length = min(sequence_length_dict.values())
    shortest_sequences_identifiers = []

    for identifer, sequence in records.items():
        if len(sequence) == shortest_sequence_length:
            shortest_sequences_identifiers.append(identifer)

    return (shortest_sequences_identifiers, shortest_sequence_length)

### Pre-conidtions: takes a dictionary contains a pattern of specific length and its number of occurences...
### Post-conditions: Returns a dictionary with removed patterns that has an occurence equals one...
def Filter_Repeat_Counts(repeat_counts):
    filtered_repeat_count = {}
    for pattern, count in repeat_counts.items():
        if count > 1:
            filtered_repeat_count[pattern] = count

    return filtered_repeat_count
### Pre-conidtions: takes the length of the wanted repeat, the sequence...
### Post-conditions: returns a dictionary of all repeats and their number of occurences...
def Find_Repeats(repeat_len, sequence):
    repeat_counts = {}
    for i in range(0, len(sequence) - repeat_len + 1):
        subsequence = sequence[i:i+repeat_len]
        if repeat_counts.get(subsequence, None):
            repeat_counts[subsequence] += 1
        else:
            repeat_counts[subsequence] = 1
        

    return Filter_Repeat_Counts(repeat_counts)

### Pre-conditions: takes a path to a Fasta file with multiple records and the length of the repeat...
### Post-conditions: Returns a dictionary that contains each repeat and its count in the whole file...
def Find_All_Repeats(file_path, repeat_len):
    records = File_Parser(file_path)
    all_repeats_counts = {}
    for _, sequenece in records.items():
        repeats_counts = Find_Repeats(repeat_len, sequenece)
        for pattern, count in repeats_counts.items():
            if all_repeats_counts.get(pattern, None):
                all_repeats_counts[pattern] += count
            else:
                all_repeats_counts[pattern] = count
    
    return all_repeats_counts

### Pre-conditions: takes a path to a Fasta file with multiple records and the length of the repeat...
### Post-conditions: Returns a dictionary that contains the most frequent repeats in the whole file...
def Find_Most_Frequent_Repeat(file_path, repeat_len):
    repeat_counts = Find_All_Repeats(file_path, repeat_len)

    most_frequent_value = max(repeat_counts.values())
    most_frequent_repeats = {}

    for pattern, count in repeat_counts.items():
        if count == most_frequent_value:
            most_frequent_repeats[pattern] = count

    return most_frequent_repeats


### Pre-conditions: takes a path to a Fasta file with multiple records and the length of the repeat...
### Post-conditions: Returns a dictionary that contains the less frequent repeats in the whole file...
def Find_Less_Frequent_Repeat(file_path, repeat_len):
    repeat_counts = Find_All_Repeats(file_path, repeat_len)

    less_frequent_value = min(repeat_counts.values())
    less_frequent_repeats = {}

    for pattern, count in repeat_counts.items():
        if count == less_frequent_value:
            less_frequent_repeats[pattern] = count

    return less_frequent_repeats


