from MultiFastaAnalyzer import *
from ORF_Methods import *

# This function is called when the user write --nor argument...
def Print_Number_Of_Records(file_path):
    print(f"The number of records in the file are: {Get_Number_Of_Records(file_path)}")

# This function is called when the user write --longestseq argument...
def Print_Longest_Sequences(file_path):
    longest_sequences = Find_Longest_Sequences(file_path)

    print(f"\nThe Following identifiers have the longest sequence in the file which equals to {longest_sequences[-1]}\n")

    for identifier in longest_sequences[0:-1]:
        print(f"{identifier[0]}\n")

# This function is called when the user write --shortestseq argument...
def Print_Shortest_Sequences(file_path):
    shortest_sequences = Find_Shortest_Sequences(file_path)

    print(f"\nThe Following identifiers have the shortest sequence in the file which equals to {shortest_sequences[-1]}\n")

    for identifier in shortest_sequences[0:-1]:
        print(f"{identifier[0]}\n")

# This function is called when the user write --longestOrflen <reading frame> argument...
def Print_Longest_ORF_Len(file_path, reading_frame):

    longest_orfs = Find_The_Longest_ORF(file_path, reading_frame)

    print(f"The length of the longest ORF in the reading frame {reading_frame} is: {longest_orfs[-1]}")


# This function is called when the user write --longestOrfid <reading frame> argument...
def Print_Longest_ORF_Identifier(file_path, reading_frame):
    
    longest_orfs = Find_The_Longest_ORF(file_path, reading_frame)

    print(f"The identifiers that has the lonest ORF in the reading frame {reading_frame} are :\n")
    for orf_identifier in longest_orfs[0:-1]:
        print(f"{orf_identifier[0]}\n")


# This function is called when the user write --getlongestOrf <identifier> argument...
def Print_Longest_ORF_To_Passed_Identifier(file_path, user_identifier):
    records = File_Parser(file_path)
    user_sequence = None
    for identifier, sequence in records.items():
        if user_identifier in identifier:
            user_sequence = sequence
    
    if user_sequence:
        max_orfs = Find_Longest_ORF_In_Sequence(user_sequence)
        print(max_orfs)
    
    else:
        print("The id is not found")


# This function is called when the user write --findrepeats <repeat len> argument...
def Print_Repeats(file_path, repeat_len):
    repeats_count = Find_All_Repeats(file_path, repeat_len)

    for repeat, count in repeats_count.items():
        print(f"The Pattern {repeat} occures in the file {count} times")


# This function is called when the user write --mostfreq <repeat len> argument...
def Print_Most_Frequent_Repeat(file_path, repeat_len):
    most_freq = Find_Most_Frequent_Repeat(file_path, repeat_len)

    for repeat, count in most_freq.items():
        print(f"The Pattern {repeat} occures in the file {count} times")


# This function is called when the user write --longestOrfstart <reading frame> argument...
def Print_Longest_ORF_Start_Position(file_path, reading_frame):
    longest_orfs = Find_The_Longest_ORF(file_path, reading_frame)

    records = File_Parser(file_path)

    for element in longest_orfs[0:-1]:
        print(f"The start position for sequence with id {element[0]} is {records[element[0]].find(element[1])}")

    