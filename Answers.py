from MultiFastaAnalyzer import *

def Answer_Q1_2(input_file_path, output_file_path):
    records = File_Parser(input_file_path)

    Analyze_Sequences_Lengths(records, output_file_path)

def Answer_Q3(input_file_path, output_file_path):
    records = File_Parser(input_file_path)

    orf_dict = Get_ORF_Dictionary(records)

    longest_orf_dict = Get_Longest_ORF(orf_dict)
    sorted_longest_orf_dict = dict(sorted(longest_orf_dict.items(), key=lambda x:x[1][1]))
    with open(output_file_path, 'a+') as outfile:
        outfile.write(f'{"*"*20}Answers to the third question{"*"*20}\n')
        for key, value in sorted_longest_orf_dict.items():
           outfile.write(f"{5*'#'}For the record with the identifier: {key}{5*'#'}\n")
           for orf in value[0]:
                outfile.write(f"\t\tThe ORF starting at position of {records[key].index(orf[0])} HAS A LENGTH OF {orf[1]}\n")


def Answer_Q4(input_file_path, output_file_path, repeat_len):
    pass


