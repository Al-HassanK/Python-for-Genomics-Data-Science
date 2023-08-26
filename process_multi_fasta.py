from OutputMethods import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--nor", help="print the number of records in the file", action="store_true")

parser.add_argument("--longestseq", help="print the identifier of the longest sequence as well to its length", action="store_true")

parser.add_argument("--shortestseq", help="print the identifier of the sortest sequence as well to its length", action="store_true")

parser.add_argument("--longestOrflen", help="print the length of the longest ORF in the file in a specified <reading frame> (1, 2, 3)",
                    type=int)

parser.add_argument("--longestOrfid", help="print the id of the longest ORF in the file in a specified <reading frame> (1, 2, 3)",
                    type=int)

parser.add_argument("--getlongestOrf", help="print the longest ORF in all reading frames to the specified sequence <identifier> as\
                    well to the starting position of each ORF", type=str)

parser.add_argument("--findrepeats", help="print all repeats of length <repeat len> in all sequences in the file", type=int)

parser.add_argument("--mostfreq", help="print the most frequent repeat in the whole file", type=int)

parser.add_argument("--longestOrfstart",help="print the start position of the longest ORF in the file in for all reading frames",
                     type=int)

parser.add_argument("infile", help="The multiFasta file should be in fasta format", type=str)


args = parser.parse_args()

if args.nor:
    Print_Number_Of_Records(args.infile)

if args.longestseq:
    Print_Longest_Sequences(args.infile)

if args.shortestseq:
    Print_Shortest_Sequences(args.infile)

if args.longestOrflen:
    Print_Longest_ORF_Len(args.infile, args.longestOrflen)

if args.longestOrfid:
    Print_Longest_ORF_Identifier(args.infile, args.longestOrfid)

if args.getlongestOrf:
    Print_Longest_ORF_To_Passed_Identifier(args.infile, args.getlongestOrf)

if args.findrepeats:
    Print_Repeats(args.infile, args.findrepeats)

if args.mostfreq:
    Print_Most_Frequent_Repeat(args.infile, args.mostfreq)

if args.longestOrfstart:
    Print_Longest_ORF_Start_Position(args.infile, args.longestOrfstart)
