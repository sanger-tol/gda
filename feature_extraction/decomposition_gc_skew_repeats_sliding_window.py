#!/usr/bin/env python3
"""
Script for finding sequence GC, GC skew, stop codon frequency and frequency of some repeat motifs in a FASTA file with sliding window
"""

import general_purpose_functions as gpf
import sys
import argparse

STOP_CODON_SEQ = ("TAG", "TAA", "TGA")
CAG_REPEAT_MOTIFS = ("TGC", "CAG", "GCT")

TELOMERIC_SEQ_DICT = {"vertebrates": "TTAGGG",
"paramecium": "TTGGGT,TTGGGG",
"apicomplexan":	"TTAGGGT,TTAGGGC",
"oxytricha": "TTTTGGGG",
"arabidopsis_thaliana":	"TTTAGGG",
"cestrum_elegans": "TTTTTTAGGG",
"allium": "CTCGGTTATGGG",
"zostera_marina": "TTAGGG",
"green_algae": "TTTTAGGG",
"insects": "TTAGG",
"roundworms": "TTAGGC",
"saccharomyces_cerevisiae":	"GGTGT,GGTGTGT,GGTGTGTGT,GGTGTGTGTGT,GGTGTGTGTGTGT,GGTGTGTGTGTGTGT,GGGTGT,GGGTGTGT,GGGTGTGTGT,GGGTGTGTGTGT,GGGTGTGTGTGTGT,GGGTGTGTGTGTGTGT",
"saccharomyces_castellii": "TCTGGGTG",
"candida_glabrata":	"GGGGTCTGGGTGCTG",
"candida_albicans":	"GGTGTACGGATGTCTAACTTCTT",
"candida_tropicalis": "GGTGTACGGATGTCACGATCATT,GGTGTAAGGATGTCACGATCATT",
"candida_maltosa": "GGTGTACGGATGCAGACTCGCTT",
"candida_guillermondii": "GGTGTAC",
"candida_pseudotropicalis": "GGTGTACGGATTTGATTAGTTATGT",
"kluyveromyces_lactis": "GGTGTACGGATTTGATTAGGTATGT",
"schizosaccharomyces_pombe": "TTACAG,TTACAGG,TTACAGGG,TTACAGGGG,TTACAGGGGG,TTACAGGGGGG,TTACAGGGGGGG,TTACAGGGGGGGG,TTACCG,TTACCGG,TTACCGGG,TTACCGGGG,TTACCGGGGG,TTACCGGGGGG,TTACCGGGGGGG,TTACCGGGGGGGG,TTACG,TTACGG,TTACGGG,TTACGGGG,TTACGGGGG,TTACGGGGGG,TTACGGGGGGG,TTACGGGGGGGG",
"dictyostelium": "AG,AGG,AGGG,AGGGG,AGGGGG,AGGGGGG,AGGGGGGG,AGGGGGGGG",
"tetrahymena": "TTGGGG"}


def get_query_seq_freq(seq, query_strings, min_repeat_units):
    """
    Finds the frequency of matches to a list of query strings in a sequence (string)
    """
    match_count = 0
    for query_string in query_strings:
        rev_comp_query_string = gpf.reverse_complement(query_string)
        match_count += seq.count(query_string * min_repeat_units)
        match_count += seq.count(rev_comp_query_string * min_repeat_units)
        match_freq = match_count / len(seq)
    return match_freq


def get_seq_gc_and_skew(seq):
    """
    Finds the GC% and skew of a DNA sequence
    """
    gc = None
    g_count = seq.count("G")
    c_count = seq.count("C")
    a_count = seq.count("A")
    t_count = seq.count("T")
    gc_count = g_count + c_count
    atgc_count = gc_count + a_count + t_count
    if atgc_count > 0:
        gc = (gc_count / atgc_count) * 100
    gc_skew = None
    at_skew = None
    if g_count + c_count > 0:
        gc_skew = (g_count - c_count) / (g_count + c_count)
    if a_count + t_count > 0:
        at_skew = (a_count - t_count) / (a_count + t_count)
    return gc, gc_skew, at_skew


def get_cpg_percentage(seq):
    """
    Finds the percentage of CpG sites (counted on one strand)
    """
    cpg_percentage = None
    seq_len = len(seq)
    if seq_len > 0:
        cpg_percentage = (seq.count("CG") / seq_len) * 100
    return cpg_percentage


def validate_custom_telomeric_seq(custom_telomeric_seq):
    """
    Function for validating the characters in the user-provided telomeric sequence
    """
    custom_telomeric_seq = custom_telomeric_seq.upper()
    allowed_chars = "ATGC,"
    chars_ok_flag = all(c in allowed_chars for c in allowed_chars)
    if chars_ok_flag == False:
        sys.stderr.write("Invalid characters found in the custom telomeric sequence provided by the user: " + custom_telomeric_seq + "\n")
        sys.exit(1)


def main(in_path, out_folder, pipeline_output_folder, chunk_size, min_repeat_unit_count, telomeric_seq_preset, custom_telomeric_seq):

    telomeric_seq = None
    if custom_telomeric_seq is not None:
        validate_custom_telomeric_seq(custom_telomeric_seq)
        telomeric_seq = custom_telomeric_seq.split(",")
    else:
        telomeric_seq = TELOMERIC_SEQ_DICT[telomeric_seq_preset].split(",")

    gpf.run_system_command("mkdir -p " + out_folder)
    fasta_filename = in_path.split("/")[-1]
    fasta_basename = fasta_filename.split(".")[0]

    out_path = out_folder + "/" + fasta_basename + "_gc_skew_repeats.txt"
    with open(out_path, "w") as f:
        fasta_data = gpf.read_fasta_in_chunks(in_path)
        f.write("scaffold\tstart_pos\tend_pos\tgc_percentage\tgc_skew\tat_skew\tstop_codon_freq\ttelomere_freq\tcag_freq\tcpg_percentage\n")
        for header, seq in fasta_data:
            seq = seq.upper()
            seq_chunks = gpf.string_to_chunks(seq, chunk_size)
            pos = 0
            for chunk in seq_chunks:
                gc, gc_skew, at_skew = get_seq_gc_and_skew(chunk)
                stop_codon_freq = get_query_seq_freq(chunk, STOP_CODON_SEQ, 1)
                telomer_freq = get_query_seq_freq(chunk, telomeric_seq, min_repeat_unit_count)
                cag_freq = get_query_seq_freq(chunk, CAG_REPEAT_MOTIFS, min_repeat_unit_count)
                cpg_percentage = get_cpg_percentage(chunk)
                chunk_end = pos + len(chunk)
                out_line = "\t".join([header, str(pos), str(chunk_end), str(gc), str(gc_skew), str(at_skew), str(stop_codon_freq), str(telomer_freq), str(cag_freq), str(cpg_percentage)])
                f.write(out_line + "\n")
                pos += chunk_size
    gpf.run_system_command("gc_skew_etc_to_bedgraph.py {} {} {}".format(out_path, pipeline_output_folder, in_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("in_path", type=str, help="Path to DNA FASTA file")
    parser.add_argument("out_folder", type=str, help="Folder for output file")
    parser.add_argument("pipeline_output_folder", type=str, help="Folder for outputting the bedgraph files")
    parser.add_argument("--chunk_size", type=int, help="Sliding window step size (bp), default: 5000", default=5000)
    parser.add_argument("--min_repeat_unit_count", type=int, help="Minimum number of repeat units in repeats", default=3)
    parser.add_argument("--telomeric_seq_preset", type=str, help="Telomeric sequence type preset (default: vertebrates)", default="vertebrates", choices=["vertebrates", "paramecium", "apicomplexan", "oxytricha", "arabidopsis_thaliana", "cestrum_elegans", "allium", "zostera_marina", "green_algae", "insects", "roundworms", "saccharomyces_cerevisiae", "saccharomyces_castellii", "candida_glabrata", "candida_albicans", "candida_tropicalis", "candida_maltosa", "candida_guillermondii", "candida_pseudotropicalis", "kluyveromyces_lactis", "schizosaccharomyces_pombe", "dictyostelium" "tetrahymena"])
    parser.add_argument("--custom_telomeric_seq", type=str, help="Optional: telomeric sequence(s) as a comma separated string, e.g. 'TTAGGGT,TTAGGGC'. If a sequence is entered here, it overrides the telomeric sequence type preset. Default: None", default=None)
    args = parser.parse_args()
    main(args.in_path, args.out_folder, args.pipeline_output_folder, args.chunk_size, args.min_repeat_unit_count, args.telomeric_seq_preset, args.custom_telomeric_seq) 

   