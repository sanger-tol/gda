#!/usr/bin/env python3
"""
Script for combining the output GFF3 files of Augustus, Barrnap and tRNAscan into one GFF3 file
"""

import os.path, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import general_purpose_functions as gpf
import os
import argparse

def main(augustus_gff_path, barrnap_gff_path, trnascan_gff_path, fasta_path, out_folder, species_id):
    fasta_basename = fasta_path.split("/")[-1]
    fasta_basename = fasta_basename.split(".")[0]

    out_path = out_folder + "/" + fasta_basename + "_annotations.gff3"

    trnascan_data = gpf.ll(trnascan_gff_path)
    barrnap_data = gpf.ll(barrnap_gff_path)
    augustus_data = gpf.ll(augustus_gff_path)

    tempfile1_path = out_folder + "/" + fasta_basename + "_temp1.gff3"
    tempfile2_path = out_folder + "/" + fasta_basename + "_temp2.gff3"

    with open(tempfile1_path, "w") as f:
        for line in trnascan_data:
            if line.startswith("#") == False:
                f.write(line + "\n")
            else:
                if line == "##gff-version 3" or line.startswith("##sequence-region"):
                    f.write(line + "\n")
        for line in barrnap_data:
            if line.startswith("#") == False:
                line = line.replace(" ", "_")
                f.write(line + "\n")
        for line in augustus_data:
            if line.startswith("#") == False:
                split_line = line.split()
                if len(split_line) >= 3:
                    if split_line[2] == "transcript":
                        split_line[2] = "mRNA"
                f.write("\t".join(split_line) + "\n")
        
    os_command = "gt gff3 -sort -retainids -tidy {} > {}".format(tempfile1_path, tempfile2_path)
    gpf.run_system_command(os_command)

    os_command2 = "add_missing_ids_to_gff3.py {} {} > {}".format(tempfile2_path, species_id, out_path)
    gpf.run_system_command(os_command2)

    os.remove(tempfile1_path)
    os.remove(tempfile2_path)
    os.remove(augustus_gff_path)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("augustus_gff_path", type=str, help="Path to GFF3 file produced by Augustus")
    parser.add_argument("barrnap_gff_path", type=str, help="Path to GFF3 file produced by Barrnap")
    parser.add_argument("trnascan_gff_path", type=str, help="Path to GFF3 file produced by tRNAscan-SE")
    parser.add_argument("fasta_path", type=str, help="Path to input FASTA file")
    parser.add_argument("out_folder", type=str, help="Folder for outputting the merged GFF3 file")
    parser.add_argument("species_id", type=str, help="Species ID to be used in the output GFF3 file")
    args = parser.parse_args()
    main(args.augustus_gff_path, args.barrnap_gff_path, args.trnascan_gff_path, args.fasta_path, args.out_folder, args.species_id)
    









