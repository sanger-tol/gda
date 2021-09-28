#!/usr/bin/env python3
"""
Script for validating the input files of the GDA feature extraction pipeline
"""

import general_purpose_functions as gpf
import mimetypes
import sys
import os
from Bio import SeqIO
import argparse


class UnicodeFileChecker:
    """
    Class for shared functions of checking FASTA and GFF files
    """
    def __init__(self, in_path):
        self.in_path = in_path

    def check_if_file_is_unicode(self):
        """
        Exits if the input file contains non-unicode symbols (is a binary file)
        """
        try:
            with open(self.in_path, "r") as f:
                for l in f:
                    len(l)
        except UnicodeDecodeError:
            sys.stderr.write("The input file " + self.in_path + " contains non-Unicode characters\n")
            sys.exit(1)


    def check_mime(self):
        """
        Uses MIME to check if the file is gzipped or an EMBL file
        """
        mime = mimetypes.guess_type(self.in_path)
        if mime[1] == "gzip":
            sys.stderr.write("The input file " + self.in_path + " appears to be gzipped. gzipped input files of this type are currently not supported\n")
            sys.exit(1)
        if mime[0] == "chemical/x-embl-dl-nucleotide":
            sys.stderr.write("The input file " + self.in_path + " appears to be in EMBL format. This format is currently not supported\n")
            sys.exit(1)


class FastaChecker(UnicodeFileChecker):
    """
    Class for validating FASTA files
    """
    def __init__(self, in_path, sequence_type, max_sequences=None, min_total_seq_len=None, min_longest_seq_len=None):
        super().__init__(in_path)
        self.sequence_type = sequence_type
        self.max_sequences = max_sequences
        self.min_total_seq_len = min_total_seq_len
        self.min_longest_seq_len = min_longest_seq_len


    def is_fasta(self):
        """
        https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
        """
        with open(self.in_path, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


    def check_fasta_contents(self):
        """
        Function for ensuring that a DNA FASTA file or protein FASTA file contains only valid characters for the FASTA type in its sequences
        """
        CHARACTERS_DICT = {"dna": ("A", "T", "G", "C", "N"), "protein": ("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T", "X", "B", "Z", "U", "J", "O", "*", "+", "#")}
        FORBIDDEN_CHARS = [":"]
        fasta_data = gpf.read_fasta_in_chunks(self.in_path)
        characters_outside_dna_set_count = 0
        headers_list = list()
        for header, seq in fasta_data:
            split_header = header.split()
            header_stem = split_header[0]
            for forbidden_char in FORBIDDEN_CHARS:
                if forbidden_char in header_stem:
                    sys.stderr.write("Header in FASTA file {} ({}) contains a character ({}) that is not allowed in the input FASTA headers\n".format(self.in_path, header_stem, forbidden_char))
                    sys.exit(1)
            if header_stem in headers_list:
                sys.stderr.write("Header in FASTA file " + self.in_path + " (" + header_stem + ") is not unique when truncated after first whitespace\n")
                sys.exit(1)
            else:
                headers_list.append(header_stem)
            seq = seq.upper()
            seq_len = len(seq)
            if seq_len == 0:
                sys.stderr.write("Sequence " + header + " in FASTA file " + self.in_path + " appears to be empty\n")
                sys.exit(1)
            counts_sum = sum([seq.count(n) for n in CHARACTERS_DICT["dna"]])
            if self.sequence_type == "protein":
                characters_outside_dna_set_count += seq_len - counts_sum
                counts_sum = sum([seq.count(n) for n in CHARACTERS_DICT["protein"]])
            if counts_sum != seq_len:
                sys.stderr.write("Unsupported characters found in " + self.sequence_type + " sequence of the DNA FASTA file " + self.in_path + ", in sequence " + header + "\n")
                sys.exit(1)
        if self.sequence_type == "protein" and characters_outside_dna_set_count == 0:
            sys.stderr.write("Warning: the file {} that was provided as a protein FASTA file might actually be a DNA FASTA file, based on its content\n".format(self.in_path))

    
    def check_fasta_sequences_count(self):
        """
        Checks whether the number of sequences in the file exceeds a limit (if a limit has been defined)
        """
        if self.max_sequences is not None:
            fasta_data = gpf.read_fasta_in_chunks(self.in_path)
            for items in enumerate(fasta_data):
                counter = items[0]
                if counter > self.max_sequences:
                    sys.stderr.write("The FASTA file {} contains more sequences than allowed. The maximum allowed number of sequences in this file is {}\n".format(self.in_path, self.max_sequences))
                    sys.exit(1)


    def check_seq_lengths(self):
        """
        Checks if the total sequence length of the FASTA file and the length of the longest sequence in the FASTA file fit minimum thresholds
        """
        scaff_sizes = list()
        fasta_data = gpf.read_fasta_in_chunks(self.in_path)
        for fasta_tuple in fasta_data:
            seq = fasta_tuple[1]
            scaff_sizes.append(len(seq))
        longest_seq_len = max(scaff_sizes)
        assembly_size = sum(scaff_sizes)
        if assembly_size < self.min_total_seq_len:
            sys.stderr.write("Input assembly size ({} bp) should not be smaller than the sliding window step size ({} bp)\n".format(str(assembly_size), str(self.min_total_seq_len)))
            sys.exit(1)
        if longest_seq_len < self.min_longest_seq_len:
            sys.stderr.write("The length of the longest sequence in the assembly ({} bp) should not be smaller than the sliding window step size ({} bp)\n".format(str(longest_seq_len), str(self.min_longest_seq_len)))
            sys.exit(1)


    def process_fasta(self):
        """
        Function for triggering multiple steps of checking a FASTA file
        """
        gpf.check_if_file_exists(self.in_path)
        FastaChecker.check_if_file_is_unicode(self)
        FastaChecker.check_mime(self)
        is_fasta_result = FastaChecker.is_fasta(self)
        if is_fasta_result == False:
            sys.stderr.write("The file {} appears to be not a FASTA file\n".format(self.in_path))
            sys.exit(1)
        FastaChecker.check_fasta_sequences_count(self)
        if self.min_longest_seq_len is not None and self.min_total_seq_len is not None:
            FastaChecker.check_seq_lengths(self)
        FastaChecker.check_fasta_contents(self)


class GFF_checker(UnicodeFileChecker):
    """
    Class for validating GFF files
    """
    def __init__(self, gff_path, fasta_path):
        super().__init__(gff_path)
        self.gff_path = gff_path
        self.fasta_path = fasta_path
    
    def gt_validate_gff(self):
        """
        Validating a GFF3 file using GenomeTools GFF3 Validator
        """
        gt_command = "gt gff3validator '" + self.gff_path + "'"
        exit_code = os.system(gt_command)
        if exit_code != 0:
            sys.stderr.write("GFF validation of the file " + self.gff_path + " using gt gff3validator failed\n")
            sys.stderr.write("The issue might get resolved by cleaning the GFF3 file with GenomeTools before running this pipeline. Example:\n")
            sys.stderr.write("gt gff3 -sort -retainids -tidy input.gff3 > input_cleaned.gff3\n")
            sys.exit(1)


    def check_gff_for_required_features(self):
        """
        Function for checking if required GFF features (gene, exon, mRNA and CDS) are present in the input GFF file
        """
        gff_data = gpf.ll(self.gff_path)
        features_dict = {"gene": False, "exon": False, "mRNA": False, "CDS": False}
        for line in gff_data:
            if line.startswith("#") == False:
                split_line = line.split()
                if len(split_line) >= 3:
                    feature_name = split_line[2]
                    if feature_name in features_dict:
                        features_dict[feature_name] = True
                    if features_dict["gene"] == True and features_dict["exon"] == True and features_dict["mRNA"] == True and features_dict["CDS"] == True:
                        break
        false_features = [n for n in features_dict.keys() if features_dict[n] == False]
        false_features = [str(n) for n in false_features]
        if len(false_features) > 0:
            sys.stderr.write("The following feature(s) are required but were not found in the input GFF3 file ({}):\n".format(self.gff_path))
            for feature in false_features:
                sys.stderr.write(feature + "\n")
            sys.stderr.write("If these types of features have other names in your GFF3 file, please change the feature names in the GFF3 file before running this pipeline\n")
            sys.exit(1)


    def check_if_gff_and_fasta_match(self):
        """
        Checking if the GFF file matches the FASTA file it has been paired with
        """

        def get_fasta_seq_lengths(fasta_path):
            """
            Extracts the lengths of sequences in a FASTA file
            """
            seq_lengths_dict = dict()
            fasta_data = gpf.read_fasta_in_chunks(fasta_path)
            for header, seq in fasta_data:
                header = header.split()[0]
                seq_lengths_dict[header] = len(seq)
            return seq_lengths_dict

        seq_lengths_dict = get_fasta_seq_lengths(self.fasta_path)
        gff_data = gpf.ll(self.gff_path)

        for counter, line in enumerate(gff_data):
            if line.startswith("#") == False and len(line) > 3:
                split_line = line.split()
                scaff = split_line[0]
                end_coord = None
                try:
                    end_coord = int(split_line[4])
                except ValueError:
                    sys.stderr.write("Failed to parse scaffold coordinates in GFF file " + self.gff_path + ", line " + str(counter + 1) + "\n")
                    sys.exit(1)
                if scaff in seq_lengths_dict:
                    if end_coord > seq_lengths_dict[scaff]:
                        sys.stderr.write("GFF entry is outside the bounds of a scaffold:\n")
                        sys.stderr.write(line + "\n")
                        sys.exit(1)
                else:
                    sys.stderr.write("Scaffold " + scaff + " appears in the GFF file (" + self.gff_path + ") but was not found in the corresponding FASTA file (" + self.fasta_path + ")\n")
                    sys.exit(1)


    def process_gff(self):
        """
        Function for triggering multiple steps of checking a GFF file
        """
        gpf.check_if_file_exists(self.gff_path)
        GFF_checker.check_if_file_is_unicode(self)
        GFF_checker.check_mime(self)
        GFF_checker.gt_validate_gff(self)
        GFF_checker.check_gff_for_required_features(self)
        GFF_checker.check_if_gff_and_fasta_match(self)



def validate_fastq_files(validatefastq_path, fq1_path, fq2_path):
    """
    Function for validating FASTQ files using https://github.com/biopet/validatefastq
    """
    gpf.check_if_file_exists(fq1_path)
    gpf.check_if_file_exists(fq2_path)
    gt_command = "java -jar {} -i {} -j {}".format(validatefastq_path, fq1_path, fq2_path)
    exit_code = os.system(gt_command)
    if exit_code != 0:
        sys.stderr.write("Validation of the FASTQ files {} and {} using validatefastq failed\n".format(fq1_path, fq2_path))
        sys.exit(1)


class GGTableChecker(UnicodeFileChecker):
    """
    Class for validating FASTA files
    """
    def __init__(self, in_path, proteomes_folder):
        super().__init__(in_path)
        self.proteomes_folder = proteomes_folder

    def validate_csv_for_gg_file(self):
        """
        Function for validating the CSV file that is used to match species IDs with proteome FASTA files for OrthoMCL
        """
        species_list = list()
        prot_files_list = list()
        csv_data = gpf.l(self.in_path)
        for line in csv_data:
            if line.count(",") == 1:
                split_line = line.split(",")
                species = split_line[0]
                prot_file_name = split_line[1]
                if species in species_list:
                    sys.stderr.write("The species name {} in the input file {} appears to be duplicated\n".format(species, self.in_path))
                    sys.exit(1)
                if prot_file_name in prot_files_list:
                    sys.stderr.write("The proteome file name {} in the input file {} appears to be duplicated\n".format(prot_file_name, self.in_path))
                    sys.exit(1)

                protein_fasta_path = None
                if "/" in prot_file_name:
                    protein_fasta_path = prot_file_name
                else:
                    protein_fasta_path = self.proteomes_folder + "/" + prot_file_name
                if os.path.isfile(protein_fasta_path) == False:
                    sys.stderr.write("The proteome file {} that is listed in the CSV file {} was not found\n".format(prot_file_name, self.in_path))
                    sys.exit(1)
                species_list.append(species)
                prot_files_list.append(prot_file_name)
            else:
                sys.stderr.write("Error parsing the CSV file {}: incorrect number of comma separated columns\n".format(self.in_path))
                sys.exit(1)
        if len(prot_files_list) == 0 or len(species_list) == 0:
            sys.stderr.write("No valid entries found in the CSV file {}\n".format(self.in_path))
            sys.exit(1)


    def process_csv_for_gg_file(self):
        """
        Function for triggering multiple steps of checking the CSV for GG file
        """
        gpf.check_if_file_exists(self.in_path)
        GGTableChecker.check_if_file_is_unicode(self)
        GGTableChecker.check_mime(self)
        GGTableChecker.validate_csv_for_gg_file(self)


def validate_orthomcl_references_folder(orthomcl_ref_folder):
    """
    Function for validating the contents of the folder of references for OrthoMCL
    """
    if orthomcl_ref_folder != "NA":
        ref_subfolders = [f.path for f in os.scandir(orthomcl_ref_folder) if f.is_dir()]
        print(ref_subfolders)
        if len(ref_subfolders) == 0:
            sys.stderr.write("No subfolders with reference files for OrthoMCL were found in the OrthoMCL references directory ({})\n".format(orthomcl_ref_folder))
            sys.exit(1)
        else:
            for ref_subfolder in ref_subfolders:
                print(ref_subfolder)
                ref_files = [f for f in os.listdir(ref_subfolder) if os.path.isfile(os.path.join(ref_subfolder, f))]
                print(ref_files)
                if "table_for_gg_file.csv" in ref_files:
                    proteome_fasta_files = [n for n in ref_files if n != "table_for_gg_file.csv"]
                    if len(proteome_fasta_files) > 0:
                        csv_for_gg_file_path = ref_subfolder + "/table_for_gg_file.csv"
                        csv_checker = GGTableChecker(csv_for_gg_file_path, ref_subfolder)
                        csv_checker.process_csv_for_gg_file()
                        print(proteome_fasta_files)
                        for proteome_fasta_file in proteome_fasta_files:
                            proteome_fasta_file_path = ref_subfolder + "/" + proteome_fasta_file
                            fc = FastaChecker(proteome_fasta_file_path, "protein")
                            fc.process_fasta()
                    else:
                        sys.stderr.write("No reference proteome files were found in OrthoMCL reference proteomes folder ({})\n".format(ref_subfolder))
                        sys.exit(1)
                        
                else:
                    sys.stderr.write("The required file named 'table_for_gg_file.csv' was not found in an OrthoMCL reference proteomes folder ({})\n".format(ref_subfolder))
                    sys.exit(1)


def check_filenames_for_bad_characters(filenames_list):
    """
    Function for checking the input file names for characters that might cause problems in downstream processing
    """
    bad_chars = ("?", " ", "\\", "*", "|", ">", "<", ":", ";", "&")
    for filename in filenames_list:
        for char in bad_chars:
            if char in filename:
                sys.stderr.write("An input file name ({}) contains a character ({}) that is not permitted in the input file names for this pipeline\n".format(filename, char))
                sys.exit(1)


def main(target_assembly_fasta_path, target_assembly_gff_path, liftoff_reference_fasta_path, liftoff_reference_gff_path, ref_mitoch_fasta_path, ref_apicoplast_fasta_path, rna_seq_fq1_path, rna_seq_fq2_path, orthomcl_ref_folder, chunk_size):
    
    input_filenames_list = [target_assembly_fasta_path, target_assembly_gff_path, liftoff_reference_fasta_path, liftoff_reference_gff_path, ref_mitoch_fasta_path, ref_apicoplast_fasta_path, rna_seq_fq1_path, rna_seq_fq2_path]
    input_filenames_list = [os.path.abspath(n) for n in input_filenames_list]
    check_filenames_for_bad_characters(input_filenames_list)
    
    script_folder = os.path.dirname(os.path.abspath(__file__))
    validatefastq_path = script_folder + "/third_party_files/validatefastq-assembly-0.1.1.jar"
    
    fc = FastaChecker(target_assembly_fasta_path, "nucleotide", max_sequences=None, min_total_seq_len=chunk_size, min_longest_seq_len=chunk_size)
    fc.process_fasta()

    if target_assembly_gff_path != "NA":
        gff_c = GFF_checker(target_assembly_gff_path, target_assembly_fasta_path)
        gff_c.process_gff()

    if liftoff_reference_fasta_path != "NA" and liftoff_reference_gff_path != "NA":
        fc = FastaChecker(liftoff_reference_fasta_path, "nucleotide")
        fc.process_fasta()
        gff_c = GFF_checker(liftoff_reference_gff_path, liftoff_reference_fasta_path)
        gff_c.process_gff()

    if ref_mitoch_fasta_path != "NA":
        fc = FastaChecker(ref_mitoch_fasta_path, "nucleotide", max_sequences=1)
        fc.process_fasta()

    if ref_apicoplast_fasta_path != "NA":
        fc = FastaChecker(ref_apicoplast_fasta_path, "nucleotide", max_sequences=1)
        fc.process_fasta()
    
    if rna_seq_fq1_path != "NA" and rna_seq_fq2_path != "NA":
        validate_fastq_files(validatefastq_path, rna_seq_fq1_path, rna_seq_fq2_path)

    if orthomcl_ref_folder != "NA":
        validate_orthomcl_references_folder(orthomcl_ref_folder)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_assembly_fasta_path", type=str, help="Path to the FASTA file of the target assembly")
    parser.add_argument("--target_assembly_gff_path", type=str, help="Path to the GFF3 annotations file for the target assembly (default: NA)", default="NA")
    parser.add_argument("--liftoff_reference_fasta_path", type=str, help="Path to the FASTA file of a reference assembly, for annotation transfer using Liftoff", default="NA")
    parser.add_argument("--liftoff_reference_gff_path", type=str, help="Path to the GFF3 annotations file of the reference assembly, for annotation transfer using Liftoff (default: NA)", default="NA")
    parser.add_argument("--ref_mitoch_fasta_path", type=str, help="Path to FASTA file with reference mitochondrion sequence (default: NA)", default="NA")
    parser.add_argument("--ref_apicoplast_fasta_path", type=str, help="Path to FASTA file with reference apicoplast sequence (default: NA)", default="NA")
    parser.add_argument("--rna_seq_fq1_path", type=str, help="Path to the fastq.gz file 1 of RNA-Seq reads (default: NA)", default="NA")
    parser.add_argument("--rna_seq_fq2_path", type=str, help="Path to the fastq.gz file 2 of RNA-Seq reads (default: NA)", default="NA")
    parser.add_argument("--orthomcl_ref_folder", type=str, help="Path to the folder with the reference proteomes for OrthoMCL (default: NA)", default="NA")
    parser.add_argument("--chunk_size", type=int, help="Sliding window step size (bp) for running GDA (default: 5000)", default=5000)
    args = parser.parse_args()
    main(args.target_assembly_fasta_path, args.target_assembly_gff_path, args.liftoff_reference_fasta_path, args.liftoff_reference_gff_path, args.ref_mitoch_fasta_path, args.ref_apicoplast_fasta_path, args.rna_seq_fq1_path, args.rna_seq_fq2_path, args.orthomcl_ref_folder, args.chunk_size)







