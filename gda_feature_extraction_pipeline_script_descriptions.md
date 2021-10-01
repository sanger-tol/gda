# Descriptions of scripts that are components of the genomic feature extraction pipeline
Most of these scripts are are triggered by the Nextflow master script for the extraction of genomic features. The users are not required to interact with these scripts directly but knowing what these scripts do might help for troubleshooting or for non-standard use cases.
<br/>
<br/>
## Main ##
**autodownsample_merged_tsv.py**<br/>
Script for automatically triggering the downsampling of the TSV table that has been derived by merging bedgraph files<br/>
<br/>
**autorun_clustering.py**<br/>
Script for automatically running the gda_parameters.py and gda_clustering.py scripts after running the genomic feature extraction<br/>
<br/>
**cluster_junctions_fisher_test.py**<br/>
Script for counting cluster junctions between windows in the GDA output 'clusters.bed' file.<br/>
Fisher test is used to determine if some types of cluster junctions occur more rarely or often that expected by chance<br/>
<br/>
**convert_einverted_output_to_gff.py**<br/>
Script for converting the output of EMBOSS einverted to GFF3 format<br/>
<br/>
**count_neighbouring_clusters_in_bed_file.py**<br/>
Script for counting neighboring clusters of each cluster in GDA output BED file<br/>
<br/>
**decomposition_gc_skew_repeats_sliding_window.py**<br/>
Script for finding sequence GC, GC skew, stop codon frequency and frequency of some repeat motifs in a FASTA file with sliding window<br/>
<br/>
**dustmasker_get_masked_seq_percentage_bedgraph.py**<br/>
Script for finding the percentages of nucleotides that got masked by Dustmasker in a FASTA file, as a sliding window<br/>
<br/>
**extract_gene_stats.py**<br/>
Script for extracting gene statistics (length, exon count, intron count, strand) from a GFF file<br/>
<br/>
**extract_gff_features.py**<br/>
Script for extracting mRNA, pseudogene, tRNA and rRNA features from GFF3 file to convert them into bedgraph format<br/>
<br/>
**filtered_ectopic_organellar_blast_hits_to_gff.py**<br/>
Script for filtering the BLAST results to detect ectopic mitochondrial and apicoplast sequences. The script outputs the ectopic mitochondrion and apicoplast BLAST hits as GFF<br/>
<br/>
**gaps_to_bedgraph.py**<br/>
Script for detecting Ns in a FASTA file with sliding window and reporting them in bedgraph format<br/>
<br/>
**gc_skew_etc_to_bedgraph.py**<br/>
Script for converting the output of decomposition_gc_skew_repeats_sliding_window.py to bedgraph format<br/>
<br/>
**gda_check_software_dependencies.py**<br/>
Script for checking if required software for the GDA feature extraction pipeline is in path<br/>
<br/>
**general_purpose_functions.py**<br/>
General purpose functions<br/>
File for functions that can be reused in many Python scripts<br/>
<br/>
**genome_decomp_pipeline_shared_functions.py**<br/>
File for functions that are shared between scripts of the decomposition pipeline<br/>
<br/>
**get_fasta_sequence_lengths.py**<br/>
Script for extracting sequence names and lengths from a FASTA file<br/>
Argument1: path to FASTA file<br/>
Output: csv table printed to STDOUT. First column: sequence names. Second column: sequence lengths.<br/>
<br/>
**gff_features_to_bedgraph.py**<br/>
Script for converting GFF features to bedgraph format. Output (STDOUT): a bedgraph file with the fractions of regions that contain the query feature in fixed length chunks of scaffolds<br/>
<br/>
**kmer_freq_sliding_window.py**<br/>
Script for counting kmer frequencies in a FASTA file using a sliding window<br/>
Output: bedgraph files for the counts of kmers in every sliding window step across the scaffolds in the input FASTA file<br/>
<br/>
**map_rna-seq_reads_and_get_coverage.py**<br/>
Script for getting RNA-Seq read coverage for decomposition analysis of genomes<br/>
<br/>
**quick_test_feature_extraction_pipeline.py**<br/>
Script for running a quick test of the genomic feature extraction pipeline of GDA with P. falciparum chromosome 1 as the input.<br/>
This script runs all the mandatory parts of the pipeline but skips most of the optional parts<br/>
<br/>
**run_blast_to_detect_ectopic_organellar_seq.py**<br/>
Script for detection of ectopic organellar sequences using BLAST<br/>
<br/>
**run_dustmasker.py**<br/>
Script for running DustMasker to detect low complexity regions in an assembly<br/>
<br/>
**run_einverted.py**<br/>
Script for running einverted to detect inverted repeats<br/>
<br/>
**run_ltrharvest_and_ltrdigest.py**<br/>
Script for running LTRharvest and LTRdigest<br/>
<br/>
**run_ltrharvest_and_ltrdigest_with_genome_chunk.py**<br/>
Script for running LTRharvest and LTRdigest<br/>
<br/>
**run_ltrharvest_ltrdigest_with_genome_chunks.py**<br/>
Script for splitting a genome assembly into chunks with fixed size and running LTRharvest + LTRdigest with the chunks<br/>
<br/>
**run_red_meshclust2.py**<br/>
Script for running Red and MeShClust2 to detect repeat families<br/>
<br/>
**run_wgsim.py**<br/>
Script for running wgsim to generate simulated reads, mapping these reads and finding their coverage<br/>
<br/>
**sam_to_sorted_indexed_bam.py**<br/>
Script conversion of .sam file with mapped reads to sorted and indexed .bam file<br/>
<br/>
**samtools_depth_to_bedgraph.py**<br/>
Script for converting coverage data (based on SAMtools depth) to bedgraph format. Output (STDOUT): a bedgraph file with mean coverage of fixed length chunks of scaffolds<br/>
<br/>
**shorten_fasta_headers.py**<br/>
Script for shortening FASTA headers, by splitting the header and keeping only the first element<br/>
<br/>
**split_assembly.py**<br/>
Script for splitting a genome assembly into chunks with fixed size<br/>
<br/>
**split_multifasta_with_fixed_nr_of_sequences.py**<br/>
Script for dividing a multiFASTA file into smaller FASTA files with a fixed number of sequences<br/>
Argument1: path to input FASTA file<br/>
Argument2: path for output folder<br/>
Argument3: number of sequences per output file<br/>
Output: original FASTA file split into individual files<br/>
<br/>
**stats_per_gene_to_bedgraph.py**<br/>
Script for converting the table of stats per each gene to bedgraph<br/>
<br/>
**validate_input_files.py**<br/>
Script for validating the input files of the GDA feature extraction pipeline<br/>
<br/>
**validate_nextflow_config.py**<br/>
Script for validating the nextflow.config file of GDA<br/>
<br/>
**validate_pipeline_run_folder.py**<br/>
Script for validating the GDA pipeline run folder before running the pipeline<br/>
<br/>
<br/>
## Genome annotation ##
**add_missing_ids_to_gff3.py**<br/>
Script for processing a GFF file to add missing IDs. The input is a GFF that has been created by combining Augustus, Barrnap and tRNAscan-SE<br/>
<br/>
**combine_annotation_gff_files.py**<br/>
Script for combining the output GFF3 files of Augustus, Barrnap and tRNAscan into one GFF3 file<br/>
<br/>
**convert_trnascan_bed_file_to_gff.py**<br/>
Script for converting tRNAscan output BED file to GFF3 format<br/>
<br/>
**gda_annotate_genes.py**<br/>
Master script for running gene annotation scripts for GDA<br/>
<br/>
**gff_to_transcripts_and_proteins.py**<br/>
Script for extracting transcripts and protein sequences from a GFF3 and a genome assembly FASTA file<br/>
<br/>
**run_augustus.py**<br/>
Script for running Augustus for genome annotation as multiple parallel jobs<br/>
This script uses snippets of code adapted from https://github.com/stephenrdoyle/generic_scripts/blob/master/random_workflows/run_augustus_split_by_contigs.sh<br/>
<br/>
**run_barrnap.py**<br/>
Script for running Barrnap for detecting rRNAs<br/>
<br/>
**run_liftoff.py**<br/>
Script for running Liftoff to transfer gene annotations<br/>
<br/>
**run_trnascan.py**<br/>
Script for running tRNAscan-SE for detecting tRNAs<br/>
<br/>
<br/>
## OrthoMCL ##
**generate_orthomcl_gg_file_from_fasta.py**<br/>
Script for generating gg_file for OrthoMCL from protein FASTA files<br/>
Argument1: path to a CSV file. First column: species identifiers (short names). Second column: names of FASTA files for each species (without folder path)<br/>
Output: gg_file for OrthoMCL<br/>
Argument2: path to folder with protein FASTA files<br/>
<br/>
**orthomcl_batch.py**<br/>
Script for running OrthoMCL as batch<br/>
<br/>
**orthomcl_conservation.py**<br/>
Script for converting OrthoMCL results into a table of paralog counts, ortholog counts and conservation ratio<br/>
<br/>
**remove_non_mrna_cds_features.py**<br/>
Script for processing a GFF3 file to remove CDS features whose parent feature is something other than 'mRNA'<br/>
<br/>
**run_orthomcl.py**<br/>
Script for running OrthoMCL (including Diamond blastp for OrthoMCL)<br/>
<br/>
<br/>
## RepeatModeler ##
**condense_simple_repeat_sequences.py**<br/>
Script for condensing a list of simple repeat sequences to remove redundant sequences.<br/>
For example: TAA and TTA are the same sequence, one is the reverse complement of the other. TAATAA repeat is the same as TAA repeat. AATAAT is the same as TAATAA but with a shifted starting point<br/>
Input: GFF with repeat locations from RepeatMasker, processed with process_repeatmasker_gffs.py to extract only simple repeat sequences<br/>
Output: simple repeats GFF with redundant sequences collapsed into one sequence<br/>
<br/>
**find_repeats_enriched_at_scaff_edges.R**<br/>
Script for using RepeatModeler's repeat families output for detecting repeats that are enriched at scaffold edges<br/>
This is not a component of the main pipeline of GDA but can be used as an extra step to get more information out of the data<br/>
<br/>
**process_repeatmasker_gffs.py**<br/>
Script for running scripts that process RepeatMasker gff files<br/>
<br/>
**reformat_repeatmasker_gff.py**<br/>
Script for splitting the simple and complex repeat lines in RepeatMasker GFF output and reformatting the GFF so that it can be used as the input for the multiple_gff_features_to_bedgraph.py script<br/>
<br/>
**repeatmasker_gff_to_bedgraph.py**<br/>
Script for converting RepeatMasker repeat coordinates from GFF to bedgraph<br/>
<br/>
**repeatmasker_simple_repeat_frequencies.py**<br/>
Script for finding simple repeat frequencies in in GFF derived from the output of RepeatModeler + RepeatMasker<br/>
<br/>
**run_repeatmasker_repeatmodeler.py**<br/>
Script for running RepeatMasker and RepeatModeler as a part of genome decomposition<br/>
<br/>
**sum_simple_or_complex_repeat_tracks.py**<br/>
Script for making bedgraph tracks that are the sum of all simple repeat tracks or sum of all complex repeat tracks<br/>
<br/>
<br/>
## Tandem Repeats Finder ##
**run_trf.py**<br/>
Script for running Tandem Repeats Finder as a part of genome decomposition<br/>
<br/>
**trf_repeat_density_sliding_window.py**<br/>
Script for finding repeat density in a genome using sliding window on genome FASTA file where repeats have been masked with Tandem Repeats Finder<br/>
Output: tab separated table. Column1: scaffold name. Column2: chunk start coordinate in the scaffold (1-based). Column3: chunk end coordinate in the scaffold.<br/>
    Column4: fraction of nucleotides in the chunk that were masked by TandemRepeatsFinder. Column5: True if the fraction of masked nucleotides exceeds a cutoff, False if not.<br/>
<br/>
**trf_repeat_density_to_bedgraph.py**<br/>
Script for converting Tandem Repeats Finder repeat density values to bedgraph format<br/>
Argument1: path to an output file of trf_repeat_density_sliding_window.py<br/>
Output (STDOUT): input file converted to bedgraph format<br/>
<br/>
**trf_repeat_density_to_gff.py**<br/>
Script for converting Tandem Repeats Finder repeat density values (that have been divided into repeat-rich and repeat-poor regions) into GFF or BED format<br/>
<br/>
