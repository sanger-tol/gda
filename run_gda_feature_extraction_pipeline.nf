#!/usr/bin/env nextflow
// Genome Decomposition Analysis pipeline feature extraction master script 

params.pipeline_output_folder = params.pipeline_run_folder + "/bedgraph_output"
params.error_stream_logs_folder = params.pipeline_run_folder + "/error_stream_logs"
params.chunk_size = params.chunk_size.toString()
gene_annotation_folder = params.pipeline_run_folder + "/gene_annotation"
gc_skew_folder = params.pipeline_run_folder + "/gc_skew_folder"
dustmasker_folder = params.pipeline_run_folder + "/dustmasker"
ectopic_organellar_seq_folder = params.pipeline_run_folder + "/ectopic_organellar_seq"
tandem_repeats_finder_folder = params.pipeline_run_folder + "/tandem_repeats_finder"
rna_seq_mapping_folder = params.pipeline_run_folder + "/rna_seq_mapping"
ltr_folder = params.pipeline_run_folder + "/ltrharvest_ltrdigest"
einverted_folder = params.pipeline_run_folder + "/einverted"
orthomcl_folder = params.pipeline_run_folder + "/orthomcl"
wgsim_folder = params.pipeline_run_folder + "/wgsim"
repeatmodeler_folder = params.pipeline_run_folder + "/repeatmodeler"
meshclust2_folder = params.pipeline_run_folder + "/red_meshclust2"
merged_bedgraph_folder = params.pipeline_run_folder + "/merged_bedgraph_table"
clustering_folder = params.pipeline_run_folder + "/gda_out"

raw_fasta_path = params.assembly_fasta_path
raw_fasta_filename = raw_fasta_path.split("/")[-1]

raw_fasta_basename = raw_fasta_filename.split("\\.")[0]

cleaned_fasta_folder = params.pipeline_run_folder + "/fasta"

assembly_path = cleaned_fasta_folder + "/" + raw_fasta_basename + ".fasta"

fasta_filename = assembly_path.split("/")[-1]
fasta_basename = fasta_filename.split("\\.")[0]

annotations_gff_path = params.gff_path
if (params.run_gene_annotation_pipeline == true) {
    annotations_gff_path = params.pipeline_run_folder + "/gene_annotation/" + fasta_basename + "_annotations.gff3"
}


process check_prerequisite_software {
    echo true
    cpus 1
    output:
    val true into prerequisite_software_check_done_ch
    script:
    """
    gda_check_software_dependencies.py
    """
}


process validate_nextflow_config {
    echo true
    cpus 1
    input:
    val prerequisite_software_check_done_flag from prerequisite_software_check_done_ch
    output:
    val true into config_validation_done_ch
    def config_file_path = params.pipeline_run_folder + "/nextflow.config"
    script:
    """
    validate_nextflow_config.py ${config_file_path}
    """
}


process set_up_pipeline_run_folder {
    cache 'deep'
    echo true
    cpus 1
    input:
    val config_validation_check_done_flag from config_validation_done_ch
    file raw_fasta_infile from file(raw_fasta_path)
    def shorten_fasta_headers_command = raw_fasta_path + " > " + assembly_path + " 2> " + params.error_stream_logs_folder + "/shorten_fasta_headers_errors.txt"
    output:
    val true into pipeline_run_folder_setup_done_ch
    script:
    """
    mkdir -p ${params.pipeline_run_folder}
    validate_pipeline_run_folder.py ${params.pipeline_run_folder}
    mkdir -p ${params.pipeline_output_folder}
    mkdir -p ${params.error_stream_logs_folder}
    mkdir -p ${cleaned_fasta_folder}
    shorten_fasta_headers.py ${shorten_fasta_headers_command} 
    """
}


process validate_input_files {
    cache 'deep'
    echo true
    cpus 1
    input:
    val pipeline_run_folder_setup_done_flag from pipeline_run_folder_setup_done_ch
    file assembly_infile from file(assembly_path)
    def gff_infile = params.gff_path == "NA" ? "NA" : file(params.gff_path)
    def liftoff_reference_fasta_file = params.reference_assembly_path == "NA" ? "NA" : file(params.reference_assembly_path)
    def liftoff_reference_gff_file = params.reference_gff_path == "NA" ? "NA" : file(params.reference_gff_path)
    def ref_mitoch_fasta_file = params.ref_mitoch_fasta_path == "NA" ? "NA" : file(params.ref_mitoch_fasta_path)
    def ref_apicoplast_fasta_file = params.ref_apicoplast_fasta_path == "NA" ? "NA" : file(params.ref_apicoplast_fasta_path)
    def rna_seq_fastq1_file = params.rna_seq_fastq_1_path == "NA" ? "NA" : file(params.rna_seq_fastq_1_path)
    def rna_seq_fastq2_file = params.rna_seq_fastq_2_path == "NA" ? "NA" : file(params.rna_seq_fastq_2_path)
    output:
    val true into input_file_validation_done_ch
    script:
    """
    validate_input_files.py ${assembly_path} --target_assembly_gff_path ${params.gff_path} --liftoff_reference_fasta_path ${params.reference_assembly_path} --liftoff_reference_gff_path ${params.reference_gff_path} --ref_mitoch_fasta_path ${params.ref_mitoch_fasta_path} --ref_apicoplast_fasta_path ${params.ref_apicoplast_fasta_path} --rna_seq_fq1_path ${params.rna_seq_fastq_1_path} --rna_seq_fq2_path ${params.rna_seq_fastq_2_path} --orthomcl_ref_folder ${params.orthomcl_references_folder} --chunk_size ${params.chunk_size} > ${params.error_stream_logs_folder}/validation_of_input_files_stdout.txt 2> ${params.error_stream_logs_folder}/validation_of_input_files_stderr.txt
    """
}


process annotate_genes {
    cache 'deep'
    echo true
    cpus params.threads
    input:
    file assembly_infile from file(assembly_path)
    def liftoff_reference_fasta_file = params.reference_assembly_path == "NA" ? "NA" : file(params.reference_assembly_path)
    def liftoff_reference_gff_file = params.reference_gff_path == "NA" ? "NA" : file(params.reference_gff_path)
    val input_files_validated_flag from input_file_validation_done_ch
    output:
    val true into gene_annotations_done_ch
    script:
    if (params.run_gene_annotation_pipeline == true) {        
        """
        gda_annotate_genes.py ${assembly_path} ${params.pipeline_run_folder} ${params.pipeline_output_folder} ${params.augustus_species} ${params.annotation_target_species_id} --augustus_chunk_overlap ${params.augustus_chunk_overlap} --augustus_chunk_size ${params.augustus_chunk_size} --reference_assembly_path ${params.reference_assembly_path} --reference_gff_path ${params.reference_gff_path} --barrnap_kingdom ${params.barrnap_kingdom} --threads ${params.threads} --chunk_size ${params.chunk_size} 2> ${params.error_stream_logs_folder}/gene_annotation_stderr.txt
        """
    } else {
        """
        echo 'Skipping the gene annotation part of the pipeline'
        """ 
    }

}


process gc_skew_etc {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into gc_skew_etc_done_ch
    """    
    mkdir -p ${gc_skew_folder}
    decomposition_gc_skew_repeats_sliding_window.py ${assembly_path} ${gc_skew_folder} ${params.pipeline_output_folder} --chunk_size ${params.chunk_size} --telomeric_seq_preset ${params.telomeric_seq_preset} 2> ${params.error_stream_logs_folder}/gc_skew_repeats_etc_stderr.txt
    """
}


process kmer_frequencies {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into kmer_frequencies_done_ch
    script:
    def kmers_size3_bedgraph_path = params.pipeline_output_folder + "/" + raw_fasta_basename + "_kmer_size3_skew.bedgraph"
    def kmers_size4_bedgraph_path = params.pipeline_output_folder + "/" + raw_fasta_basename + "_kmer_size4_skew.bedgraph"
    """
    kmer_freq_sliding_window.py ${assembly_path} --chunk_size ${params.chunk_size} --kmer_size 3 > ${kmers_size3_bedgraph_path} 2> ${params.error_stream_logs_folder}/kmer_frequencies_size3_stderr.txt
    kmer_freq_sliding_window.py ${assembly_path} --chunk_size ${params.chunk_size} --kmer_size 4 > ${kmers_size4_bedgraph_path} 2> ${params.error_stream_logs_folder}/kmer_frequencies_size4_stderr.txt
    """
}


process dustmasker_low_complexity_percentage {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into dustmasker_done_ch
    """
    mkdir -p ${dustmasker_folder}
    run_dustmasker.py ${assembly_path} ${dustmasker_folder} ${params.pipeline_output_folder} --chunk_size ${params.chunk_size} 2> ${params.error_stream_logs_folder}/dustmasker_stderr.txt
    """
}


process ectopic_organellar_sequences {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    def ref_mitoch_fasta_file = params.ref_mitoch_fasta_path == "NA" ? "NA" : file(params.ref_mitoch_fasta_path)
    def ref_apicoplast_fasta_file = params.ref_apicoplast_fasta_path == "NA" ? "NA" : file(params.ref_apicoplast_fasta_path)
    output:
    val true into ectopic_organellar_seq_done_ch
    """
    mkdir -p ${ectopic_organellar_seq_folder}
    run_blast_to_detect_ectopic_organellar_seq.py ${assembly_path} ${params.ref_mitoch_fasta_path} ${ectopic_organellar_seq_folder} ${params.pipeline_output_folder} --ref_apicoplast_fasta_path ${params.ref_apicoplast_fasta_path} --chunk_size ${params.chunk_size} 2> ${params.error_stream_logs_folder}/ectopic_organellar_seq_stderr.txt
    """
}


process tandem_repeats_finder {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into tandem_repeats_finder_done_ch
    """
    mkdir -p ${tandem_repeats_finder_folder}
    run_trf.py ${assembly_path} ${tandem_repeats_finder_folder} ${params.pipeline_output_folder} --chunk_size ${params.chunk_size} > ${params.error_stream_logs_folder}/tandem_repeats_finder_stdout.txt 2> ${params.error_stream_logs_folder}/tandem_repeats_finder_stderr.txt 
    """
}


process process_gff_annotations {
    cache 'deep'
    echo true
    cpus 1
    def denovo_annotation_path = gene_annotation_folder + "/" + raw_fasta_basename + "_annotations.gff3"
    def gene_stats_csv_path = gene_annotation_folder + "/" + fasta_basename + "_gene_stats.csv"
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    val gene_annotations_done_flag from gene_annotations_done_ch
    file assembly_infile from file(assembly_path)
    file gene_stats_csv_file from file(gene_stats_csv_path)
    def gff_infile = params.gff_path == "NA" ? file(denovo_annotation_path) : file(params.gff_path)
    output:
    val true into processing_gff_annotations_done_ch 
    script:
    def gff_feature_extraction_command = annotations_gff_path + " " + assembly_path + " " + params.pipeline_output_folder + " --chunk_size " + params.chunk_size
    
    def stats_per_gene_to_bedgraph_command = assembly_path + " " + gene_stats_csv_path + " gene_properties " + params.pipeline_output_folder + " " + params.assembly_title + " --chunk_size " + params.chunk_size
    if (params.custom_gff_tags != "" && params.custom_gff_tags != "NA") {
        gff_feature_extraction_command += " --custom_gff_tags " + params.custom_gff_tags
    }
    if (params.include_gene_strand_bias == true) {
        stats_per_gene_to_bedgraph_command += " --include_gene_strand_bias"
    }
    if ((params.run_gene_annotation_pipeline == false && annotations_gff_path != "NA") || (params.run_gene_annotation_pipeline == true)) {
        """
        mkdir -p ${gene_annotation_folder}
        extract_gff_features.py ${gff_feature_extraction_command} 2> ${params.error_stream_logs_folder}/gff_annotations_stderr.txt
        extract_gene_stats.py ${annotations_gff_path} > ${gene_stats_csv_path} 2> ${params.error_stream_logs_folder}/gene_stats_extraction_stderr.txt
        stats_per_gene_to_bedgraph.py ${stats_per_gene_to_bedgraph_command} 2> ${params.error_stream_logs_folder}/gene_stats_to_bedgraph_stderr.txt
        """
    } else {
        """
        echo 'Skipping the processing of GFF annotations'
        """
    }    
}


process rna_seq_coverage {
    cache 'deep'
    echo true
    cpus params.threads
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    def rna_seq_fastq1_file = params.rna_seq_fastq_1_path == "NA" ? "NA" : file(params.rna_seq_fastq_1_path)
    def rna_seq_fastq2_file = params.rna_seq_fastq_2_path == "NA" ? "NA" : file(params.rna_seq_fastq_2_path)
    output:
    val true into rna_seq_coverage_done_ch 
    script:
    def rna_seq_mapping_command = assembly_path + " " + params.rna_seq_fastq_1_path + " " + params.rna_seq_fastq_2_path + " " + rna_seq_mapping_folder + " " + params.pipeline_output_folder + " --chunk_size " + params.chunk_size + " --threads " + params.threads
    if (params.rf_stranded == true) {
        rna_seq_mapping_command += " --rf_stranded"
    }
    if (params.rna_seq_fastq_1_path != "NA" && params.rna_seq_fastq_2_path != "NA") {
        """
        mkdir -p ${rna_seq_mapping_folder}
        map_rna-seq_reads_and_get_coverage.py ${rna_seq_mapping_command} > ${params.error_stream_logs_folder}/rna_seq_coverage_stdout.txt 2> ${params.error_stream_logs_folder}/rna_seq_coverage_stderr.txt
        """
    } else {
        """
        echo 'Skipping RNA-Seq coverage'
        """
    }
}


process ltrharvest_and_ltrdigest {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into ltrharvest_ltrdigest_done_ch 
    script:
    def ltrharvest_ltrdigest_command = assembly_path + " " + ltr_folder + " " + params.pipeline_output_folder + " --chunk_size " + params.chunk_size
    def ltrharvest_ltrdigest_stdout_path = ltr_folder + "/" + raw_fasta_basename + "_ltrharvest_ltrdigest_stdout.txt"
    """
    mkdir -p ${ltr_folder}
    mkdir -p ${ltr_folder}/temp_files
    run_ltrharvest_and_ltrdigest.py ${ltrharvest_ltrdigest_command} > ${ltrharvest_ltrdigest_stdout_path} 2> ${params.error_stream_logs_folder}/ltrharvest_ltrdigest_stderr.txt
    """ 
}


process einverted {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into einverted_done_ch
    script:
    def einverted_command = assembly_path + " " + einverted_folder + " " + params.pipeline_output_folder + " --chunk_size " + params.chunk_size
    """
    mkdir -p ${einverted_folder}
    run_einverted.py ${einverted_command} 2> ${params.error_stream_logs_folder}/einverted_stderr.txt
    """ 
}


process orthomcl {
    cache 'deep'
    echo true
    cpus params.threads
    def denovo_annotation_path = gene_annotation_folder + "/" + raw_fasta_basename + "_annotations.gff3"
    def target_proteome_path = orthomcl_folder + "/" + raw_fasta_basename + "_proteome.faa"
    def orthomcl_command = assembly_path + " " + annotations_gff_path + " " + orthomcl_folder + " " + params.pipeline_output_folder + " " + params.annotation_target_species_id + " " + params.orthomcl_references_folder + " --threads " + params.threads + " --memory_limit " + params.diamond_memory_limit + " --chunk_size " + params.chunk_size 
    def orthomcl_stdout_path = orthomcl_folder + "/" + raw_fasta_basename + "_orthomcl_stdout.txt"
    input:
    def gff_infile = params.gff_path == "NA" ? file(denovo_annotation_path) : file(params.gff_path)
    def target_proteome_file = file(target_proteome_path)
    file assembly_infile from file(assembly_path)
    val input_files_validated_flag from input_file_validation_done_ch
    val gene_annotations_done_flag from gene_annotations_done_ch
    output:
    val true into orthomcl_done_ch
    script:
    if (params.orthomcl_references_folder != "NA" && (params.gff_path != "NA" || params.run_gene_annotation_pipeline == true)) {
        if (params.diamond_sensitive == false) {
            orthomcl_command += " --diamond_nonsensitive"
        }        
        """
        mkdir -p ${orthomcl_folder}
        orthomcl_batch.py ${orthomcl_command} > ${orthomcl_stdout_path} 2> ${params.error_stream_logs_folder}/orthomcl_stderr.txt
        """
    } else {
        """
        echo 'Skipping the OrthoMCL part of the pipeline'
        """ 
    }

}



process wgsim_coverage {
    cache 'deep'
    echo true
    cpus params.threads
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into wgsim_done_ch
    script:
    def wgsim_command = assembly_path + " " + wgsim_folder + " " + params.pipeline_output_folder + " --chunk_size " + params.chunk_size + " --threads " + params.threads  + " --target_coverage " + params.wgsim_target_coverage
    """
    mkdir -p ${wgsim_folder}
    run_wgsim.py ${wgsim_command} 2> ${params.error_stream_logs_folder}/wgsim_stderr.txt
    """ 
} 



process detect_repeat_families {
    cache 'deep'
    cpus params.threads
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into repeat_family_detection_done_ch
    script:
    def repeatmodeler_command = raw_fasta_basename + " " + assembly_path + " " + repeatmodeler_folder + " " + params.pipeline_output_folder + " --threads " + params.threads + " --chunk_size " + params.chunk_size
    def red_meshclust2_command = assembly_path + " " + meshclust2_folder + " " + params.pipeline_output_folder + " --chunk_size " + params.chunk_size + " --threads " + params.threads
    if (params.run_repeat_family_detection == true) {
        if (params.repeat_family_detection_engine == "repeatmodeler") {
            """
            mkdir -p ${repeatmodeler_folder}
            run_repeatmasker_repeatmodeler.py ${repeatmodeler_command} 2> ${params.error_stream_logs_folder}/repeatmodeler_repeatmasker_stderr.txt
            """
        } else if (params.repeat_family_detection_engine == "meshclust2") {
            """
            mkdir -p ${meshclust2_folder}
            run_red_meshclust2.py ${red_meshclust2_command} 2> ${params.error_stream_logs_folder}/red_meshclust2_stderr.txt
            """
        }  
    } else {
        """
        echo 'Skipping the detection of repeat families'
        """
    }
     
}


process gaps {
    cache 'deep'
    echo true
    cpus 1
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    file assembly_infile from file(assembly_path)
    output:
    val true into gaps_done_ch
    script:
    def gaps_bedgraph_path = params.pipeline_output_folder + "/" + raw_fasta_basename + "_gaps.bedgraph"
    """
    gaps_to_bedgraph.py ${assembly_path} --chunk_size ${params.chunk_size} > ${gaps_bedgraph_path} 2> ${params.error_stream_logs_folder}/gaps_stderr.txt 
    """
}


process merge_bedgraph_files {
    cache 'deep'
    echo true
    cpus 1
    def bedgraph_tsv_path = merged_bedgraph_folder + "/" + raw_fasta_basename + "_merged_bedgraph.tsv"
    def merge_bedgraph_command = assembly_path + " " + params.pipeline_output_folder + " " + params.chunk_size + " '" + params.assembly_title + "'" + " true " + " > " + bedgraph_tsv_path
    def autodownsample_tsv_command = assembly_path + " " + bedgraph_tsv_path + " --chunk_size " + params.chunk_size
    
    input:
    val input_files_validated_flag from input_file_validation_done_ch
    val gene_annotations_done from gene_annotations_done_ch
    val gc_skew_etc_done_ch
    val dustmasker_flag from dustmasker_done_ch
    val kmers_flag from kmer_frequencies_done_ch
    val ectopic_organellar_seq_flag from ectopic_organellar_seq_done_ch
    val tandem_repeats_finder_flag from tandem_repeats_finder_done_ch
    val gff_annotations_flag from processing_gff_annotations_done_ch
    val rna_seq_coverage_flag from rna_seq_coverage_done_ch
    val ltrharvest_ltrdigest_flag from ltrharvest_ltrdigest_done_ch
    val einverted_flag from einverted_done_ch
    val orthomcl_flag from orthomcl_done_ch
    val wgsim_flag from wgsim_done_ch
    val repeat_families_flag from repeat_family_detection_done_ch
    val gaps_flag from gaps_done_ch
    file assembly_infile from file(assembly_path)
    script:
    """
    mkdir -p ${merged_bedgraph_folder}
    merge_bedgraph_files ${merge_bedgraph_command} 2> ${params.error_stream_logs_folder}/merge_bedgraph_stderr.txt
    autodownsample_merged_tsv.py ${autodownsample_tsv_command} 2> ${params.error_stream_logs_folder}/autodownsample_tsv_stderr.txt
    """
}
