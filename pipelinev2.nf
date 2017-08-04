#!/usr/bin/env nextflow

/*if( params.help ) {
    help()
    return
}*/

if( !nextflow.version.matches('0.25+') ) {
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    exit 1
}

if( params.index ) { 
    index = Channel.fromPath(params.index).toSortedList() 
    if( !index.exists() ) exit 1, "Genome index files could not be found: ${params.index}"    
}

if( params.reference ) {     
    reference = file(params.reference)
    if( !reference.exists() ) exit 1, "Reference genome file could not be found:${params.reference}"          
}

threads = params.threads
/*
 * Aligns fastq files to hg38
 */

Channel
    .fromPath( "${params.input_path}/*fastq")
    .ifEmpty { exit 1, "Fastq files could not be found in: ${params.input_path}" }
    .set { gene_trap_insertions }

/*
 * Aligns fastqs to genome via bwa
 */

process AlignToGenome {

    publishDir "${params.output_dir}/Alignment", mode: "copy"

    input:
    file fastq from gene_trap_insertions
    file idx from index.first()

    output:
    file alignment into aligned_sams

    """
    bwa mem index fastq -t ${threads} > sam
    samtools view -bS sam | samtools sort -o alignment
    """
}

/*
 * Counts the gene trap insertions in genes, and compiles metrics
 */

process MakeInsertionTables {
    input:
    file bed from aligned_sams

    output:
    file insertion_table into insertion_table_csvs

    """
    intron_exon_counts.py -i $params.intron_bed -e $params.exon_bed -b bed -o insertion_table
    """
}

