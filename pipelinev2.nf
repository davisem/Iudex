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

threads = params.threads
index = Channel.fromPath(params.index).toSortedList()
intron_bed = file(params.intron_bed)
exon_bed = file(params.exon_bed)

/*
 * Aligns fastq files to hg38
 */

Channel
    .fromPath( "${params.input_path}/*fastq")
    .ifEmpty { exit 1, "Fastq files could not be found in: ${params.input_path}" }
    .set { gene_trap_insertions }


process AlignToGenome {

    publishDir "${params.output_dir}/Alignment", mode: "copy"

    input:
    file fastq from gene_trap_insertions
    file idx from index

    output:
    file "${fastq.baseName}.aln" into aligned_bams

    """
    bwa mem hg38.fa ${fastq} -t ${threads} > sam
    samtools view -bS sam | samtools sort -o "${fastq.baseName}.aln"
    """
}

/*
 * Counts the gene trap insertions in genes, and compiles metrics
 */

process MakeInsertionTables {
    input:
    
    file bam from aligned_bams
    file intron_bed
    file exon_bed


    output:
    file "${bam.baseName}.table" into insertion_table_csvs

    """

    python /src/intron_exon_counts.py -i ${intron_bed} -e ${exon_bed} -b ${bam} -o "${bam.baseName}.table"
    """
}
