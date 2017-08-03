#!/usr/bin/env nextflow
params.help = false
params.threads = 1

params.genome_index = "$baseDir/venv/hg38"
params.input_path = "$PWD"
params.output_dir = "$PWD"
params.intron_bed = ''
params.exon_bed = ''


/*
 * Aligns fastq files to hg19
 */

Channel
    .fromPath( "${params.input_path}/*fastq")
    .ifEmpty { exit 1, "Fastq files could not be found in: ${params.input_path}" }
    .set { gene_trap_insertions }


process AlignToGenome {

    publishDir "${params.output_dir}/Alignment", mode: "copy"

    input:
    file reads from gene_trap_insertions

    output:
    file alignment.bam into aligned_sams

    """
    bowtie2 -x ${params.genome_index} -p ${params.threads} -U $reads -S alignment.sam
    samtools view -bS alignment.sam | samtools sort -@ ${params.threads} -o alignment.bam
    """
}


process MakeInsertionTables {
    input:
    file bed from aligned_sams

    output:
    file insertion_table into insertion_table_csvs

    """
    intron_exon_counts.py -i $params.intron_bed -e $params.exon_bed -b bed -o insertion_table
    """
}

