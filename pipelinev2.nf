#!/usr/bin/env nextflow

params.ref_genome = '/Users/ericdavis/bin/bowtie-0.12.8/indexes/hg19.ebwt/hg19'
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
    file alignment into aligned_sams

    """
    bowtie -m1 -v1 --best -3 14 $params.ref_genome $reads -S alignment
    """
}

process samtToBam {
    
    input:
    file alignment from aligned_sams

    output:
    file output into aligned_bams

    """
    samtools view -bS -o output $alignment
    """
}

process bamToBed {
    
    input: 
    file bam from aligned_bams

    output:
    file bed into beds

    """
    bedtools bamtobed -i $bam > bed
    """
}

process MakeInsertionTables {
    input:
    file bed from beds

    output:
    file insertion_table into insertion_table_csvs

    """
    intron_exon_counts.py -i $params.intron_bed -e $params.exon_bed -b bed -o insertion_table
    """
}