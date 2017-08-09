#!/usr/bin/env nextflow


/*
vim: syntax=groovy
-*- mode: groovy;-*-
================================================================================
=             S C R E E N I N G | H A P L O I D | P I P E L I N E              =
================================================================================
@Author
Eric Davis <emdavis48@gmail.com>
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/davisem/iudex
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/davisem/iudex/blob/master/README.md
--------------------------------------------------------------------------------
 @Licence
 https://github.com/davisem/iudex/blob/master/LICENSE
--------------------------------------------------------------------------------
 Processes overview
 - Run the Iudex analysis pipeline
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

if ( params.help ) {
    help()
    return
}

if( !nextflow.version.matches('0.25+') ) {
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    exit 1
}

if ( params.index ) {
    index = Channel.fromPath(params.index).toSortedList()
}

threads = params.threads
intron_bed = file(params.intron_bed)
exon_bed = file(params.exon_bed)
false_positive_probability = params.false_positive_probability

Channel
    .fromPath( "${params.input_path}/*fastq" )
    .ifEmpty { exit 1, "Fastq files could not be found in: ${params.input_path}" }
    .set { gene_trap_insertions }

/*
 * Pre-filters duplicate reads from the input fastqs.
 */

process FilterFastq {

    publishDir "${params.output_dir}/FilterFastq", mode: "copy"
    
    input:
    file fastq from gene_trap_insertions

    output:
    file "${fastq.baseName}.filtered" into filtered_fastqs

    """
    /src/fastq_filterer ${fastq} "${fastq.baseName}.filtered" ${false_positive_probability}
    """

}

/*
 * Aligns fastq files to hg38
 */

process AlignToGenome {

    publishDir "${params.output_dir}/Alignment", mode: "copy"

    input:
    file fastq from filtered_fastqs
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

    publishDir "${params.output_dir}/Insertions", mode: "copy"

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
