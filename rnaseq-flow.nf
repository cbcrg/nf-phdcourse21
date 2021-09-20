/* 
 * include required tasks
 */
include { index; quantification; fastqc; multiqc  } from './rnaseq-modules.nf'

/* 
 * define the data analysis workflow 
 */
workflow rnaseq_flow {
    // required inputs
    take:
    transcriptome
    read_files

    // workflow implementation
    main:
    index(transcriptome)
    quantification(index.out, read_files)
    fastqc(read_files)
    multiqc( quantification.out.mix(fastqc.out).collect() )

    // Optionally, outputs can be named using the emit option
    // emit:
    // multiqc_report = multiqc.out.multiqc_report
}