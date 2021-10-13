nextflow.enable.dsl=2

/* 
 * include required modules
 */
include { index; quantification; fastqc; multiqc  } from './rnaseq-modules.nf'

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcript = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcript}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )

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

    // Optionally, outputs can be named using the emit option see https://www.nextflow.io/docs/latest/dsl2.html#process-named-output
    // emit:
    // multiqc_report = multiqc.out.multiqc_report
}

/* 
 * Invocation of the data analysis workflow 
 */
workflow {
    rnaseq_flow( params.transcript, read_pairs_ch )
}