/* 
 * Example showing pipeline modularizaion 
 * Using Nextfloq DSL2
 */
nextflow.enable.dsl=2

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

include { rnaseq_flow } from './rnaseq-flow-implicit.nf'

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )

workflow {
    rnaseq_flow( params.transcript, read_pairs_ch )
}

// You can also declare again here implicitely the workflow
// workflow rnaseq_implicit {
//     rnaseq_flow( params.transcript, read_pairs_ch )
// }

// workflow {
//     rnaseq_implicit()
// }