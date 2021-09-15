/* 
 * This code enables the new dsl of Nextflow. 
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


/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {
    
    input:
    path transcriptome
    
    output:
    path 'index'

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

/*
 * A workflow consists of a number of invocations of processes which are fed
 * with the expected input channels as if they were custom functions. 
 * You can only invoke once a process per workflow.
 */

workflow {
    result = index(params.transcript)
}