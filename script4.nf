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

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )

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
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process quantification {
    
    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path(pair_id)

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

workflow {
    index( params.transcript )
    quantification( index.out, read_pairs_ch )
}