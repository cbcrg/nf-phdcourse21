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

/*
 * Run fastQC to check quality of reads files
 */
process fastqc {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}

/*
 * Create a report using multiQC for the quantification
 * and fastqc processes
 */
process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    path('*')

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """
}

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )

workflow {
    index( params.transcript )
    quantification( index.out, read_pairs_ch )
    fastqc( read_pairs_ch )
    multiqc( quantification.out.mix(fastqc.out).collect() )
}

workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
