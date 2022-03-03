#!/usr/bin/env nextflow
/*
========================================================================================
                         assembly-Pipeline
========================================================================================
Nanopore Genome Analysis Pipeline.
#### Homepage / Documentation
#https://github.com/XXXXXXX
#### Authors
Dr. Hiren Ghosh  <hiren.ghosh@gmail.com>
UniKlink Freiburg
Date: June 2021
========================================================================================
*/

def helpMessage() {
  
log.info """

Usage:

The typical command for running the pipeline is as follows:

nextflow run main.mf --reads /path/to/fastq --assembler <name of the assembler> --outdir /path/to/outdir
       
Mandatory arguments:

Assembly workflow:

--reads				reads fastq file
--assembler			name of the assembler
--outdir 			name of the output dir

Optional arguments:

--threads			Number of CPUs to use during the job
--help				This usage statement
-with-report			html report file
--fast5				fast5 path
        """
}

println """\
         Nananopore Assembly - Pipeline
         ===================================
         assembler    : ${params.assembler}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


/*
------------------------------------------------------------------------------
                       C O N F I G U R A T I O N
------------------------------------------------------------------------------
*/

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '1.0dev'

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
// pipeline parameter

params.fast5 = ''
params.reads = ''
params.genome_size = ''

params.output = ''
params.cpus = ''
params.mem = ''

// choose the assembler
params.assembler = ''
if (params.assembler != 'flye' && params.assembler != 'canu' \
    && params.assembler != 'unicycler') {
        exit 1, "--assembler: ${params.assembler}. \
        Assembly parameter missing : Should be 'flye', 'canu' or 'unicycler'"
}
// requires --genome_size for canu
if (params.assembler == 'canu' && params.genomeSize == '20000')
// requires --output
if (params.output == '') {
    exit 1, "--output is a required parameter"
}
//requires --reads
if (params.reads == '') {
    exit 1, "--reads is a required parameter"
}

//requires --output
if (params.outdir == '') {
    exit 1, "--outdir is a required parameter"
}

process adapter_trimming {
    input:
  file(reads) from file(params.reads)

    output:
  file('trimmed.fastq') into trimmed_reads

	script:
        """
    porechop -i "${reads}" -t "${task.cpus}" -o trimmed.fastq
    
    
        """
}
// Trimmed reads will be use by both assembly and concensus
trimmed_reads.into {trimmed_for_assembly; trimmed_for_concensus }

process filter_long {
    input:
  file reads from trimmed_for_assembly

    output:
  file('trimmed.fastq.gz') into filtlong_reads

	script:
        """
    filtlong --min_length 1000 ${reads} | gzip > trimmed.fastq.gz
    
        """
}

process assembly {

publishDir "Assembly_out", mode: 'copy', pattern: 'assembly.fasta'

    input:
  	file reads from filtlong_reads
  
    output:
  	file 'assembly.fasta' into assembly

	script:
	if(params.assembler == 'flye')
      
        """
        mkdir fly_out
    	flye --nano-raw ${reads} --genome-size 1m --out-dir  fly_out
    	mv fly_out/assembly.fasta assembly.fasta
    	
    	"""
    
    else if(params.assembler == 'canu')
       
        """
        canu -p assembly -d canu_out genomeSize=2m -nanopore-raw "${reads}" maxThreads=8 gnuplotTested=true useGrid=false 
        mv canu_out/assembly.contigs.fasta assembly.fasta
       
        """
    else if(params.assembler == 'unicycler')
        """
        unicycler -l "${reads}" -o unicycler_out
        mv unicycler_out/assembly.fasta assembly.fasta
        """        
}
