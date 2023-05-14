nextflow.enable.dsl=2

params {
    threads = 16
    reads = "$baseDir/data/*_R{1,2}.f*q.gz"
    outdir = 'results'
}

println "BINNING   PIPELINE"
println "================================="
println "reads              : ${params.reads}"

process fastp {

    container 'docker://dgg32/fastp:0.20.1'

    input:
    tuple val(sample_id), file(reads)

    output:
    tuple val(sample_id), file('*.fastq.gz')

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
    -o ${sample_id}_trim_R1.fastq.gz -O ${sample_id}_trim_R2.fastq.gz \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -g --detect_adapter_for_pe -l 100 -w ${params.threads}
    """
}

process megahit {

    container 'docker://dgg32/megahit:v1.2.9'

    input:
    tuple val(sample_id), file(reads) from fastp

    output:
    tuple val(sample_id), path("${sample_id}_megahitout"), file(reads)

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} -o ${sample_id}_megahitout -t ${params.threads}
    """
}

process maxbin {

    container 'docker://nanozoo/maxbin2:v2.2.7'

    input:
    tuple val(sample_id), path(megahitout), file(reads) from megahit

    output:
    tuple val(sample_id), path("${sample_id}_maxbinout")

    script:
    """
    mkdir ${sample_id}_maxbinout
    run_MaxBin.pl -contig ${megahitout}/final.contigs.fa -out ${sample_id}_maxbinout/maxbin -reads ${reads[0]} -reads2 ${reads[1]} -thread ${params.threads}
    """
}

process checkm {

    container 'docker://dgg32/checkm:v1.1.3'

    input:
    tuple val(sample_id), path(maxbinout) from maxbin

    output:
    tuple val(sample_id), path('*_checkmout')

    script:
    """
    checkm lineage_wf -t 32 -x fasta $maxbinout ${sample_id}_checkmout
    """
}

help {
    """
    Usage: nextflow run main.nf --reads 'reads_pattern'

    Pipeline Options:
    --reads     Input read files pattern (required)
    --outdir    Output directory (default: 'results')

    Other Options:
    --threads   Number of threads to use (default: 16)

    Example:
    nextflow run main.nf --reads '/path/to/reads/*_R{1,2}.fastq.gz'
    """
}

workflow {

    Channel
        .fromFilePairs(params.reads)
        .ifEmpty { error "Cannot find any reads matching" }
        .set { read_pairs }

    Channel
        .from(read_pairs)
        .map { sample_id, reads -> tuple(sample_id, file(reads)) }
        .into(fastp)

    polished_reads = fastp.out

    assembly = megahit.out

    bins = maxbin.out
