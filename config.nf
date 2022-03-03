docker.enabled = flase

profiles {
  standard {
    process.container = 'nextflow/rnaseq-nf'
    docker.enabled = true
  }
  
  cluster {
    process.executor = 'slurm'
    process.queue = 'short'
    process.memory = '60 GB' 
    process.time = '60 min'
    process.cpus = 8     
  }

  batch {
    process.container = 'nextflow/rnaseq-nf:latest' 
    process.executor = 'awsbatch'
    process.queue = 'nextflow-ci'
    workDir = 's3://nextflow-ci/work'
    aws.region = 'eu-west-1'
    aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'
  }
} 
