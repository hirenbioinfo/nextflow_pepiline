docker.enabled = true

process {
    memory = 60.GB
}


profiles {

    local {
        process.executor = 'local'
        params.threads = 16
    }

    cluster {
        process.executor = 'sge'
        process.memory = '64GB'
        params.threads = 32
    }
