version 1.0

task downloadTask {
    input {
      String sampleName
      File url
      Int cpu
      Int memory
      Int maxRetries = 3
      Int runtime_preemptible = 100 
      Int runtime_disk_gb = 25
    }

    command {
        wget --no-check-certificate "${url}"        
    }

    output {
        File wgs_reads="${sampleName}.tar.bz2"
    }

    runtime {
        disks: 'local-disk ${runtime_disk_gb} GB SSD'
        memory: '${memory} GB'
        docker: 'alpine:3.17.0'
        preemptible: "${runtime_preemptible}"
        maxRetries: "${maxRetries}"
    }

    meta {
        author: "Joe Boktor"
    }
}

workflow downloadHMPWorkflow {
    input {
        String sampleName
        File url
        Int cpu
        Int memory
    }

    call downloadTask {
        input:
            sampleName=sampleName,
            url=url,
            cpu=cpu,
            memory=memory
            }

    output {
        File wgs_reads = downloadTask.wgs_reads
    }
}