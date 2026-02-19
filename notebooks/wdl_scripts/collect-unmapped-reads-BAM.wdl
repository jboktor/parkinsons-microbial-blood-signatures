version 1.0

task filterUnmappedReadsTask {
    input {
      String sampleName
      File bam
      Int cpu
      Int memory
      Int maxRetries = 2
      Int runtime_preemptible = 25 
      Float runtime_disk_multiplier = 3
      Int runtime_disk_gb = ceil(size(bam, 'GB') * runtime_disk_multiplier)
    }

    command {
        # collect samtools flagstat on bam
        samtools flagstat -@ "${cpu}" -O tsv "${bam}" > "${sampleName}_flagstat.tsv"

        # filter unmapped reads from bam
        samtools view -@ "${cpu}" -h -b -f 4 "${bam}" > "${sampleName}_unmapped.bam"

        # collect samtools flagstat on filtered bam
        samtools flagstat -@ "${cpu}" -O tsv "${sampleName}_unmapped.bam" > "${sampleName}_unmapped-f4_flagstat.tsv"
    }

    output {
        File bam_flagstat="${sampleName}_flagstat.tsv"
        File unmapped_bam="${sampleName}_unmapped.bam"
        File unmapped_bam_flagstat="${sampleName}_unmapped-f4_flagstat.tsv"
    }

    runtime {
        disks: 'local-disk ${runtime_disk_gb} SSD'
        memory: '${memory} GB'
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        preemptible: "${runtime_preemptible}"
        maxRetries: "${maxRetries}"
    }

    meta {
        author: "Joe Boktor"
    }
}

workflow filterUnmappedReadsWorkflow {
    input {
        String sampleName
        File bam
        Int cpu
        Int memory
    }

    call filterUnmappedReadsTask {
        input:
            sampleName=sampleName,
            bam=bam,
            cpu=cpu,
            memory=memory
            }

    output {
        File bam_flagstat = filterUnmappedReadsTask.bam_flagstat
        File unmapped_bam = filterUnmappedReadsTask.unmapped_bam
        File unmapped_bam_flagstat = filterUnmappedReadsTask.unmapped_bam_flagstat
    }
}