version 1.0

task filterUnmappedReadsTask {
    input {
      String sampleName
      File cram
      File refFasta
      Int cpu
      Int memory
      Int maxRetries = 2
      Int runtime_preemptible = 10 
      Float runtime_disk_multiplier = 3
      Int runtime_disk_gb = ceil(size(cram, 'GB') * runtime_disk_multiplier)
    }

    command {
        # collect samtools flagstat on human-genome aligned cram file
        samtools flagstat -@ "${cpu}" -O tsv "${cram}" > "${sampleName}_flagstat.tsv"

        # convert to bam
        samtools view -@ "${cpu}" -h -b -f 4 -T "${refFasta}" "${cram}" > "${sampleName}_unmapped.bam"

        # collect unfiltered bam file stats
        samtools flagstat -@ "${cpu}" -O tsv "${sampleName}_unmapped.bam" > "${sampleName}_unmapped-f4_flagstat.tsv"
    }

    output {
        File cram_flagstat="${sampleName}_flagstat.tsv"
        File bam_flagstat="${sampleName}_unmapped-f4_flagstat.tsv"
        File unmapped_bam="${sampleName}_unmapped.bam"
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
        File cram
        File refFasta
        Int cpu
        Int memory
    }

    call filterUnmappedReadsTask {
        input:
            sampleName=sampleName,
            cram=cram,
            refFasta=refFasta,
            cpu=cpu,
            memory=memory
            }

    output {
        File cram_flagstat = filterUnmappedReadsTask.cram_flagstat
        File bam_flagstat = filterUnmappedReadsTask.bam_flagstat
        File unmapped_bam = filterUnmappedReadsTask.unmapped_bam
    }
}