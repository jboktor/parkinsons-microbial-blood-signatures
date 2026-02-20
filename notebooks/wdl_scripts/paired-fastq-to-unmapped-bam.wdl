version 1.0

## Modified version of : https://github.com/amp-pd/amp-pd-workflows/blob/master/wgs/paired-fastq-to-unmapped-bam/paired-fastq-to-unmapped-bam.wdl
## 
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - One or more read groups, one per pair of FASTQ files 
##
## Outputs :
## - Set of unmapped BAMs, one per read group


# WORKFLOW DEFINITION
workflow ConvertPairedFastQsToUnmappedBamWf {
  Array[String] readgroup_list
  Map[String, Array[File]] fastq_pairs
  Map[String, Array[String]] metadata

  Int preemptible_tries

  # Convert multiple pairs of input fastqs in parallel
  scatter (readgroup in readgroup_list) {

    # Convert pair of FASTQs to uBAM
    call PairedFastQsToUnmappedBAM {
      input:
        fastq_1 = fastq_pairs[readgroup][0],
        fastq_2 = fastq_pairs[readgroup][1],
        readgroup_name = readgroup,
        sample_name = metadata[readgroup][0],
        library_name = metadata[readgroup][1],
        platform_unit = metadata[readgroup][2],
        run_date = metadata[readgroup][3],
        platform_name = metadata[readgroup][4],
        sequencing_center = metadata[readgroup][5],
        platform_model = metadata[readgroup][6],

        preemptible_tries = preemptible_tries
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bams = PairedFastQsToUnmappedBAM.output_bam
  }
}

# TASK DEFINITIONS

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  File fastq_1
  File fastq_2
  String readgroup_name
  String sample_name
  Int disk_size
  String mem_size
  String docker
  String gatk_path
  Int preemptible_tries

  command {
    ${gatk_path} --java-options "-Xmx3000m" \
      FastqToSam \
      --FASTQ ${fastq_1} \
      --FASTQ2 ${fastq_2} \
      --OUTPUT ${readgroup_name}.unmapped.bam \
      --READ_GROUP_NAME ${readgroup_name} \
      --SAMPLE_NAME ${sample_name} \
      --LIBRARY_NAME ${library_name} \
      --PLATFORM_UNIT ${platform_unit} \
      --RUN_DATE ${run_date} \
      --PLATFORM ${platform_name} \
      --SEQUENCING_CENTER ${sequencing_center} \
      --PLATFORM_MODEL ${platform_model} 
  }
  runtime {
    docker: docker
    memory: mem_size
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${readgroup_name}.unmapped.bam"
  }
}
