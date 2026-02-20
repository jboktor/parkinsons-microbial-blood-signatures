version 1.0

task arcasHLAanalysis {
    input {
        File bam
        File bai
        String sampleName
        Int threads
        Int memoryGB
        Float diskMultiplier = 3
        Int runtimeDiskGB = ceil((size(bam, 'GB') * diskMultiplier) + 50)
        String dockerImage = "jboktor/arcashla_base:latest"
        }

    command {
        cd /cromwell_root && \
        git clone https://github.com/RabadanLab/arcasHLA.git && \
        conda env create -f /cromwell_root/arcasHLA/environment.yml && \
        echo "source activate arcas-hla" > ~/.bashrc && \
        PATH=/opt/conda/envs/arcas-hla/bin:$PATH && \
        mkdir -p "/cromwell_root/tmp/arcasHLA_${sampleName}" && \
        /cromwell_root/arcasHLA/arcasHLA extract \
            --threads ${threads} \
            --outdir "/cromwell_root" \
            --log "/cromwell_root/${sampleName}.extract.log" \
            --temp "/cromwell_root/tmp/arcasHLA_${sampleName}" \
            --verbose \
            ${bam} && \
        /cromwell_root/arcasHLA/arcasHLA genotype \
            --threads ${threads} \
            --outdir "/cromwell_root" \
            --log "/cromwell_root/${sampleName}.genotype.log" \
            --temp "/cromwell_root/tmp/arcasHLA_${sampleName}" \
            --verbose \
            "/cromwell_root/${sampleName}.star.extracted.1.fq.gz" \
            "/cromwell_root/${sampleName}.star.extracted.2.fq.gz"
    }

    output {
        File extract_f1 = "/cromwell_root/${sampleName}.star.extracted.1.fq.gz"
        File extract_f2 = "/cromwell_root/${sampleName}.star.extracted.2.fq.gz"
        File extract_log = "/cromwell_root/${sampleName}.extract.log"
        File alignment = "/cromwell_root/${sampleName}.alignment.p"
        File em_json = "/cromwell_root/${sampleName}.em.json"
        File genotype_json = "/cromwell_root/${sampleName}.genotype.json"
        File genotype_log = "/cromwell_root/${sampleName}.genotype.log"
    }

    runtime {
        docker: "${dockerImage}"
        memory: "${memoryGB} GB"
        cpu: threads
        disks: "local-disk ${runtimeDiskGB} HDD"
    }
}

workflow arcasHLAworkflow {
    input {
        File bam
        File bai
        String sampleName
        Int threads
        Int memoryGB
    }

    call arcasHLAanalysis {
        input:
        bam = bam,
        bai = bai,
        sampleName = sampleName,
        threads = threads,
        memoryGB = memoryGB
    }

    output {
        File extract_f1 = arcasHLAanalysis.extract_f1
        File extract_f2 = arcasHLAanalysis.extract_f2
        File extract_log = arcasHLAanalysis.extract_log
        File alignment = arcasHLAanalysis.alignment
        File em_json = arcasHLAanalysis.em_json
        File genotype_json = arcasHLAanalysis.genotype_json
        File genotype_log = arcasHLAanalysis.genotype_log
    }
}
