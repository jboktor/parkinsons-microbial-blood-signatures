import os
import sys
import time
import subprocess as sp
import time


wkdir = '/central/groups/MazmanianLab/joeB/PDMBS/workflow/WGS/'
scriptsdir = '/central/groups/MazmanianLab/joeB/PDMBS/parkinsons-microbial-blood-signatures/notebooks/shell_scripts/'
fastq_raw = os.path.join(wkdir, "fastqs/")
fastq_clean = os.path.join(wkdir, "clean_fastqs/")
fastq_clean_stats = os.path.join(wkdir, "clean_fastqs_stats/")

for fastq_file in os.listdir(fastq_raw):
    outputloc = os.path.join(fastq_clean, fastq_file)
    if (not os.path.isfile(outputloc)):
        if (sp.getoutput('squeue -u jboktor | wc -l') > 5000): 
            time.sleep(5)
        command = 'sbatch '+ scriptsdir + 'bbduk_readQC.sh -i ' + fastq_raw + ' -o ' + fastq_clean + ' -s ' + fastq_file
        print(command)
        process = os.popen(command, 'r')
        print(process.read())
    else:
        print("File exists: " + output)