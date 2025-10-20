#!/bin/bash

# Add 10x to path
export PATH=/homevol/apattison/Code/bin/cellranger/cellranger-7.2.0:$PATH

# https://kb.10xgenomics.com/hc/en-us/articles/360035250172-Why-do-I-not-see-any-FASTQs-or-incomplete-FASTQs-for-the-SRA-of-my-interest-
# wget the bam file
#wget https://sra-pub-src-1.s3.amazonaws.com/SRR7159843/mouse_DSS_1.bam.1

# Use 10x bam to fastq to get the reads back out for reprocessing
#cd /pvol/andrew/shani_single_cell/raw/

# cellranger bamtofastq --nthreads=22 ./mouse_DSS_1.bam.1 ./mouse_DSS_1
# cellranger bamtofastq --nthreads=22 ./Mouse_DSS_2.bam.1 ./mouse_DSS_2
# cellranger bamtofastq --nthreads=8 ./Mouse_DSS_3.bam.1 ./mouse_DSS_3
# cellranger bamtofastq --nthreads=3 ./Mouse_HC_1.bam.1 ./mouse_HC_1
# cellranger bamtofastq --nthreads=3 ./Mouse_HC_2.bam.1 ./mouse_HC_2
# cellranger bamtofastq --nthreads=3 ./Mouse_HC_3.bam.1 ./mouse_HC_3

# Reading the outputs
# https://kb.10xgenomics.com/hc/en-us/articles/360058600992-How-do-I-find-out-which-FASTQ-files-belong-to-which-library-in-10x-Genomics-bamtofastq-output-folders-

#source /homevol/apattison/Code/bin/cellranger/cellranger-7.2.0/sourceme.bash

cd /pvol/andrew/shani_single_cell/results

#samtools view -H mouse_DSS_1.bam.1
# https://github.com/10XGenomics/bamtofastq/issues/23 need to do a comma separated list is bam2fastq 
# outs multiple files

f1=/pvol/andrew/shani_single_cell/raw/mouse_DSS_1/DSS1_3L_MissingLibrary_1_HKJCWBBXX
f2=/pvol/andrew/shani_single_cell/raw/mouse_DSS_1/DSS1_3L_MissingLibrary_1_HKJNNBBXX

cellranger count --id=SRR7159843 \
                   --transcriptome=/homevol/apattison/Code/bin/cellranger/refdata-gex-mm10-2020-A \
                   --fastqs=$f1,$f2 \
                   --localcores=30 \
                   --localmem=100
                   
f1=/pvol/andrew/shani_single_cell/raw/mouse_DSS_2/DSS2_3L_MissingLibrary_1_HKJCWBBXX
f2=/pvol/andrew/shani_single_cell/raw/mouse_DSS_2/DSS2_3L_MissingLibrary_1_HKJNNBBXX

cellranger count --id=SRR7159844 \
                   --transcriptome=/homevol/apattison/Code/bin/cellranger/refdata-gex-mm10-2020-A \
                   --fastqs=$f1,$f2 \
                   --localcores=30 \
                   --localmem=100                  

f1=/pvol/andrew/shani_single_cell/raw/mouse_DSS_3/DSS3_3L_MissingLibrary_1_HKJCWBBXX
f2=/pvol/andrew/shani_single_cell/raw/mouse_DSS_3/DSS3_3L_MissingLibrary_1_HKJNNBBXX

cellranger count --id=SRR7159845 \
                   --transcriptome=/homevol/apattison/Code/bin/cellranger/refdata-gex-mm10-2020-A \
                   --fastqs=$f1,$f2 \
                   --localcores=30 \
                   --localmem=100                      
                   
f1=/pvol/andrew/shani_single_cell/raw/mouse_HC_1/HC1_3L_MissingLibrary_1_HKJCWBBXX
f2=/pvol/andrew/shani_single_cell/raw/mouse_HC_1/HC1_3L_MissingLibrary_1_HKJNNBBXX

cellranger count --id=SRR7159840 \
                   --transcriptome=/homevol/apattison/Code/bin/cellranger/refdata-gex-mm10-2020-A \
                   --fastqs=$f1,$f2 \
                   --localcores=30 \
                   --localmem=100  
                   
f1=/pvol/andrew/shani_single_cell/raw/mouse_HC_2/HC2_3L_MissingLibrary_1_HKJCWBBXX
f2=/pvol/andrew/shani_single_cell/raw/mouse_HC_2/HC2_3L_MissingLibrary_1_HKJNNBBXX

cellranger count --id=SRR7159841 \
                   --transcriptome=/homevol/apattison/Code/bin/cellranger/refdata-gex-mm10-2020-A \
                   --fastqs=$f1,$f2 \
                   --localcores=30 \
                   --localmem=100   
                   
f1=/pvol/andrew/shani_single_cell/raw/mouse_HC_3/HC3_3L_MissingLibrary_1_HKJCWBBXX
f2=/pvol/andrew/shani_single_cell/raw/mouse_HC_3/HC3_3L_MissingLibrary_1_HKJNNBBXX

cellranger count --id=SRR7159842 \
                   --transcriptome=/homevol/apattison/Code/bin/cellranger/refdata-gex-mm10-2020-A \
                   --fastqs=$f1,$f2 \
                   --localcores=30 \
                   --localmem=100   