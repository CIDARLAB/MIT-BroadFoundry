
echo Total reads > ~/cir_read_data.txt

samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T1_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T1_2.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T2_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T2_2.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T3_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T3_2.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T4_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T4_2.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T5_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T5_2.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T6_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T6_2.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T7_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T7_2.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T8_1.bam.bam >> ~/cir_read_data.txt
samtools view -c /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T8_2.bam.bam >> ~/cir_read_data.txt

echo Mapped reads >> ~/cir_read_data.txt

samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T1_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T1_2.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T2_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T2_2.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T3_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T3_2.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T4_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T4_2.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T5_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T5_2.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T6_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T6_2.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T7_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T7_2.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T8_1.bam.bam >> ~/cir_read_data.txt
samtools view -c -F 4 /btl/projects/seq/foundry/H9FULADXX/pkgs/SN0025683/SplitFastq/BAMS/Circuit_sample_T8_2.bam.bam >> ~/cir_read_data.txt
