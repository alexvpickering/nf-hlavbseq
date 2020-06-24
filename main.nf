 // Define the default parameters
params.reads = ""
params.outdir = "./results"
params.samplename = ""
params.readgroup = "${params.samplename}"
params.reference = "~/hg19/Homo_sapiens_assembly19.fasta"

/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 


process map_and_sort {
	publishDir "${params.outdir}/mapped_reads"

	input:
	tuple val(pair_id), path(reads) from read_pairs_ch

	output:
	file "${params.samplename}.bam" into bamfile_ch1, bamfile_ch2, bamfile_ch3, bamfile_ch4
	
	"""
	bwa mem -t 6 -M -R '@RG\\tID:${params.readgroup}\\tSM:${params.samplename}\\tPL:ILLUMINA' $params.reference $reads \
	| samtools sort -@6 -o ${params.samplename}.bam -
	"""	
}

process index_bam {
	publishDir "${params.outdir}/mapped_reads"

	input:
	file bamfile from bamfile_ch1

	output:
	file "${bamfile}.bai" into bamidx_ch

	script:
	"""
	samtools index $bamfile
	"""
}

process extract_hla_readnames {
	publishDir "${params.outdir}/hlavbseq"

	input:
	file bamfile from bamfile_ch2
	file bamidx from bamidx_ch

	output:
	file "${bamfile.baseName}_partial_reads.txt" into hla_readnames_ch

	script:
	"""
	samtools view $bamfile 6:29907037-29915661 6:31319649-31326989 6:31234526-31241863 \
	6:32914391-32922899 6:32900406-32910847 6:32969960-32979389 6:32778540-32786825 \
	6:33030346-33050555 6:33041703-33059473 6:32603183-32613429 6:32707163-32716664 \
	6:32625241-32636466 6:32721875-32733330 6:32405619-32414826 6:32544547-32559613 \
	6:32518778-32554154 6:32483154-32559613 6:30455183-30463982 6:29689117-29699106 \
	6:29792756-29800899 6:29793613-29978954 6:29855105-29979733 6:29892236-29899009 \
	6:30225339-30236728 6:31369356-31385092 6:31460658-31480901 6:29766192-29772202 \
	6:32810986-32823755 6:32779544-32808599 6:29756731-29767588 \
	| awk '{print \$1}' | sort | uniq > ${bamfile.baseName}_partial_reads.txt
	"""
}

process search_hla_readnames {
	publishDir "${params.outdir}/hlavbseq"

	input:
	file bamfile from bamfile_ch3
	file hla_readnames from hla_readnames_ch

	output:
	file "${bamfile}.bai.idx" into readnames_idx_ch
	file "${bamfile.baseName}_partial.sam" into partial_sam_ch
	file "${bamfile.baseName}_partial_1.fastq" into partial_fastq1_ch
	file "${bamfile.baseName}_partial_2.fastq" into partial_fastq2_ch

	script:
	"""
	java -jar ~/HLAVBSeq/bamNameIndex.jar index $bamfile
	java -jar ~/HLAVBSeq/bamNameIndex.jar search $bamfile --name $hla_readnames --output ${bamfile.baseName}_partial.sam
	java -jar /usr/bin/picard.jar SamToFastq I=${bamfile.baseName}_partial.sam F=${bamfile.baseName}_partial_1.fastq F2=${bamfile.baseName}_partial_2.fastq
	"""
}

process extract_unmapped_reads {
	publishDir "${params.outdir}/hlavbseq"

	input:
	file bamfile from bamfile_ch4

	output:
	file ${bamfile.baseName}_unmapped.bam into unmapped_bam_ch
	file "${bamfile.baseName}_unmapped_1.fastq" into unmapped_fastq1_ch
	file "${bamfile.baseName}_unmapped_2.fastq" into unmapped_fastq2_ch

	script:
	"""
	samtools view -bh -f 12 $bamfile > ${bamfile.baseName}_unmapped.bam
	java -jar /usr/bin/picard.jar SamToFastq I=${bamfile.baseName}_unmapped.bam F=${bamfile.baseName}_unmapped_1.fastq F2=${bamfile.baseName}_unmapped_2.fastq
	"""
}

process combine_reads {
	publishDir "${params.outdir}/hlavbseq"

	input:
	file partial_fastq1 from partial_fastq1_ch
	file partial_fastq2 from partial_fastq2_ch
	file unmapped_fastq1 from unmapped_fastq1_ch
	file unmapped_fastq2 from unmapped_fastq2_ch

	output:
	file "${params.samplename}_part_1.fastq" into combined_fastq1_ch
	file "${params.samplename}_part_2.fastq" into combined_fastq2_ch

	script:
	"""
	cat $partial_fastq1 $unmapped_fastq1 > ${params.samplename}_part_1.fastq
	cat $partial_fastq2 $unmapped_fastq2 > ${params.samplename}_part_2.fastq
	"""
}

process align_hla {
	publishDir "${params.outdir}/hlavbseq"

	input:
	file combined_fastq1 from combined_fastq1_ch
	file combined_fastq1 from combined_fastq2_ch

	output:
	file "${params.samplename}_part.sam" into combined_sam_ch

	script:
	"""
	bwa index ~/HLAVBSeq/hla_all_v2.fasta
	bwa mem -t 8 -P -L 10000 -a ~/HLAVBSeq/hla_all_v2.fasta $combined_fastq1 $combined_fastq2 > ${params.samplename}_part.sam
	"""
}

process estimate_hla_types {
	publishDir "${params.outdir}/hlavbseq"

	input:
	combined_sam from combined_sam_ch

	output:
	file "${params.samplename}_result.txt" into result_ch

	script:
	"""
	java -jar ~/HLAVBSeq/HLAVBSeq.jar ~/HLAVBSeq/hla_all_v2.fasta $combined_sam ${params.samplename}_result.txt --alpha_zero 0.01 --is_paired
	"""
}

process call_hla_digits {

	input:
	file result from result_ch

	output:
	file "${params.samplename}_report.txt" into report_ch

	script:
	"""
	python ~/HLAVBSeq/call_hla_digits.py -v $result \
	-a ~/HLAVBSeq/Allelelist_v2.txt -r 90 -d 4 --ispaired \
	> ${params.samplename}_report.txt
	"""
}