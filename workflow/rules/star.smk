rule star_index:
    input:
        genome=get_ref_file(config["ref"]["dna_fasta"]),
        gtf=get_ref_file(config["ref"]["gtf"]),
    output:
        chrNameLength="results/star_index/chrNameLength.txt",
    log:
        "logs/star/star_index.log",
    params:
        index_dir=lambda w, output: os.path.dirname(output[0]),
        tmp="results/star_index/STARtmp",
        read_length=lambda _: int(config["read_length"]) - 1,
    conda:
        "../envs/star.yaml"
    shell:
        "STAR"
        " --runMode genomeGenerate"
        " --runThreadN {resources.threads}"
        " --genomeDir {params.index_dir}"
        " --genomeFastaFiles {input.genome}"
        " --sjdbGTFfile {input.gtf}"
        " --sjdbOverhang {params.read_length}"
        " --outTmpDir {params.tmp}"
        " &> {log}"


rule star_align:
    input:
        fq1="results/fastp/{sample}_trimmed_R1.fastq.gz",
        fq2="results/fastp/{sample}_trimmed_R2.fastq.gz",
        idx="results/star_index/chrNameLength.txt",
    output:
        bam="results/star_align/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    log:
        "logs/star/align/{sample}.log",
    params:
        prefix=lambda w, output: output.bam.replace("Aligned.sortedByCoord.out.bam", ""),
        index_dir=lambda w, input: os.path.dirname(input.idx),
    conda:
        "../envs/star.yaml"
    shell:
        "STAR --runMode alignReads"
        " --genomeDir {params.index_dir}"
        " --readFilesIn {input.fq1} {input.fq2}"
        " --runThreadN {resources.threads}"
        " --outFileNamePrefix {params.prefix}"
        " --readFilesCommand zcat"
        " --outReadsUnmapped Fastx"
        " --twopassMode Basic"
        " --outSAMtype BAM SortedByCoordinate"
        " --outSAMattributes NH HI AS NM MD XS"
        " &> {log}"
        #" --outFilterScoreMinOverLread 0.3" 
        #" --outFilterMatchNminOverLread 0.3" 
        #" --outFilterMultimapNmax 100" 
        #" --winAnchorMultimapNmax 100"
        #" --alignEndsProtrude 5 ConcordantPair"
