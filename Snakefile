import os
configfile: "config.yaml"

data_d = 'data'
star_d ='stardata'
read = ['1', '2']
fastqc_d = 'fastqc'
multiqc_dir = 'multiqc'
trimmed_dir = 'trimmed'
idx_g= 'idx_genome'
align_out= 'aligned'
sam_d = 'sam_sorted'
fcounts_d='fcounts'
deseq_d= 'deseq'

rule all:
    input:
        os.path.join(multiqc_dir, 'multiqc_report.html'),
        os.path.join(multiqc_dir, trimmed_dir, 'multiqc_report.html'),
        os.path.join(idx_g, 'genomeParameters.txt'),
        os.path.join(idx_g, 'Genome'),
        expand(os.path.join(align_out, '{sample}', 'Aligned.sortedByCoord.out.bam'), sample=config['samples']),
        expand(os.path.join(sam_d, '{sample}.bam'), sample=config['samples']),
        os.path.join(fcounts_d, 'fcount_result.txt'),   
        expand(os.path.join(deseq_d, '{prep}_DE_volcano.pdf'), prep=config['prep']),
        expand(os.path.join(deseq_d, '{prep}_DE_results.csv'), prep=config['prep'])        
rule fastqc:
    input:
       os.path.join(data_d, '{sample}_R{read}_001.fastq')
    output:
        os.path.join(fastqc_d, '{sample}_R{read}_001_fastqc.html')
    shell:
        "fastqc {input} -o {fastqc_d}"

rule multiqc: 
    input: 
        expand(os.path.join(fastqc_d, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=read)
    output:
        os.path.join(multiqc_dir, 'multiqc_report.html')
    shell: 
        "multiqc {fastqc_d} -o {multiqc_dir}"
        
        
rule trimm:
    input:
        r1 = os.path.join(data_d, '{sample}_R1_001.fastq'),
        r2 = os.path.join(data_d, '{sample}_R2_001.fastq')
    output:
        out1 = os.path.join(trimmed_dir, '{sample}_R1_001.fastq'),
        out2 = os.path.join(trimmed_dir, '{sample}_R2_001.fastq')
    shell: 
        "bbmap/bbduk.sh in1={input.r1} out1={output.out1} in2={input.r2} out2={output.out2} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"
        
        
rule fastqc_trimmed:
    input:
       os.path.join(trimmed_dir, '{sample}_R{read}_001.fastq')
    output:
        os.path.join(fastqc_d, trimmed_dir, '{sample}_R{read}_001_fastqc.html')
    shell:
        "fastqc {input} --outdir={fastqc_d}/{trimmed_dir}"


rule multiqc_trimmed: 
    input: 
        expand(os.path.join(fastqc_d, trimmed_dir, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=read)
    output:
        os.path.join(multiqc_dir, trimmed_dir, 'multiqc_report.html')
    shell: 
        "multiqc {fastqc_d}/{trimmed_dir} -o {multiqc_dir}/{trimmed_dir}"

rule indexing:
     input:
        gtf=os.path.join(star_d, 'chr19_20Mb.gtf'),
        fa=os.path.join(star_d, 'chr19_20Mb.fa')
     output:
        os.path.join(idx_g, 'genomeParameters.txt'),
        os.path.join(idx_g, 'Genome')
     shell:
        "STAR --genomeDir {idx_g} --runMode genomeGenerate --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeSAindexNbases 11 --outFileNamePrefix star_files"

rule aligning:
      input:
        os.path.join(idx_g, 'genomeParameters.txt'),
        os.path.join(idx_g, 'Genome'),
        a1 = os.path.join(trimmed_dir, '{sample}_R1_001.fastq'),
        a2 = os.path.join(trimmed_dir, '{sample}_R2_001.fastq'),
        gtf = os.path.join(star_d, 'chr19_20Mb.gtf')
      output:
        os.path.join(align_out, '{sample}', 'Aligned.sortedByCoord.out.bam')
      shell:
        "STAR --readFilesIn {input.a1} {input.a2} --genomeDir {idx_g} --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {input.gtf} --outFileNamePrefix {align_out}/{wildcards.sample}/"

rule sam_sort:
    input:
        os.path.join(align_out, '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        os.path.join(sam_d, '{sample}.bam')
    shell:
        "samtools sort -n -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

def fcounts_inputs(wildcards):
    files = expand(os.path.join(sam_d, '{sample}.bam'), sample=config['samples'])
    return files

rule f_counts:
    input:
        fcounts_inputs
    output:
        os.path.join(fcounts_d, 'fcount_result.txt')
    shell:
        "featureCounts {input} -p -t exon -g gene_id -a {star_d}/chr19_20Mb.gtf -o {output} -s 1"

rule deseq:
    input:
        os.path.join(fcounts_d, 'fcount_result.txt'),
    params:
        deseq_d,
        config['prep']
    output:
        os.path.join(deseq_d, '{prep}_DE_volcano.pdf'),
        os.path.join(deseq_d, '{prep}_DE_results.csv')
    script:
        "deseq2.R"
