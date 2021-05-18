import os
configfile: "config.yaml"

data_d = 'data'
read = ['1', '2']

fastqc_d = 'fastqc'
multiqc_dir = 'multiqc'
trimmed_dir = 'trimmed'

rule all:
    input:
        os.path.join(multiqc_dir, 'multiqc_report.html'),
        os.path.join(multiqc_dir, trimmed_dir, 'multiqc_report.html')   
                
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
