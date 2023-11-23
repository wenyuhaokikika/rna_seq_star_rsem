'''
@File    :   Snakefile
@Time    :   2023/11/16 21:40:10
@Author  :   wenyuhao 
@Version :   1.0
@Contact :   925201392@qq.com
@Desc    :   RNA-seq Pipeline (STAR+RSEM)
'''

# here put the import lib

import os
configfile: "config/config.yaml"

# here put the import lib
import os
from functools import partial
from snakemake.shell import shell
import pandas as pd

units = pd.read_table(config["units"], dtype=str,comment='#')
samples = units['sample'].to_list()
libraryD = units.set_index('sample')['library'].to_dict()

def get_fastq(wildcards):
    """
    Get fastq files of given sample-unit.
    need columns ["sample",'unit',"fq1","fq2"].
    """
    # fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    fastqs = units.loc[units["sample"]==wildcards.sample,:].iloc[0]#["fq1", "fq2"].dropna()
    #print(fastqs,dict(fastqs))
    if not pd.isna(fastqs.fq2):#len(fastqs) == 2:
        return {"fq1": fastqs.fq1, "fq2": fastqs.fq2}
    return {"fq1": fastqs.fq1}

def get_key(key,wildcards):
    return units[units['sample'] == wildcards.sample].iloc[0].to_dict()[key]

# BASE_DIR = os.path.dirname(workflow.snakefile)
rule all:
    input:
        transcriptAGG = config['workspace']+'/isoforms.results',
        geneAGG = config['workspace']+'/genes.results'

rule trim:
    output: #config['workspace']+)
        fq1=config['workspace']+'/samples/{sample}/preprocess/r1_trimed.fq.gz',
        fq2=config['workspace']+'/samples/{sample}/preprocess/r2_trimed.fq.gz',
        json=config['workspace']+'/samples/{sample}/qc/fastp.json',
        html=config['workspace']+'/samples/{sample}/qc/fastp.html'
    input:
        unpack(get_fastq)
    log:
        config['workspace']+'/log/preprocess/{sample}/fastp.log'
    threads:
        workflow.cores
    shell:
        'fastp -w {threads} -i {input.fq1} -I {input.fq2}'
        ' -o {output.fq1} -O {output.fq2}'
        ' -j {output.json} -h {output.html} 2>{log}'

rule mapping:
    output:
        config['workspace']+'/samples/{sample}/mapping/1_Aligned.toTranscriptome.out.bam'
    input:
        fq1=rules.trim.output.fq1,
        fq2=rules.trim.output.fq2
    log:
        config['workspace']+'/log/mapping/{sample}/star.log'
    params:
        library = partial(get_key,'library'),
        platform = partial(get_key,'platform'),
        genome_dir= os.path.join(config['reference']['refDir'],'STARindex'), #rules.STARindex.output,
        prefix=config['workspace']+'/samples/{sample}/mapping/1_'
    threads:
        #workflow.cores #20 if workflow.cores > 20 else workflow.cores
        config['STAR']['threads']
    priority:
        10
    shell:
        'STAR  --readFilesIn {input.fq1} {input.fq2}'
        ' --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:{params.library} PL:{params.platform}'
        ' --alignIntronMax 1000000'
        ' --alignIntronMin 20'
        ' --alignMatesGapMax 1000000'
        ' --alignSJDBoverhangMin 1'
        ' --alignSJoverhangMin 8'
        ' --alignSoftClipAtReferenceEnds Yes'
        ' --chimJunctionOverhangMin 15'
        ' --chimMainSegmentMultNmax 1'
        ' --chimOutType Junctions SeparateSAMold WithinBAM SoftClip'
        ' --chimSegmentMin 15'
        ' --genomeDir {params.genome_dir}'
        ' --genomeLoad NoSharedMemory'
        ' --limitSjdbInsertNsj 1200000'
        ' --outFileNamePrefix {params.prefix}'
        ' --outFilterIntronMotifs None'
        ' --outFilterMatchNminOverLread 0.33'
        ' --outFilterMismatchNmax 999'
        ' --outFilterMismatchNoverLmax 0.1'
        ' --outFilterMultimapNmax 20'
        ' --outFilterScoreMinOverLread 0.33'
        ' --outFilterType BySJout'
        ' --outSAMattributes NH HI AS nM NM ch'
        ' --outSAMstrandField intronMotif'
        ' --outSAMtype BAM Unsorted'
        ' --outSAMunmapped Within'
        ' --quantMode TranscriptomeSAM GeneCounts'
        ' --readFilesCommand zcat'
        ' --runThreadN {threads}'
        ' --twopassMode Basic >{log} 2>&1'
    
rule run_rsem:
    output:
        iso = config['workspace']+'/samples/{sample}/rsem/rsem.genes.results',
        gene = config['workspace']+'/samples/{sample}/rsem/rsem.isoforms.results'
    input:
        bam = config['workspace']+'/samples/{sample}/mapping/1_Aligned.toTranscriptome.out.bam' #rules.mapping.output,
    threads:
        10
    params:
        seed = config['seed'],
        ref = os.path.join(config['reference']['refDir'],'rsem_ref/rsemRef'),
        prefix = config['workspace']+'/samples/{sample}/rsem/rsem'
    log:
        config['workspace']+'/log/rsem/{sample}/rsem-calculate-expression.log'
    shell:
        'rsem-calculate-expression'
        ' -p {threads}'
        ' --alignments'
        ' --paired-end'
        ' --bam {input}'
        ' --seed {params.seed}'
        ' {params.ref}'
        ' {params.prefix}  >{log} 2>&1'

rule aggregate:
    output:
        transcript = config['workspace']+'/isoforms.results',
        gene = config['workspace']+'/genes.results'
    input:
        transcript = expand(config['workspace']+'/samples/{sample}/rsem/rsem.isoforms.results',sample=samples),
        gene = expand(config['workspace']+'/samples/{sample}/rsem/rsem.genes.results',sample=samples)
    # params:
    #     names=units['sample'].to_list()
    run:
        import pandas as pd
        from functools import reduce
        isoforms = dict(zip(samples,[pd.read_table(i) for i in input.transcript]))
        gene = dict(zip(samples,[pd.read_table(i) for i in input.gene]))
        
        isoforms = [isoforms[k].assign(**{f'{cols}_{k}':lambda x:x[cols] for cols in ['expected_count','TPM','FPKM','IsoPct']})\
                                                .drop(['gene_id','length','effective_length','expected_count','TPM','FPKM','IsoPct'],axis=1)\
                                                    for k in isoforms.keys()]
        isoforms = reduce(lambda x,y:x.merge(y,on = 'transcript_id',how='outer'),isoforms)

        gene = [gene[k].assign(**{f'{cols}_{k}':lambda x:x[cols] for cols in ['expected_count','TPM','FPKM']})\
                                                .drop(['transcript_id(s)','length','effective_length','expected_count','TPM','FPKM'],axis=1)\
                                                    for k in gene.keys()]
        gene = reduce(lambda x,y:x.merge(y,on = 'gene_id',how='outer'),gene)
        isoforms.to_csv(output.transcript,sep='\t',index=False)
        gene.to_csv(output.gene,sep='\t',index=False)
