
SAMPLES, = glob_wildcards(config["data_dir"]+"CD630_{sample}.fasta")


rule all:
  input:   
    expand("{sample}/intarna_out_{sample}.txt", sample=SAMPLES),
    expand("{sample}/passFilter_Complete_{sample}.txt", sample=SAMPLES),
    expand("{sample}/passFilter_list_{sample}.tsv", sample=SAMPLES),
    expand("{sample}/passFilter_OnlyLink_{sample}.txt", sample=SAMPLES),
    expand("{sample}/linkPosition_{sample}.gff", sample=SAMPLES),
    expand("{sample}/linkFiles/{sample}_links_inta.txt", sample=SAMPLES),
    expand("{sample}/linkFiles/{sample}_links_rcd_inta.txt", sample=SAMPLES),
    expand("{sample}/linkFiles/{sample}_label_inta.txt", sample=SAMPLES),
    expand("plots/{sample}_circosR_inta.png", sample=SAMPLES)
      


rule intarna:
  output:
    "{sample}/intarna_out_{sample}.txt"
  input:
    query=config["data_dir"]+"CD630_{sample}.fasta",
    allfasta=config["data_dir"]+config["allseq_file"]
  conda:
    config["env_dir"]+"condaEnv_intaRNA.yml"
  log:
    out="Logs/inta_{sample}.stdout",
    err="Logs/inta_{sample}.stderr"
  shell:
    "IntaRNA -q {input.query} -t {input.allfasta} --outMode D --out {wildcards.sample}/intarna_out_{wildcards.sample}.txt"

rule parseIntaRes:
  output:
    "{sample}/passFilter_Complete_{sample}.txt",
    "{sample}/passFilter_list_{sample}.tsv",
    "{sample}/passFilter_OnlyLink_{sample}.txt"
  input:
    "{sample}/intarna_out_{sample}.txt"
  params:
    scriptPath=config["work_dir"]+config["scripts_dir"],
    engThd=config["eng_threshold"],
  log:
    out="Logs/intaParse_{sample}.stdout",
    err="Logs/intaParse_{sample}.stderr"
  shell:
    "python3 {params.scriptPath}parseIntaOut.py {input} {params.engThd} {wildcards.sample}"

rule intaGff:
  output:
    "{sample}/linkPosition_{sample}.gff"
  input:
    "{sample}/passFilter_Complete_{sample}.txt"
  params:
    scriptPath=config["work_dir"]+config["scripts_dir"],
    tlList=config["path2TlList"],
    pPos=config["data_dir"]+config["path2PPos"],
    mPos=config["data_dir"]+config["path2MPos"]
  log:
    out="Logs/intaGff_{sample}.stdout",
    err="Logs/intaGff_{sample}.stderr"
  shell:
    "python3 {params.scriptPath}gffFromInta.py {params.tlList} {input} {wildcards.sample} {params.pPos} {params.mPos}"


rule intaLinkFile:
  output:
    "{sample}/linkFiles/{sample}_links_inta.txt",
    "{sample}/linkFiles/{sample}_links_rcd_inta.txt",
    "{sample}/linkFiles/{sample}_label_inta.txt"
  input:
    "{sample}/passFilter_list_{sample}.tsv"
  params:
    scriptPath=config["work_dir"]+config["scripts_dir"],
    gffFile=config["data_dir"]+config["path2gff"],
    labelFile=config["data_dir"]+config["path2label"],
    symbolFile=config["data_dir"]+config["path2symbol"]
  log:
    out="Logs/intaLink_{sample}.stdout",
    err="Logs/intaLink_{sample}.stderr"
  shell:
    "python3 {params.scriptPath}create_linkFile_inta.py {input} {params.gffFile} {wildcards.sample} {params.labelFile} {params.symbolFile}"

rule intaCircos:
  output:
    "plots/{sample}_circosR_inta.png"
  input:
    "{sample}/linkFiles/{sample}_links_inta.txt",
    "{sample}/linkFiles/{sample}_links_rcd_inta.txt",
    "{sample}/linkFiles/{sample}_label_inta.txt"
  params:
    scriptPath=config["work_dir"]+config["scripts_dir"],
    engThd=config["eng_threshold"]
  log:
    out="Logs/intaCircos_{sample}.stdout",
    err="Logs/intaCircos_{sample}.stderr"
  shell:
    "Rscript {params.scriptPath}intacircos.r {params.engThd} $(pwd) {wildcards.sample}"
