# Advanced Snakemake Usage

---

# Threads

Use `threads` directive to define threads

The number of threads the jobs need is considered by the Snakemake scheduler
  * the sum of the threads of all running jobs does not exceed a given number of available CPUs cores
  * cores can be given with the `--cores` option

Threads directive in a rule is interpreted as a maximum
  * when less cores than threads are provided, the number of threads a rule uses will be reduced to the number of given cores


```
rule bwa_map:
        input:
                fa="data/genome.fa",
                fq="data/samples/{sample}.fastq"
        output:
                "mapped_reads/{sample}.bam"
        threads: 8
        shell:
                "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

---

# Config Files

Uses either JSON or YAML
`configfile` directive

Typically added to the top of Snakefiles. Store contents into a global dictionary named `config`

Config can store samples

`config.yaml`

```
samples:
        A: data/samples/A.fastq
        B: data.samples/B.fastq

```

in `Snakefile`


```
configfile: "config.yaml"

rule bcftools_call:
        input:
                fa="data/genome.fa",
                bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        ...

```

---

# Input Functions

can use lambda expressions or normal functions

Input functions take as a single argument a `wildcards` object

Required to return a string or a list of strings
  * typically paths to input files

```
rule bwa_map:
        input:
                "data/genome.fa",
                lambda wildcards: config["samples"][wildcards.sample]
        ...
```

