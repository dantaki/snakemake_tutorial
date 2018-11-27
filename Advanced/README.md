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
---

# Rule parameters

define arbitrary parameters

```
rule bwa_map:
        input:
                ...
        params:
                rg=r"@RG\tID:{sample}\tSM:{sample}"
        shell:
                "bwa mem -R '{params.rg}' ... "
```

---

# Logging

Define log files with the `log` directive.

not subject to rule matching and not cleaned up when a job fails


modify shell command to collect `STDERR` output of both `bwa` and `samtools` and
pipe it into the file referred by `{log}`

** must contain the same wildcards as the output file to avoid collisions **

```
rule bwa_map:
  ...

        log:
                "logs/bwa_mem/{sample}.log"
        shell:
                "(bwa mem -R '{params.rg}' ... | "
                "samtools view -Sb -> {output}) 2> {log}" 
```

---

# Temporary and Protected Files

Snakemake can mark output files as temporary
  * deleted once every consuming job has been executed

```
rule bwa_map:
        ...
        
        output:
                temp("mapped_reads/{sample}.bam")
        ...

```

Protect files from accidental deletion or modification

```
rule samtools_sort:
        ...
        
        output:
                protected("sorted_reads/{sample}.bam")
        ...
```

---

# Benchmarking

Use the `benchmark` directive to measure the wall clock time and memory usage of a job

** must contain the same wildcards as output to avoid collisions **

can repeat benchmarks `benchmark("...", 3)`


```
rule bwa_map:
        ...

        benchmark:
                "benchmarks/{sample}.bwa.benchmark.txt"
        ...
```

---

# Cluster Execution

```
$ snakemake --cluster qsub --jobs 100
```

Each job will be compiled into a shell script that is submitted with the given command (`qsub`)

`--jobs` limits the number of concurrent jobs

Snakemake supports a separate configuration file for execution on a cluster

Parameters in the cluster config file are accessed by the `cluster.*` wildcard

`cluster.yaml`

```
__default__:
        n: 1
        queue: condo
        e: '{log}.err'
        o: '{log}.out'
        time: "48:00:00"

```

`$ snakemake --cluster-config cluster.yaml --cluster "qsub -q {cluster.queue} -n {cluster.n} -W {cluster.time}" `

---


