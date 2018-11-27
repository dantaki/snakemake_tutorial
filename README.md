# Skeleton for the Snakemake tutorial

This repository hosts the skeleton code needed for the [Snakemake tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html).

---

# Notes

* A Snakemake workflow is defined by specifying rules in a Snakefile
* **Rules** decompose the workflow into small steps
  * rules contain directives, which are followed by lists of files used or created by rules
  * example directives: `input`, `output`, `shell`, `script`

* `shell` directive is followed by a python string containing a shell command
  * elements of the rule are refer by braces notation (similar to `.format()` )

* Snakemake tries to generate target files: `$ snakemake {TARGET_FILE}`

* Snakemake only re-runs jobs if one of the input files is newer than one of the output files
  * or if one of the input files will be updated by another job

* Rules are generalized by named wildcards
  * `{sample}`
  * all output files of a rule have to contain exactly the same wildcards

* compose `snakemake` for multiple targets
  * `$ snakemake -np mapped_reads/{A,B}.bam`

* Snakemake automatically creates missing directories

* To let the **shell** access defined **wildcards**, use the `wildcards` object ( `{wildcards.sample}` )  for each wildcard

* Snakemake can aggregate files with `expand(<dir/{sample}.bam>, sample=SAMPLES)`
  * a list of files is obtained given the pattern of `"dir/{sample}.bam"` and a given list of samples `SAMPLES`
  * `expand("fastq/{sample}.{pair}.bam", sample=SAMPLES, pair=[1,2])`
  * Since snakemake is **top-down** define the list of samples at the **top** of the file
    * `SAMPLES = ["A", "B"]`

* Names can be specified for input or output
  * `fa="data/genome.fa"`
  * `samtools mpileup -g -f {input.fa}`

* Rules can be targets, provided the rule has no wildcards
  * Best practice: have a rule `all` at the top of the workflow
    * contains all typically desired target files as input

---

# Tutorial

## Step 1

* Mapping reads

```
rule bwa_map:
        input:
                "data/genome.fa",
                "data/samples/A.fastq"
        output:
                "mapped_reads/A.bam"
        shell:
                "bwa mem {input} | samtools view -Sb - > {output}"
```

Since the rule `bwa_map` has 2 `input` files, they will be **concatenated by whitespace in the `shell` directive** 

## Step 2

* Generalizing Rules

Using wildcards to generalize rules

```
rule bwa_map:
        input:
                "data/genome.fa",
                "data/samples/{sample}.fastq"
        output:
                "mapped_reads/{sample}.bam"
        shell:
                "bwa mem {input} | samtools view -Sb - > {output}"
```

## Step 3

* Sorting reads

This will make a new rule that takes as input the output from the `bwa_map` rule

By running `$ snakemake -np sorted_reads/B.bam` it will first run the rule `bwa_map` and then the rule `samtools_sort`

```
rule samtools_sort:
        input:
                "mapped_reads/{sample}.bam"
        output:
                "sorted_reads/{sample}.bam"
        shell:
                "samtools sort -T sorted_reads/{wildcards.sample} "
                "-O bam {input} > {output}"
```

Snakemake allows to access wildcards in the shell command via the `wildcards` object that has an attribute with the value for each wildcard

## Step 4

* Visualizing DAG of jobs

Add an index step

```
rule samtools_index:
        input:
                "sorted_reads/{sample}.bam"
        output:
                "sorted_reads/{sample}.bam.bai"
        shell:
                "samtools index {input}"
```

Take a look at the DAG of jobs with `$ snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg`

## Step 5

* Expanding arguments

`expand()` is a helper function that collects input files

`expand()` can take multiple wildcards

```
expand("fastq/{sample}.{pair}.fastq", sample=SAMPLES, pair=[1,2])
```

The list `SAMPLES` would need to be defined at the top of the file.

```
SAMPLES = ["A", "B"]
```

Add a rule for joint variant calling

Names can be specified for input or output


```
rule bcftools_call:
        input:
                fa="data/genome.fa",
                bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
                bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES),
        output:
                "calls/all.vcf"
        shell:
                "samtools mpileup -g -f {input.fa} {input.bam} | "
                "bcftools call -mv - > {output}"

```

## Step 6

* Adding custom scripts

Use the `script` directive to use custom scripts

Script paths are always relative to the referring Snakefile.

In the script, all properties of the rule like `input`, `output`, `wildcards` are available as attributes of a global `snakemake` object.

Script logic is separate from workflow logic. **Can be shared between workflows**

Best practice: use `script` whenever an inline code block would have more than a few lines of code

R scripts can be used. S4 object named `snakemake`
  * `snakemake@input[[1]]`


```
rule plot_quals:
        input:
                "calls/all.vcf"
        output:
                "plots/quals.svg"
        script:
                "scripts/plot-quals.py"
```

## Step 7

* Target Rule



```
rule all:
        input:
                "plots/quals.svg"
```
                                

