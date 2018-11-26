# Skeleton for the Snakemake tutorial

This repository hosts the skeleton code needed for the [Snakemake tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html).

---

# Notes

* A Snakemake workflow is defined by specifying rules in a Snakefile
* **Rules** decompose the workflow into small steps
  * rules contain directives, which are followed by lists of files used or created by rules
  * example directives: `input`, `output`, `shell`

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
