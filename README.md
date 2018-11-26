# Skeleton for the Snakemake tutorial

This repository hosts the skeleton code needed for the [Snakemake tutorial](http://snakemake.readthedocs.io/en/latest/tutorial/welcome.html).

---

# Notes

* A Snakemake workflow is defined by specifying rules in a Snakefile
* **Rules** decompose the workflow into small steps
  * rules contain directives, which are followed by lists of files used or created by rules
  * example directives: `input`, `output`, `shell`

* `shell` directive is followed by a python string containing a shell command
  * elements of the rule are refer by braces notation (similar to `.format()`


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

