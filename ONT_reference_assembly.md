# Reference assembly using ONT Nanopore data

This tutorial is a basic introduction to reference assembly commands for ONT data. 


## Quality filtering and trimming

The first step in processing ONT data is to trim low quality 
reads and check the overall quality of the data.

```
/software/Porechop-v0.2.4/porechop-runner.py -i /path/to/directory/with/reads/ -o runA.fastq --discard_middle
```

You can change `runA.fastq` to a name that is more descriptive to your project.

The general statistics of the quality filtered data can then be obtained with:

```
Nanostats --fastq runA.fastq
```

The output from this will give you lots of useful information about the quality of the
 data. It being ONT data, you expect to see some nice long reads, a mean >5000 is not 
 uncommon and you expect a good proportion of your bases to be >Q15 in terms of 
 quality. Having a large number of long reads will make it easier for detecting 
 recombination or scaffolding genomes.
 
## Mapping reads to a reference

Assuming you have a reference in a fasta format in a file called `reference.fa`, 
you can use `minimap2` combined with `samtools` by piping one command into the next: 

```
minimap2 -ax map-ont reference.fa runA.fastq | samtools view -bS -F 4 | samtools sort -o runAtoRef.bam
```

You should name the output BAM file to something more memorable for your project. By using the 
command `-F 4` we are only keeping in the BAM file reads that map to the reference.


We can use `NanoStat` again to find out the stats of the reads that have mapped:

```
NanoStat --bam runAtoRef.bam
```

If mapping to a large genome like the human genome you need to use an additional parameter `split-prefix`. 
We have the human genome already stored on the server at `/db/wgs/Human/ensembl99/bwa/Human.fa` 

```
minimap2 -ax map-ont --split-prefix /tmp/runAtohuman /db/wgs/Human/ensembl99/bwa/Human.fa runA.fastq | samtools view -bS -F 4 | samtools sort -o runAtoHuman.bam
```

It is then a good idea to download the BAM file and view it in Tablet to get a feel 
for how the reads are spread across the genome.

## Extracting specific reads

It can be useful to extract reads mapping in specific region of your reference genome. 
For example, we may want to extract all the identifiers of reads that start in the 
first 1578 bases of the reference genome:

```
samtools view -S runAtoRef.bam | awk '$4 > 0 && ($4 < 1578){print $1}' > ListOfInteresting_readIDs.txt
``` 

Likewise we can look for reads mapping on specific chromosomes of the human genome:

```
samtools view -S runAtoHuman.bam 9:138331210 | awk '{print $1 "\t" $4 "\t" $10 "\t" $6}' > Read_info.txt
```

