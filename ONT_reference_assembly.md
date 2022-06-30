# Reference assembly using ONT Nanopore data



## Quality filtering and trimming

The first step in porcessing ONT data is to trim low quality 
reads and check the overall quality of the data.

```
/software/Porechop-v0.2.4/porechop-runner.py -i /path/to/directory/with/reads/ -o output..fastq --discard_middle
```

