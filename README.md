# Iso-seq data pipeline

## PacBio concepts
PacBio SMRT sequencing operates within a silicon chip (**SMRTcell**) fabricated to contain a large number of microscopic holes (**ZMWs**, or zero-mode waveguides), each assigned a hole number.

Within a ZMW, PacBio SMRT sequencing is performed on a circularized molecule called a **SMRTbell**. The SMRTbell, depicted below, consists of:

* the customer’s double-stranded DNA insert (I)
* (optional) double-stranded DNA barcodes (sequences BL, BR) used for multiplexing DNA samples. While the barcodes are optional, they must be present at both ends if present at all.
* SMRTbell adapters (sequences AL, AR), each consisting of a double stranded stem and a single-stranded hairpin loop. 
  
![image](https://user-images.githubusercontent.com/31697487/135611262-770421a4-0538-4738-865b-750f25bae351.png)
  
Specifically the raw data or each SMRT-cell will be in files named  ```.subreads.bam```,  ```.subreads.xml```, and  ```.subreads.pbi```.

The ```.bam data``` can be converted to fastq or fasta files with bamtools.


### What is Circular Consensus Sequence or ccs?
ccs combines multiple subreads of the same SMRTbell molecule using a statistical model to produce one highly accurate consensus sequence, also called a HiFi read.
Single Molecule, Real-Time (SMRT) Sequencing technology has evolved to a different type of long read, known as highly accurate long reads, or **HiFi reads**. PacBio is the only sequencing technology to offer HiFi reads that provide accuracy of >99.9%, on par with short reads and Sanger sequencing.
 
 
![alt text](https://ccs.how/img/generate-hifi.png)


## Step 0: Input
For each SMRT cell a duck.subreads.bam is needed for processing.

## Step 1: Circular Consensus Sequence (ccs) calling 

Each sequencing run is processed by ccs to generate one representative circular consensus sequence (CCS) for each ZMW.
* **Input**: Subreads from a single movie in PacBio BAM format (```.subreads.bam```).
* **Output**: Consensus reads in a format inferred from the file extension: unaligned BAM (```.bam```); bgzipped FASTQ (```.fastq.gz```); 
or SMRT Link XML (```.consensusreadset.xml```) which also generates a corresponding unaligned BAM file.

more info at [how does ccs work](https://ccs.how/how-does-ccs-work.html)
```bash
ccs duck.subreads.bam duck.ccs.bam --min-rq 0.9
```
### Parallelize by chunk

For this, the ```.subreads.bam``` file must accompanied by a ```.pbi``` file. To generate the index ```subreads.bam.pbi```, use ```pbindex```, which can be installed with ```conda install pbbam```.
```bash
pbindex duck.subreads.bam
```
An example workflow, all ccs invocations can run simultaneously:
```bash
ccs duck.subreads.bam duck.ccs.1.bam --chunk 1/10 -j <THREADS>
ccs duck.subreads.bam duck.ccs.2.bam --chunk 2/10 -j <THREADS>
...
ccs duck.subreads.bam duck.ccs.10.bam --chunk 10/10 -j <THREADS>
```
Merge chunks with pbmerge and index with pbindex
```bash
pbmerge -o duck.ccs.bam duck.ccs.*.bam
pbindex duck.ccs.bam
 ```
or use samtools
```bash
samtools merge -@8 duck.ccs.bam duck.ccs.*.bam
```
## Step 2: Primer removal and demultiplexing 
Removal of primers and identification of barcodes is performed using lima that offers a specialized ```--isoseq``` mode. 
If there are more than two sequences in your primer.fasta file or better said more than one pair of 5' and 3' primers.
Use lima with ```--peek-guess``` to remove spurious false positive signal. Lima will remove unwanted combinations and orient sequences to 5' → 3' orientation.
```bash
lima duck.ccs.bam barcoded_primers.fasta duck.fl.bam --isoseq --peek-guess
```
Output files will be called according to their primer pair. Example for single sample libraries: ```duck.fl.NEB_5p--NEB_Clontech_3p.bam```

## Step 3a: Refine
Data now contains full-length (FL) reads, but still needs to be refined by:
* Trimming of poly(A) tails.
* Rapid concatemer identification and removal.

So, input and output for this step should be:
* **Input**: The input file for refine is one demultiplexed CCS file with full-length reads and the primer fasta file.
* **Output**: The output files of refine contain *Full Length Non-Chimeric (FLNC)* reads.
```bash
isoseq refine duck.NEB_5p--NEB_Clontech_3p.fl.bam primers.fasta duck.flnc.bam
```

If your sample has poly(A) tails, use ```--require-polya```. This filters for FL reads that have a poly(A) tail with at least 20 base pairs (```--min-polya-length```) and removes identified tail:
```bash
isoseq refine duck.NEB_5p--NEB_Clontech_3p.fl.bamduck.flnc.bam --require-polya
```

## Step 3b: Merge SMRT Cells (optional)
If you used more than one SMRT cells, list all of your ```<duck>.flnc.bam``` in one ```flnc.fofn```, a file of filenames:
```bash
ls duck*.flnc.bam duck*.flnc.bam duck*.flnc.bam > flnc.fofn
 ```

## Step 4: Clustering (optional)
This step takes the FLNC reads and clusters the transcript sequences by similarity. It then makes a multiple alignment of each cluster and performs error correction using this alignment. This step is optional and if you are interested in allele specific expression/transcript phasing, you should not run this step as it can removed allele specific sequence variation. Give this step as many cores as possible

* **Input** The input file for cluster is one FLNC file:
```<duck>.flnc.bam``` or ```flnc.fofn```

* **Output** The following output files of cluster contain polished isoforms:

   *  ```duck.bam ```
   *  ```duck.hq.fasta.gz ```  with predicted accuracy ≥ 0.99
   *  ```duck.lq.fasta.gz ```  with predicted accuracy < 0.99
   *  ```duck.bam.pbi ```
   *  ```duck.transcriptset.xml ```
 
```bash
isoseq cluster flnc.fofn clustered.bam --verbose --use-qvs
  ```
   
more info of [Isoseq3 pipeline](https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md)

## Step 4: Mapping reads against genome with Minimap2
This step can take inputs from either the FLNC reads or the Polished reads. Minimap2 maps the transcript sequences onto a genome assembly.
```bash
minimap2 -ax splice -uf -C5 genome_duck_ref.fa duck_iso-seq.fq > duck_aln.sam
```
Option ```-C5``` reduces the penalty on non-canonical splicing sites. It helps to align such sites correctly for data with low error rate such as Iso-seq reads and traditional cDNAs. On this example, minimap2 makes one junction error. Applying ```--splice-flank=no``` fixes this alignment error.
splice:hq?
  
more info in [Minimap2](https://github.com/lh3/minimap2)

Note that the command line above is optimized for the final Iso-seq reads. PacBio's Iso-seq pipeline produces intermediate sequences at varying quality. For example, some intermediate reads are not stranded. For these reads, option -uf will lead to more errors. Please revise the minimap2 command line accordingly.

## Step 5: Collapse 
TAMA Collapse is a tool that allows you to collapse redundant transcript models in your Iso-Seq data.
This step takes a bam/sam file from the transcript mapping and collapses redundant transcripts based on genomic location.

Also provides the following information: 
* Source information for each predicted feature
* Variation calling
* Genomic poly-A detection 
* Strand ambiguity

```bash
python tama_collapse.py -s duck_aln.sam -f genome_duck_ref.fa -p prefix -x capped
```
more info of [TAMA Collapse](https://github.com/GenomeRIK/tama/wiki/Tama-Collapse)

## Step 5: Merge
Merging is the process of combining multiple transcriptomes. For instance, if you have Iso-Seq data for different tissue types you might want to process them separately and then combine them at the end to use as a transcriptome for downstream analysis. However, the act of merging transcriptomes is non-trivial with respect to choosing what criteria to use to merge transcript models. You probably would also like to keep a record of which models from your merged transcriptome came from which source. 
```bash
python tama_merge.py -f filelist.txt -p merged_annos
```
NOTE: If you do not see "TAMA Merge has completed successfully!" as the last line of the terminal output, then TAMA Merge has not completed and there are likely issues with the input files.

more info of [TAMA Merge](https://github.com/GenomeRIK/tama/wiki/Tama-Merge)
