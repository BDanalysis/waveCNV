# waveCNV

## Installation
The following software must be installed on your machine:

Linux LST 14 or 16 <br/>
Python : tested with version 3.7<br/>
Python >= 3.7<br/>

### Python dependencies
* numpy == 1.21.2
* seaborn == 0.11.2
* pysam == 0.18.0
* matplotlib == 3.5.1

You can install the above package using the following command：

pip install numpy sklearn  pysam matplotlib


## Running
HBOS-CNV requires two input files, a bam file after sort and a reference folder,
the folder contains the reference sequence for each chromosome of the bam file.

Note: at present, the approach only supports single chromosome samples, like chr5, chrX...
The whole genomes can be detected by batch file
### runnig command
python HBOS-CNV.py [reference] [bamfile] [binSize] [chrom]

[reference]: the reference folder path

[bamfile]: a bam file after sort

[binSize]: the window size ('1000'by default)

[chrom]: Chromosome number.

### run the default example
python process2.py ./reference/ ./test.sort.bam 1000 ./output/ 21

Finally, the result files will be generated in the folder called './output/'.


