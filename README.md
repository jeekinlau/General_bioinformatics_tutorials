# General_bioinformatics_tutorials
The general idea of this tutorial is to get you familar with general bioinformatics. Much of the examples are pertinent to the Texas A&M University (TAMU) High Performance Research Compututing (HPRC) resources.

## Table of contents
1. [Introduction](#Basics)
2. [Downloading Data](#downloading-data)
3. [Data QC](#quality-control)

## Basics <a name="Basics"></a>
### Logging in 
There are a few easy ways to access the TAMU HPRC.    
**Important:You must be physically on campus or logged into TAMU VPN  to access the clusters**

1. Access using SSH. Open terminal on linux or macOS, cmd or powershell on Windows.     
```
# terminal command
ssh NetID@terra.tamu.edu
ssh NetID@grace.tamu.edu
ssh NetID@faster.hprc.edu
```

2. You can use MobaXterm on windows to access the terminal easier just reduces having to type in password and it has a file manager like sidebar that allows for FTP.

3. You can also use the HPRC Open OnDemand portal https://portal.hprc.tamu.edu/

Side note: I like to use a combination of all three depending on what you are doing and ease of procedure.


Let's login to GRACE because I do most of my work on GRACE.

Moving around in linux
```
pwd    # this tells you your current directory
cd     # change directory
ls     # list what is in current directory add -l to see list view
mv     # move
mkdir  # make directory
rm     # remove BE CAREFUL with this one

cd ..  # go one folder up
ls ..  # look at contents one folder structure up

```




In order to load a program we need for analysis or any step we should run the following lines.

```
module spider R 
ml spider R

module spider R/4.2.0
```
## Downloading data <a name="downloading-data"></a>
First we need to download all the files. Now the tool to do this is called fastq-dump which is part of the sra-tools from ncbi. What this does is downloads the accession numbers that you tell it to from NCBI. This code below will download SRR7060177 then compress the fasta file to .gz format.
```
fastq-dump SRR7060177 --gzip
```

Now you dont want to have to type all this in so we can utilize many cores and many instance of the HPRC to do this relatively quickly. However we cannot do this from the SSH login window or you will get a nasty email and you could be put in HPRC timeout for a week or two. So we need to create what is called a slurm job. All a slurm Job is file that asks the job scheduler for certain resources you want (i.e. cores, ram, runtime, etc.). Otherwise the rest of the slurm file is just a bash script. The specs you can ask for are located here https://hprc.tamu.edu/wiki/Grace:Batch


Below we are going to download the dataset from the paper https://acsess.onlinelibrary.wiley.com/doi/full/10.3835/plantgenome2018.02.0010

```
#!/bin/bash
#SBATCH --export=NONE
#SBATCH --job-name=download_NCBI_sequences
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=350G
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL

ml GCC/10.2.0  OpenMPI/4.0.5 SRA-Toolkit/2.10.9 parallel/20210322 WebProxy
#mkdir $SCRATCH/tomatoes
#mkdir $SCRATCH/tomatoes/rawdata
cd $SCRATCH/tomatoes/rawdata

parallel -j 45 fastq-dump {} --gzip ::: SRR7060177	SRR7060178	SRR7060179  SRR7060180	SRR7060181	SRR7060182	SRR7060183	SRR7060184	SRR7060185	SRR7060186	SRR7060187	SRR7060188	SRR7060189	SRR7060190	SRR7060191	SRR7060192	SRR7060193	SRR7060194	SRR7060195	SRR7060196	SRR7060197	SRR7060198	SRR7060199	SRR7060200	SRR7060201	SRR7060202	SRR7060203	SRR7060204	SRR7060205	SRR7060206	SRR7060207	SRR7060208	SRR7060209	SRR7060210	SRR7060211	SRR7060212	SRR7060213	SRR7060214	SRR7060215	SRR7060216	SRR7060217	SRR7060218	SRR7060219	SRR7060220	SRR7060221	SRR7060222	SRR7060223	SRR7060224	SRR7060225	SRR7060226	SRR7060227	SRR7060228	SRR7060229	SRR7060230	SRR7060231	SRR7060232	SRR7060233	SRR7060234	SRR7060235	SRR7060236	SRR7060237	SRR7060238	SRR7060239	SRR7060240	SRR7060241	SRR7060242	SRR7060243	SRR7060244	SRR7060245	SRR7060246	SRR7060247	SRR7060248	SRR7060249	SRR7060250	SRR7060251	SRR7060252	SRR7060253	SRR7060254	SRR7060255	SRR7060256	SRR7060257	SRR7060258	SRR7060259	SRR7060260	SRR7060261	SRR7060262	SRR7060263	SRR7060264	SRR7060265	SRR7060266	SRR7060267	SRR7060268	SRR7060269	SRR7060270	SRR7060271	SRR7060272	SRR7060273	SRR7060274	SRR7060275	SRR7060276	SRR7060277	SRR7060278	SRR7060279	SRR7060280	SRR7060281	SRR7060282	SRR7060283	SRR7060284	SRR7060285	SRR7060286	SRR7060287	SRR7060288	SRR7060289	SRR7060290	SRR7060291	SRR7060292	SRR7060293	SRR7060294	SRR7060295	SRR7060296	SRR7060297	SRR7060298	SRR7060299	SRR7060300	SRR7060301	SRR7060302	SRR7060303	SRR7060304	SRR7060305	SRR7060306	SRR7060307	SRR7060308	SRR7060309	SRR7060310	SRR7060311	SRR7060312	SRR7060313	SRR7060314	SRR7060315	SRR7060316	SRR7060317	SRR7060318	SRR7060319	SRR7060320	SRR7060321	SRR7060322	SRR7060323	SRR7060324	SRR7060325	SRR7060326

```

Copy this text into a file and name it `download.SLURM` then execute by going to the folder where the slurm file is and `sbatch download.SLURM`


THIS FAILED because the compute nodes are not connected to the internet. IN the module load section we will add WebProxy and it should work. Also comment out the two mkdir lines since the folders already exist.

## Data QC <a name="quality-control"></a>
Next we will use two software, fastcq and multqc to look at the quality of all the samples we just downloaded.

If we are just doing one sample 

```
fastqc SRR7060177.fasta.gz
```

This will produce an [HTML file](../resources/SRR7060177_fastqc.html)
