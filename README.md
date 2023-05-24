# WORKFLOW FOR GENERATING GENE ALIGNMENTS & RUNNING THE BRANCH-SITE TEST IN CODEML 
This repository contains scripts for generating single copy alignments from orthogroups, and running the branch-site test in codeml. 
The scripts were made for study systems that have very few single-copy orthogroups (not uncommon in plants). To generate enough single-copy alignments,  orthogroups that contain multiple gene copies/transcripts per species are divided into subsets based on smallest genetic distance (an approximation). 

The approach for alignment generation and set-up for the branch-site test (2a in the Figure below) is described in Birkeland et al. 2020, Mol Biol Evol + Supplementary, with some minor changes  
https://academic.oup.com/mbe/article/37/7/2052/5804990?login=false

![Pipelinechart from Birkeland et al. 2020](https://github.com/siribi/CODEML_HYPHY_PIPELINE/blob/main/images/FigureS2_Pipeline_chart.png)

The scripts are written in bash, python and R for a slurm based computing cluster like this one: https://documentation.sigma2.no/jobs/job_scripts/saga_job_scripts.html

**Get to know the branch-site test in CodeML:** <br />
CodeML is known to have a steep learning curve, but there are a lot of useful resources out there which can help you get a good idea of how the program works: <br />
- Check out the section on CodeML in the [PAML Frequently Asked Questions](http://abacus.gene.ucl.ac.uk/software/pamlFAQs.pdf). Here you can get advice on how many species are required, divergence levels etc. <br />
- [This tutorial](https://github.com/romainstuder/evosite3d/blob/update_tutorials/tutorials/branchsite.md) on CodeML by Romain Studer gives a good and short overview of the basics of the branch-site test. I discovered that I had the original tutorial with images in pdf version - you can find it in the [resources folder](https://github.com/siribi/CODEML_WORKFLOW/blob/main/resources/RomainStuder_tutorial_codeml_Evosite3D.pdf). <br />
- Álvarez-Carretero, Paschalia Kapli and Ziheng Yang (the developer of PAML) recently published [A Beginner’s Guide on the Use of PAML to Detect Positive Selection](https://academic.oup.com/mbe/article/40/4/msad041/7140562?login=true) in Molecular Biology and Evolution. <br />
- Jeffares et al. gives a nice introduction on how to use CodeML to calculate ω and identify positive selection using all genes in a genome. This tutorial also gives a good introduction to the different models of adaptive evolution implemented in CodeML: [A Beginners Guide to Estimating the Non-synonymous to Synonymous Rate Ratio of all Protein-Coding Genes in a Genome](https://www.researchgate.net/publication/268231340_A_Beginners_Guide_to_Estimating_the_Non-synonymous_to_Synonymous_Rate_Ratio_of_all_Protein-Coding_Genes_in_a_Genome) <br />
- Join the [PAML discussion group](https://groups.google.com/g/pamlsoftware) to ask questions, search for answers and follow discussions on CodeML/PAML <br />
- The original paper describing the branch-site test in codeml by Zhang et al. 2005 [Evaluation of an Improved Branch-Site Likelihood Method for Detecting Positive Selection at the Molecular Level](https://academic.oup.com/mbe/article/22/12/2472/1009544?login=true) <br />
- The branch-site test is based on the dn/ds ratio. Need a refresher on how dn/ds usually is calculated? You can find a [handout](https://github.com/siribi/CODEML_WORKFLOW/blob/main/resources/How_to_calculate_dNdS_ratios.pdf) from a precursor of this course: [Searching for Natural Selection: dN/dS](https://www.coursera.org/lecture/genetics-evolution/searching-for-natural-selection-dn-ds-s-6VocF) in the resources folder. <br />

**Note on convergence issues:** In this workflow, I run CodeML from four to six times per model with different initial parameter values to overcome convergence problems. The run with the highest likelihood score is used in subsequent analyzes. You can read more about this issue in [Wong et al. 2004](https://academic.oup.com/genetics/article/168/2/1041/6059620?login=true).  <br />

Before you can run the branch-site test in CodeML you need to generate gene alignments and a species tree or a gene tree for each alignment. 

# PART 1. MAKING ALIGNMENTS FROM ORTHOGROUPS
#################################################################################### <br />

I use the following programs for making gene alignments: <br />
[ORTHOFINDER](https://github.com/davidemms/OrthoFinder) <br />
[MAFFT](https://mafft.cbrc.jp/alignment/software/) <br />
[EMBOSS](https://emboss.sourceforge.net/download/#Gettingstarted) (distmat algorithm) <br />
[GUIDANCE2](https://github.com/anzaika/guidance)  <br />
[PRANK](http://wasabiapp.org/software/prank/) <br />
[PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) <br />
You will also need Python 2 (I have not updated the scripts to newer Python versions)

NOTE ON GUIDANCE2: 
- GUIDANCE2 does not seem to have been updated since 2016: http://guidance.tau.ac.il/source 
- I have experienced problems with the MSA_parser.pm in GUIDANCE2 and contacted the developers regarding this. 
They sent me a new version, but I don't think it is implemented in the program (I had the exact same problem when I downloaded it years later). 
With the new MSA_parser.pm, the final alignment file is called: MSA.PRANK.aln.Sorted.With_Names. 
You can replace the file 
```
/cluster/home/siribi/nobackup/programs/guidance.v2.02/www/Guidance/../bioSequence_scripts_and_constants/MSA_parser.pm 
```
with the MSA_parser.pm found in the [GUIDANCE2_bugfix directory](https://github.com/siribi/CODEML_WORKFLOW/tree/main/GUIDANCE2_bugfix). There is no need for new compilation/installation. <br />

To make alignments with this set up, you first need to: 
1) have a set of cds and peptide files for the species you work with, and 
2) run OrthoFinder with the peptide files from your species set. You have the possibility to do multiple sequence alignment with MAFFT within Orthofinder, or you could use the files in "Orthogroup_Sequences" and align them with MAFFT afterwards.     

NB: Note that this workflow may be sensitive to gene sequences that are not complete or weirdly formatted (e.g., they are lacking start codons, contain internal stop codons or are not in frame). It can be a good idea to do a sanity check first, especially with some new non-model genomes. In certain cases I had to run cds files through Transdecoder first.

Formatting of fastas: <br />
It can be smart to add a species abbreviation to all fasta entries, so that they can be searched and filtered on the basis of this name. My scripts do some sorting and filtering on the basis of such species abbreviations. Scripts to edit fasta headers can be found in the [fasta_prep directory](https://github.com/siribi/CODEML_WORKFLOW/tree/main/scripts/fasta_prep).

In the scripts I use the following terminology: <br />
PARENT = path to your working directory <br />
SCRIPTS = path to the directory with the relevant scripts <br />

#################################################################################### <br />
**Part 1a. Preparation of input files and running the distmat algorithm in EMBOSS** 

Introduction: First part of the pipeline involves calculating protein distance between all genes in an orthogroup using the distmat algorithm 
in EMBOSS. To do this we need to prepare the input files using the scripts: The main script called [distmat_prep.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distmat_prep.sbatch) which needs 
the four accessory scripts: 1) [if_else_mafft.sh](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/if_else_mafft.sh) (optional), 2) [sort_gene_of_interest_v8.py](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/sort_gene_of_interest_v8.py), 3) [make_100files_directories.sh](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/make_100files_directories.sh) and 4) [distance_worker.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distance_worker.sbatch) 

BEFORE starting [distmat_prep.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distmat_prep.sbatch)
1. Either i) copy unaligned sequences from **OrthoFinder/Results_XXXXX/Orthogroup_Sequences/** to a directory called **orthofasta/** for alignment with MAFFT (should end with .fa), or ii) copy all alignments from the **OrthoFinder/Results_XXXXX/MultipleSequenceAlignments/** directory in OrthoFinder to a directory called **orthofasta/** in the PARENT directory (should end with .aligned.fa). If you already have aligned sequences you can hashtag out the run of [if_else_mafft.sh](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/if_else_mafft.sh) (protein alignment) in [distmat_prep.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distmat_prep.sbatch). <br />

2. Change A) "path" and B) "pattern" in [sort_gene_of_interest_v8.py](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/sort_gene_of_interest_v8.py) <br />

	**path = Path to orthofasta/ (containing the .aligned.fa you will get after running mafft)** <br />
	**pattern = Species abbreviation to search for. Protein distances will be calculated in relation to this species (see below)** <br />
	
	This script starts the process of dividing the orthogroup fastas into 
	orthogroup subsets - ultimately resulting in alignments with just one gene copy 
	per species. The process is slightly complicated (and just an approximation), 
	but it starts with defining a species of interest in our dataset (for instance 
	"Draniv" for Draba nivalis) and then dividing each orthogroup into subsets 
	based on the gene copies of that species. For instance, the script will make 
	three orthogroup subsets called OG0000334_Draniv_000.fasta 
	OG0000334_Draniv_001.fasta and OG0000334_Draniv_002.fasta if there are three 
	gene copies from "Draniv" in orthogroup OG0000334. To find which gene copies 
	from the other species that have the smallest genetic distance to each of these 
	gene copies, we place each gene copy together with the rest of the orthogroup 
	sequences and calculate protein distance. (Note that the species of interest 
	is always put last in the file to make file processing easier downstream).

3. Change A) **$PARENT** and B) **$SCRIPTS** in [distmat_prep.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distmat_prep.sbatch) 

4. Submit [distmat_prep.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distmat_prep.sbatch) 

NB: distmat_prep.sbatch will divide all *sorted.aligned files into directories with 100 in each using "make_100files_directories.sh". However, it differs from system to system how you should number these directories. The script gives each directory the name dir.1, dir.2, dir.3 etc., but some systems only reads dir.001, dir.002, dir.003)

After [distmat_prep.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distmat_prep.sbatch) is finished: <br />
Now we should be ready to do an array run of distmat with this script: [distance_worker.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distance_worker.sbatch), which should have been copied to **$PARENT/orthofasta/sorted_aligned/**. <br />

1. Change **$DISTMAT_OUT** in [distance_worker.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distance_worker.sbatch). This will be the main output directory for the array run. It should already have been made by the previous script [distmat_prep.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distmat_prep.sbatch) <br />

2. Submit an array run of [distance_worker.sbatch](https://github.com/siribi/CODEML_WORKFLOW/blob/main/scripts/alignments/distance_worker.sbatch) (Differs from system to system, so you might have to figure out how these work on your cluster. Some need both a submit and worker script, while some only use a worker script. On my current cluster I use the following command for e.g. 118 dir.*: 
```
sbatch --array=1-118 distance_worker.sbatch 
```

#################################################################################### <br />
**Part 1b. Make files for GUIDANCE2**

Based on the results from the distmat algorithm from EMBOSS, next step is to build new files with only one gene copy for each species. The "make_files_for_guidance.sbatch" will run several helper scripts to process the output from distmat. Importantly it will pick out one gene copy per species based on the smalles genetic distance to the gene copies from the species of interest. Finally, it will assemble new single copy orthogroup fasta with all nucleotide sequences instead of protein sequences.

1. Concatenate all the cds files so that we can search them for nucleotide sequences (we want codon alignments), e.g.

```
cat *.cds > CladeB_Draba_concat.cds
```

2. In  the "make_files_for_guidance.sbatch" script, change directory paths A) "PARENT" B) "SCRIPTS" and C) "DISTMAT". 

3. Edit the the grep_and_transpose_*.sh script to fit your dataset. This is the file that pick out one gene per species, so make sure all  your species    are there. 

4. In GUIDANCE_prep.py Change A) "path", B) "paths", C) "fdna" and D) "pahti" (fnda is the concatenated cds file). This script assembles the nucleotide fasta files. 

#################################################################################### <br />
**Part 1c. Run Guidance** <br />
Previous studies have shown that a phylogeny aware aligner such as PRANK in combination with alignment filtering with GUIDANCE improves positive selection inference on simulated data (Jordan and Goldman 2012; Privman et al. 2012). 
With the help of GUIDANCE2, all alignments with a sequence score <0.6, as well as all columns scoring <0.93 (default), are removed from the data set. 
I have chosen to remove the entire alignment if one of the sequences fall below 0.6 just for simplicity/smoother running of CODEML (this happens in the next step - 1d).

1. Copy new fasta files from single_copy_orthofasta to new directory. Divide into directories (50 files in each). Run make_50_directories.sh or execute this command: 

```
i=0; 
for f in *; 
do 
d=dir.$(printf %1d $((i/50+1))); 
mkdir -p $d; 
mv "$f" $d; 
let i++; 
done
```

2. Make the output directories

```
mkdir FASTAS_Guidance_Edits MSA_Scores_Guidance Seq_Scores_Guidance
```

3. Run worker_guidance_fi.sbatch
- This run requires a lot of space and also takes a lot of computing hours <br />
- Notice that the main output of GUIDANCE2 is renamed: 
```
mv MSA.PRANK.Without_low_SP_Col.With_Names $FILEBASE.guidance.edit.fasta
```

#################################################################################### <br />
**Part 1d. Removing alignments with bad sequence scores**

Run detect_lower.py in the Seq_Scores_Guidance, to get a list (badseqscores.txt) of alignments containing sequences with scores less than 0.6. 

Remember to load the right python module. Here are examples of python2 modules that I have used on two different clusters (named Abel and Saga). 
```
module load python2/2.7.10.gnu (Abel) 
module load Python/2.7.15-fosscuda-2018b (Saga)
```

adding file extentions to badseqscores.txt
```
sed -i 's/$/\.guidance.edit.fasta/' ./badseqscores.txt 
```

In directory with FASTAS_Guidance_Edits, move list of bad seq scores to new directory
```
cat badseqscores.txt | xargs mv -t badseqscores/
```

If needed, remove alignments containing gaps and internal stop codons (this should actually be done beforehand though...)
Here I made a list of sequences in Arabis alpina and Cardamine hirsuta that contained Ns (and also one for internal STOP codons)
Use this file for the next steps (Remember to remove e.g. > with  sed -i 's/>//g')
Grepping for files with gaps and stop-codons and feeding them into a file:
```
grep -f Aalp_sequences_with_gaps_EDIT2.file *.fasta >> Guidance_edit_files_with_Aalp.file
```
Delete everything after : in Guidance_edit_files_with_Aalp.file
```
sed 's/:.*//g' Guidance_edit_files_with_Aalp.file >> List_of_files_with_Aalp_gaps.file
```
```
cat List_of_files_with_Aalp_gaps.file | xargs mv -t Aalp_GAPS/
```


# PART 2. RUNNING THE BRANCH-SITE TEST IN CODEML
#################################################################################### <br />

**Part 2a. Preparing to run codeml**

In directory where you want to run Codeml, copy over FASTAS_Guidance_Edits and the scripts directory of interest
All the scripts directories should be ok now, but if needed it's easy to transform codeml_prep.sbatch for Cochlearia with a series of sed -i 's///g' commands! e.g. 
```
sed -i 's/Bulk_run_Draba.tree/Bulk_run_Cochlearia.tree/g' codeml_prep_Cochlearia.sbatch
```
Remember to change PARENT and SCRIPTS directories in codeml_prep.sbatch
The script codeml_prep.sbatch makes all necessary files and directories to run all models of CODEML.
To run the controls, you only need to change the tree file :)  

#################################################################################### <br />
**Part 2b. Running the branch-site test in codeml**

Run each model! 
You need to cd into the model-run directory and change the worker and submit script. 
In the worker script copy the path of the mlc-output directory (one directory up), and insert the path to the model directory (where you have all the directories).
In the submit script you need to change the number of directories. 
