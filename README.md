# CODEML & HYPHY PIPELINE
This repository contains scripts for generating single copy alignments from orthogroups, and running the branch-site test in codeml or RELAX in hyphy. 
The scripts were made for study systems that have very few single-copy orthogroups (not uncommon in plants), and orthogroups that contain multiple gene copies/transcripts per species are divided into subsets based on smallest genetic distance (an approximation). 

The approach for alignment generation and set-up for the branch-site test is described in Birkeland et al. 2020, Mol Biol Evol + Supplementary, with some minor changes  
https://academic.oup.com/mbe/article/37/7/2052/5804990?login=false

[embed]https://github.com/siribi/CODEML_HYPHY_PIPELINE/blob/main/images/Birkeland_et_al_2020_MBE_PipelineChart.pdf[/embed]

The scripts are written in bash, python and R for a slurm based computing cluster like this one: https://documentation.sigma2.no/jobs/job_scripts/saga_job_scripts.html

PARENT = path to your working directory <br />
SCRIPTS = path to the directory with the relevant scripts <br />

# PART 1. MAKING ALIGNMENTS FROM ORTHOGROUPS
#################################################################################### <br />

I use the following programs for making gene alignments: <br />
ORTHOFINDER <br />
EMBOSS (distmat algorithm) <br />
GUIDANCE2  <br />
PRANK <br />
You will also need Python 2 (I have not updated the scripts to newer Python versions)

NOTE ON GUIDANCE2: 
- GUIDANCE2 does not seem to have been updated since 2016: http://guidance.tau.ac.il/source 
- I have experienced problems with the MSA_parser.pm in GUIDANCE2 and contacted the developers regarding this. 
They sent me a new version, but I don't think it is implemented in the program (I had the exact same problem when I downloaded it years later). 
With the new MSA_parser.pm, the final alignment file is called: MSA.PRANK.aln.Sorted.With_Names. 
You can replace the file /cluster/home/siribi/nobackup/programs/guidance.v2.02/www/Guidance/../bioSequence_scripts_and_constants//MSA_parser.pm 
with the MSA_parser.pm found in the GUIDANCE2_bugfix directory. There is no need for new compilation/installation. <br />
- GUIDANCE2 sometimes errors out on predicted genes from genome cds, so I have run such files through Transdecoder first.

To make alignments with this set up, you first need to: 
1) have a set of cds and peptide files for the species you work with, and 
2) run OrthoFinder with the peptide files from your species set. Remember to include the multiple sequence alignment option (https://github.com/davidemms/OrthoFinder).    
NB: Note that GUIDANCE2, codeml and RELAX may be sensitive to gene sequences that are not complete or weirdly formatted (e.g., they are lacking start codons, contain internal stop codons etc). This is an issue that may apply to some non-model genomes. 

Formatting of fastas:
It can be smart to add a species abbreviation to all fasta entries, so that they can be searched and filtered on the basis of this name.
My scripts do some sorting and filtering on the basis of such abbreviations. 
Scripts to edit fasta headers can be found in the fasta_prep directory.

####################################################################################
**Part 1a. Preparation of input files and running the distmat algorithm in EMBOSS** 

Introduction: First part of the pipeline involves calculating protein distance between all genes in an orthogroup using the distmat algorithm 
in EMBOSS. To do this we need to prepare the input files using the scripts: The main script called "distmat_prep.sbatch" which needs 
the two accessory scripts: sort_gene_of_interest_v8.py and make_directories.sh

BEFORE starting ditstmat_prep.sbatch
1. Copy all alignments from the MultipleSequenceAlignments in the OrthoFinder results directory to a directory called e.g. orthofasta in the PARENT directory
2. Concatenate all the cds files so that we can search them for nucleotide sequences.
	cat *.cds > CladeB_Draba_concat.cds 
3. Change A) "path" and B) "pattern" in sort_gene_of_interest_v8.py (explanation in the file)
	This script starts the process of dividing the orthogroup fastas into 
	orthogroup subsets - ultimately resulting in alignments with just one gene copy 
	per species. The process is slightly complicated (and just an approximation), 
	but it starts with defining a species of interest in our dataset (here 
	"Poptre" for Populus tremula) and then dividing each orthogroup into subsets 
	based on the gene copies of that species. For instance, the script will make 
	three orthogroup subsets called OG0000334_Poptre_001.fasta 
	OG0000334_Poptre_002.fasta and OG0000334_Poptre_003.fasta if there are three 
	gene copies from "Poptre" in orthogroup OG0000334. To find which gene copies 
	from the other species that have the smallest genetic distance to each of these 
	gene copies, we place each gene copy together with the rest of the orthogroup 
	sequences and calculate protein distance. (Note that the species of interest 
	is always put last in the file to make file processing easier downstream).
3. Change A) $PARENT and B) $SCRIPTS in distmat_prep.sbatch 
4. Submit distmat_prep.sbatch

NB: distmat_prep.sbatch will divide all *sorted.aligned files into directories with 100 in each using "make_directories.sh". However, it differs from system to system how you should number these directories. The script gives each directory the name dir.1, dir.2, dir.3 etc., but some systems only reads dir.001, dir.002, dir.003)

4. After the distmat_prep.sbatch is finished, we are ready to do an array run of distmat! distance_worker.sbatch is copied to the relevant directory.
	Change $DISTMAT_OUT in the worker script called "distance_worker.sbatch" (i.e. the distmat directory)
	submit an array run of "distance_worker.sbatch" (Differs from system to system, so you might have to figure out how these work on your cluster. Some need both a submit and worker script, while some only use a worker script. On my current cluster I use the following command for e.g. 118 dir.*: 
	sbatch --array=1-118 distance_worker.sbatch 

#################################################################################### <br />
**Part 1b. Make files for GUIDANCE2**

Based on the results from the distmat algorithm from EMBOSS, next step is to build new files with only one gene copy for each species. The "make_files_for_guidance.sbatch" will run several helper scripts to process the output from distmat. Importantly it will pick out one gene copy per species based on the smalles genetic distance to the gene copies from the species of interest. Finally, it will assemble new single copy orthogroup fasta with all nucleotide sequences instead of protein sequences.

1. In  the "make_files_for_guidance.sbatch" script, change directory paths A) "PARENT" B) "SCRIPTS" and C) "DISTMAT". 
2. Edit the the grep_and_transpose_*.sh script to fit your dataset. This is the file that pick out one gene per species, so make sure all  your species are there. 
3. In GUIDANCE_prep.py Change A) "path", B) "paths", C) "fdna" and D) "pahti" (fnda is the concatenated cds file). This script assembles the nucleotide fasta files. 

#################################################################################### <br />
**Part 1c. Run Guidance** <br />
Previous studies have shown that a phylogeny aware aligner such as PRANK in combination with alignment filtering with GUIDANCE improves positive selection inference on simulated data (Jordan and Goldman 2012; Privman et al. 2012). 
With the help of GUIDANCE2, all alignments with a sequence score <0.6, as well as all columns scoring <0.93 (default), are removed from the data set. 
I have chosen to remove the entire alignment if one of the sequences fall below 0.6 just for simplicity/smoother running of CODEML (this happens in the next step - 1d).

1. Copy new fasta files from single_copy_orthofasta to new directory. Divide into directories (50 files in each). Run make_50_directories.sh or execute this command: 

	i=0; 
	for f in *; 
	do 
    		d=dir.$(printf %1d $((i/50+1))); 
    		mkdir -p $d; 
    		mv "$f" $d; 
    		let i++; 
	done

2. Make the output directories 
	mkdir FASTAS_Guidance_Edits MSA_Scores_Guidance Seq_Scores_Guidance

3. Run worker_guidance_fi.sbatch
	Note 1: This run requires a lot of space. Consider running it in the work directory. 
	Note 2: Notice that the main output of GUIDANCE2 gets this name: mv MSA.PRANK.Without_low_SP_Col.With_Names $FILEBASE.guidance.edit.fasta
	Note 3: See note about bug in GUIDANCE2 above

#################################################################################### <br />
**Part 1d. Removing alignments with bad sequence scores**

Run detect_lower.py in the Seq_Scores_Guidance, to get a list (badseqscores.txt) of alignments containing sequences with scores less than 0.6. 
Remember to load the right python module: 
	module load python2/2.7.10.gnu (Abel) 
	module load Python/2.7.15-fosscuda-2018b (Saga)
adding file extentions to badseqscores.txt
	sed -i 's/$/\.guidance.edit.fasta/' ./badseqscores.txt 
In directory with FASTAS_Guidance_Edits, move list of bad seq scores to new directory
	cat badseqscores.txt | xargs mv -t badseqscores/

If needed, remove alignments containing gaps and internal stop codons (this should actually be done beforehand though...)
Here I made a list of sequences in Arabis alpina and Cardamine hirsuta that contained Ns (and also one for internal STOP codons)
Use this file for the next steps (Remember to remove e.g. > with  sed -i 's/>//g')
Grepping for files with gaps and stop-codons and feeding them into a file:
	grep -f Aalp_sequences_with_gaps_EDIT2.file *.fasta >> Guidance_edit_files_with_Aalp.file
Delete everything after : in Guidance_edit_files_with_Aalp.file
	sed 's/:.*//g' Guidance_edit_files_with_Aalp.file >> List_of_files_with_Aalp_gaps.file
	cat List_of_files_with_Aalp_gaps.file | xargs mv -t Aalp_GAPS/



# PART 2. RUNNING THE BRANCH-SITE TEST IN CODEML
#################################################################################### <br />
**Part 2a. Preparing to run codeml**

In directory where you want to run Codeml, copy over FASTAS_Guidance_Edits and the scripts directory of interest
All the scripts directories should be ok now, but if needed it's easy to transform codeml_prep.sbatch for Cochlearia with a series of sed -i 's///g' commands!
e.g. sed -i 's/Bulk_run_Draba.tree/Bulk_run_Cochlearia.tree/g' codeml_prep_Cochlearia.sbatch
Remember to change PARENT and SCRIPTS directories in codeml_prep.sbatch
The script codeml_prep.sbatch makes all necessary files and directories to run all models of CODEML.
To run the controls, you only need to change the tree file :)  

#################################################################################### <br />
**Part 2b. Running the branch-site test in codeml**

Run each model! 
You need to cd into the model-run directory and change the worker and submit script. 
In the worker script copy the path of the mlc-output directory (one directory up), and insert the path to the model directory (where you have all the directories).
In the submit script you need to change the number of directories. 
