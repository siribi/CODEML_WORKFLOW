# CODEML_HYPHY_PIPELINE
Scripts for generating alignments from orthogroups and running the branch-site test in codeml and RELAX in hyphy

The scripts are written in bash, python and R for a slurm based computing cluster like this one: https://documentation.sigma2.no/jobs/job_scripts/saga_job_scripts.html

PART 1. MAKING ALIGNMENTS FROM ORTHOGROUPS
#################################################################################### <br />

I use the following programs for making gene alignments: <br />
EMBOSS (distmat algorithm) <br />


To make alignments with this set up, you first need to: 
1) have a set of cds and peptide files for the species you work with, and 
2) run OrthoFinder with the peptide files from your species set (https://github.com/davidemms/OrthoFinder).    

SB: WRITE SOMETHING ABOUT FILE FORMATTING OF FASTAS: 

Note: I have experienced that Guidance errors out on predicted genes from genome cds, so I have run such files through Transdecoder first. 

####################################################################################
Part 1a: Preparation of input files and running the distmat algorithm in EMBOSS: 

Introduction: First part of the pipeline involves calculating protein distance between all genes in an orthogroup using the distmat algorithm in EMBOSS. To do this we need to prepare the input files using three scripts: 1) The main script called "distmat_prep.sbatch" which needs the following python script: sort_gene_of_interest_v8.py starts the process of dividing the orthogroup fastas into orthogroup subsets - ultimately resulting in alignments with just one gene copy per species. The process is slightly complicated (and just an approximation), but it starts with defining a species of interest in our dataset (here "Poptre" for Populus tremula) and then dividing each orthogroup into subsets based on the gene copies of that species. For instance, the script will make three orthogroup subsets called OG0000334_Poptre_001.fasta OG0000334_Poptre_002.fasta and OG0000334_Poptre_003.fasta if there are three gene copies from "Poptre" in orthogroup OG0000334. To find which gene copies from the other species that have the smallest genetic distance to each of these gene copies, we place each gene copy together with the rest of the orthogroup sequences and calculate protein distance. (Note that the species of interest is always put last in the file to make file processing easier downstream).

## sort_gene_of_interest_v8.py ##
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

Starting ditstmat_prep.sbatch: 
First we need to cocatenated the peptide fastas used by orthofinder and a set of files with gene identifiers made from the Orthogroups.tsv file, and then we run the three scripts "distmat_prep.sbatch", "MAFFTprep_allOrthogroups.py" and "sort_gene_of_interest_v8.py".  
1. Making each line of the Orthogroups.txt into a separate file called an index file:
	Copy Orthofasta.tsv and the appropriate scripts_2a* library for the species to the PARENT directory
	Change $PARENT path in make_indexfiles_from_orthogroups.sh and run it  

2. Concatenate fastas (both amino acid and cds files) so that we can search them for gene sequences:
	make one directory with all the cds and one directory with all the aa files. Concatenate them into one file. 
	cat *.fasta > CladeB_Draba_concat.pep 

3. Running distmat_prep.sbatch (in the scripts directory) ###NOTE TO SELF: CHECK THAT NUMBERING CORRESPONDS TO SCRIPTS AND CHANGE TO SMALL LETTERS OR i) ii) iii) etc. 
	distmat_prep.sbatch: Change	A) $PARENT	and B) $SCRIPTS 
	MAFFTprep_allOrthogroups.py: Change A) "path", B) "paths", C) "fdna" and D) "pahti" (explanation in the file)
	sort_gene_of_interest_v8.py: Change A) "path" and B) "pattern" in sort_gene_of_interest_v8.py (explanation in the file)
	NB: distmat_prep.sbatch will divide all *sorted.aligned files into directories with 100 in each using "make_directories.sh". However, it differs from system to system how you should number these directories. The script gives each directory the name dir.1, dir.2, dir.3 etc., but some systems only reads dir.001, dir.002, dir.003)

4. After the distmat_prep.sbatch is finished, we are ready to do an array run of distmat! 
	Change $DISTMAT_OUT in the worker script called "distance_worker.sbatch" (i.e. the distmat directory)
	submit an array run of "distance_worker.sbatch" (Differs from system to system, so you might have to figure out how these work on your cluster. Some need both a submit and worker script, while some only use a worker script. On my current cluster I use the following command for e.g. 118 dir.*: 
	sbatch --array=1-118 distance_worker.sbatch 

######################################################################################################################

## B: Make files for GUIDANCE2

Based on the results from the distmat algorithm from EMBOSS, next step is to build new files with only one gene copy for each species. The "make_files_for_guidance.sbatch" will run several helper scripts to process the output from distmat. Importantly it will pick out one gene copy per species based on the smalles genetic distance to the gene copies from the species of interest. Finally, it will assemble new single copy orthogroup fasta with all nucleotide sequences instead of protein sequences.

1. In  the "make_files_for_guidance.sbatch" script, change directory paths A) "PARENT" B) "SCRIPTS" and C) "DISTMAT". 
2. Edit the the grep_and_transpose_*.sh script to fit your dataset. This is the file that pick out one gene per species, so make sure all  your species are there. 
3. In GUIDANCE_prep.py Change A) "path", B) "paths", C) "fdna" and D) "pahti" (fnda is the concatenated cds file). This script assembles the nucleotide fasta files. 

######################################################################################################################

## C: Run Guidance

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
	Note: This run requires a lot of space. Consider running it in the work directory. 
	NB: I had problems with the MSA_parser.pm in GUIDANCE2 and contacted the developers regarding this. They sent me a new version, but I don't think it is implemented in the program yet (I had the exact same problem when I downloaded it years later). With the new MSA_parser.pm, the final alignment file is called: MSA.PRANK.aln.Sorted.With_Names. 

######################################################################################################################

## D: Prepping for codeml after guidance run: 

Run detect_lower.py in the Seq_Scores_Guidance, to get a list (badseqscores.txt) of alignments containing sequences with scores less than 0.6. #module load python2/2.7.10.gnu (Abel) or Python/2.7.15-fosscuda-2018b (Saga)
#adding file extentions to badseqscores.txt
	sed -i 's/$/\.guidance.edit.fasta/' ./badseqscores.txt 
#In directory with FASTAS_Guidance_Edits, move list of bad seq scores to new directory
cat badseqscores.txt | xargs mv -t badseqscores/

### If needed, remove alignment containing gaps and internal stop codons (this should actually be done beforehand though...)
#Here I made a list of sequences in Arabis alpina and Cardamine hirsuta that contained Ns (and also one for internal STOP codons)
#Use this file for the next steps (Remember to remove e.g. > with  sed -i 's/>//g')
#Grepping for files with gaps and stop-codons and feeding them into a file:
grep -f Aalp_sequences_with_gaps_EDIT2.file *.fasta >> Guidance_edit_files_with_Aalp.file
# Delete everything after : in Guidance_edit_files_with_Aalp.file
sed 's/:.*//g' Guidance_edit_files_with_Aalp.file >> List_of_files_with_Aalp_gaps.file
cat List_of_files_with_Aalp_gaps.file | xargs mv -t Aalp_GAPS/

### In directory where you want to run Codeml, copy over FASTAS_Guidance_Edits and the scripts directory of interest
### All the scripts directories should be ok now, but if needed it's easy to transform codeml_prep.sbatch for Cochlearia with a series of sed -i 's///g' commands!
#e.g. sed -i 's/Bulk_run_Draba.tree/Bulk_run_Cochlearia.tree/g' codeml_prep_Cochlearia.sbatch
### Remember to change PARENT and SCRIPTS directories in codeml_prep.sbatch
### The script codeml_prep.sbatch makes all necessary files and directories to run all models of CODEML.
### To run the controls, you only need to change the tree file :)  

######################################################################################################################

## E: Running codeml

###Run each model! 
###You need to cd into the model-run directory and change the worker and submit script. 
#In the worker script copy the path of the mlc-output directory (one directory up), and insert the path to the model directory (where you have all the directories).
#In the submit script you need to change the number of directories. 

