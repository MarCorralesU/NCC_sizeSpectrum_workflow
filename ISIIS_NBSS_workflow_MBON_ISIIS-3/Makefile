# Makefile
# Date: January 14th, 2025
#
# Needed to run shell commands correctly
.ONESHELL:

#
# Setup the conda environment since its lost when using the Makefile
# system for some reason.
#

# LINES 13-29 added on 10/7/2024, to see if shell version is detected automatically and the appropiate conda env is chosen
# Detect node architecture
ARCH := $(shell uname -m)

# Select the right conda env prefix and partition
#ifeq ($(ARCH),ppc64le)
#  CONDA_ENV_PREFIX := /home/zoo/corralma/miniconda3_ppc64le/envs/plankton
#  PARTITION := cqls_ppc64le
#  CONDA_EXE := /home/zoo/corralma/miniconda3_ppc64le/bin/conda
#else
#  CONDA_ENV_PREFIX := /home/zoo/corralma/miniconda3_x86/envs/plankton
#  PARTITION := all.q  # <-- set your x86 partition here
#  CONDA_EXE := /fs1/local/cqls/software/x86_64/miniforge3-24.7.1-0/bin/conda
#endif

# Point at a known conda executable (the site miniforge is fine)

# Use conda run instead of "conda activate"
PY := $(CONDA_EXE) run -p $(CONDA_ENV_PREFIX) python3



SHELL = /bin/bash
#CONDA_ENV=plankton
CONDA_ENV=;conda activate plankton
CONDA_ACTIVATE=source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
CONDA_ACTIVATE_array=eval "$(command conda shell.bash hook)" && conda activate plankton
#CONDA_ACTIVATE_ARRAY := eval "$$(command conda shell.bash hook)" && conda activate $(CONDA_ENV)
#CONDA_ACTIVATE_le = source /home/zoo/corralma/miniconda3_ppc64le/etc/profile.d/conda.sh 


# $(CONDA_ACTIVATE) $(CONDA_ENV); SLURM_Batch ....
# $(CONDA_ACTIVATE) different_env; SLURM_Batch ....

#
# Number of processors to use for all jobs
#
PROCS=1

#
# CPU/GPU device to use
#
DEVICE=0

#
# DIRPATH is the default path for everything and is set to the current folder. 
# Change the DIRPATH if needed variable below to match the folder structure.
# You should only have to provide the top level "some_experiment" 
# path and the READS and other folders are included. Read files can
# have any base file name (ie Read1a or anything you need) they must
# end in fastq. 
DIRPATH:=$(shell pwd)
SEG=`which segment`
SEG_NEW=`which segment_new`
STYLE=basic
PARTITION="cqls_ppc64le"
#PARTITION="cqls_ppc64le" #"ceoas" or "all.q"
#NODELIST="-N beech[0-1]"
#NODELIST="-N aspen[5-9]"
#NODELIST="-N kawashiro[01-08]"
#NODELIST="-N yatagarasu[01-12]"
QUEUE_OPTIONS=-q $(PARTITION) -P $(PROCS)
#QUEUE_OPTIONS=-q $(PARTITION) -P $(PROCS) $(NODELIST)
QUEUE_OPTIONS_batch=--partition=$(PARTITION_sbatch) --cpus-per-task=$(PROCS) $(NODELIST)
SEG_OPTIONS=-l 66 -r 23 -f -m 50
THRESHOLD=240

#
# You are welcome to change the resources below if needed but most
# should work automatically. 
#

#
# Make list of files in different directory
#
#FILELIST=`cat inputDirs.txt  | awk -F/ '{print $$9}'`
#FILELIST=`cat inputDirs.txt`
FILELIST=`cat inputTars.txt`
FILELIST_SUBSET= ` cat ./ancillary/tar_lists/inputTars_1.txt` #`cat ./ancillary/tar_lists/inputTars_1.txt`
TARFILEPATH= `sed '5100!d' ./ancillary/input_tars_MBON.txt`
FILELIST_merged=`cat merged_files.txt`
SEG_FILES=`ls -1 SEGMENTATION/ | head -10`
SEG_FILES_COUNT=`ls -1 SEGMENTATION/ | wc -l`

CSV_FILES=`ls -1 $(DIRPATH)/CSV/`


#
# Directory List For Checking or Creating
#
DIRECTORIES=DOWNLOADS CSV SEGMENTATION SCRIPTS data ancillary figures TMP


#
# Main stanza used to print when user types "make" and nothing else. 
#
all: mkdirs
	@echo ""
	@echo -e "\e[36;1;202m   This folder is use to process Video or Image\e[0m"
	@echo -e "\e[36;1;202m   generate the base folders needed for an analysis by typing \"make\"\e[0m"
	@echo ""
	@echo -e "    The SLURM_Batch opt : \e[32;1;202m$(QUEUE_OPTIONS)\e[0m"
	@echo ""
	@echo "    The TAR list:    "
	@for file in $(FILELIST); do \
	echo -e "                                        \e[32;1;202m$${file}\e[0m"; \
	done ;
	@echo ""
	@echo "    The SEG list:    "
	@for file in $(SEG_FILES); do \
	echo -e "                                        \e[32;1;202m$${file}\e[0m"; \
	done ;
	@echo ""
	@echo ""
	@echo "    The SEG Count:    "
	@echo -e "                                        \e[32;1;202m$(SEG_FILES_COUNT)\e[0m"
	@echo ""
	@echo ""
	@echo "    Please type one of the following....."
	@echo "    	make stats"
	@echo "    	make status"
	@echo "    	make mkdirs"
	@echo "         make get_tar_list"
	@echo "    	make process_tar"
	@echo "         make process_subset_tar"
	@echo "    	make process_tars"
	@echo "    	make process_tars_array"
	@echo "    	make process_tars_array_combineCSVs"
	@echo "    	make process_tars_array_combineCSVs_test_le"
	@echo "    	make process_dirs"
	@echo "    	make process_file"
	@echo "    	make process_generate_csv_file_list"
	@echo "         make process_combine_merged_size"
	@echo "         make process_bin_size_data"
	@echo "         make process_calculate_size_spectra"
	@echo "         make process_make_diagnostic_figs"
	@echo "         make figures"
	@echo "    	make vars"
	@echo "    	make docs"
	@echo "    	make clean"
	@echo "    	make realclean"
	@echo "    	make superclean"
	@echo ""

docs:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Documentation for Using this Makefile ++==>\e[0m"
	@echo ""
	@echo -e "\tThe basics are listed below and by default the"
	@echo -e "\t\"all\" recipe is run when a user types \"make\""
	@echo -e "\tby default. The all recipe will show some basic"
	@echo -e "\tstats on the data from the settings at the top of"
	@echo -e "\tthe file. Those setting should be edited by the"
	@echo -e "\tend user and changed to match the experiment or"
	@echo -e "\tinfratstructure."
	@echo -e ""
	@echo -e "\tSettings:"
	@echo -e "\t\t"
	@echo -e "\t\t"
	@echo ""
	@echo ""
	@echo ""
	@echo -e "\tCommand: make mkdirs"
	@echo -e "\t\tThis command will autogenerate the directories needed to do the analysis."
	@echo -e "\t\tThis command runs automatically when you type the \"make\" command with nothing else."
	@echo ""
	@echo ""
	@echo -e "\tCommand: make vars"
	@echo -e "\t\tThis command will show all the variables used to setup a run"
	@echo -e "\t\tas well as show you the current values for those variables. "
	@echo -e "\t\t"
	@echo -e "\t\tFull set of Variables:"
	@echo -e "\t\t\tPARTITION"
	@echo -e "\t\t\tNODELIST"
	@echo -e "\t\t\tQUEUE_OPTIONS"
	@echo -e "\t\t\tDEVICE"
	@echo -e "\t\t\tPROCS"
	@echo -e "\t\t\tDIRPATH"
	@echo ""
	@echo -e "\t\tAny of these can be changed at the command line by typing VARIABLE_NAME=VALUE"
	@echo -e "\t\t\tExample:"
	@echo -e "\t\t\t\tmake PARTITION=ceoas vars"
	@echo ""
	@echo ""
	@echo -e "\tCommand: make status"
	@echo -e "\t\tShow the status of the squeue command if your user has jobs running as well as"
	@echo -e "\t\tusing srun to login to the GPU node and output nvidia-smi."
	@echo -e "\t\t\tExample:"
	@echo -e "\t\t\t\tmake status"
	@echo -e "\t\t\t\tmake PARTITION=all.q NODELIST=cqls-gpu1 status"
	@echo ""
	@echo ""
	@echo -e "\tCommand: make clean"
	@echo -e "\t\tThis command will clean out the folders generated in the main folder with StdOut on them."
	@echo -e "\t\tThese files come from running the SLURM_Batch command and hold the run files for the runs"
	@echo -e "\t\tdone on the cluser including the batch script used to submit the job."
	@echo ""
	@echo ""
	@echo -e "\tCommand: make realclean"
	@echo -e "\t\tThis command will clean out the StdOut folders and anything under TRAINING/images,"
	@echo -e "\t\tTRAINING/labels, TRAINING/modules for the experiment to be reset to the begining."
	@echo -e "\t\tIt will not remove any backed up runs under TRAINING/run[#]."
	@echo ""
	@echo ""
	@echo -e "\tCommand: make superclean"
	@echo -e "\t\tThis command will take the folder back to the original Makefile and that's it folks...."
	@echo ""
	@echo ""

mkdirs:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Checking Directories: Making sure we have all the needed input and output directories ++==>\e[0m"
	@for dir in $(DIRECTORIES); do \
	! test -d $(DIRPATH)/$${dir} && echo -e "\e[32;3;202m          Making: $${dir}\e[0m";\
	test -d $(DIRPATH)/$${dir} || mkdir -p $(DIRPATH)/$${dir};\
	done;
	@echo ""

get_tar_list:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Let's make a list of the tarfiles to get the size data ++==>\e[0m"
	@test -d $(DIRPATH)/ancillary || mkdir -p $(DIRPATH)/ancillary
	@test -d $(DIRPATH)/ancillary/tar_lists || mkdir -p $(DIRPATH)/ancillary/tar_lists
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/0.1.0_get_tarfile_list.R -d /jet/home/corrales/shared/SPECTRA -s noSplit -P MEZCAL" -r get_tars

process_tar:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Processing a Tar File ++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -t $(TARFILEPATH) -v" -r ProcessFile-StdOut $(QUEUE_OPTIONS)

process_subset_tar:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files ++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@c=1
	@for tar in $(FILELIST_SUBSET); do \
		$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -t $${tar}" -r ProcessDirs-$$c-StdOut $(QUEUE_OPTIONS) ;\
        c=$$(( $${c:-$(COUNT)}+1 )); \
	done ;
#	@/bin/rm -rf $(DIRPATH)/SEGMENTATION/*
	

process_tars_array:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files ++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@echo "Submitting first array job..."
	jobid1=$$($(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON/$$FILE" -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON -R "6,10p" -r ProcessDirs-test-StdOut_1 -N aspen8 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}')
	echo "Submitting second array job..."
	jobid2=$$($(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON/$$FILE" -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON -R "11,14p" -r ProcessDirs-test-StdOut_2 -N aspen9 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}')
	echo "Submitted job1: $$jobid1, job2: $$jobid2"
	test -d CSV_combining-StdOut || mkdir -p CSV_combining-StdOut
	sbatch --dependency=afterok:$$jobid1,$$jobid2 --job-name=CSV_combining --output=CSV_combining-StdOut/CSV_combining_%j.out --error=CSV_combining-StdOut/CSV_combining_%j.err $(QUEUE_OPTIONS_batch) --wrap="$(DIRPATH)/SCRIPTS/Combine_size_csvs.py --base-dir $(DIRPATH)/CSV" 
	#$(DIRPATH)/SCRIPTS/SLURM_Batch -c '$(DIRPATH)/SCRIPTS/Combine_size_csvs.py --base-dir $(DIRPATH)/CSV' -r CSV_combining-StdOut --dependency=afterok:$$jobid1,$$jobid2 -N beech1 $(QUEUE_OPTIONS) 
	#@$(DIRPATH)/SCRIPTS/SLURM_Batch -c '$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON/$$FILE' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON -R "8001,9000p" -r ProcessDirs-StdOut_3 -N beech0 $(QUEUE_OPTIONS)
	#@$(DIRPATH)/SCRIPTS/SLURM_Batch -c '$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON/$$FILE' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON -R "9001,10000p" -r ProcessDirs-StdOut_4 -N beech1 $(QUEUE_OPTIONS)
	#@$(DIRPATH)/SCRIPTS/SLURM_Batch -c '$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON/$$FILE' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON -R "10001,11000p" -r ProcessDirs-StdOut_5 -N beech0 $(QUEUE_OPTIONS)
	#@$(DIRPATH)/SCRIPTS/SLURM_Batch -c '$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON/$$FILE' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/input_tars_MBON -R "11001,12000p" -r ProcessDirs-StdOut_7 -N beech1 $(QUEUE_OPTIONS)
#more than 5000 tars can't be done in bridges. the Stdout_6 was not submitted



#Note 06/02/2025, processing of inputTars is running up to 22500, continue from there 

process_tars_array_le:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files and combining the size CSVs afterwards++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@test -d CSV_combining-StdOut || mkdir -p CSV_combining-StdOut
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "1, 1500p" -r ProcessDirs-MBON-le-StdOut_1 $(QUEUE_OPTIONS)
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "1501, 3000p" -r ProcessDirs-MBON-le-StdOut_2 $(QUEUE_OPTIONS)
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "3001, 4500p" -r ProcessDirs-MBON-le-StdOut_3 $(QUEUE_OPTIONS)
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "4501, 6000p" -r ProcessDirs-MBON-le-StdOut_4 $(QUEUE_OPTIONS)
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "6001, 7700p" -r ProcessDirs-MBON-le-StdOut_5 $(QUEUE_OPTIONS)


process_trainingLib_sizes:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files and combining the size CSVs afterwards++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "which $(CONDA_EXE); $(PY) $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars_trainingLib/\$$FILE" -D $(DIRPATH)/inputTars_trainingLib -R "1, 2p" -r ProcessDirs-trainingLib-StdOut $(QUEUE_OPTIONS)


process_tars_array_86:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files and combining the size CSVs afterwards++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@test -d CSV_combining-StdOut || mkdir -p CSV_combining-StdOut
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "23035, 23500p" -r ProcessDirs-SPECTRA-86-StdOut_6 -q all.q -P 1 -N aspen12 
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "23501, 24500p" -r ProcessDirs-SPECTRA-86-StdOut_2 -q all.q -P 1 -N aspen13 
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "24501, 25500p" -r ProcessDirs-SPECTRA-86-StdOut_3 -q all.q -P 1 -N aspen12 
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "25501, 26500p" -r ProcessDirs-SPECTRA-86-StdOut_4 -q all.q -P 1 -N aspen13
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -c TRUE -t $(DIRPATH)/inputTars/\$$FILE" -D $(DIRPATH)/inputTars -R "26501, 27500p" -r ProcessDirs-SPECTRA-86-StdOut_5 -q all.q -P 1 -N aspen13


combine_size_csvs_le:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining the size CSVs afterwards++==>\e[0m"
	@test -d CSV_latest_combining-StdOut || mkdir -p CSV_latest_combining-StdOut
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE_le) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Combine_size_csvs.py --base-dir $(DIRPATH)/CSV_latest" -r CSV_latest_combining-StdOut $(QUEUE_OPTIONS) 


combine_size_csvs_le_test:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining the size CSVs afterwards++==>\e[0m"
	@test -d CSV_le_combining-StdOut || mkdir -p CSV_le_combining-StdOut
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(CONDA_ACTIVATE) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Combine_size_csvs.py --base-dir $(DIRPATH)/CSV_test_le" -r CSV_le_combining-StdOut $(QUEUE_OPTIONS) 

process_tars_array_combineCSVs_le-old: 
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files and combining the size CSVs afterwards++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@test -d CSV_combining-StdOut || mkdir -p CSV_combining-StdOut
	@echo "$(CONDA_ACTIVATE) $(CONDA_ENV)"; 
	@echo "Submitting first array job..."; \
        jobid1=$$($(CONDA_ACTIVATE) $(CONDA_ENV); $(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5301,5310p" -r ProcessDirs-StdOut_1 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting second array job..."; \
        jobid2=$$($(CONDA_ACTIVATE) $(CONDA_ENV); $(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5311,5320p" -r ProcessDirs-StdOut_2 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting third array job..."; \
        jobid3=$$($(CONDA_ACTIVATE) $(CONDA_ENV); $(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5321,5330p" -r ProcessDirs-StdOut_3 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting fourth array job..."; \
        jobid4=$$($(CONDA_ACTIVATE) $(CONDA_ENV); $(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5331,5340p" -r ProcessDirs-StdOut_4 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting fifth array job..."; \
        jobid5=$$($(CONDA_ACTIVATE) $(CONDA_ENV); $(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5341,5350p" -r ProcessDirs-StdOut_5 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitted job1: $$jobid1, job2: $$jobid2, job3: $$jobid3, job4: $$jobid4, job5: $$jobid5"; \
        sbatch --dependency=afterok:$$jobid1,$$jobid2,$$jobid3,$$jobid4,$$jobid5 --job-name=CSV_combining_MBON \
          --output=CSV_combining-StdOut/CSV_combining_%j.out \
          --error=CSV_combining-StdOut/CSV_combining_%j.err \
          $(QUEUE_OPTIONS_batch) \
          --wrap="$(CONDA_ACTIVATE) $(CONDA_ENV); $(DIRPATH)/SCRIPTS/Combine_size_csvs.py --base-dir $(DIRPATH)/CSV" --mem=16000


process_tars_array_combineCSVs_test_le: 
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files and combining the size CSVs afterwards++==>\e[0m"
	@test -d $(DIRPATH)/CSV_test_le || mkdir -p $(DIRPATH)/CSV_test_le
	@test -d CSV_test_le_combining-StdOut || mkdir -p CSV_test_le_combining-StdOut
	@echo "Submitting first array job..."; \
        jobid1=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5301,5310p" -r ProcessDirs-StdOut_1 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting second array job..."; \
        jobid2=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5311,5320p" -r ProcessDirs-StdOut_2 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting third array job..."; \
        jobid3=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5321,5330p" -r ProcessDirs-StdOut_3 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting fourth array job..."; \
        jobid4=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5331,5340p" -r ProcessDirs-StdOut_4 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting fifth array job..."; \
        jobid5=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5341,5350p" -r ProcessDirs-StdOut_5 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitted job1: $$jobid1, job2: $$jobid2, job3: $$jobid3, job4: $$jobid4, job5: $$jobid5"; \
        sbatch --dependency=afterok:$$jobid1,$$jobid2,$$jobid3,$$jobid4,$$jobid5 --job-name=CSV_test_le_combining_MBON \
          --output=CSV_test_le_combining-StdOut/CSV_test_le_combining_%j.out \
          --error=CSV_test_le_combining-StdOut/CSV_test_le_combining_%j.err \
          $(QUEUE_OPTIONS_batch) \
          --wrap="$(CONDA_ACTIVATE) $(CONDA_ENV); which python3; $(DIRPATH)/SCRIPTS/Combine_size_csvs.py --base-dir $(DIRPATH)/CSV_test_le" --mem=16000




process_tars_array_combineCSVs: 
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files and combining the size CSVs afterwards++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@test -d CSV_combining-StdOut || mkdir -p CSV_combining-StdOut
	@echo "Submitting first array job..."; \
	jobid1=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
	  -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
	  -D $(DIRPATH)/input_tars_MBON \
	  -R "5301,5310p" -r ProcessDirs-StdOut_1 -N aspen10 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
	echo "Submitting second array job..."; \
	jobid2=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
	  -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
	  -D $(DIRPATH)/input_tars_MBON \
	  -R "5311,5320p" -r ProcessDirs-StdOut_2 -N aspen11 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
	echo "Submitting third array job..."; \
        jobid3=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5321,5330p" -r ProcessDirs-StdOut_3 -N aspen10 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
	echo "Submitting fourth array job..."; \
        jobid4=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5331,5340p" -r ProcessDirs-StdOut_4 -N aspen11 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
        echo "Submitting fifth array job..."; \
        jobid5=$$($(DIRPATH)/SCRIPTS/SLURM_Batch \
          -c "$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -s /scratch -t $(DIRPATH)/input_tars_MBON/\$$FILE" \
          -D $(DIRPATH)/input_tars_MBON \
          -R "5341,5350p" -r ProcessDirs-StdOut_5 -N aspen10 $(QUEUE_OPTIONS) | awk '/Submitted batch job/ {print $$NF}'); \
	echo "Submitted job1: $$jobid1, job2: $$jobid2, job3: $$jobid3, job4: $$jobid4, job5: $$jobid5"; \
	sbatch --dependency=afterok:$$jobid1,$$jobid2,$$jobid3,$$jobid4,$$jobid5 --job-name=CSV_combining_MBON \
	  --output=CSV_combining-StdOut/CSV_combining_%j.out \
	  --error=CSV_combining-StdOut/CSV_combining_%j.err \
	  $(QUEUE_OPTIONS_batch) \
	  --wrap="$(CONDA_ACTIVATE_array); which python3; $(DIRPATH)/SCRIPTS/Combine_size_csvs.py --base-dir $(DIRPATH)/CSV" --mem=16000




process_tars:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Lots of Tar Files ++==>\e[0m"
	@test -d $(DIRPATH)/CSV || mkdir -p $(DIRPATH)/CSV
	@c=1
	@for tar in $(FILELIST); do \
	$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(DIRPATH)/SCRIPTS/Process_TarBall_PixelCount.py -t $${tar}" -r ProcessDirs-$$c-StdOut $(QUEUE_OPTIONS) ;\
	c=$$(( $${c:-$(COUNT)}+1 )); \
	done ;


process_file:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go doing SEGMENTATION using Threshold-MSER ++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(DIRPATH)/SCRIPTS/calculate_dark_pixels.py -i inputfile.jpg -c output.csv -t $(THRESHOLD)" -r ProcessFile-StdOut $(QUEUE_OPTIONS)
	./SCRIPTS/Process_TarBall_PixelCount.py -t /home/zoo/corralma/SPECTRA_proj/data/SPECTRA_tars/Camera_Top-03-09-22-01-20-17.898.tar.gz -v #/nfs3/Cowen_Lab/Mezcal/2018/Winter/BigCam/101A/segmentation/20180215210008.264.tar.gz -v

process_dirs:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Process Directories ++==>\e[0m"
	@for dir in $(FILELIST); do \
	$(DIRPATH)/SCRIPTS/SLURM_Batch -c "$(DIRPATH)/SCRIPTS/calculate_dark_pixels.py -i /home/cqls/sullichr/Labs/Cowen/marco/testData/$${dir}/seg -c $(DIRPATH)/CSV/$${dir}.csv -t $(THRESHOLD)" -r ProcessDirs-$${dir}-StdOut $(QUEUE_OPTIONS) ;\
	done ;

process_generate_csv_file_list:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining merged files and size data files++==>\e[0m"
	@test -d $(DIRPATH)/data || mkdir -p $(DIRPATH)/data
	@test -d $(DIRPATH)/data/combined_data || mkdir -p $(DIRPATH)/data/combined_data
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/1_combine_merge_size.R -m /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/merged_files/$$FILE -l YES -c CSV' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/merged_files -R "1p" -r GenerateCSV_sizeList-StdOut -q all.q -P 1 


NODE := aspen5

process_combine_merged_size:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining merged files and size data files++==>\e[0m"
	@test -d $(DIRPATH)/data || mkdir -p $(DIRPATH)/data
	@test -d $(DIRPATH)/data/combined_data || mkdir -p $(DIRPATH)/data/combined_data
	#@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'R --max-ppsize=500000 --no-save --slave --file=$(DIRPATH)/SCRIPTS/1_combine_merge_size.R --args -m /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files/$$FILE -l NO -c CSV' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files -R "1,20p" -r CombineMergedSize_MBON_ISIIS-3-StdOut -q all.q -P 1 -m 32g
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/1_combine_merge_size.R -m /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files/$$FILE -l NO -c CSV -t ancillary/transect_lats.csv' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files -R "1,5p" -r CombineMergedSize_MBON_ISIIS-3-_1-5_StdOut -q all.q -P 1 -m 32g -N aspen3
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/1_combine_merge_size.R -m /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files/$$FILE -l NO -c CSV -t ancillary/transect_lats.csv' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files -R "6,10p" -r CombineMergedSize_MBON_ISIIS-3-_6-10_StdOut -q all.q -P 1 -m 32g -N aspen5
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/1_combine_merge_size.R -m /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files/$$FILE -l NO -c CSV -t ancillary/transect_lats.csv' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files -R "11,15p" -r CombineMergedSize_MBON_ISIIS-3_11-15_StdOut -q all.q -P 1 -m 32g -N aspen11
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/1_combine_merge_size.R -m /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files/$$FILE -l NO -c CSV -t ancillary/transect_lats.csv' -D /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/merged_files -R "16,20p" -r CombineMergedSize_MBON_ISIIS-3_16-20_StdOut -q all.q -P 1 -m 32g -N aspen12


CRUISE := w22
TRANSECT := NH
TOW_ID := D
TAXON := sebastes

process_subset_by_taxon:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go subsetting the combined size data by taxon++==>\e[0m"
	@test -d $(DIRPATH)/data || mkdir -p $(DIRPATH)/data
	@test -d $(DIRPATH)/data/binned_data || mkdir -p $(DIRPATH)/data/binned_data
	@for combined_dir in $$(find $(DIRPATH)/data/combined_data -mindepth 1 -maxdepth 1 -type d -name "*$(TRANSECT)*"); do \
	echo "Processing directory: $$combined_dir"; \
	$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/1.1_subset_size_data_by_taxon.R -c $${combined_dir}/\$$FILE -p $(TAXON) -l YES -v SPECTRA -t ancillary/transect_lats.csv -b ancillary/shelf_break_200.csv" -D $${combined_dir} -R "1,200p" -r SubsetByTaxonSizeData-$(TAXON)-$$(basename $$combined_dir)-StdOut -q all.q -P 1 -N aspen11 ;\
        done ;

process_combine_subsetted_taxon_data:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go  combining all the files for $(TAXON) ++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/1.1.2_combined_subsetted_taxa.R -p $(TAXON) -t ancillary/transect_lats.csv' -r CombineTaxonFiles-$(TAXON)-StdOut -q all.q -P 1

process_check_depth_from_merged:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go checking Depth from merged files for $(CRUISE) $(TRANSECT) ++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/supp_exploreDepthStamps.R -c $(CRUISE) -t $(TRANSECT)" -r CheckDepthFromMerged-$(CRUISE)-$(TRANSECT)-StdOut -q all.q -P 1

process_MBON_combine_cams:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining the camera data by $(CRUISE) and $(TRANSECT)++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/2.0.1_combineCam_merged_data.R -c $(CRUISE) -t $(TRANSECT) -w $(TOW_ID)' -r CombineCams-$(CRUISE)-$(TRANSECT)-$(TOW_ID)-StdOut -q all.q -P 1 -N $(NODE) 

process_MBON_get_non_overlapping_Cam_data:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining the camera data by $(CRUISE) and $(TRANSECT)++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/2.0.2_get_singleCam_data.R -c $(CRUISE) -t $(TRANSECT) -w $(TOW_ID) -m 5' -r Non_overlap_data-$(CRUISE)-$(TRANSECT)-$(TOW_ID)-StdOut -q all.q -P 1 -m 128g -N $(NODE) 


process_SPECTRA_remove_NAN_2Cams:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining the camera data by $(CRUISE) and $(TRANSECT)++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/supp_remove_NA_from_combinedCams.R -c $(CRUISE) -t $(TRANSECT)' -r remove_NAN_fromCombined2Cams-$(CRUISE)-$(TRANSECT)-StdOut -q all.q -P 1 -m 128g 


OPTS_ABUNDANCE := -f

process_SPECTRA_binned_abundance_7hrCorr:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining the camera data by $(CRUISE) and $(TRANSECT)++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/2.1_ISIIS_binned_abundance_data_7hrCorr.R -c $(DIRPATH)/data/combined_data_2Cams/$(CRUISE)/$(TRANSECT)/$(TOW_ID) -i 3 -v MBON_ISIIS-3 -t ancillary/transect_lats.csv -b ancillary/shelf_break_200.csv -d ancillary/transect_shelf_land.csv $(OPTS_ABUNDANCE)' -r Binned_abundance-$(CRUISE)-$(TRANSECT)-$(TOW_ID)-StdOut_7hrCorr -q all.q -P 1 -m 32g -N $(NODE)

process_SPECTRA_get_cast_depth_grouping:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go getting the files correspodnign to each cast id in  $(CRUISE) and $(TRANSECT)++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/2_NEW_create_file_groups.R -c $(DIRPATH)/data/combined_data_2Cams/$(CRUISE)/$(TRANSECT)/$(TOW_ID)' -r FileList_cast-id-$(CRUISE)-$(TRANSECT)-$(TOW_ID)-StdOut -q all.q -P 1 -m 32g -N $(NODE)


process_get_meanDateTime_locations:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go getting locations for each date time stamp for  $(CRUISE) ++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/2.1.1_get_dateTime_locations.R -c $(CRUISE)' -r Locations_dateTime-$(CRUISE)-StdOut -q all.q -P 1 -m 32g

process_SPECTRA_binned_abundance:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go combining the camera data by $(CRUISE) and $(TRANSECT)++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/2.1_ISIIS_binned_abundance_data.R -c $(DIRPATH)/data/combined_data_2Cams/$(CRUISE)/$(TRANSECT) -i 3 -v SPECTRA -t ancillary/transect_lats.csv -b ancillary/shelf_break_200.csv -d ancillary/transect_shelf_land.csv -a FALSE -s TRUE -f FALSE' -r Binned_abundance-$(CRUISE)-$(TRANSECT)-StdOut -q all.q -P 1 -m 32g


THRESHOLD:= TRUE

process_calculate_size_spectra_NEW:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go estimating Size spectra following Dugenne et al. 2024 ++==>\e[0m"
	@test -d $(DIRPATH)/data/size_spectra || mkdir -p $(DIRPATH)/data/size_spectra
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/3_NEW_ISIIS_bin_and_NBSS.R \
                                                -c /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON_ISIIS-3/data/combined_data_2Cams/$(CRUISE)/$(TRANSECT)/$(TOW_ID) \
                                                -i 3 \
                                                -v SPECTRA \
						-s contour \
                                                -e time_Pacific,Latitude,Longitude,Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l \
                                                -g 5 \
						-x 5 \
						-p 0.3 \
   						-t ancillary/transect_lats.csv \
						-b ancillary/shelf_break_200.csv \
						-d ancillary/transect_shelf_land.csv \
						-r $(THRESHOLD)" \
						-r CalcSizeSpectra-plankton-$(THRESHOLD)-$(CRUISE)-$(TRANSECT)-$(TOW_ID)-StdOut -q all.q -P 1 -m 128g -N $(NODE) 


process_SPECTRA_krige_plankton_abundance:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go kriging $(CRUISE) and $(TRANSECT) binned abundance data++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/4_Kriging_new_abundance.R -c $(CRUISE) \
												 -t $(TRANSECT) \
												 -e Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l \
										    		 -d ancillary/transect_shelf_land.csv' -r Kriging_abundance-$(CRUISE)-$(TRANSECT)-StdOut -q all.q -P 1 -m 32g

process_LISST_size_spectra:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go calculating LISST size spectra++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/supp_LISST_size_spectra_processing.R -d 5' -r LISST_sizeSpectra-Calc-thresholded-StdOut -q all.q -P 1 -m 32g



process_merge_LISST-ISIIS_NBSS:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go estimating Size spectra following Dugenne et al. 2024 ++==>\e[0m"
	@test -d $(DIRPATH)/data/ISIIS-LISST_mergedNBSS_5m_allParticles || mkdir -p $(DIRPATH)/data/ISIIS-LISST_mergedNBSS_5m_allParticles
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/5_LISST_ISIIS_NBSS.R \
                                                -c $(CRUISE) \
                                                -t $(TRANSECT) \
                                                -s contour \
                                                -e time_Pacific,Latitude,Longitude,Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l \
                                                -l ancillary/transect_lats.csv \
                                                -b ancillary/shelf_break_200.csv \
                                                -d ancillary/transect_shelf_land.csv" -r LISST-ISIIS_SizeSpectra-allParticles-$(CRUISE)-$(TRANSECT)-StdOut -q all.q -P 1 -m 128g 






ENV_VAR := Altitude
process_add_enviro_vars:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go calculating LISST size spectra++==>\e[0m"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c 'Rscript $(DIRPATH)/SCRIPTS/supp_add_enviro_vars.R -v $(ENV_VAR) \
											      -s sizeSpectra \
											      -p contour \
											      -d 5
											      -t 10' -r addEnvVar_$(ENV_VAR)-StdOut -q all.q -P 1 -m 32g



process_bin_size_data:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go binning ISIIS size data++==>\e[0m"
	@test -d $(DIRPATH)/data || mkdir -p $(DIRPATH)/data
	@test -d $(DIRPATH)/data/binned_data || mkdir -p $(DIRPATH)/data/binned_data
	@echo "Processing directory: data/combined_data_2Cams_forBinning/$(CRUISE)/$(TRANSECT)"
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/2_ISIIS_bin_size_data.R -c $(DIRPATH)/data/combined_data_2Cams_forBinning/$(CRUISE)/$(TRANSECT) -i 3 -v SPECTRA -s contour -t ancillary/transect_lats.csv -b ancillary/shelf_break_200.csv" -r BinCombinedSizeData-$(CRUISE)-$(TRANSECT)-StdOut -N $(NODE) 

process_bin_size_data_old:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go binning ISIIS size data++==>\e[0m"
	@test -d $(DIRPATH)/data || mkdir -p $(DIRPATH)/data
	@test -d $(DIRPATH)/data/binned_data || mkdir -p $(DIRPATH)/data/binned_data
	@for combined_dir in $$(find $(DIRPATH)/data/combined_data -mindepth 1 -maxdepth 1 -type d); do \
        echo "➡️  Processing directory: $$combined_dir"; \
	$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/2_ISIIS_bin_size_data.R -c $${combined_dir}/\$$FILE -i 2 -v MEZCAL -s contour -t ancillary/transect_lats.csv -b ancillary/shelf_break_200.csv" -D $${combined_dir} -R "1,200p" -r BinCombinedSizeData-$$(basename $$combined_dir)-StdOut -N aspen11 $(QUEUE_OPTIONS) ;\
        done ;

process_calculate_size_spectra:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go estimating Size spectra following Dugenne et al. 2024 ++==>\e[0m"
	@test -d $(DIRPATH)/data/size_spectra || mkdir -p $(DIRPATH)/data/size_spectra
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/3_ISIIS_calculate_size_spectra.R \
						-b /nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/data/binned_data \
						-s contour \
						-i 2 \
						-g season,Transect_name,shelf_loc_200,depth_mid,lat_mid,lon_mid \
						-e Temp,Salinity,Fluoro,Oxygen,Irrandiance \
                                                -t FALSE" -r CalcSizeSpectra-StdOut -N aspen12 -m 16g $(QUEUE_OPTIONS)

process_make_diagnostic_figs:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here we go Processing a Tar File ++==>\e[0m"
	@test -d $(DIRPATH)/figures || mkdir -p $(DIRPATH)/figures
	@$(DIRPATH)/SCRIPTS/SLURM_Batch -c "Rscript $(DIRPATH)/SCRIPTS/4_ISIIS_PSS_diagnostic_plots.R" -r ProcessFigs-StdOut -N aspen12 $(QUEUE_OPTIONS)


status:
	@$(eval PART=`echo $(PARTITION) | sed 's/-q//g' | sed 's/ //g'`)
	@$(eval NODE=`echo $(NODELIST) | sed 's/-N//g' | sed 's/ //g'`)
	@$(eval SQ=`squeue | grep ${USER} | wc -l`)
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Current Status of Jobs and Process on GPU ++==>\e[0m"
	@echo ""
	@echo -e "\e[32;1;202m   SQUEUE for user ${USER}:\e[0m"
	@if [ $(SQ) != "0" ]; then \
	squeue | head -1 ; \
	squeue | grep ${USER} ; \
	echo "" ; \
	echo -e "\e[32;1;202m   NVIDIA-SMI for $(PART)@$(NODE):\e[0m" ; \
	srun -n 1 -N 1 -w $(NODE) -p $(PART) --propagate=NONE nvidia-smi ; \
	else \
	echo -e "\tYou have no jobs running....\n" ; \
	fi
	@echo ""

vars:
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Here the set of default variables ++==>\e[0m"
	@echo ""
	@echo -e "\tFull set of Default Variables and Current Values:"
	@echo -e "\t\tPARTITION     =  $(PARTITION)"
	@echo -e "\t\tNODELIST      =  $(NODELIST)"
	@echo -e "\t\tQUEUE_OPTIONS =  $(QUEUE_OPTIONS)"
	@echo -e "\t\tDEVICE        =  $(DEVICE)"
	@echo -e "\t\tPROCS         =  $(PROCS)"
	@echo -e "\t\tDIRPATH       =  $(DIRPATH)"
	@echo -e "\t\tSEG_OPTIONS   =  $(SEG_OPTIONS)"
	@echo ""
	@echo -e "\t* Any of these can be changed at the command line by typing VARIABLE_NAME=VALUE"
	@echo -e "\t\tExamples:"
	@echo -e "\t\t\tmake PARTITION=ceoas vars"
	@echo -e "\t\t\tmake PARTITION=ceoas-arm NODELIST=youmu vars"
	@echo -e "\t\t\tmake PARTITION=cqls_gpu NODELIST=cqls-gpu1 DEVICE=1 vars"
	@echo -e "\t\t\tmake PARTITION=cqls_ppc64le DEVICE=1 vars"
	@echo ""


clean: 
	@echo ""
	@echo -e "\e[36;1;202m   <==++ Removing StdOut folders  ++==>\e[0m"
	@while [ -z "$$CONTINUE" ]; do \
	    read -r -p "      This will get rid of all StdOut Folders that has been generated... [y/n] " CONTINUE; \
	done ; \
	if [ ! $$CONTINUE == "y" ]; then \
	if [ ! $$CONTINUE == "Y" ]; then \
	    echo ""; \
	    echo -e "\e[31;5;202m   OK if you want to keep them! Exiting Now!\e[0m" ; echo ""; exit 1 ; \
	fi \
	fi
	@echo ""
	@echo -e "      \e[32;1;202mRemoving all the StdOut folders.....\e[0m"
	@/bin/rm -r $(DIRPATH)/*StdOut*
	@echo ""

realclean:
	@echo ""
	@echo "Removing Everything Created:"
	@while [ -z "$$CONTINUE" ]; do \
	    read -r -p "      This will get rid of everything that has been generated... [y/n] " CONTINUE; \
	done ; \
	if [ ! $$CONTINUE == "y" ]; then \
	if [ ! $$CONTINUE == "Y" ]; then \
	    echo ""; \
	    echo -e "\e[31;5;202m   Yeah probably a good idea not to do that! Exiting Now!\e[0m" ; echo ""; exit 1 ; \
	fi \
	fi
	@echo ""
	@echo -e "      \e[32;1;202mRemoving basically everything that was generated....\e[0m"
	@/bin/rm -r $(DIRPATH)/CSV/*  $(DIRPATH)/*StdOut* $(DIRPATH)/SEGMENTATION/*
	@echo ""

superclean:
	@echo ""
	@echo "Removing Everything:"
	@while [ -z "$$CONTINUE" ]; do \
	    read -r -p "      This will get rid of everything that has been generated... [y/n] " CONTINUE; \
	done ; \
	if [ ! $$CONTINUE == "y" ]; then \
	if [ ! $$CONTINUE == "Y" ]; then \
	    echo ""; \
	    echo -e "\e[31;5;202m   Yeah probably a good idea not to do that! Exiting Now!\e[0m" ; echo ""; exit 1 ; \
	fi \
	fi
	@echo ""
	@echo -e "      \e[32;1;202mRemoving basically everything that was every created taking it back to just a Makefile.....\e[0m"
	@/bin/rm -r $(DIRPATH)/CSV $(DIRPATH)/DOWNLOADS $(DIRPATH)/*StdOut* $(DIRPATH)/SEGMENTATION
	@echo ""
