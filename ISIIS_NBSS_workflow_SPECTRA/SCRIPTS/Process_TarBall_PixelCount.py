#!/usr/bin/env python3
#####################################################################
#
#  OSU <putsomethinghere> Pipeline
#                                                                   
#  Copyright 2024
#
#  Put YOUR NAME HERE
#
#  College of Earth, Ocean, and Atmospheric Sciences
#  Oregon State University
#  Corvallis, OR 97331
#
#  email: youronid@oregonstate.edu
#
# This program is not free software; you can not redistribute it
# and/or modify it at all.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#                                                                   
#                                                                   
#####################################################################
from datetime import datetime
import time
import glob
import os
import sys
import getopt
import re
import subprocess
from pathlib import Path
import pandas as pd
import shutil
import sys
import time
import random
#####################################################################
#                                                                   #
#                         Start of Program                          #
#                                                                   #
#####################################################################

#####################################################################
# Setup verbose mode as no to we has a silent run by default        #
#####################################################################
verbose = "No"
total_processed = 0

#####################################################################
# Set the paths                                                     #
#####################################################################
path = os.getcwd()
script_path = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(script_path, os.pardir))
processing_dir = ""
scratch_dir = path + "/TMP"
tarball_file = ""
threshold = 240

#####################################################################
# datetime object containing current date and time                  #
#####################################################################
current_time = datetime.now()
current_unix = (time.mktime(current_time.timetuple()))

opts, args = getopt.getopt(sys.argv[1:],"hvd:t:o:s:c:",["processing_dir=","tarball_file=", "scratch_dir=", "check_CSV_exists="])
check_CSV_exists = False  # Default

for opt, arg in opts:
    if opt == '-h':
        print ("\n\tYou did not provide all the correct command line arguments\n")
        print ("\tUsage:\n")
        print ("\t" + str(sys.argv[0]) + " -d <processing_dir> -t <tarball_file> -h -v\n")
        print ("\t-d The data processing directory multiple tarball files (current: " + processing_dir + ")")
        print ("\t-t Single tarball file (current: " + tarball_file + ")")
        print ("\t-s The scratch directory for processing (current: " + scratch_dir + ")")
        print ("\t-v Verbose mode / Print dialog")
        print ("\t-h Help message")
        print ("\n")
        sys.exit()
    elif opt == '-v':
        verbose = "Yes"
    elif opt in ("-d", "--processing_dir"):
        processing_dir = arg
    elif opt in ("-t", "--tarball_file"):
        tarball_file = arg
    elif opt in ("-s", "--scratch_dir"):
        scratch_dir = arg
    elif opt in ("-c", "--check_CSV_exists"):
        print(f"DEBUG: received arg for -c was: {arg}", file=sys.stderr)
        check_CSV_exists = str(arg).strip().upper() in ["TRUE", "YES", "1"]
#####################################################################
# Start the run                                                     #
#####################################################################
print("\n\tStart of Run!", file=sys.stdout, flush=True)
print(f"sys.argv: {sys.argv}", file=sys.stderr)

#####################################################################
# Make a quick list of the variables                                #
#####################################################################
if verbose == "Yes":
    print ('\n\t  Current Paths:')
    print ("\t\t Main Directory:", path)
    if tarball_file:
        print ('\t\t   TarBall Path: ', tarball_file)
    elif processing_dir:
        print ('\t\tInput Directory: ', processing_dir)
    if scratch_dir:
        print ('\t\tScratch Directory: ', scratch_dir)
    print("")

print ("\tProcessing Data: ", file=sys.stdout, flush=True)


processing_dir_tmp = scratch_dir

first_iteration = True

if tarball_file:
    if verbose == "Yes":
        print ('\t\t   Path To File: ', os.path.join(path, tarball_file))
    if tarball_file.endswith('.tar') or tarball_file.endswith('.tar.gz'):
        base_name = os.path.basename(tarball_file)
        root_name, _ = os.path.splitext(base_name)
        root_name2, _ = os.path.splitext(root_name)
        if check_CSV_exists:
            target_csv = os.path.join(path, "CSV" ,root_name2, root_name2 + ".csv")
            print(f"checking for the existence of {target_csv}", file=sys.stderr)
            if os.path.exists(target_csv):
                try:
                    with open(target_csv, 'r') as f:
                        line_count = sum(1 for _ in f)
                        print(f"line count is {line_count}", file=sys.stderr)
                    if line_count > 1:
                        print(f"\n\t🚫 Process terminated: {target_csv} already exists and contains data", file=sys.stdout, flush=True)
                        sys.exit()
                except Exception as e:
                    print(f"\n\t⚠️ Error reading {target_csv}: {e}", file=sys.stderr, flush=True)
        if tarball_file.endswith('.tar.gz'):
            if not os.path.exists( str(processing_dir_tmp) + str('/') + str(root_name2)): # this line and next was added 8/12/2024, to match the tar structure, check this maybe with bing to ensure i am getting the right tarballs also, on 9/12/2024 Marco removed str(path) + str('/') + from the if not statement, i think it was producing incorrect directories, same with the makedirs below
                os.makedirs(str(processing_dir_tmp) + str('/') + str(root_name2))# str(path) + str('/') + 
            if not os.path.exists(str(path) + str("/CSV/") + str(root_name2)):
                os.makedirs(str(path) + str("/CSV/") + str(root_name2))
            # Do something with "tarball_file" variable
            
            # Create a dedicated extraction directory under scratch, 06/09/2025, lines 134-142 try indefinitely until successful
            extraction_dir = os.path.join(scratch_dir, f"{root_name2}")
            target_dir = os.path.join(processing_dir_tmp, root_name2)
            time.sleep(random.uniform(0, 2))
            # Retry loop until the directory is created successfully
            while True:
                try:
                    os.makedirs(extraction_dir, exist_ok=True)
                    break
                except PermissionError as e:
                    print(f"⚠️ Retrying mkdir for {extraction_dir}: {e}", file=sys.stderr)
                    time.sleep(1)

            # Extract tarball into that directory safely
            tar_cmd = f"tar --no-same-permissions --no-same-owner -xzf {tarball_file} -C {extraction_dir}"
            os.system(tar_cmd)

            # Move extracted content into SEGMENTATION directory
            
            #loop below only needed for ISIIS 2's nested dirs inside tarballs

            # Cleanup extraction directory
        else:
            if not os.path.exists(str(path) + str("/CSV/") + str(root_name)):
                os.makedirs(str(path) + str("/CSV/") + str(root_name))
            # Do something with "tarball_file" variable
            if not os.path.exists(str(scratch_dir) + str('/') + str(root_name2)):
                os.makedirs(str(scratch_dir) + str('/') + str(root_name2))
            os.system("tar --no-overwrite-dir --no-same-permissions --no-same-owner -xf " + str(tarball_file) + " -C " + str(scratch_dir))

        # Process the images through the pixel counter
        if verbose == "Yes":
            print ("\t\t Processing Tar:")
        for directory in [d for d in os.listdir(processing_dir_tmp) if root_name2 in d]:
            print (tarball_file)
            seg_dir = processing_dir_tmp + str("/") + directory + "/corrected_crop" # only applicable to tarballs from SPECTRA
            if verbose == "Yes":
                print ("\t\t\tSeg Folder: " + str(seg_dir))
            os.system(path + str("/SCRIPTS/calculate_dark_pixels.py -i ") + str(seg_dir) + str(" -t ") + str(threshold) + str(" -c ./CSV/") + str(directory) + str("/") + str(directory) + str(".csv"))# line modified on 12/11/20204, it was str(root_name2) + str("/") + str(directory) + str(".csv")
            print(path + str("/SCRIPTS/calculate_dark_pixels.py -i ") + str(seg_dir) + str(" -t ") + str(threshold) + str(" -    c ./CSV/") + str(directory) + str("/") + str(directory) + str(".csv"))
            total_processed = total_processed + 1 

            if verbose == "Yes":# modified on 2/4/2025 by M. Corrales. The rm -rf was included in the loop to make a targeted deletion of the tar directory in SEGMENTATION, else the remove commands might conflict with each other. Take these and the next two lines out of the loop to go back to previous version
                print ("\t\t    Cleaning Up:")
            print('deleting as rm -rf' +str(processing_dir_tmp) + str("/") + directory)
            os.system(str("rm -rf ") + str(processing_dir_tmp) + str("/") + directory)# this is the part that is failing, files under SEGMENTATION are not being deleted. removed the /bin/ for debugging in bridges infrastructure
    else:
        if verbose == "Yes":
            print ('\t\tExiting - Not and TarBall File: ', tarball_file)



elif processing_dir:
    if os.path.exists(processing_dir):
        for filename in os.listdir(processing_dir):
            if filename.endswith('.tar') or filename.endswith('.tar.gz'):
                tarball_file = filename
                total_processed = total_processed + 1
                if verbose == "Yes":
                    print ('\t\t  File: ', tarball_file)
                # Do something with "tarball_file" variable
            else:
                if verbose == "Yes":
                    print ('\t\tExiting - Not a tar or tar.gz file: ', tarball_file)
    else:
        print ('\t\tExiting - processing directory does not exists: ', processing_dir)

print ("\n\tTotal Processed Files: " + str(total_processed), file=sys.stdout, flush=True)

#
#
#
#
#
#
#
## === Call block inserted before sys.exit() ===
#    # Detect current tar ID from context (assume it is passed or inferred)
#    
#    
#        
#    
#    
#    
#
#        # Clean up only relevant subdirs

print ("\tEnd of Run!\n", file=sys.stdout, flush=True)

sys.exit()

#####################################################################
#                                                                   #
#                           End of Program                          #
#                                                                   #
#####################################################################
