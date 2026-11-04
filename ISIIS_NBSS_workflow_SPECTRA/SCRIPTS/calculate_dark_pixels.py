#!/usr/bin/env python3
import cv2
import numpy as np
import os
import sys
import csv
import math as m 

def help_msg():
    print ("\n\tYou did not provide all the correct command line arguments")
    print ("\tUsage:")
    print ("\t" + str(sys.argv[0]) + " -i <processing_dir|input_crop> -t <threshold_value> -c <csv_output_file> -h -v")
    print("\t-i: this is input for a crop directory or individual img. EX -i ./crops or -i ./img.jpg")
    print("\t-t: this is threshold for what is considered light and dark(anything over this threshold will be consdiered dark). EX -t 240 (uses RGB values so in reality it would be 240,240,240)")
    print("\t-c: this is the input for what you want to name the output csv. EX -c test.csv")
    print ("\t-h Help message")
    print ("\n")
    sys.exit()


# Directory containing images
# OR single img
image_dir = ''  
lower = 240
upper = 255
csv_name = 'default'

for i in range(1,len(sys.argv)):
    if sys.argv[i] == '-i':
        image_dir = sys.argv[i+1]
    elif sys.argv[i] == '-t':
        lower = int(sys.argv[i + 1])
    elif sys.argv[i] == '-c':
        csv_name = sys.argv[i+1]
    elif sys.argv[i] == '-h':
        help_msg()

with open(csv_name,'w',newline='') as csvfile:
    fieldnames = ['crop','total_pixels','light','dark','percent_dark', 'contour', 'major_ellipse_axis', 'minor_ellipse_axis', 'ellipse_area']
    writer = csv.DictWriter(csvfile,fieldnames=fieldnames)

    if os.stat(csv_name).st_size == 0:
        writer.writeheader()


    if os.path.isdir(image_dir):
        # Loop through all images in the directory
        for filename in os.listdir(image_dir):
            if filename.endswith(('.png', '.jpg', '.jpeg')):  # Add more extensions if needed
                img_path = os.path.join(image_dir, filename)
                try: 
                    # Read the image
                    img = cv2.imread(img_path)

                    # Get the dimensions of the image
                    height, width = img.shape[:2]
                    total_pixels = height * width

                    # Extract the ROI from the original image (using the entire image)
                    imgCrop = img  # In this case, we use the entire image

                    # Define the color range you want to count (e.g., for grayscale)
                    lowerbound = np.array([lower, lower, lower])
                    upperbound = np.array([upper, upper, upper])

                    # Create a mask for the pixel range
                    mask = cv2.inRange(imgCrop, lowerbound, upperbound)

                    # Count non-zero pixels in the mask
                    countInRange = cv2.countNonZero(mask)

                    percent_in_range = (countInRange / total_pixels) * 100
                    dark = total_pixels - countInRange

                    # added 2/12/2024 by M. Corrales-Ugalde: process to get ellipse and contour 
                    # find contour area, gets the maximum particle contour from the image. image needs to be converted to grayscale and be inverted
                    countours,hierarchy=cv2.findContours(cv2.bitwise_not(cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)),cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
                    max_contour_area = cv2.contourArea(max(countours, key = cv2.contourArea))
                    try:
                        ellipse = cv2.fitEllipse(max(countours, key = cv2.contourArea)) # this function returns a tuple, second element are the major and minor axes, next line extracts these values directly from the fucntion
                        (x, y), (MA, ma), angle = cv2.fitEllipse(max(countours, key = cv2.contourArea)) # in addition, this allows to get the axis of the ellipse!
                        major_axis = max(MA,ma)
                        minor_axis = min(MA, ma)
                        area_ellipse=m.pi * (MA/2) * (ma/2)
                    except:
                        major_axis = float('NaN')
                        minor_axis = float('NaN')
                        area_ellipse = float('NaN')

                    writer.writerow({'crop': filename, 'total_pixels': total_pixels,'light': countInRange ,'dark': dark, 'percent_dark': percent_in_range, 'contour':max_contour_area, 'major_ellipse_axis':major_axis, 'minor_ellipse_axis':minor_axis, 'ellipse_area':area_ellipse})
                except:
                    pass
                # Output or store the count as needed
    else:
        if image_dir.endswith(('.png','.jpg','.jpeg','.tiff')):
            img_path = os.path.join(image_dir)
            try:
                # Read the image
                img = cv2.imread(img_path)

                # Get the dimensions of the image
                height, width = img.shape[:2]
                total_pixels = height * width

                # Extract the ROI from the original image (using the entire image)
                imgCrop = img  # In this case, we use the entire image

                # Define the color range you want to count (e.g., for grayscale)
                lowerbound = np.array([lower, lower, lower])
                upperbound = np.array([upper, upper, upper])

                # Create a mask for the pixel range
                mask = cv2.inRange(imgCrop, lowerbound, upperbound)

                # Count non-zero pixels in the mask
                countInRange = cv2.countNonZero(mask)

                percent_in_range = (countInRange / total_pixels) * 100
                dark = total_pixels - countInRange

                # added 2/12/2024 by M. Corrales-Ugalde: process to get ellipse and contour 
                # find contour area, gets the maximum particle contour from the image. image needs to be converted to grayscale and be inverted
                countours,hierarchy=cv2.findContours(cv2.bitwise_not(cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)),cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
                max_contour_area = cv2.contourArea(max(countours, key = cv2.contourArea))
                try:
                    ellipse = cv2.fitEllipse(max(countours, key = cv2.contourArea)) # this function returns a tuple, second element are the major and minor axes, next line extracts these values directly from the fucntion
                    (x, y), (MA, ma), angle = cv2.fitEllipse(max(countours, key = cv2.contourArea)) # in addition, this allows to get the axis of the ellipse!
                    major_axis = max(MA,ma)
                    minor_axis = min(MA, ma)
                    area_ellipse=m.pi * (MA/2) * (ma/2)
                except:
                    major_axis=float('NaN')
                    minor_axis=float('NaN')
                    area_ellipse=float('NaN')


                writer.writerow({'crop': filename, 'total_pixels': total_pixels,'light': light ,'dark': countInRange, 'percent_dark': percent_in_range, 'contour':max_contour_area, 'major_ellipse_axis':major_axis, 'minor_ellipse_axis':minor_axis, 'ellipse_area':area_ellipse})
            except:
                pass
            # Output or store the count as needed


