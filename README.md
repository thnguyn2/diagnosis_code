#Prostate cancer diagnosis source code
Quantitative Light Imaging Lab, University of Illinois at Urbana-champaign
#Introduction
This is the source code for the prostate cancer diagnosis from QPI images. It it written in MATLAB by Tan H. Nguyen (huutan86@gmail.com). If you use our code, please consider citing it since it takes a lot of time and resources to develop.

#Installation
Create a folder where you want to install the code, cd in to that folder. In git, type

	git clone https://github.com/thnguyn2/diagnosis_code.git

Now, check out the master branch

	git checkout master

If you want to get updated with the latest version of the code

	git pull origin master

After cloning the source code, open iccv_paper/demo.m to begin.

#Source code organization
The project tree is as follows:

 * Folders:

	iccv_paper/: a folder containing working version of the source code
	non-relevant_old_code/: a folder with old, depricated functions
	support/: a folder with supporting function

 * Files:
	
The following files in the repo. folders are just for summarzing the results and new ideas. They are not very useful though. I will clean them up in a near future
	
	EM for diagnosis.jnt: a journal file containing Tan's derivation on how to train a hierachical tree to learn the structure of the multi-class diagnosis problem.
	generate_gt_from_code_diagnosis.m: [TBA]
	truncate_segmented_map_with_from_code_label.m: [TBA]
	Fuzzy map for feature review.pptx: a summary of the class likelihood produced by the classifier
	Note on the TMA2B data.txt: what i learnt from a discussion with Dr. Sridharan at UC Davis
	README.md: this file	

The following files are the actual code	

	iccv_paper/demo.m: the main source code to process all the data.

This file will call the following functions

	support/findFileNameFromRois.m: find the name of all tiff files available. The data need to be organized in a specific structure as described below
	
The following files are just for summarzing the results and new ideas. They are not very useful though. I will clean them up in a near future
	
	EM for diagnosis.jnt: a journal file containing Tan's derivation on how to train a hierachical tree to learn the structure of the multi-class diagnosis problem.
	generate_gt_from_code_diagnosis.m: [TBA]
	truncate_segmented_map_with_from_code_label.m: [TBA]
	Fuzzy map for feature review.pptx: a summary of the class likelihood produced by the classifier
	Note on the TMA2B data.txt: what i learnt from a discussion with Dr. Sridharan at UC Davis
	README.md: this file	


#Data organization
The QPI data is required to run the code. Please contact our lab director, Dr. Gabriel Popescu at gpopescu@illinois.edu if you want to obtain the data and run the code with it,

#How to use the code

#Acknowlegdement
