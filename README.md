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
The QPI data is required to run the code. Please contact our lab director, Dr. Gabriel Popescu at gpopescu@illinois.edu if you want to obtain the data and run the code with it. On the ozymandias server, the data can be found at raid6\Tan\Prostate_cancer_diagnosis_data\TMA_code_and_diagnosis\. The structure of this folder is as follows
  
 * d[x]+[y]\: where x, y in {2,3,4,5}. A folder contains all the data for the diagnosis group with Gleason score [x]+[y]. Inside this folder, you can find the raw QPI image, e.g., A19.tif, 3x downsampled QPI image, e.g., A19_small.tif, manual segmentation results, e.g., A19_seg_gt.tif, A19_seg_gt_3072.tif, A19_seg_gt_3072_hf.tif, final segmentation result, A19_seg_multi_res_3072_strict_hf.tif. Note that "hf" stands for hole-filling where the segmentation results have been post-process to remove holes inside the gland and stroma likelihood, e.g., A19_seg_ws[zz]_fuzzy.tif. Herezz is the size of the window in which the histogram of texton is computed. 

 * dBPH\, dHGPIN\, dNormal\: similar to d[x]+[y]. They store data for the BPH, HGPIN and benign group.

 * TMA2D_all_cores\: all SLIM images of all cores in the TMA.

 * texdir\: contains [corename]_lm_fr.mat: Leung-Malik filter response, [corename]_texton_index_[zz]_map.tif: texton index map where zz is the number of texton_map.tif: texton index map where zz is the number of texton. If [zz]='', the number of texton is 50. [corename]_texton_hist_[ww].mat: histogram of texton index data where ww is the window size for computing. kmeans_res_[zz]_clusters.mat: coordinates of textons where zz is the number of textons.
 
 * Verified_diagnosis_results: manual segmentaion with biopsy diagnosis, verified by Shamira.

#How to use the code

#Acknowlegdement
