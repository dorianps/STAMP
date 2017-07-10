# STAMP
Stacked Multimodal Predictions  
*****  
  
This is a collection of scripts used to achieve the prediction of post-stroke aphasia scores from the integration of multiple modalities. Each modality is used to create one or more predictions of aphasia based on some measures (i.e., graph theory measures, connectivity matrix, lesion load, etc.). These preliminary predictions are then combined in a final multimodal prediction. The hierarchical nature of the predictions gives it the term "stacked" in the name. All predictions use random forests.

#### IMPORTANT
This is by no means a finished product. I am making the scripts public as an example on how STAMP can work. The code was built and expanded as part of a research study, and is dirty and far from ready to use. Even for me is difficult to read at time. It takes the results of the "recursive feature elimination" pipeline (caret package) to use only the necessary variables from each source of information. I am not placing the RFE results on github, they are too many files anyway. It also relies on a number of functions that may or may not be present in the github repository (cannot go back to find what calls what from the beginning of the study). As I said, this is just a collection of scripts that might be useful to researchers that want to apply the same concepts. If you see swear words in the code, they are just part of the standard researcher's frustration.  
  
  
#### File purpose:  
* `densitySearch_allThresholds_HBMrevision.R`: contains a full run of STAMP framework, produces also plots. Most of the code is same as in `densitysearch.R`
* `IITcomputeConnectomes.R`: compute virtual lesioned connectomes for each patient
* `STAMPboxplots4paper_gray.R`: Creates SMAP plots in black and white, for a journal that (in 2017) still doesn't accept online figures with colors.
* `strokeRest.R`: This is a fully functional function to run resting state analysis on stroke data. The inside documentation is not great yet, but the function can run everything from preprocessing to the production of multiple connectomes from different atlasses. It eliminates the lesioned area when computing nuisance variables during preprocessing.
* `pennStrokeRest.R/sh`: Files used to compute resting state connectomes from Penn patients. One is a simple R file, the other is the shell script used to call the R file in a cluster environment. These are used to run things in a computing cluster.
* `rfeCV.R`: Script used to compute full-split cross validation.
* `restAnalysis.R`: Some code I used to understand the potential of resting state or dti data.
* `selectRFE.singlejob.R/sh`: Main script that runs a single recursive feature elimination in the cluster. The concept was to prepare a small Rdata file with all what was needed for RFE, than call this script. It would load the Rdata file, run RFE and save the results on the same file. This way, all I need to keep track is Rdata filenames, instead of the data themselves.
* `selectVarsRF*`: a bunch of files probably similar to selectRFE, don't remember exactly at what stage did I use these scripts.

