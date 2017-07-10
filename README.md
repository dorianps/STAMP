# STAMP
Stacked Multimodal Predictions  
*****  
  
This is a collection of scripts used to achieve the prediction of post-stroke aphasia scores from the integration of multiple modalities. Each modality is used to create one or more predictions of aphasia based on some measures (i.e., graph theory measures, connectivity matrix, lesion load, etc.). These preliminary predictions are then combined in a final multimodal prediction. The hierarchical nature of the predictions gives it the term "stacked" in the name. All predictions use random forests.

#### IMPORTANT
This is by no means a finished product. I am making the scripts public as an example on how STAMP can work. The code was built and expanded as part of a research study, and is dirty and far from ready to use. Even for me is difficult to read at time. It takes the results of the "recursive feature elimination" pipeline (caret package) to use only the necessary variables from each source of information. I am not placing the RFE results on github, they are too many files anyway. It also relies on a number of functions that may or may not be present in the github repository (cannot go back to find what calls what from the beginning of the study). As I said, this is just a collection of scripts that might be useful to researchers that want to apply the same concepts.

