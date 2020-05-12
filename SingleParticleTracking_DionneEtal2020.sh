# Instructions for the installation of `htsimaging`: https://github.com/rraadd88/htsimaging#installation
# Instructions for the usage of `htsimaging`: run `python endocytosis.py run-trials -h`
# Inputs:
# $project_dir : contains the subdirectories containing time-lapse images for individual trials.
# $bright_fn_marker: In each subdirectory containing images for individual trials the filename of bright-field image contains a unique label. 

# Activate the Anaconda environment 
conda activate htsimaging
# Execute the command for the analysis of images 
python endocytosis.py run-trials $project_dir $bright_fn_marker
