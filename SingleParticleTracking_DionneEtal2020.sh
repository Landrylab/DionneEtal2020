# Instructions for the installation of `htsimaging`: https://github.com/rraadd88/htsimaging#installation
# Instructions for the usage of `htsimaging`: `python endocytosis.py run-trials -h`
# Inputs:
# $project_dir : contains the subdirectories containing time-lapse images for individual trials. In each subdirectory containing images for individual trials the filename of bright-field image contains an unique label.
# $bright_fn_marker: unique label in the filename of bright-field images.

# Activate the Anaconda environment 
conda activate htsimaging
# Execute the command for the analysis of time-lapse images 
python endocytosis.py run-trials $project_dir $bright_fn_marker
