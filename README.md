## Predicting IFN Response as a Function of ADAR1 Editing, L1 RT, and dsRNA-Producing SINE Expression
There are two Python scripts designed to train and optimize the predicted IFN response by comparing it to the defined 20 RIG-I IFN gene expressions from the Sun et al. (2024) Immunity paper.

- helper_funcs.py contains several data-wrangling functions to generate the immunity_model.tsv table. An example table is provided in the data folder.

- fit_model.py iterates over different L1 and ADAR1 weights to minimize the differences between predicted and observed IFN responses from RNA-seq data.

To run the code, you need to clone this repository and then simply execute the following:

```bash
# Activate the conda environment
conda activate Immunity_env![GraphicalAbstract](https://github.com/user-attachments/assets/bf122e72-f8f7-480f-afb7-a90f9500f9bd)


# Run the Python script
python fit_model.py
