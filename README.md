# adjusted-mf
This repository contains the code of the experiments described in "Misalignment problem in matrix decomposition with missing values" by Sofia Fernandes, Mário Antunes, Diogo Gomes and Rui L. Aguiar

## Installation
The code is written in Python 3.6. To install the required libraries run:
```
pip install -r requeriments.txt
```

## Usage
### Datasets
The datasets used are publicly available, and should be downloaded and renamed as follows:
- the gas sensor dataset is available [here](https://archive.ics.uci.edu/ml/machine-learning-databases/00309/dataset_two_gas_sources_.zip) -  the file should be unzipped into folder experiments_on_gas_sensor_dataset/data/gas_sensor
- daily average temperature datasets should be downloaded using the script get_noaa_temp_data.R at experiments_on_temperature_datasets/data, along with a token that should be requested at [NOAA's  National Climatic Data Center](https://www.ncdc.noaa.gov/cdo-web/token)
- [Seattle traffic dataset](https://github.com/zhiyongc/Seattle-Loop-Data) can be obtained from [here](https://raw.github.com/xinychen/transdim/master/datasets/Seattle-data-set/mat.csv) - the file should renamed to 'seattle.csv' and placed at experiments_on_traffic_datasets/data
- [Guangzhou traffic dataset](https://doi.org/10.5281/zenodo.1205229) can be obtained from [here](https://raw.github.com/xinychen/transdim/master/datasets/Guangzhou-data-set/tensor.mat) - the file should be renamed to 'guangzhou.mat' and placed at experiments_on_traffic_datasets/data

The daily average temperature and traffic datasets require further preprocessing that is carried out by running the preprocessing scripts (found in the corresponding 'data' folders).

### Run
The experiments scripts should be run in the corresponding experiment folder, for example to run the experiments on the gas sensor data, run:
```
cd experiments_on_gas_sensor
mkdir results
python gas_NMF_adjustment_exps.py
```
the script outputs the PRESS for the multiple settings and stores all the results in the results folder.
