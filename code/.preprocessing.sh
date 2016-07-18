#find ../data/output/ -name "*.rds" -delete
#find ../data/output/ -name "*.csv" -delete

python data_processing/preprocessing.py
Rscript data_processing/preprocessing.R
