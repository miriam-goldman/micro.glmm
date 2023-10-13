# micro-glmm
## Purpose to Make an deploy Micro-glmm R package and files to run it
<img width="985" alt="image" src="https://github.com/miriam-goldman/micro-glmm/assets/46382670/7f1fde4b-7439-4c0c-a99d-bcba26fbe5f3">

# Instructions on how to run
### Install Conda
```
conda env create -f microglmm.yml
conda activate microglmm
```
### Move into R folder

```
cd src
```
### run prep snps
```
Rscript prep_snps_command_line.R  --species_id exp_species -v TRUE  -o "../example_output"   -i "../example_data/exp_species.snps_info.tsv"   -f "../example_data/exp_species.snps_freqs.tsv"  -d "../example_data/exp_species.snps_depth.tsv"   -a 3 -m 5  -n 5   
 ```
### run prep genes

```
Rscript prep_genes_command_line.R --species_id exp_species -v TRUE -o "../example_output" -c "../example_data/exp_species.genes_copynum.tsv" -d "../example_data/exp_species.genes_depth.tsv" -m 5 -n 5 --GRM ../example_output/exp_species.GRM.tsv --write_csv TRUE
 ```
### run population structure test

```
Rscript pop_structure_test_command_line.R --species_id exp_species -v TRUE -o ../example_output --GRM ../example_output/exp_species.GRM.tsv --metadata ../example_data/IBD_metadata.tsv --tau 1 --n_tau 10

```

### Run Maker test

```
Rscript marker_test_command_line.R --species_id exp_species -v TRUE -o ../example_output --Rdata ../example_output/exp_species.model_obj.Rdata --copy_number ../example_output/exp_species.copynumber_data.csv

```




### 2023-09-27: Made Repository
### 2023-10-02: updated function prep-snps; example of how to run with simplest options
```
Rscript prep_snps_command_line.R  --species_id exp_species -v TRUE  -o "../example_output"   -i "../example_data/exp_species.snps_info.tsv"   -f "../example_data/exp_species.snps_freqs.tsv"  -d "../example_data/exp_species.snps_depth.tsv"   -a 3 -m 5  -n 5   
 ```
 still to do: make sure qp and pangeomne run (waiting on new pangeome)
 
 ### 2023-10-03: updated function prep-genes; example of how to run with simplest options
```
Rscript prep_genes_command_line.R --species_id exp_species -v TRUE -o "../example_output" -c "../example_data/exp_species.genes_copynum.tsv" -d "../example_data/exp_species.genes_depth.tsv" -m 5 -n 5 --GRM ../example_output/exp_species.GRM.tsv --write_csv TRUE
 ```

## 2023-10-04: prep_genes split into two functions with a command line call
## 2023-10-05: prep_snps split into two functions with a command line call
## 2023-10-05: also added log file output for both functions
## 2023-10-10: updated function and added output file for structure test example of pop structre test 
```
Rscript pop_structure_test_command_line.R --species_id exp_species -v FALSE -o ../example_output --GRM ../example_output/exp_species.GRM.tsv --metadata ../example_data/IBD_metadata.tsv --tau 1 --n_tau 10

```

## 2023-10-11: added marker_test_command_line example of function below:
```
Rscript marker_test_command_line.R --species_id exp_species -v FALSE -o ../example_output --Rdata ../example_output/exp_species.model_obj.Rdata --copy_number ../example_output/exp_species.copynumber_data.csv

```
## 2023-10-13: added conda install and reconfigured folder to work better with r library
