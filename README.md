![micro_glmm_logo_11072023](https://github.com/miriam-goldman/micro.glmm/assets/46382670/88ea180f-13d9-4787-8e2f-5655fa53fad0)
## Purpose to Make an deploy Micro.glmm R package and files to run it
<img width="950" alt="image" src="https://github.com/miriam-goldman/micro.glmm/assets/46382670/1d814a2d-ae54-4690-b70f-8353bae1457f">

# Instructions on how to run
### Install Conda

```
conda env create -f microglmm.yml
conda activate microglmm
```
### install package
From R
```
R
library(devtools)
install_github('miriam-goldman/micro.glmm')
```

### Move into src folder
```
cd src
```
### run prep snps
```
Rscript prep_snps_command_line.R  --species_id exp_species -v TRUE  -o "../example_output"
-i "../example_data/exp_species.snps_info.tsv"
-f "../example_data/exp_species.snps_freqs.tsv"
-d "../example_data/exp_species.snps_depth.tsv"
-a 3 -m 5  -n 5   
 ```
### run prep genes

```
Rscript prep_genes_command_line.R --species_id exp_species -v TRUE -o "../example_output"
-c "../example_data/exp_species.genes_copynum.tsv"
-d "../example_data/exp_species.genes_depth.tsv"
-m 5 -n 5 --GRM ../example_output/exp_species.GRM.tsv --write_csv TRUE
 ```
### run population structure test

```
Rscript pop_structure_test_command_line.R --species_id exp_species -v TRUE -o ../example_output
--GRM ../example_output/exp_species.GRM.tsv
--metadata ../example_data/IBD_metadata.tsv
--tau 1 --n_tau 10
```

### Run Maker test

```
Rscript marker_test_command_line.R --species_id exp_species -v TRUE
-o ../example_output
--Rdata ../example_output/exp_species.model_obj.Rdata
--copy_number ../example_output/exp_species.copynumber_data.csv
```
