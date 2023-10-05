# micro-glmm
## Purpose to Make an deploy Micro-glmm R package and files to run it
<img width="985" alt="image" src="https://github.com/miriam-goldman/micro-glmm/assets/46382670/7f1fde4b-7439-4c0c-a99d-bcba26fbe5f3">

### 2023-09-27: Made Repository
### 2023-10-02: updated function prep-snps; example of how to run
```
 Rscript prep_snps.R 
 --species_id 100077
  -v TRUE 
  -o "~/Projects/micro-glmm/test_output" 
  -i "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/snps/100077/100077.snps_info.tsv" 
  -f "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/snps/100077/100077.snps_freqs.tsv" 
  -d "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/snps/100077/100077.snps_depth.tsv" 
  -r FALSE 
  --centroid_prevalence_cutoff .8 
  --run_qp FALSE 
  -u 5 
  -l .2 
  -a 3
   -m 5 
   -n 5
 ```
 still to do: make sure qp and pangeomne run (waiting on new pangeome)
 
 ### 2023-10-03: updated function prep-genes; example of how to run
```
 Rscript prep_genes.R
--species_id 104158
-v TRUE
 -o "~/Projects/micro-glmm/test_output"
-c "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/genes/104158/104158.genes_copynum.tsv"
 -d "/Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/genes/104158/104158.genes_depth.tsv"
-m 5 -n 5 --genes_summary /Users/miriamgoldman/Projects/IBD_folder/midas2_data/merge/genes/genes_summary.tsv
--GRM /Users/miriamgoldman/Projects/micro-glmm/test_output/104158.GRM.tsv
--metadata ~/Projects/micro-glmm/IBD_metadata.tsv
 ```

## 2023-10-04: prep_genes split into two functions with a command line call
## 2023-10-05: prep_snps split into two functions with a command line call
## 2023-10-05: also added log file output for both functions
Today+ tomorrow add structure test function
