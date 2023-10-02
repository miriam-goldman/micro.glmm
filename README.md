# micro-glmm
## Purpose to Make an deploy Micro-glmm R package and files to run it
2023-09-27: Made Repository
2023-10-02: updated function prep-snps; example of how to run
'''
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
 '''
 still to do: make sure qp and pangeomne run (waiting on new pangeome)
 
 Next update: prep genes