

profile_dir="/mnt/cz/ip"
bam_dir="/mnt/cz/fb"

~/bin/time-1.9/time -v inStrain compare -o /mnt/cz/ibd/instrain_compare_output \
  --breadth 0.4 -s /mnt/cz/ibd/bowtie2_indexes/repgenomes.stb --database_mode \
  --skip_plot_generation --skip_popANI -p 32 --group_length 1000000 \
  -i $(cat i_string) \
  --bams $(cat b_string) &> /mnt/cz/ibd/instrain_log/compare_32.log

#https://instrain.readthedocs.io/en/latest/user_manual.html
