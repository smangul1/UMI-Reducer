#. /u/local/Modules/default/init/modules.sh
#module load python/2.7.3
pip install pysam
python /u/home/s/serghei/code2/rop/gprofile.py example/afterCollapsing_keepingMultiMappedSort//example_afterCollapsing_keepingMultiMappedSort.bam test/example_afterCollapsing_keepingMultiMappedSort.csv --perCategory --mouse --multi
python /u/home/s/serghei/code2/rop/gprofile.py test/example_afterCollapsing_keepingMultiMappedSort_perCategory/ --mouse
