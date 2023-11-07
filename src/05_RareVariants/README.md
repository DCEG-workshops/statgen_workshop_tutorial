## Preprocessing steps to generate 1000G agds files for rare variants analysis

Step 1: Download genotype files in VCF format and pedigree file.
1. download_ped.r: 
Download 1000G pedigree file.
2. download_vcf.r:
Download 1000G WGS Data in VCF format.

Step 2: Generate GDS files from VCF files
3. vcf2gds.r and vcf2gds.sh:
generate GDS files from VCF files.

Step 3: Annotate GDS files using FAVORannotator
4. download_FAVORDB.r
download FAVOR database.
5. copy_GDS.r
generate a copy of GDS files to a new aGDS folder. 
Then annotate the GDS files in the aGDS folder.
6. Varinfo_gds.R, Varinfo_gds.sh
7. Annotate.R, Annotate.sh
8. vcf2gds.r, vcf2gds.sh
9. Add_QC_label.R (add a variant-QC label to be all "PASS", see Reference folder)
Details of 6-8 see the instruction from https://github.com/xihaoli/STAARpipeline-Tutorial. Note that this step requires xsv.

Step 4: Calculate sparse GRM using FastSparseGRM
1. generate bed files from gds files:
gds2bed.r, gds2bed.sh  
2. merge 22 chromosomes to one file:
mergebed.sh
3. LD pruning:
LDpruned.sh
4. generate sparse GRM: 
Install FastSparseGRM
FastSparseGRM.sh (Option 2. Run the entire pipeline with one wrapper function)

Details see the instruction from https://github.com/rounakdey/FastSparseGRM. Note that this step requires KING and PLINK.







 
