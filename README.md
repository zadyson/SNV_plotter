# Spatial_SNV_plotter
Plots the distribution of heterozygous (het) and homozygous (hom) SNV calls from RedDog mapping pipeline output (see https://github.com/katholt/reddog for details).  Generally we use this for assessing the spatial distribution of het SNVs and if these might be associated with either (a) a mixed sample or (b) recombination.  For (a) these would be randomly distributed throughtout the genome sequences, and for (B) these would be clustered in specific regions in the sequence.  We generally would assess this for sequences that have high levels of het SNVs compared to hom SNVs. &nbsp;

## What is in this repo?
**SNV_plotter.R** - R script for creating plots and summary tables suitable for sequence QC analysis using both high quality (q30.vcf) and heterozygous (het.vcf) files generated usign RedDog.  &nbsp;

**SNV_plotter.R** - Same as the above script but it is much faster as it only creates plots &nbsp;

#### Example output (2 pages of figures per sample):

![image](https://user-images.githubusercontent.com/8507671/144292933-4795321d-68d1-43dd-a8d5-ddf5bd431e18.png)

![image](https://user-images.githubusercontent.com/8507671/144292977-129cdfbf-9788-46f3-a525-8b9f42385793.png)
