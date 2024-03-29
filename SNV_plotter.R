# QC Script re: RedDog Het/Hom ratios
#
# A QC script for summarising and plotting
# het and hom SNPs location in excluded
# and included regions

# Import required packages
library(vcfR)
library(ape)
library(dplyr)
library(stringr)
library(tidyr)
library(reshape2)
library(data.table)
library(ggplot2)


#------------------------------------------------
# Last updated Aug 20 2017
#
# TO DO:
#       - reshape the summary table to be counts for
#         each type and location and add in ratios
#         for each
#       - work out why the second plot
#         (commented out) crashes the pdf, or change
#         to ggplot
#
#------------------------------------------------

#----------------------Example usage:------------
#
# After sourcing the script provide the parameters below and run.
#
#       ref_id is what is listed in the RedDog output (the reference genome),
#       refernce genome must be supplied as fasta. Excluded regions should be
#       supplied as a csv with one region per line and no headers.
#
#-----Example parameters:---------
#
#    working_dir_path <- "/Users/zdyson/Dropbox/Folder/15052017_QC_Script"
#
#    vcf_path <- "vcf"
#
#    ids <- c("20087_8#59","20087_8#60","20087_8#63",
#         "20087_8#64","20087_8#68","20087_8#69",
#         "20087_8#70","20087_8#72")
#
#    ref_id <- "AL513382"
#
#    ref_genome_path <- "test_data_reference/Typhi_CT18.gbk.fasta"
#
#    excluded_regions_path <- "test_data_excluded_regions/CT18_final_repeats_phage_forParseSNPTable.csv"
#
#--------To run:------------------
#
#    plot_hets(working_dir_path, vcf_path, ids, ref_id, ref_genome_path, excluded_regions_path)
#
#--------------------------------



plot_hets <- function(working_dir_path, vcf_path, ids, ref_id, ref_genome_path, excluded_regions_path){

  # Set working directory
  setwd(working_dir_path)

  # Mapped isolates to be examined
  list_of_ids <- ids

  # Set chromosome that reads have been
  # mapped to ino order to produce the
  # vcf files
  reference_id <- ref_id

  # Parse reference and annotations
  ref_genome <- ape::read.dna(ref_genome_path,
                              format = "fasta")

  excluded_regions <- read.table(
    excluded_regions_path,
    sep=",")
  colnames(excluded_regions) <- c("region_start","region_end")

  # Initialise graphical devices
  pdf("compiled_snp_distributions_in_excluded_and_non_excluded_regions.pdf", width=15, height=20)

  # Analyse the het/hom SNPs for each isolate
  for (name in 1:length(list_of_ids)){ # should multithread this with for each

    # Parse het and hom vcf files
    vcf_het <- read.vcfR(
      paste0(vcf_path, "/",
             list_of_ids[name],"_",
             reference_id,"_het.vcf"),
      verbose = TRUE )

    vcf_hom <- read.vcfR(
      paste0(vcf_path, "/",
             list_of_ids[name],"_",
             reference_id,"_q30.vcf"),
      verbose = TRUE )

    # Create tables of het and hom snps from vcf
    het_vcf_df <- tbl_df(vcf_het@fix)
    hom_vcf_df <- tbl_df(vcf_hom@fix)


    # create new columns for snp type + code and location in (un)excluded region(s)
    hom_vcf_df <- mutate(hom_vcf_df, SNP_type="homozygous",
                         SNP_code="1",SNP_region="undetermined")
    het_vcf_df <- mutate(het_vcf_df, SNP_type="heterozygous",
                         SNP_code="2",SNP_region="undetermined")

    # merge het and hom snps into a single table
    het_hom_vcf_df <- rbind(het_vcf_df, hom_vcf_df)

    # remove indels
    het_hom_vcf_df <- het_hom_vcf_df %>%
      filter(!grepl('INDEL', INFO))

    # Get details of isolate mapping
    isolate_name <- list_of_ids[name]
    chromosome_id <- het_hom_vcf_df$CHROM[1]

    # Pull out DP4 and calculate prop reads mapped to reference and altenrate
    # and output a tidy summary table of this information
    het_hom_vcf_df <- het_hom_vcf_df  %>%
      # Extract and create a new column for DP4 from INFO
      mutate(DP4=sub(".*?DP4=(.*?);.*", "\\1",
                     het_hom_vcf_df$INFO)) %>%
      # Create a new column with the name of the isolate
      mutate(Isolate_ID=isolate_name) %>%
      # Split DP4 across 4 columns
      separate(DP4, c("DP4_fwd_ref","DP4_rev_ref",
                      "DP4_fwd_alt","DP4_rev_alt"),
               sep=",") %>%
      # Calculate the total proportion of reads mapped
      # to each SNP
      mutate(total_DP4=(as.numeric(DP4_fwd_ref)+
                          as.numeric(DP4_rev_ref)+
                          as.numeric(DP4_fwd_alt)+
                          as.numeric(DP4_rev_alt))) %>%
      # Determine proportion of reads mapped to ref and alt
      mutate(prop_ref=
               as.numeric(((as.numeric(DP4_fwd_ref)+
                              as.numeric(DP4_rev_ref))/
                             as.numeric(total_DP4))),
             prop_alt=
               as.numeric(((as.numeric(DP4_fwd_alt)+
                              as.numeric(DP4_rev_alt))/
                             as.numeric(total_DP4)))) %>%
      # Remove commas from ALT column
      mutate(ALT=gsub(",", ".", ALT)) %>%
      # Select required columns for tidy data frame
      select(Isolate_ID, CHROM, POS, REF, ALT,
             prop_ref, prop_alt, SNP_region,
             SNP_type, SNP_code) %>%
      # Reorder data frame by SNP position
      arrange(as.numeric(POS))

    # build vector of excluded regions
    all_excuded_regions <- NULL
    for (excluded_region in 1:nrow(excluded_regions)){
      all_excuded_regions <- c(all_excuded_regions, c( excluded_regions[excluded_region,1]: excluded_regions[excluded_region,2]))
    }

    # Annotate SNPs from excluded regions and create tables
    for (snp in 1:nrow(het_hom_vcf_df)) {
      if (as.numeric(het_hom_vcf_df[snp,3]) %in% all_excuded_regions){
        het_hom_vcf_df[snp,8] <- "excluded"
      }else{
        het_hom_vcf_df[snp,8] <- "included"
      }
    }

    # Create an included regions only summary table
    het_hom_vcf_included_only_df <- het_hom_vcf_df[het_hom_vcf_df$SNP_region == "included",]

    # Create an included regions only het SNP summary table
    het_hom_vcf_included_het_only_df <-
      het_hom_vcf_df[het_hom_vcf_df$SNP_region == "included" &
                       het_hom_vcf_df$SNP_type=="heterozygous",]

    # Create an included regions only hom SNP summary table
    het_hom_vcf_included_hom_only_df <-
      het_hom_vcf_df[het_hom_vcf_df$SNP_region == "included" &
                       het_hom_vcf_df$SNP_type=="homozygous",]


    #--------Pannel figure of all SNPs------------------
    par(mar=c(4,4,4,4))
    par(oma=c(2,2,2,2))
    par(mai=c(1,1,1,1))
    par(mfrow=c(2,2))

    # Plot all SNPs
    plot(het_hom_vcf_df$prop_ref~het_hom_vcf_df$POS,
         xlim=c(0,length(ref_genome)), ylim=c(0,1),
         col="black",
         xlab="Position in genome",
         ylab="Percentage of reads",
         main=paste("Distribution of all SNPs in",
                    isolate_name))
    points(het_hom_vcf_df$prop_alt~het_hom_vcf_df$POS, col="red")


    # Plot all included SNPs
    plot(het_hom_vcf_included_only_df$prop_ref~het_hom_vcf_included_only_df$POS,
         xlim=c(0,length(ref_genome)), ylim=c(0,1),
         col="black",
         xlab="Position in genome",
         ylab="Percentage of reads",
         main=paste("Distribution of all non-excluded SNPs in",
                    isolate_name))
    points(het_hom_vcf_included_only_df$prop_alt~het_hom_vcf_included_only_df$POS,
           col="red")

    # Plot only included het SNPs
    plot(
      het_hom_vcf_included_het_only_df$prop_ref~het_hom_vcf_included_het_only_df$POS,
      xlim=c(0,length(ref_genome)), ylim=c(0,1),
      col="black",
      xlab="Position in genome",
      ylab="Percentage of reads",
      main=paste("Distribution of all non-exclued heterozygous SNPs in",
                 isolate_name))
    points(
      het_hom_vcf_included_het_only_df$prop_alt~het_hom_vcf_included_het_only_df$POS,
      col="red")

    # Plot only included hom SNPs
    plot(het_hom_vcf_included_hom_only_df$prop_ref~het_hom_vcf_included_hom_only_df$POS,
         xlim=c(0,length(ref_genome)), ylim=c(0,1),
         col="black",
         xlab="Position in genome",
         ylab="Percentage of reads",
         main=paste("Distribution of all non-excluded homozygous SNPs in",
                    isolate_name))
    points(het_hom_vcf_included_hom_only_df$prop_alt~het_hom_vcf_included_hom_only_df$POS,
           col="red")

    #-----------------Figure of SNP density--------------------------
    # Set up plot layout
    layout(matrix(c(1:4), nrow=2, byrow=TRUE),
           widths=c(0.80,0.20), heights=c(0.20,0.80))
    par(mar=c(2,2,0,0))
    par(oma=c(0,0,4,0))

    # Plot top histogram
    par(mai=c(0,0.4,0,0))
    hom_density <- NULL
    het_density <- NULL
    if (length(hom_vcf_df$POS) > 1){
      hom_density <- density(as.numeric(hom_vcf_df$POS))
    }
    if (length(het_vcf_df$POS) > 1){
      het_density <- density(as.numeric(het_vcf_df$POS))
    }
    plot(hom_density, col="green",
         xlim=c(0,length(ref_genome)),
         ylim=(c(0,max(c(hom_density$y, het_density$y)))),
         type="l",
         lwd=2,
         main="",
         xlab="",
         ylab="",
         bty="n", xaxt="n", yaxt="n")
    lines(het_density, col="orange",
          xlim=c(0,length(ref_genome)), type="l",
          lwd=2)
    text("top", "hets", col="orange")
    text("top", "homs", col="green")

    # empty
    plot(1:10,1:10,bty="n", xaxt="n", yaxt="n", xlab="", ylab="",col="white")

    # Plot all SNPs
    par(mai=c(0.4,0.4,0,0))
    plot(het_hom_vcf_df$prop_ref~het_hom_vcf_df$POS,
         xlim=c(0,length(ref_genome)), ylim=c(0,1),
         col="black",
         xlab="Position in genome",
         ylab="Percentage of reads",
         main="")
    points(het_hom_vcf_df$prop_alt~het_hom_vcf_df$POS, col="red")

    # Plot right histogram
    par(mai=c(0.4,0,0,0))
    prop_ref_density <- density(het_hom_vcf_df$prop_ref)
    plot(prop_ref_density$y, prop_ref_density$x, col="black",
         type="l", ylim=c(0,1), lwd=2,xlab="",ylab="",bty="n",
         xaxt="n", yaxt="n")
    prop_alt_density <- density(het_hom_vcf_df$prop_alt)
    lines(prop_alt_density$y, prop_alt_density$x, col="red",
          lwd=2)

    title(paste("Distribution of all SNPs in",isolate_name), outer=T)

    # Write the summary tables to file (fwrite does a better job with large tables)
    fwrite(het_hom_vcf_df,
           paste0(isolate_name,"_processed_vcf_files_long_format.csv"),
           quote=F, row.names=F)

    # Merge loci without SNPs into data tables
    all_loci <- tbl_df(1:length(ref_genome))

    colnames(all_loci) <- "POS"

    het_hom_vcf_df <- merge(het_hom_vcf_df, all_loci, by="POS", all=T)

    het_hom_vcf_df$SNP_type[is.na(
      het_hom_vcf_df$SNP_type)] <- "no_snp"

    het_hom_vcf_df$SNP_code[is.na(
      het_hom_vcf_df$SNP_code)] <- "0"

    het_hom_vcf_df$Isolate_ID[is.na(
      het_hom_vcf_df$Isolate_ID)] <- isolate_name

    het_hom_vcf_df$CHROM[is.na(
      het_hom_vcf_df$CHROM)] <- chromosome_id

    # Sort by SNP position
    het_hom_vcf_df <- arrange(het_hom_vcf_df, as.numeric(POS))

    # Create a matrix of SNPs
    het_hom_vcf_df$POS <- factor(het_hom_vcf_df$POS, c(1:length(ref_genome)))
    het_hom_matrix <- dcast(data=het_hom_vcf_df, formula=POS~Isolate_ID)

    # Write the summary tables to file (fwrite does a better job with large tables)
    fwrite(het_hom_matrix,
           paste0(isolate_name,"_processed_vcf_matrix.csv"),
           quote=F)
  }

  # Switch off graphical device and write plots to file
  dev.off()

  # Compile data files
  compiled_data <- do.call("rbind",
                           lapply(list.files(pattern="*long_format.csv"),
                                  read.csv, header=T))

  # Compile matrix by lines (as strings) to save memory
  compiled_matrix <- do.call("cbind",
                             lapply(list.files(pattern="*vcf_matrix.csv"),
                                    read.csv, header=T, row.names=1, check.names=F))

  # Summarise SNP counts by region and type
  summary_table <- summarise(group_by(compiled_data, Isolate_ID, SNP_type, SNP_region), count=n())

  # Write compiled and summarised data tables to file
  fwrite(compiled_data, "compiled_data.csv", quote=F)
  fwrite(compiled_matrix, "compiled_matrix.csv", quote =F,row.names=T)
  fwrite(summary_table, "compiled_summary_table.csv", quote=F)

  print("Finished")
}
