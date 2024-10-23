#---
#title: "Distribution of N50 Contig Abundances for Co-assembly and Individual Assembly (exclusivley)"
#author: "Kate Mortensen"
#date: "5/10/2023"
#---

#---------Library Setup

library(optparse)
library(ggplot2)
library(latex2exp)

# Arguments------------------------

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="/path/to/assembly2sp_count.tsv", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="path/to/output_directory", metavar="character"),
  make_option(c("-n", "--basename"), type="character", default=NULL, 
              help="basename of output file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Variables------------------------------
input = opt$input
outdir = opt$outdir
basename = opt$basename
fig_height = 5
fig_width = 7

dir.create(outdir, showWarnings = FALSE)

# Test Variables--------------------------------

# input = '/home/kate/ADMB/DASTOOL_MAGS/results/dist_N50contig_abund.tsv'
# outdir = "/home/kate/ADMB/DASTOOL_MAGS/results/R_figs"
# basename = 'dist_N50contig_abund'

# input = '/home/kate/ADMB/METAWRAP_MAGS/results/dist_N50contig_abund.tsv'
# outdir = "/home/kate/ADMB/METAWRAP_MAGS/results/R_figs"
# basename = 'dist_N50contig_abund'

# Distribution of N50 Contig Abundance ----------------------

df <- read.csv(input, header=TRUE, sep = "\t")
colnames(df)[colnames(df) == "log10_arrcpm_10e8"] ="abund"

# Base Plot -------------------------------------------

gg <- ggplot(df, aes(fill=assembly_type, x=abund)) + 
  geom_histogram(alpha=0.8, position="identity") +
  labs(x='Abundance as log10(arrcpm*10e8)', y='Frequency',
  title='Distribution of N50 Contig Abundance')
  
# Modify theme components -------------------------------------------

gg <- gg + theme(plot.title=element_text(size=20, 
                                    lineheight=1.2),  # title
            axis.title.x=element_text(size=15),  # X axis title
            axis.title.y=element_text(size=15),  # Y axis title
            axis.text.x=element_text(size=10),  # X axis text
            axis.text.y=element_text(size=10)) +  # Y axis text
            theme_minimal() +
            theme(aspect.ratio=1)

# Modify Legend -------------------------------------------

gg <- gg + scale_fill_viridis_d(name="Assembly Type",
  labels=unname(TeX(c("co-assembly", "individual", "MAG")))) 

# gg <- gg + scale_fill_grey(name="Assembly Type",
#   labels=unname(TeX(c("co-assembly", "individual", "MAG")))) 

gg <- gg + theme(legend.title=element_text(size=12),
            legend.justification=c("right", "top"), 
            legend.background = element_blank())
gg

ggsave(paste(basename, ".pdf", sep=""), gg, height = fig_height, width = fig_width, path = outdir, device = "pdf")
ggsave(paste(basename, ".png", sep=""), gg, height = fig_height, width = fig_width, path = outdir, device = "png")


# References:
# http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
# https://r-graph-gallery.com/239-custom-layout-legend-ggplot2.html
