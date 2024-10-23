#---
#title: "assembly spacer counts by assembly type and assembly"
#author: "Kate Mortensen"
#date: "5/10/2023"
#---

#---------Library Setup

library(forcats) # fct_inorder
library(optparse)
library(ggplot2)

# Arguments------------------------

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="/path/to/assembly2sp_count.tsv", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="path/to/output_directory", metavar="character"),
  make_option(c("-n", "--basename"), type="character", default=NULL, 
              help="basename of output file", metavar="character"),
  make_option(c("-a", "--accessions"), type="character", default=NULL, 
              help="/path/to/SRR_list.txt", metavar="character"),
  make_option(c("-c", "--topmost_bin_count"), type="integer", default=NULL, 
              help="choose the number of most abundant bins to plot", metavar="integer")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Variables------------------------------

input = opt$input
outdir = opt$outdir
accessions = scan(opt$accessions, what="", sep="\n")
bin_count = as.integer(opt$topmost_bin_count)
basename = opt$basename
fig_height = 4
fig_width = 8

dir.create(outdir, showWarnings = FALSE)

# Test Variables--------------------------------

# input = "/home/kmorten/ADMB/DASTOOL_MAGS/results/assembly2sp_summary.tsv"
# outdir = "/home/kmorten/ADMB/DASTOOL_MAGS/results/R_figs"
# accessions_path = "/home/kmorten/ADMB/SRR_list.txt"
# accessions = scan(accessions_path, what="", sep="\n")
# bin_count = 5
# basename = 'sp_count'
# fig_height = 4
# fig_width = 8

# dir.create(outdir, showWarnings = FALSE)


# Spacer Counts (rNA) -------------------------------

df <- read.csv(paste(input), header=TRUE, sep = "\t")

# subset by repeat guide
df1 = df[ df$repeat_guide %in% c('rNA'), ]
# exclude certain assemblies
df1 = df1[ !(df1$assembly %in% c('coassembly', 'unbinned', accessions)), ]
# order df by spacer count
df1 = df1[order(df1$spacer_count, decreasing = TRUE), ]
# reset indicies 
rownames(df1) <- NULL
# subset topmost bins
bins <- df1$assembly[1:bin_count]
# assemblies to include in plot 
include = c(accessions)
for (i in bins) { include <- c(include,i)}
# subset by repeat guide (again)
df1 = df[ df$repeat_guide %in% c('rNA'), ]
# subset by assemblies
df1 = df1[ df1$assembly %in% include, ]
# order df by spacer count
df1 = df1[order(df1$spacer_count, decreasing = TRUE), ]
# reset indicies 
rownames(df1) <- NULL

graph <- ggplot(df1, aes(x=reorder(assembly, spacer_count), y=spacer_count, fill = assembly_type))+
  geom_col(color='black')+
  theme_minimal() + 
  xlab('Assembly') + 
  ylab('Spacer Count') +
  ggtitle('Assembly Spacer Count', subtitle = "(no repeat guide)") +
  labs(fill = "Assembly Type") + 
  scale_fill_grey() + 
  #scale_fill_manual('Assembly Type', values=c('coral2', 'steelblue', 'pink')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  #facet_wrap( ~type, scales="free", shrink=TRUE) 
graph 

ggsave(paste(basename, "_rNA.pdf", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "pdf")
ggsave(paste(basename, "_rNA.png", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "png")

# Spacer Counts (rALL) -------------------------------

df <- read.csv(paste(input), header=TRUE, sep = "\t")

# subset by repeat guide
df1 = df[ df$repeat_guide %in% c('rALL'), ]
# exclude certain assemblies
df1 = df1[ !(df1$assembly %in% c('coassembly', 'unbinned', accessions)), ]
# order df by spacer count
df1 = df1[order(df1$spacer_count, decreasing = TRUE), ]
# reset indicies 
rownames(df1) <- NULL
# subset topmost bins
bins <- df1$assembly[1:bin_count]
# assemblies to include in plot 
include = c(accessions)
for (i in bins) { include <- c(include,i)}
# subset by repeat guide (again)
df1 = df[ df$repeat_guide %in% c('rALL'), ]
# subset by assemblies
df1 = df1[ df1$assembly %in% include, ]
# order df by spacer count
df1 = df1[order(df1$spacer_count, decreasing = TRUE), ]
# reset indicies 
rownames(df1) <- NULL

graph <- ggplot(df1, aes(x=reorder(assembly, spacer_count), y=spacer_count, fill = assembly_type))+
  geom_col(color='black')+
  theme_minimal() + 
  xlab('Assembly') + 
  ylab('Spacer Count') +
  ggtitle('Assembly Spacer Count', subtitle = "(repeat guides)") +
  labs(fill = "Assembly Type") + 
  scale_fill_grey() + 
  #scale_fill_manual('Assembly Type', values=c('coral2', 'steelblue', 'pink')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  #facet_wrap( ~type, scales="free", shrink=TRUE) 
graph 

ggsave(paste(basename, "_rALL.pdf", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "pdf")
ggsave(paste(basename, "_rALL.png", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "png")
