#---
#title: "assembly spacer counts by assembly method and repeat guide"
#author: "Kate Mortensen"
#date: "5/9/2023"
#---

#---------Library Setup

library(optparse)
library(ggplot2)

# Arguments------------------------

option_list = list(
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="path/to/output_directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="/path/to/rep2div2spcount.tsv", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Variables------------------------------

input = opt$input
outdir = opt$outdir
fig_height = 4
fig_width = 8

dir.create(outdir, showWarnings = FALSE)

# Test Variables--------------------------------

# input = "/home/kmorten/ADMB/DASTOOL_MAGS/results/rep2div2spcount.tsv"
# outdir = "/home/kmorten/ADMB/DASTOOL_MAGS/results/R_figs"
# fig_height = 4
# fig_width = 8

dir.create(outdir, showWarnings = FALSE)

# Repeat-guided Spacer Discovery----------------------

df <- read.csv(input, header=TRUE, sep = "\t")

df = df[ !(df$repeat_guide %in% c('rNA', 'rALL')), ]

graph <- ggplot(df, aes(fill=assembly_type_diversity, y=spacer_count, x=reorder(repeat_guide, -spacer_count))) + 
  geom_bar(position='stack', stat='identity') +
  theme_minimal() + 
  labs(x='Repeat Guide', y='Spacer Count', title='Repeat-guided Spacer Discovery') +
  #labs(x='Repeat Guide', y='Spacer Count') +
  #theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) +
  scale_fill_grey(name="Assembly Type",labels=c("Co-assembly", "Individual Assembly", "In Common")) + 
  #scale_fill_manual('Assembly Type', values=c('coral2', 'darkmagenta','cyan3')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

graph

ggsave("rep2div2spcount.pdf", graph, height = fig_height, width = fig_width, path = outdir, device = "pdf")
ggsave("rep2div2spcount.png", graph, height = fig_height, width = fig_width, path = outdir, device = "png")