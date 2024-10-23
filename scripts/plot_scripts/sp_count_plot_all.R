#---
#title: "assembly spacer counts by assembly type and assembly"
#author: "Kate Mortensen"
#last update: "9/4/2023"
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

input = "/home/kmorten/ADMB/results/assembly2sp_summary_all.tsv"
outdir = "/home/kmorten/ADMB/results/R_figs"
accessions_path = "/home/kmorten/ADMB/SRR_list.txt"
accessions = scan(accessions_path, what="", sep="\n")
bin_count = 5
basename = 'sp_count_all'
fig_height = 4
fig_width = 8

dir.create(outdir, showWarnings = FALSE)

# Spacer Counts (rNA) -------------------------------

df <- read.csv(paste(input), header=TRUE, sep = "\t")
# define assembly tools 
tools = unique(df$tool)
# subset by repeat guide
df1 = df[ df$repeat_guide %in% c('rNA'), ]
# exclude certain assemblies
df1 = df1[ !(df1$assembly %in% c('coassembly', 'unbinned', accessions)), ]
# order df by spacer count
df1 = df1[order(df1$spacer_count, decreasing = TRUE), ]
# start an empty dataframe
df2 = data.frame()
# subet by MAG assembly toools
for (t in tools) {
    # subet by MAG assembly toools
    df_tmp = df1[df1$tool %in% c(t), ]
    if (nrow(df_tmp) > 1) {
        # order df by spacer count
        df_tmp = df_tmp[order(df_tmp$spacer_count, decreasing = TRUE), ]
        # reset indicies 
        rownames(df_tmp) <- NULL
        # subset topmost bins
        bins <- df_tmp$assembly[1:bin_count]
        # subset by bin
        df_new = df_tmp[1:bin_count,]
        # concat new df
        df2 = rbind(df2, df_new)
    }
}
# subset by repeat guide (again)
df3 = df[ df$repeat_guide %in% c('rNA'), ]
# subset by assemblies
df3 = df3[ df3$assembly %in% c(accessions), ]
# order df by spacer count
df3 = df3[order(df3$spacer_count, decreasing = TRUE), ]
# reset indicies 
rownames(df3) <- NULL
# concat individual assemblies 
df_final = rbind(df2, df3)

level_order = c()
for (i in unique(df_final$tool)) { level_order <- c(level_order,i)}

graph <- ggplot(df_final, aes(x=reorder(assembly, spacer_count), y=spacer_count, fill = assembly_type))+
  geom_col(color='black')+
  theme_minimal() + 
  xlab('Assembly') + 
  ylab('Spacer Count') +
  ggtitle('Assembly Spacer Count', subtitle = "(no repeat guide)") +
  labs(fill = "Assembly Type") + 
  scale_fill_grey() + 
  #scale_fill_manual('Assembly Type', values=c('coral2', 'steelblue', 'pink')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   facet_wrap( ~tool, scales="free", shrink=TRUE) 
  facet_grid(~factor(tool, level_order), 
            scales = "free_x", # Let the x axis vary across facets.
            space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
            switch = "x")      # Move the facet labels to the bottom.
graph 

ggsave(paste(basename, "_rNA_all.pdf", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "pdf")
ggsave(paste(basename, "_rNA_all.png", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "png")

# Spacer Counts (rALL) -------------------------------

df <- read.csv(paste(input), header=TRUE, sep = "\t")
# define assembly tools 
tools = unique(df$tool)
# subset by repeat guide
df1 = df[ df$repeat_guide %in% c('rNA'), ]
# exclude certain assemblies
df1 = df1[ !(df1$assembly %in% c('coassembly', 'unbinned', accessions)), ]
# order df by spacer count
df1 = df1[order(df1$spacer_count, decreasing = TRUE), ]
# start an empty dataframe
df2 = data.frame()
# subet by MAG assembly toools
for (t in tools) {
    # subet by MAG assembly toools
    df_tmp = df1[df1$tool %in% c(t), ]
    if (nrow(df_tmp) > 1) {
        # order df by spacer count
        df_tmp = df_tmp[order(df_tmp$spacer_count, decreasing = TRUE), ]
        # reset indicies 
        rownames(df_tmp) <- NULL
        # subset topmost bins
        bins <- df_tmp$assembly[1:bin_count]
        # subset by bin
        df_new = df_tmp[1:bin_count,]
        # concat new df
        df2 = rbind(df2, df_new)
    }
}
# subset by repeat guide (again)
df3 = df[ df$repeat_guide %in% c('rALL'), ]
# subset by assemblies
df3 = df3[ df3$assembly %in% c(accessions), ]
# order df by spacer count
df3 = df3[order(df3$spacer_count, decreasing = TRUE), ]
# reset indicies 
rownames(df3) <- NULL
# concat individual assemblies 
df_final = rbind(df2, df3)

level_order = c()
for (i in unique(df_final$tool)) { level_order <- c(level_order,i)}

graph <- ggplot(df_final, aes(x=reorder(assembly, spacer_count), y=spacer_count, fill = assembly_type))+
  geom_col(color='black')+
  theme_minimal() + 
  xlab('Assembly') + 
  ylab('Spacer Count') +
  ggtitle('Assembly Spacer Count', subtitle = "(repeat guides)")  +
  labs(fill = "Assembly Type") + 
  scale_fill_grey() + 
  #scale_fill_manual('Assembly Type', values=c('coral2', 'steelblue', 'pink')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   facet_wrap( ~tool, scales="free", shrink=TRUE) 
  facet_grid(~factor(tool, level_order), 
            scales = "free_x", # Let the x axis vary across facets.
            space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
            switch = "x")      # Move the facet labels to the bottom.
graph 

ggsave(paste(basename, "_rALL_all.pdf", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "pdf")
ggsave(paste(basename, "_rALL_all.png", sep=""), graph, height = fig_height, width = fig_width, path = outdir, device = "png")
