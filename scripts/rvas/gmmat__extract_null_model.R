#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
    optparse::make_option(
      "--null_model_rds",
      type = "character",
      default = '',
      help = "RDS of null model"
    ),
    
    optparse::make_option(
      "--metadata",
      type = "character",
      default = '',
      help = "Metadata file"
    ),
    
    optparse::make_option(
      "--trait",
      type = "character",
      default = '',
      help = "Trait. Should be at least one column in `metadata`"
    ),
    
    optparse::make_option(
      "--out_file",
      type = "character",
      default = "results.tsv.gz",
      help = "Output file"
    ),

    optparse::make_option(
      "--verbose",
      action = "store_true",
      default = FALSE,
      help = ""
    )
)

parser <- optparse::OptionParser(
    usage = "%prog",
    option_list = optionList,
    description = paste0(
        "Plots results from differential gene expression."
    )
)

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
    optionStrings <- character()
    for (item in parserObj@options) {
        optionStrings <- append(optionStrings,
                                c(item@short_flag, item@long_flag))
    }
    optionStrings
}

optStrings <- getOptionStrings(parser)
arguments <- optparse::parse_args(parser, positional_arguments = TRUE)
################################################################################

######################## Required Packages #####################################=
suppressPackageStartupMessages(library(GMMAT))
################################################################################

################################ Functions #####################################

################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose

null_mod <- readRDS(arguments$options$null_model_rds)
metadata <- read.csv(
  arguments$options$metadata,
  sep='\t',
  header=T
)
trait <- arguments$options$trait

# filter metadata to samples and phenotype
cols_retain <- c(trait, sprintf('%s__untransformed', trait))
metadata <- metadata[c(
  'sample_id',
  colnames(metadata)[colnames(metadata) %in% cols_retain]
)]

# now get residuals
ids <- null_mod$id_include

resids <- null_mod$residuals
scale_resids <- null_mod$scaled.residuals
names(resids) <- ids
names(scale_resids) <- ids

# add to metadata 
metadata[sprintf('%s__residuals', trait)] <- resids[metadata$sample_id]
metadata[sprintf('%s__scaled_residuals', trait)] <- scale_resids[
  metadata$sample_id
]

gz_file <- gzfile(arguments$options$out_file, "w", compression = 9)
write.table(
  metadata,
  file = gz_file,
  sep = '\t',
  row.names = F,
  col.names = T
)
close(gz_file)

################################################################################

if (verbose) {
    cat("Done.\n")
}