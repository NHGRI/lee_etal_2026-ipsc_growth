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
      "--gds_file",
      type = "character",
      default = '',
      help = "Genotype GDS file"
    ),
    
    optparse::make_option(
      "--mask_file",
      type = "character",
      default = NULL,
      help = "Mask file"
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

# Run GWAS #####################################################################
rez <- GMMAT::SMMAT(
  null_mod,
  group.file = arguments$options$mask_file,
  geno.file = arguments$options$gds_file,
  MAF.range = c(0, 0.5), # we take care of MAF in group file
  miss.cutoff = 1,
  method = "davies",
  tests = c("O", "E") # default to run everything
)

saveRDS(rez, file = sprintf('%s.rds', arguments$options$out_file))

write.table(
  rez,
  file = arguments$options$out_file,
  sep = '\t',
  row.names = F,
  col.names = T
)

################################################################################

if (verbose) {
    cat("Done.\n")
}