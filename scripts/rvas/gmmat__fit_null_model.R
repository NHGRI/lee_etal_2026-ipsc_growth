#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(1234)

##################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
    optparse::make_option(
      "--metadata",
      type = "character",
      default = '',
      help = "Tab-delimited metadata file"
    ),
    
    optparse::make_option(
      "--trait",
      type = "character",
      default = '',
      help = "Trait to model"
    ),
    
    optparse::make_option(
      "--covariates",
      type = "character",
      default = '',
      help = "
        Comma-separated list of covariates. Covariates must appear in metadata.
        Any column not in `covariates_discrete` will be cast to continuous.
      "
    ),
    
    optparse::make_option(
      "--covariates_discrete",
      type = "character",
      default = '',
      help = "Comma-separated list of discrete covariates."
    ),
    
    optparse::make_option(
      "--grm_matrix",
      type = "character",
      default = 'grm.mtx',
      help = "GRM matrix to model"
    ),
    
    optparse::make_option(
      "--grm_id",
      type = "character",
      default = 'grm.grm.id',
      help = "GRM IDs to model"
    ),
    
    optparse::make_option(
      "--additional_matrix_covariates",
      type = "character",
      default = NULL,
      help = "Comma-separated list of additional matrix covariates to include"
    ),
    
    optparse::make_option(
      "--trait_type",
      type = "character",
      default = NULL,
      help = "Trait type. Currently only binary and quantitative are supported."
    ),

    optparse::make_option(
      "--out_file",
      type = "character",
      default = "null_model.rds",
      help = "Output R object with the fit null model."
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

######################## Required Packages #####################################
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GMMAT))
suppressPackageStartupMessages(library(Matrix))
################################################################################

################################ Functions #####################################
# build_grm <- function(grm_file, id_file) {
#   grm <- data.table::fread(
#     grm_file,
#     sep='\t',
#     header=F,
#     col.names = c('sample1', 'sample2', 'n_snps', 'relate_score')
#   )
#   ids <- read.csv(id_file, sep='\t', header = F)$V1 
#   
#   # Build matrx
#   grm <- grm %>%
#     dplyr::mutate(
#       sample1 = ids[as.integer(sample1)],
#       sample2 = ids[as.integer(sample2)]
#     ) %>%
#     tidyr::pivot_wider(
#       id_cols = 'sample1',
#       names_from = 'sample2',
#       values_from = 'relate_score'
#     ) %>%
#     tibble::column_to_rownames('sample1') %>%
#     as.matrix(.)
#   
#   n <- ncol(grm)
#   for (i in 1:n) {
#     k <- 1
#     while (k < i) {
#       grm[k, i] <- grm[i, k]
#       k <- k + 1 
#     }
#   }
#   return(grm)
# }

get_grm_mtx <- function(mtx_file, id_file) {
  grm <- Matrix::readMM(mtx_file)
  
  # read IDs and assign
  ids <- read.csv(id_file, sep='\t', header = F)$V1
  rownames(grm) <- ids
  colnames(grm) <- ids
  return(grm)
}

cast_covariates <- function(
    df,
    cols,
    cast_func,
    cast_func_description,
    verbose = T
) {
  if (verbose) {
    print(sprintf("Casting columns to be %s...", cast_func_description))
  }
  for (col in cols) {
    if (!(col %in% colnames(df))) {
      print(sprintf("Column `%s` not in dataframe. Skipping...", col))
      next()
    }
    df[col] <- cast_func(df[[col]])
  }
  return(df)
}
################################################################################

######################## Read Data & Manipulate ################################
verbose <- arguments$options$verbose
out_file <- arguments$options$out_file

# Prepare metadata #############################################################
# Get covariates -- cast any not marked as discrete as continuous
covs <- strsplit(
  x = arguments$options$covariates,
  split = ',',
  fixed = T
)[[1]]

covs_disc <- strsplit(
  x = arguments$options$covariates_discrete,
  split = ',',
  fixed = T
)[[1]]

metadata <- read.csv(
  arguments$options$metadata,
  sep='\t',
  header=T
)

# continuous
metadata <- cast_covariates(
  metadata,
  setdiff(covs, covs_disc),
  as.numeric,
  'Numeric',
  verbose = T
)

# discrete
metadata <- cast_covariates(
  metadata,
  covs_disc,
  as.character,
  'Character',
  verbose = T
) 

# Prepare matrices #############################################################
grm <- get_grm_mtx(arguments$options$grm_matrix, arguments$options$grm_id)
matrices <- list(grm)
if (!is.null(arguments$options$additional_matrix_covariates)) {
  files <- strsplit(
    x = arguments$options$additional_matrix_covariates,
    split = ',',
    fixed = T
  )[[1]]
  
  for (file in files) {
    mtx <- as.matrix(read.csv(
      file = file,
      sep = '\t',
      header = T,
      row.names = T
    ))
    matrices <- c(matrices, mtx)
  }
}

# Get trait type ###############################################################
if (arguments$options$trait_type == 'binary') {
  dist_family <- binomial(link = "logit")
} else if (arguments$options$trait_type == 'quantitative') {
  dist_family <- gaussian(link = "identity")
} else {
  stop('Trait type not supported...')
}

# Fit null model ###############################################################
form <- formula(sprintf(
  '%s ~ %s',
  arguments$options$trait,
  paste0(covs, collapse = ' + ')
))

null_model <- GMMAT::glmmkin(
  form,
  data = metadata,
  kins = matrices,
  id = "sample_id",
  family = dist_family
)

saveRDS(object = null_model, file = out_file)

if (verbose) {
    cat("Done.\n")
}