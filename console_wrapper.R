
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(arrow))

parser <- ArgumentParser()

parser$add_argument("-c", "--counts",
                      help="Path to a .csv file containing the counts")
parser$add_argument("-m", "--meta",
                    help="Path to a file containing patient meta data")
parser$add_argument("-n", "--network",
                    help="Path to a reference network")
parser$add_argument("-o", "--output",
                    help="Path to the output file")

parser$add_argument( "--tf_col",
                    help="Name of the transcription factor column in the reference network")

parser$add_argument( "--tg_col",
                     help="Name of the target gene column in the reference network")

args <- parser$parse_args()

count_path <- args$counts
network_path <- args$network
meta_path <- args$meta
output_path <- args$output
tf_col <- args$tf_col
tg_col <- args$tg_col

#count_path <- "/home/johannes/Resources/DysRegNet/PANCAN/BRCA/expected_counts.csv"
#network_path <- "/home/johannes/Resources/DysRegNet/reference_networks/HTRIdb_data.csv" 
#meta_path <- "/home/johannes/Resources/DysRegNet/PANCAN/BRCA/meta.csv"
#output_path <- "/home/johannes/Resources/DysRegNet/PANCAN/BRCA/nb_edges_exp"
#tf_col <- "SYMBOL_TF"
#tg_col <- "SYMBOL_TG"




suppressPackageStartupMessages(source("DysRegNet.R"))


# read and format data
raw_counts <- as.matrix(fread(count_path), rownames = 1)

network <- fread(network_path)
network$tf <- network[,..tf_col]
network$tg <- network[,..tg_col]
network <- network[,.(tf,tg)]

meta <- fread(meta_path)
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample
meta$sample <- NULL
meta <- meta[colnames(raw_counts),]

#run DysRegNetR

drnr <- DysRegNetR(counts=raw_counts,
                   network=network,
                   meta=meta,
                   design = tg ~ log(tf),
                   con_col = "condition",
                   norm_method = "mrn",
                   model_type = "negative-binomial",
                   cor_filter = F,
                   gof_filter = F,
                   norm_filter = F,
                   residual_type = "deviance",
                   cooks_filter = T,
                   case_for_disp = F)

drnr <- run(drnr)

# print summary

print(drnr)


# get results
result <- getResults(drnr)

# format results for feather output

edge_names <- sapply(strsplit(rownames(result),":"), function(pair){
  concat <- paste0("('", pair[1], "', '", pair[2], "')") 
})

rownames(result) <- edge_names
result <- t(result)
result<- as.data.table(result, keep.rownames = "patient id")
arrow::write_feather(result, output_path)

