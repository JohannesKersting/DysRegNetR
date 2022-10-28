
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("-c", "--counts",
                      help="Path to a .csv file containing the counts")
parser$add_argument("-m", "--meta",
                    help="Path to a file containing patient meta data")
parser$add_argument("-n", "--network",
                    help="Path to a reference network")
parser$add_argument("-o", "--output",
                    help="Path to the output file")

args <- parser$parse_args()


suppressPackageStartupMessages(source("DysRegNet.R"))

count_path <- args$counts
network_file <- args$network
meta_file <- args$meta


raw_counts <- as.matrix(fread(count_path), rownames = 1)
network <- fread(network_file)[,V1:=NULL][,.(tf,tg)]

meta <- fread(meta_file)
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample
meta$sample <- NULL
meta <- meta[colnames(raw_counts),]

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

print(drnr)


