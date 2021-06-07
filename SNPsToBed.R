
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, Xmisc, cowplot, ashr)

parser <- ArgumentParser$new()
parser$add_description("Script convert a list of SNPs and loadings into a bed file for merging. Contains additional functionality to filter by some metric, such as SDs from the mean or percentage")
parser$add_argument("--mapping", type = 'character', help = "Matching the RSIDs to the chr:pos ids")
parser$add_argument("--snp_list", type = 'character', help = "List of SNPs we actually want- accepting in RSID format for now.", default = "")
parser$add_argument("--loading", type = 'character', help = "Specify the loadings file", default = "")
parser$add_argument("--loading_num", type = 'numeric', help = "Specify the loading number", default = 1)
parser$add_argument("--percent", type = "numeric", help = "Specify the top quantile of loadings to get, as a proportion (between 0 and 1)", default = 1)
parser$add_argument("--sd", type = "numeric", help = "Specify the number of standard deviations above mean to get", default = 0)
parser$add_argument("--outdir", type  = "character", help = "Specify an output directory")
parser$add_argument("--rownames", type  = "logical", help = "Specify if there are rownames on loadings file", default = F, action = "store_true")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--test_run',type='logical',action='store_true',help='Do a test run of the program.', default = F)
parser$helpme()
args <- parser$get_args()

list_in <- args$mapping #maps RSIDs to reuluar ids
optional_data <- args$loading #loadings if you want it
add_col <- args$loading_num
snp_ids <- args$snp_list
if(args$test_run)
{
  list_in <- "/Users/ashton/Documents/JHU/Research/LocalData/snp_network/gsea_scripts/seed2_thresh0.9_h2-0.1_vars1e-5.pruned_rsids.txt"
  optional_data <- "/Users/ashton/Documents/JHU/Research/LocalData/snp_network/seed2_thresh0.9_h2-0.1_vars1e-5/loadings.flashr.2.txt"
  snp_ids <- "/Users/ashton/Documents/JHU/Research/LocalData/snp_network/seed2_thresh0.9_h2-0.1_vars1e-5/snp.ids..txt"
  add_col <- 3
}

#read in list of SNPs
id_list <- fread(list_in, header = FALSE)
colnames(id_list) <- c("ID", "RSID")
snp_list <- fread(snp_ids)
if(optional_data != "")
{
  opt_dat <- fread(optional_data)
  if(args$rownames)
  {
    opt_dat <-  opt_dat[,-1]
  }
    add_col <- unlist(opt_dat[, ..add_col])
  join_list <- data.frame("RSID" = snp_list$ids, "loadings" = add_col)
  #if we are filtering
  if(args$percent != 1){
    keep_count <- floor(args$percent * nrow(join_list))
    sel_list <-  (join_list %>% arrange(-abs(loadings)))[1:keep_count, ]
    join_list <- filter(join_list, RSID %in% sel_list$RSID)
    
  }  else if(args$sd != 0){
    sd <- sd(join_list$loadings)
    mean <- mean(join_list$loadings)
    join_list <- filter(join_list, abs(loadings) > mean + (args$sd) * sd)
  }
  else{
    join_list <- join_list #no change
  }
  
}
#assuming order the same
#PROBLEM: we don't actually have the list of snps that made the cutoff.....
#z_scores <- fread("/Users/ashton/Documents/JHU/Research/LocalData/snp_network/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.z.tsv") %>% drop_na()


#id_list <-  
#id_list <-  
if(optional_data != "")
{
  bed_out <- id_list %>% separate(ID, into = c("chr", "pos"), sep = ":", remove = F) %>% filter(RSID %in% snp_list$ids) %>%
    merge(., join_list, by = "RSID") %>% arrange(chr, pos) %>% mutate("beside" = as.numeric(pos) + 1) %>% mutate("chr" = paste0("chr", chr)) %>%
    select(chr, pos, beside, RSID, loadings)
} else{
  bed_out <- id_list %>% separate(ID, into = c("chr", "pos"), sep = ":", remove = F) %>% filter(RSID %in% snp_list$ids) %>% 
    arrange(chr, pos) %>% mutate("beside" = as.numeric(pos) + 1) %>% mutate("chr" = paste0("chr", chr)) %>%
    select(chr, pos, beside, RSID)
}

write_tsv(bed_out, args$outdir, col_names = FALSE)
