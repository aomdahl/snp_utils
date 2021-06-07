
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, Xmisc, cowplot, ashr)

parser <- ArgumentParser$new()
parser$add_description("Script convert a list of SNPs and loadings into a bed file for merging. Contains additional functionality to filter by some metric, such as SDs from the mean or percentage")
parser$add_argument("--snp_list", type = 'character', help = "List of SNPs we actually want- Please submit as a 2 column list mapping RSIDs to the chr:pos ids. Col1: SNP ID, Col2: RSID. tab-delimited file, no header.", default = "")
parser$add_argument("--loading", type = 'character', help = "Specify the loadings file", default = "")
parser$add_argument("--loading_num", type = 'numeric', help = "Specify the column (loading) number", default = 1)
parser$add_argument("--percent", type = "numeric", help = "Specify the top quantile of loadings to get, as a proportion (between 0 and 1)", default = 1)
parser$add_argument("--sd", type = "numeric", help = "Specify the number of standard deviations above mean to get", default = 0)
parser$add_argument("--outdir", type  = "character", help = "Specify an output directory")
parser$add_argument("--rownames", type  = "logical", help = "Specify if there are rownames on loadings file. Default is no, assumes SNPs ordered same as given in list.", default = F, action = "store_true")
parser$add_argument("--tsv_OC", type  = "logical", help = "Specify this if you want to make output in the Open Cravat TSV form, basically the same thing.", default = F, action = "store_true")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$add_argument('--test_run',type='logical',action='store_true',help='Do a test run of the program.', default = F)
parser$helpme()
args <- parser$get_args()

#read in the SNPs
list_in <- args$snp_list #maps RSIDs to reuluar ids
id_list <- fread(list_in, header = FALSE)
colnames(id_list) <- c("ID", "RSID")
#read in mapping, if given.
if(list_in == "")
{
 print("For the time being becasue I am lazy, this script requires RSID mapping inthe input. Please provide, its easy to get with your scripts...")
  quit()
}

optional_data <- args$loading #loadings if you want it
add_col <- args$loading_num
if(args$test_run)
{
  list_in <- "/Users/ashton/Documents/JHU/Research/LocalData/snp_network/gsea_scripts/seed2_thresh0.9_h2-0.1_vars1e-5.pruned_rsids.txt"
  optional_data <- "/Users/ashton/Documents/JHU/Research/LocalData/snp_network/seed2_thresh0.9_h2-0.1_vars1e-5/loadings.flashr.2.txt"
  add_col <- 3
}



if(optional_data != "")
{
  opt_dat <- fread(optional_data)
  if(args$rownames)
  {
    opt_dat <-  opt_dat[,-1]
  }
    add_col <- unlist(opt_dat[, ..add_col])
  join_list <- data.frame("RSID" = id_list$ID, "loadings" = add_col)
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

if(optional_data != "")
{
  
  bed_out <- id_list %>% separate(ID, into = c("chr", "pos"), sep = ":", remove = F) %>%
    merge(., join_list, by = "RSID") %>% arrange(chr, pos) %>% mutate("beside" = as.numeric(pos) + 1) %>% mutate("chr" = paste0("chr", chr)) %>%
    select(chr, pos, beside, RSID, loadings)
  if(args$tsv_OC){
    #chr10	8097619	+	A	T	s3	var001
    bed_out <- id_list %>% separate(ID, into = c("chr", "pos"), sep = ":", remove = F) %>%
      merge(., join_list, by = "RSID") %>% arrange(chr, pos) %>% mutate("strand" = "-", "ref" = "-", "alt" = "-") %>% mutate("chr" = paste0("chr", chr)) %>%
      select(chr, pos, strand, ref, alt, RSID, loadings)
  }
  
} else{
  bed_out <- id_list %>% separate(ID, into = c("chr", "pos"), sep = ":", remove = F) %>% 
    arrange(chr, pos) %>% mutate("beside" = as.numeric(pos) + 1) %>% mutate("chr" = paste0("chr", chr)) %>%
    select(chr, pos, beside, RSID)
  
  if(args$tsv_OC){
    #chr10	8097619	+	A	T	s3	var001
    bed_out <- id_list %>% separate(ID, into = c("chr", "pos"), sep = ":", remove = F) %>%
      merge(., join_list, by = "RSID") %>% arrange(chr, pos) %>% mutate("strand" = "-", "ref" = "-", "alt" = "-") %>% mutate("chr" = paste0("chr", chr), "loadings" = "sigtest") %>%
      select(chr, pos, strand, ref, alt, RSID, loadings)
  }
}

write_tsv(bed_out, args$outdir, col_names = FALSE)
