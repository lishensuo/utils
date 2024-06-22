library(tidyverse)
library(org.Mm.eg.db)

aa = data.table::fread("Gene_basic_annotation_gencode_v23.csv")
table(is.na(aa$Symbol.MGI))

all_genes <- keys(org.Mm.eg.db, keytype = "SYMBOL")
# all_genes[1213:1215]
## [1] "Del(17)18H"       "Del(17)1t<w18>-d" "Del(17)1t<w18>-p"
all_genes_2 = grep("^[^<>]*$", all_genes, value = T)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "dec2021.archive.ensembl.org")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "dec2021.archive.ensembl.org")


genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
               values = all_genes_2, 
               mart = mouse, 
               attributesL = c("hgnc_symbol"), 
               martL = human, uniqueRows=T)

genes.list = split(genes$MGI.symbol, genes$HGNC.symbol)

genes_2 = lapply(genes.list, function(x){
  # x = genes.list[[1]]
  paste(x, collapse = ',')
}) %>% do.call(rbind, .) %>% as.data.frame() %>% 
  tibble::rownames_to_column("HGNC.symbol") %>% 
  dplyr::rename("MGI.symbol"="V1")

aa2 = aa %>% 
  dplyr::select(!Symbol.MGI) %>% 
  dplyr::left_join(genes_2, by = c("Symbol"="HGNC.symbol"))

write.csv(aa2, file = "Gene_basic_annotation_gencode_v23.csv", row.names = FALSE)
