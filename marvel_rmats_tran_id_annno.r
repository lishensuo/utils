## https://wenweixiong.github.io/Splicing_Nomenclature

## SE
# library(tidyverse)
# rmats_se = data.table::fread("./rMATS2/fromGTF.SE.txt", data.table=F)
# head(rmats_se)
# subset(rmats_se, strand=="+")[1,6:length(rmats_se)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = 1:6) 

# subset(rmats_se, strand=="-")[1,6:length(rmats_se)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = c(5,6,3,4,1,2)) 

se_func = function(rmats_se){
	trans = apply(rmats_se, 1,function(x){
		if(x["strand"]=="+"){
			tran = paste0(x["chr"],":",x["upstreamES"],":",x["upstreamEE"],":+@",
						  x["chr"],":",x["exonStart_0base"],":",x["exonEnd"],":+@",
						  x["chr"],":",x["downstreamES"],":",x["downstreamEE"])
		} else {
			tran = paste0(x["chr"],":",x["downstreamES"],":",x["downstreamEE"],":-@",
						  x["chr"],":",x["exonStart_0base"],":",x["exonEnd"],":-@",
						  x["chr"],":",x["upstreamES"],":",x["upstreamEE"])
		}
		return(tran)
	})
	trans = gsub(" ","",trans)
	return(trans)
}
# rmats_se$tran_id = se_func(rmats_se)


## MXE
# library(tidyverse)
# rmats_mxe = data.table::fread("./rMATS2/fromGTF.MXE.txt", data.table=F)
# head(rmats_mxe)

# subset(rmats_mxe, strand=="+")[1,6:length(rmats_mxe)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = 1:8) 

# subset(rmats_mxe, strand=="-")[1,6:length(rmats_mxe)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = c(7,8,5,6,3,4,1,2)) 

mxe_func = function(rmats_mxe){
	trans = apply(rmats_mxe, 1,function(x){
		if(x["strand"]=="+"){
			tran = paste0(x["chr"],":",x["upstreamES"],":",x["upstreamEE"],":+@",
						  x["chr"],":",x["1stExonStart_0base"],":",x["1stExonEnd"],":+@",
						  x["chr"],":",x["2ndExonStart_0base"],":",x["2ndExonEnd"],":+@",
						  x["chr"],":",x["downstreamES"],":",x["downstreamEE"])
		} else {
			tran = paste0(x["chr"],":",x["downstreamES"],":",x["downstreamEE"],":-@",
						  x["chr"],":",x["2ndExonStart_0base"],":",x["2ndExonEnd"],":-@",
						  x["chr"],":",x["1stExonStart_0base"],":",x["1stExonEnd"],":-@",
						  x["chr"],":",x["upstreamES"],":",x["upstreamEE"])
		}
		return(tran)
	})
	trans = gsub(" ","",trans)
	return(trans)
}
# rmats_mxe$tran_id = mxe_func(rmats_mxe)


## RI
# library(tidyverse)
# rmats_ri = data.table::fread("./rMATS2/fromGTF.RI.txt", data.table=F)
# head(rmats_ri)

# subset(rmats_ri, strand=="+")[1,6:length(rmats_ri)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::distinct(.[1]) %>% 
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = 1:4) 

# subset(rmats_ri, strand=="-")[1,6:length(rmats_ri)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::distinct(.[1]) %>% 
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = c(4,3,2,1)) 

ri_func = function(rmats_ri){
	trans = apply(rmats_ri, 1,function(x){
		if(x["strand"]=="+"){
			tran = paste0(x["chr"],":",x["riExonStart_0base"],":",x["upstreamEE"],":+@",
						  x["chr"],":",x["downstreamES"],":",x["riExonEnd"])
		} else {
			tran = paste0(x["chr"],":",x["riExonEnd"],":",x["downstreamES"],":-@",
						  x["chr"],":",x["upstreamEE"],":",x["riExonStart_0base"])
		}
		return(tran)
	})
	trans = gsub(" ","",trans)
	return(trans)
}

# rmats_ri$tran_id = ri_func(rmats_ri)

## A5SS
# library(tidyverse)
# rmats_a5ss = data.table::fread("./rMATS2/fromGTF.A5SS.txt", data.table=F)
# head(rmats_a5ss)

# subset(rmats_a5ss, strand=="+")[1,6:length(rmats_a5ss)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::distinct(.[1]) %>% 
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = 1:5) 

# subset(rmats_a5ss, strand=="-")[1,6:length(rmats_a5ss)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::distinct(.[1]) %>% 
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = c(4,5,2,3,1)) 

a5ss_func = function(rmats_a5ss){
	trans = apply(rmats_a5ss, 1,function(x){
		if(x["strand"]=="+"){
			tran = paste0(x["chr"],":",x["longExonStart_0base"],":",x["shortEE"],"|",
						  x["longExonEnd"],":+@",
						  x["chr"],":",x["flankingES"],":",x["flankingEE"])
		} else {
			tran = paste0(x["chr"],":",x["longExonEnd"],":",x["longExonStart_0base"],"|",
						  x["shortES"],":+@",
						  x["chr"],":",x["flankingES"],":",x["flankingEE"])
		}
		return(tran)
	})
	trans = gsub(" ","",trans)
	return(trans)
}

# rmats_a5ss$tran_id = a5ss_func(rmats_a5ss)



## A3SS
# library(tidyverse)
# rmats_a3ss = data.table::fread("./rMATS2/fromGTF.A3SS.txt", data.table=F)
# head(rmats_a3ss)

# subset(rmats_a3ss, strand=="+")[1,6:length(rmats_a3ss)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::distinct(.[1]) %>% 
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = 1:5) 

# subset(rmats_a3ss, strand=="-")[1,6:length(rmats_a3ss)] %>% 
# 	t() %>% as.data.frame() %>%
# 	dplyr::distinct(.[1]) %>% 
# 	dplyr::arrange(.[1]) %>% 
# 	dplyr::mutate(idx = c(5,3,4,1,2)) 

a3ss_func = function(rmats_a3ss){
	trans = apply(rmats_a3ss, 1,function(x){
		if(x["strand"]=="+"){
			tran = paste0(x["chr"],":",x["flankingES"],":",x["flankingEE"],
						  ":+@",x["longExonStart_0base"],
						  x["chr"],"|",x["shortES"],":",x["longExonEnd"])
		} else {
			tran = paste0(x["chr"],":",x["flankingES"],":",x["flankingEE"],
						  ":+@",x["shortEE"],
						  x["chr"],"|",x["longExonEnd"],":",x["longExonStart_0base"])
		}
		return(tran)
	})
	trans = gsub(" ","",trans)
	return(trans)
}

# rmats_a3ss$tran_id = a5ss_func(rmats_a3ss)
