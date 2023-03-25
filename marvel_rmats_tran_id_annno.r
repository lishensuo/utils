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
			tran = paste0(x["chr"],":",as.numeric(x["upstreamES"])+1,":",x["upstreamEE"],":+@",
						  x["chr"],":",as.numeric(x["exonStart_0base"])+1,":",x["exonEnd"],":+@",
						  x["chr"],":",as.numeric(x["downstreamES"])+1,":",x["downstreamEE"])
		} else {
			tran = paste0(x["chr"],":",as.numeric(x["downstreamES"])+1,":",x["downstreamEE"],":-@",
						  x["chr"],":",as.numeric(x["exonStart_0base"])+1,":",x["exonEnd"],":-@",
						  x["chr"],":",as.numeric(x["upstreamES"])+1,":",x["upstreamEE"])
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
			tran = paste0(x["chr"],":",as.numeric(x["upstreamES"])+1,":",x["upstreamEE"],":+@",
						  x["chr"],":",as.numeric(x["1stExonStart_0base"])+1,":",x["1stExonEnd"],":+@",
						  x["chr"],":",as.numeric(x["2ndExonStart_0base"])+1,":",x["2ndExonEnd"],":+@",
						  x["chr"],":",as.numeric(x["downstreamES"])+1,":",x["downstreamEE"])
		} else {
			tran = paste0(x["chr"],":",as.numeric(x["downstreamES"])+1,":",x["downstreamEE"],":-@",
						  x["chr"],":",as.numeric(x["2ndExonStart_0base"])+1,":",x["2ndExonEnd"],":-@",
						  x["chr"],":",as.numeric(x["1stExonStart_0base"])+1,":",x["1stExonEnd"],":-@",
						  x["chr"],":",as.numeric(x["upstreamES"])+1,":",x["upstreamEE"])
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
			tran = paste0(x["chr"],":",as.numeric(x["riExonStart_0base"])+1,":",x["upstreamEE"],":+@",
						  x["chr"],":",as.numeric(x["downstreamES"])+1,":",x["riExonEnd"])
		} else {
			tran = paste0(x["chr"],":",as.numeric(x["riExonEnd"])+1,":",x["downstreamES"],":-@",
						  x["chr"],":",as.numeric(x["upstreamEE"])+1,":",x["riExonStart_0base"])
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
			tran = paste0(x["chr"],":",as.numeric(x["longExonStart_0base"])+1,":",x["shortEE"],"|",
						  x["longExonEnd"],":+@",
						  x["chr"],":",as.numeric(x["flankingES"])+1,":",x["flankingEE"])
		} else {
			tran = paste0(x["chr"],":",as.numeric(x["longExonEnd"])+1,":",x["longExonStart_0base"],"|",
						  x["shortES"],":-@",
						  x["chr"],":",as.numeric(x["flankingES"])+1,":",x["flankingEE"])
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
			tran = paste0(x["chr"],":",as.numeric(x["flankingES"])+1,":",x["flankingEE"],
						  ":+@",as.numeric(x["longExonStart_0base"])+1,
						  x["chr"],"|",x["shortES"],":",x["longExonEnd"])
		} else {
			tran = paste0(x["chr"],":",as.numeric(x["flankingES"])+1,":",x["flankingEE"],
						  ":-@",as.numeric(x["shortEE"])+1,
						  x["chr"],"|",x["longExonEnd"],":",x["longExonStart_0base"])
		}
		return(tran)
	})
	trans = gsub(" ","",trans)
	return(trans)
}

# rmats_a3ss$tran_id = a5ss_func(rmats_a3ss)
