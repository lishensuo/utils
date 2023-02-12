# custom receiver celltype's different genes list when NicheNet  analysis
nichenet_seuratobj_aggregate2 = function(receiver, seurat_obj, condition_colname, condition_oi, condition_reference, sender = "all",ligand_target_matrix,lr_network,weighted_networks,
                                        expression_pct = 0.10, geneset_oi, filter_top_ligands = TRUE ,top_n_ligands = 20,
                                        top_n_targets = 200, cutoff_visualization = 0.33,
                                        organism = "human",verbose = TRUE, assay_oi = NULL)
{
  requireNamespace("Seurat")
  requireNamespace("dplyr")

  # input check
  if(! "RNA" %in% names(seurat_obj@assays)){
    if ("Spatial" %in% names(seurat_obj@assays)){
      warning("You are going to apply NicheNet on a spatial seurat object. Be sure it's ok to use NicheNet the way you are planning to do it. So this means: you should have changes in gene expression in receiver cells caused by cell-cell interactions. Note that in the case of spatial transcriptomics, you are not dealing with single cells but with 'spots' containing multiple cells of the same of different cell types.")

      if (class(seurat_obj@assays$Spatial@data) != "matrix" & class(seurat_obj@assays$Spatial@data) != "dgCMatrix") {
        warning("Spatial Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$Spatial@data' for default or 'seurat_obj@assays$SCT@data' for when the single-cell transform pipeline was applied")
      }
      if (sum(dim(seurat_obj@assays$Spatial@data)) == 0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$Spatial@data'")
      }
    }} else {
      if (class(seurat_obj@assays$RNA@data) != "matrix" &
          class(seurat_obj@assays$RNA@data) != "dgCMatrix") {
        warning("Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data or seurat_obj@assays$SCT@data for when the single-cell transform pipeline was applied")
      }

      if ("integrated" %in% names(seurat_obj@assays)) {
        if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$integrated@data)) ==
            0)
          stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data")
      }
      else if ("SCT" %in% names(seurat_obj@assays)) {
        if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$SCT@data)) ==
            0) {
          stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$SCT@data' for data corrected via SCT")
        }
      }
      else {
        if (sum(dim(seurat_obj@assays$RNA@data)) == 0) {
          stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data'")
        }
      }
    }

  if(!condition_colname %in% colnames(seurat_obj@meta.data))
    stop("Your column indicating the conditions/samples of interest should be in the metadata dataframe")
  if(sum(condition_oi %in% c(seurat_obj[[condition_colname]] %>% unlist() %>% as.character() %>% unique())) != length(condition_oi))
    stop("condition_oi should be in the condition-indicating column")
  if(sum(condition_reference %in% c(seurat_obj[[condition_colname]] %>% unlist() %>% as.character() %>% unique())) != length(condition_reference))
    stop("condition_reference should be in the condition-indicating column")
  if(sum(receiver %in% unique(Idents(seurat_obj))) != length(receiver))
    stop("The defined receiver cell type should be an identity class of your seurat object")
  if(length(sender) == 1){
    if(sender != "all" & sender != "undefined"){
      if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
        stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
      }
    }
  } else {
    if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
      stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
    }
  }
  if(organism != "mouse" & organism != "human")
    stop("Organism should be 'mouse' or 'human'")
  # if(geneset != "DE" & geneset != "up" & geneset != "down")
  #   stop("geneset should be 'DE', 'up' or 'down'")
  if("integrated" %in% names(seurat_obj@assays)){
    warning("Seurat object is result from the Seurat integration workflow. Make sure that the way of defining expressed and differentially expressed genes in this wrapper is appropriate for your integrated data.")
  }
  # Read in and process NicheNet networks, define ligands and receptors
  if (verbose == TRUE){print("Read in and process NicheNet's networks")}
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

  if (organism == "mouse"){
    lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
    colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
    rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
    ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
    weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
  }

  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

  if (verbose == TRUE){print("Define expressed ligands and receptors in receiver and sender cells")}

  # step1 nichenet analysis: get expressed genes in sender and receiver cells

  ## receiver
  list_expressed_genes_receiver = receiver %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
  names(list_expressed_genes_receiver) = receiver %>% unique()
  expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

  ## sender
  if (length(sender) == 1){
    if (sender == "all"){
      sender_celltypes = Idents(seurat_obj) %>% levels()
      list_expressed_genes_sender = sender_celltypes %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

    } else if (sender == "undefined") {
      if("integrated" %in% names(seurat_obj@assays)){
        expressed_genes_sender = union(seurat_obj@assays$integrated@data %>% rownames(),rownames(ligand_target_matrix)) %>% union(colnames(ligand_target_matrix))
      } else {
        expressed_genes_sender = union(seurat_obj@assays$RNA@data %>% rownames(),rownames(ligand_target_matrix)) %>% union(colnames(ligand_target_matrix))
        }
    } else if (sender != "all" & sender != "undefined") {
      sender_celltypes = sender
      list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes %>% unique()
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
    }
  } else {
    sender_celltypes = sender
    list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
    names(list_expressed_genes_sender) = sender_celltypes %>% unique()
    expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  }

  # step2 nichenet analysis: define background and gene list of interest: here differential expression between two conditions of cell type of interest
  # if (verbose == TRUE){print("Perform DE analysis in receiver cell")}

  # seurat_obj_receiver= subset(seurat_obj, idents = receiver)
  # seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[[condition_colname]])
  # DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = expression_pct) %>% rownames_to_column("gene")

  # SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_receiver)

  # if(SeuratV4 == TRUE){
  #   if (geneset == "DE"){
  #     geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= lfc_cutoff) %>% pull(gene)
  #   } else if (geneset == "up") {
  #     geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= lfc_cutoff) %>% pull(gene)
  #   } else if (geneset == "down") {
  #     geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC <= lfc_cutoff) %>% pull(gene)
  #   }
  # } else {
  #   if (geneset == "DE"){
  #     geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= lfc_cutoff) %>% pull(gene)
  #   } else if (geneset == "up") {
  #     geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC >= lfc_cutoff) %>% pull(gene)
  #   } else if (geneset == "down") {
  #     geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC <= lfc_cutoff) %>% pull(gene)
  #   }
  # }


  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  if (length(geneset_oi) == 0){
    stop("No genes were differentially expressed")
  }
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  # step3 nichenet analysis: define potential ligands
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  if (length(expressed_ligands) == 0){
    stop("No ligands expressed in sender cell")
  }
  if (length(expressed_receptors) == 0){
    stop("No receptors expressed in receiver cell")
  }
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  if (length(potential_ligands) == 0){
    stop("No potentially active ligands")
  }


  if (verbose == TRUE){print("Perform NicheNet ligand activity analysis")}

  # step4 perform NicheNet's ligand activity analysis
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>%
    arrange(-pearson) %>%
    mutate(rank = rank(desc(pearson)),
           bona_fide_ligand = test_ligand %in% ligands_bona_fide)

  if(filter_top_ligands == TRUE){
    best_upstream_ligands = ligand_activities %>% top_n(top_n_ligands, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  } else {
    best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  }

  if (verbose == TRUE){print("Infer active target genes of the prioritized ligands")}

  # step5 infer target genes of the top-ranked ligands
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_targets) %>% bind_rows() %>% drop_na()
  if(nrow(active_ligand_target_links_df) > 0){
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff_visualization)
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

    order_targets = order_targets %>% intersect(rownames(active_ligand_target_links))
    order_ligands = order_ligands %>% intersect(colnames(active_ligand_target_links))

    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands,drop=FALSE] %>% t()
    p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) #+ scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))
  } else {
    vis_ligand_target = NULL
    p_ligand_target_network = NULL
    print("no highly likely active targets found for top ligands")
  }
  # combined heatmap: overlay ligand activities
  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

  rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
  p_ligand_pearson

  figures_without_legend = cowplot::plot_grid(
    p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
    p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
    align = "hv",
    nrow = 1,
    rel_widths = c(ncol(vis_ligand_pearson)+10, ncol(vis_ligand_target)))
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h")

  combined_plot = cowplot::plot_grid(figures_without_legend,
                                     legends,
                                     rel_heights = c(10,2), nrow = 2, align = "hv")

  # ligand-receptor plot
  # get the ligand-receptor network of the top-ranked ligands
  if (verbose == TRUE){print("Infer receptors of the prioritized ligands")}

  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  if (nrow(lr_network_top_matrix) > 1){
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
  } else {
    order_receptors = rownames(lr_network_top_matrix)
  }
  if (ncol(lr_network_top_matrix) > 1) {
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  } else {
    order_ligands_receptor = colnames(lr_network_top_matrix)
  }

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  dim(vis_ligand_receptor_network) = c(length(order_receptors), length(order_ligands_receptor))
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

  # bona fide ligand-receptor
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

  if (nrow(lr_network_top_df_large_strict) == 0){
    print("Remark: no bona fide receptors of top ligands")
    vis_ligand_receptor_network_strict = NULL
    p_ligand_receptor_network_strict = NULL
    lr_network_top_df_large_strict =  NULL

  } else {
    if (nrow(lr_network_top_matrix_strict) > 1){
      dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
      hclust_receptors = hclust(dist_receptors, method = "ward.D2")
      order_receptors = hclust_receptors$labels[hclust_receptors$order]
    } else {
      order_receptors = rownames(lr_network_top_matrix)
    }
    if (ncol(lr_network_top_matrix_strict) > 1) {
      dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
      hclust_ligands = hclust(dist_ligands, method = "ward.D2")
      order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
    } else {
      order_ligands_receptor = colnames(lr_network_top_matrix_strict)
    }
    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

    vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
    dim(vis_ligand_receptor_network_strict) = c(length(order_receptors), length(order_ligands_receptor))

    rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

    p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")

    lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% rename(ligand = from, receptor = to)
  }

  # DE analysis for each sender cell type -- of course only possible when having sender cell types
  if (length(sender) > 1){
    are_there_senders = TRUE
  }
  if(length(sender) == 1){
    if(sender != "undefined"){
      are_there_senders = TRUE
    } else {
      are_there_senders = FALSE
    }
  }

  if (are_there_senders == TRUE){
    if (verbose == TRUE){print("Perform DE analysis in sender cells")}
    seurat_obj = subset(seurat_obj, features= potential_ligands)

    DE_table_all = Idents(seurat_obj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seurat_obj, condition_colname = condition_colname, condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = expression_pct, celltype_col = NULL) %>% reduce(full_join, by = "gene") # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
    DE_table_all[is.na(DE_table_all)] = 0

    # Combine ligand activities with DE information
    ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene), by = "ligand")
    ligand_activities_de[is.na(ligand_activities_de)] = 0

    # make LFC heatmap
    lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
    rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

    order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
    vis_ligand_lfc = lfc_matrix[order_ligands,]
    vis_ligand_lfc = vis_ligand_lfc %>% as.matrix(ncol = length(Idents(seurat_obj) %>% levels() %>% intersect(sender_celltypes)))
    colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

    p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))

    # ligand expression Seurat dotplot
    real_makenames_conversion = lr_network$from %>% unique() %>% magrittr::set_names(lr_network$from %>% unique() %>% make.names())
    order_ligands_adapted = real_makenames_conversion[order_ligands]
    names(order_ligands_adapted) = NULL

    seurat_obj_subset = seurat_obj %>% subset(idents = sender_celltypes)
    seurat_obj_subset = SetIdent(seurat_obj_subset, value = seurat_obj_subset[[condition_colname]]) %>% subset(idents = condition_oi) ## only shows cells of the condition of interest
    rotated_dotplot = DotPlot(seurat_obj %>% subset(cells = Cells(seurat_obj_subset)), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
    rm(seurat_obj_subset)

    # combined plot
    figures_without_legend = cowplot::plot_grid(
      p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
      rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
      p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
      p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
      align = "hv",
      nrow = 1,
      rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

    legends = cowplot::plot_grid(
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
      ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
      nrow = 1,
      align = "h", rel_widths = c(1.5, 1, 1, 1))

    combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
    combined_plot

  } else {
    rotated_dotplot = NULL
    p_ligand_lfc = NULL
  }

  return(list(
    ligand_activities = ligand_activities,
    top_ligands = best_upstream_ligands,
    top_targets = active_ligand_target_links_df$target %>% unique(),
    top_receptors = lr_network_top_df_large$to %>% unique(),
    ligand_target_matrix = vis_ligand_target,
    ligand_target_heatmap = p_ligand_target_network,
    ligand_target_df = active_ligand_target_links_df,
    ligand_expression_dotplot = rotated_dotplot,
    ligand_differential_expression_heatmap = p_ligand_lfc,
    ligand_activity_target_heatmap = combined_plot,
    ligand_receptor_matrix = vis_ligand_receptor_network,
    ligand_receptor_heatmap = p_ligand_receptor_network,
    ligand_receptor_df = lr_network_top_df_large %>% rename(ligand = from, receptor = to),
    ligand_receptor_matrix_bonafide = vis_ligand_receptor_network_strict,
    ligand_receptor_heatmap_bonafide = p_ligand_receptor_network_strict,
    ligand_receptor_df_bonafide = lr_network_top_df_large_strict,
    geneset_oi = geneset_oi,
    background_expressed_genes = background_expressed_genes
  ))
}
