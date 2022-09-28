# This script was borrowwed from https://github.com/mckellardw/scMuscle/tree/main/R_scripts

# Helper functions used in each of the other scripts

# Calculate the number of PCs that contain some proportion (95%) of the variance
#     Output: number of dimensions that contain the desired percentage of variance
npcs <- function(
  seu, 
  var.toal=0.95, 
  reduction="pca"
){
  if(is.null(seu@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  
  tmp.var <- (seu@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }
  
  return(n.pcs)
}


# Function to generate mean expression profiles for individual cell types
#     Output: genes-by-celltype matrix of mean expression profiles
get_cell_type_model <- 
  function(
    seur, #seurat object
    assay='RNA',
    slot='data',
    cluster.names, #cell type model; pass in metadata column name for desired cell types
    ignore.clusters=NULL, # vector of chars; cluster IDs to ignore in model generation
    cells.use=NULL, # vector of cell barcodes to use; works in combination with ignore.clusters, doublet.classifications
    genes.ignore=NULL, # genes to remove from model
    doublet.classifications=NULL, # name of doublet classification metadata column; default is to ignore classifications
    # nCores=1, #TODO: parallelize with dofor
    verbose=T
  ){
    require(Seurat)
    
    #build singlet_cell type_gene expression matrix
    ref.mat <- list() #initialize reference expression matrix
    exp.mat <- GetAssayData(seur, assay=assay, slot=slot) # Pull out expression data
    
    # subset based on desired cells
    if(!is.null(cells.use)){ 
      exp.mat <- exp.mat[,cells.use]
    }
    
    # remove ignored genes
    if(!is.null(genes.ignore)){
      exp.mat = exp.mat[!rownames(seur)%in%genes.ignore,]
    }
    
    #Grab meta.data, filtered for singlets
    if(is.null(cells.use)){ #TODO- subset based on cells.use
      meta <- seur@meta.data  #ignore sinlget/doublet classifications
    }else{
      meta <- seur@meta.data[cells.use,] 
    }
    
    cell.types <- sort(unlist(unique(meta[[cluster.names]]))) #Get cell type names
    if(!is.null(ignore.clusters)){
      cell.types <- cell.types[(!cell.types%in%ignore.clusters)]  
    }
    print(cell.types)
    
    # Remove NAs from exp.mat and meta
    exp.mat <- exp.mat[,!is.na(meta[[cluster.names]])]
    meta <- meta[!is.na(meta[[cluster.names]]),]
    
    # Generate cell type references, one cell type at a time
    if(verbose){cat("Building reference profiles... \n")}
    if(verbose & !is.null(cells.use)){cat("   Using the provided cells... \n")}
    
    for(i in 1:length(cell.types)){
      tmp.cells = rownames(meta)[ as.vector(meta[[cluster.names]])==cell.types[i] ]
      
      if(verbose){cat(cell.types[i]," (",i,"/",length(cell.types),")","... ", sep="")}
      
      if(length(tmp.cells)==0){
        if(verbose){cat(" NO CELLS", sep="")}
      }
      
      if(length(tmp.cells)==1){
        if(verbose){cat(" Only 1 cell!", sep="")}
        ref.mat[[i]] <- exp.mat[,rownames(meta)%in%tmp.cells]  # grab expression data
      }else{
        if(verbose){cat(length(tmp.cells)," cells.", sep="")}
        ref.mat[[i]] <- 
          Matrix::rowMeans( #TODO: should these expression profiles be built on row means?
            exp.mat[,rownames(meta)%in%tmp.cells]  # grab expression data
            
          )
      }
      
      if(verbose){cat("\n")}
    }
    ref.mat <- do.call(cbind, ref.mat)
    colnames(ref.mat) <- cell.types
    
    if(verbose){cat("Done!\n")}
    
    return(ref.mat)
  }


# DoubletFinder helper functions ####

# Simple linear regression to spit out a list of estimated doublet rates, based on 10x's published rates
estimateDoubletRate.DWM <- function(seur.list, doublet.dist=NULL){
  if(is.null(doublet.dist)){
    doublet.dist <- data.frame(cells.recovered=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
                               multi.rate=     c(0.4,  0.8,  1.6,  2.3,  3.1,  3.9,  4.6,  5.4,  6.1,  6.9,   7.6)
    )
  }
  
  fit <- lm(multi.rate~cells.recovered, data = doublet.dist)
  fit.func <- function(x){
    return(as.numeric(fit$coefficients['cells.recovered']*x + fit$coefficients['(Intercept)']))
  }
  
  ncells.list <- lapply(seur.list, ncol)
  out.list <- lapply(ncells.list, fit.func)
  return(unlist(out.list))
}

# Added in ability to set meta data column name (eases merging many datasets)
doubletFinder_V3.DWM_v2 <- 
  function(
    seu, PCs, pN = 0.25, pK, nExp='auto', 
    pANN.cutoff=NULL, # hard limit for distinguishing singlets/doublets
    assay='RNA', reduction='pca',slot='counts', 
    classification.name=NULL, pANN.name=NULL, # meta data naming parameters
    reuse.pANN = FALSE, #changed in my implementation 
    sct = FALSE 
  ) {
    require(Seurat); require(fields); require(KernSmooth)
    
    #TODO: check passed seurat object for assay, reduction, etc.
    #TODO: add parallelization with future to the Seurat section
    
    ## Generate new list of doublet classificatons from existing pANN vector to save time
    if(is.character(reuse.pANN)){ #reuse pANN values, but rename classifications
      if(is.null(seu[[reuse.pANN]])){
        cat("pANN.name '", reuse.pANN, " does not exist... Please try again. \n")
        return(seu)
      }else{
        #Doublet Classification
        cat("Classifying doublets based on previous pANN values (", reuse.pANN, ")...\n",sep = '')
        
        pANN.old <- seu@meta.data[ , reuse.pANN]
        classifications <- rep("Singlet", length(pANN.old))
        
        if(is.null(pANN.cutoff)){
          classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
        }else{
          classifications[pANN.old>=pANN.cutoff] <- "Doublet"
        }
        
        # Add in metadata columns with classifications 
        if(is.null(classification.name)){
          seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
        }else if(is.character(classification.name)){
          seu@meta.data[, classification.name] <- classifications
        }else{
          cat("Doublet classifications labeled as 'DF.classifications'...\n")
          seu@meta.data[, "DF.classifications"] <- classifications
        }
        
        #Don't need to add pANN values again
        
        return(seu)
      }
    }
    if(reuse.pANN){
      if(is.null(seu[[pANN.name]])){
        cat("pANN.name '", pANN.name, " does not exist... Please try again. \n")
        return(seu)
      }else{
        #Doublet Classification
        cat("Classifying doublets based on previous pANN values...\n")
        
        pANN.old <- seu[[pANN.name]]
        classifications <- rep("Singlet", length(pANN.old))
        
        if(is.null(pANN.cutoff)){
          classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
        }else{
          classifications[pANN.old>=pANN.cutoff] <- "Doublet"
        }
        
        # Add in metadata columns with classifications 
        if(is.null(classification.name)){
          seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
        }else if(is.character(classification.name)){
          seu@meta.data[, classification.name] <- classifications
        }else{
          cat("Doublet classifications labeled as 'DF.classifications'...\n")
          seu@meta.data[, "DF.classifications"] <- classifications
        }
        
        #Don't need to add pANN values again
        
        return(seu)
      }
    }else{
      
      ## Make merged real-artifical data
      real.cells <- rownames(seu@meta.data)
      data <- GetAssayData(seu, assay=assay, slot=slot)
      n_real.cells <- length(real.cells)
      n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
      
      cat("Creating ", n_doublets, " artificial doublets from ", n_real.cells, " cells...\n")
      
      real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
      real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
      
      doublets <- (data[, real.cells1] + data[, real.cells2])/2 # Literally taking the average of two cells...
      colnames(doublets) <- paste0("X", 1:n_doublets)
      data_wdoublets <- cbind(data, doublets)
      
      ## Store important pre-processing information
      orig.commands <- seu@commands
      
      ## Initialize Seurat object
      if(slot=='counts'){
        cat("Creating Seurat object with artificial doublets...\n")
        seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
        if(assay!='SCT'){
          seu_wdoublets <- NormalizeData(seu_wdoublets,
                                         normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                         scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                         margin = orig.commands$NormalizeData.RNA@params$margin)
        }
      }else if(slot=='data'){
        cat("Creating Seurat object with artificial doublets...\n")
        seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets) #don't renormalize if normalized data is provided
      }
      
      ## Preprocess Seurat object
      if(assay=='SCT'){
        require(sctransform)
        
        if(slot=='data'){
          cat('WARNING: running SCTransform on normalized data!\n')
        }
        
        cat("Running SCTransform & PCA...\n")
        seu_wdoublets <- SCTransform(
          seu_wdoublets
        ) %>% RunPCA(
          npcs = length(PCs),
          reduction.name='DOUBLETFINDER_PCA'
        )
        
        pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
        cell.names <- rownames(seu_wdoublets@meta.data)
        nCells <- length(cell.names)
        
        rm(seu_wdoublets); gc()
        
      }else if(assay=='RNA'){
        cat("     Piping FindVariableFeatures(), ScaleData(), and RunPCA()...\n")
        seu_wdoublets <- FindVariableFeatures(
          seu_wdoublets,
          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
        ) %>% ScaleData(
          features = orig.commands$ScaleData.RNA$features,
          model.use = orig.commands$ScaleData.RNA$model.use,
          do.scale = orig.commands$ScaleData.RNA$do.scale,
          do.center = orig.commands$ScaleData.RNA$do.center,
          scale.max = orig.commands$ScaleData.RNA$scale.max,
          block.size = orig.commands$ScaleData.RNA$block.size,
          min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block
        )%>% RunPCA(
          features = orig.commands$ScaleData.RNA$features,
          npcs = length(PCs),
          rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
          weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
          reduction.name='DOUBLETFINDER_PCA',
          verbose=FALSE
        )
        
        pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
        cell.names <- rownames(seu_wdoublets@meta.data)
        nCells <- length(cell.names)
        rm(seu_wdoublets); gc() # Free up memory
      }else{
        cat("     Piping FindVariableFeatures(), ScaleData(), and RunPCA()...\n")
        seu_wdoublets <- FindVariableFeatures(
          seu_wdoublets,
          selection.method = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$selection.method,
          loess.span = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$loess.span,
          clip.max = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$clip.max,
          mean.function = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$mean.function,
          dispersion.function = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$dispersion.function,
          num.bin = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$num.bin,
          binning.method = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$binning.method,
          nfeatures = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$nfeatures,
          mean.cutoff = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$mean.cutoff,
          dispersion.cutoff = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$dispersion.cutoff
        ) %>% ScaleData(
          features = orig.commands[[paste('ScaleData', assay, sep='.')]]$features,
          model.use = orig.commands[[paste('ScaleData', assay, sep='.')]]$model.use,
          do.scale = orig.commands[[paste('ScaleData', assay, sep='.')]]$do.scale,
          do.center = orig.commands[[paste('ScaleData', assay, sep='.')]]$do.center,
          scale.max = orig.commands[[paste('ScaleData', assay, sep='.')]]$scale.max,
          block.size = orig.commands[[paste('ScaleData', assay, sep='.')]]$block.size,
          min.cells.to.block = orig.commands[[paste('ScaleData', assay, sep='.')]]$min.cells.to.block
        )%>% RunPCA(
          features = orig.commands[[paste('ScaleData', assay, sep='.')]]$features,
          npcs = length(PCs),
          rev.pca =  orig.commands[[paste('RunPCA', assay, sep='.')]]$rev.pca,
          weight.by.var = orig.commands[[paste('RunPCA', assay, sep='.')]]$weight.by.var,
          reduction.name='DOUBLETFINDER_PCA',
          verbose=FALSE
        )
        
        pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
        cell.names <- rownames(seu_wdoublets@meta.data)
        nCells <- length(cell.names)
        
        rm(seu_wdoublets); gc() # Free up memory
      }
      
      ## Compute PC distance matrix
      cat("Calculating PC distance matrix...\n")
      dist.mat <- fields::rdist(pca.coord)
      
      ## Compute pANN
      cat("Computing pANN...\n")
      pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
      rownames(pANN) <- real.cells
      colnames(pANN) <- "pANN"
      cat('   nCells = ', nCells,' \n')
      k <- round(nCells * pK)
      cat('   k = ', k,' \n')
      for (i in 1:n_real.cells) {
        neighbors <- order(dist.mat[, i])
        neighbors <- neighbors[2:(k + 1)] 
        neighbor.names <- rownames(dist.mat)[neighbors]
        pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      }
      
      #Smooth pANN values, compute cutoff for doublet identification - DWM
      #TODO: find local minima in pANN values, and calculate the tiers of multiplets
      if(nExp=='auto'){
        #Smooth pANN
        
        #Find the two maxes and one local min
        
        # Set cutoff to the local min
      }
      
      
      
      #Doublet Classification
      cat("Classifying doublets...\n")
      classifications <- rep("Singlet",n_real.cells)
      classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
      
      # Add in metadata columns with classifications and pANN values
      if(is.null(classification.name)){
        seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
      }else if(is.character(classification.name)){
        seu@meta.data[, classification.name] <- classifications
      }else{
        seu@meta.data[, "DF.classifications"] <- classifications
      }
      
      if(is.null(pANN.name)){
        seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
      }else if(is.character(pANN.name)){
        seu@meta.data[, pANN.name] <- pANN[rownames(seu@meta.data), 1]
      }else{
        seu@meta.data[, "DF.pANN"] <- pANN[rownames(seu@meta.data), 1]
      }
      
      return(seu)
    }
  }


# https://github.com/mckellardw/scMuscle/blob/main/R_scripts/scMuscle_github_v1.R
# Preprocess seurat objects
seuPreProcess <- function(seu, assay='RNA', n.pcs=50, res=0.8){
  # NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)
  
  seu = NormalizeData(
    seu
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = F
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = F,
    npcs = n.pcs
  )
  
  #find pcs to use (custom function, see helper functions script for details)
  n.pcs.use = npcs(seu, reduction = pca.name, var.toal = 0.95)
  
  # FindNeighbors %>% RunUMAP, FindClusters
  seu <- FindNeighbors(
    seu,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = FALSE
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    reduction.name=umap.name,
      verbose = F
  )
  
  seu@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use
  
  seu <- FindClusters(object = seu,resolution = res, verbose = F)
  seu[[paste0('RNA_res.',res)]] <- as.numeric(seu@active.ident)
  
  return(seu)
}