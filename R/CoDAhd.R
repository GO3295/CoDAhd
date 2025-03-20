#' Run various CoDA LR transformations on scRNA-seq data matrix with different count addition
#'
#' This function transforms the raw scRNA-seq matrix to different CoDA log-ratio-based matrix. To avoid zeros, different count addition scheme could be used.
#' @param dat The raw count scRNA-seq matrix (gene(row) x cell(column)).
#' @param countAdd The 'count' added to the whole matrix to avoid zero. Options included: s/gm (sum/geometric mean(sum))(default); s/max (sum/maximum(sum)); s/10000(equivalent to LogNorm-CLR); user-defined number(e.g.,1,0.1). s refers to the total count of each cell.
#' @param method The CoDA log-ratio method used. Options included: CLR(default); IQLR; LVHA; mdCLR; ILR; manual(should provide gene list as denominator); groupIQLR(should provide group information); groupLVHA(should provide group information).
#' @param gene The gene list used as denominator. Should be used with method="manual". Default is NULL.
#' @param group The group information used for groupIQLR or groupLVHA. Should be used with method="groupIQLR" or "groupLVHA". Default is NULL.
#' @return The output log-ratio matrix.
#' @examples
#' CLR_dat <- runCoDAhd_CountAdd(raw_dat,countAdd="s/gm",method="CLR");
#' IQLR1_dat <- runCoDAhd_CountAdd(raw_dat,countAdd=1,method="IQLR");
#' groupIQLR_dat <- runCoDAhd_CountAdd(raw_dat,countAdd="s/10000",method="groupIQLR",group=metadata$celltype);
#' HKGLR_dat <- runCoDAhd_CountAdd(raw_dat,countAdd="s/10000",method="manual",gene=c("GAPDH","UBC"));
#' @export
runCoDAhd_CountAdd <- function(dat,countAdd="s/gm",method="CLR",gene=NULL,group=NULL){

    rawdat <- dat[,]
    dat <- as.matrix(dat)

    if(countAdd=="s/gm"){
        gm_s <- exp(mean(log(apply(dat,2,sum))))
        dat <- apply(dat,2,function(x) x+(sum(x)/gm_s) )
        dat <- t(CLOSE(t(dat)))
    }else if(countAdd=="s/max"){
        max_s <- max(apply(dat,2,sum))
        dat <- apply(dat,2,function(x) x+(sum(x)/max_s) )
        dat <- t(CLOSE(t(dat)))
    }else if(countAdd=="s/10000"){
        dat <- apply(dat,2,function(x) x+(sum(x)/10000) )
        dat <- t(CLOSE(t(dat)))
    }else{
        countAdd <- as.numeric(countAdd)
        dat <- dat+countAdd
        dat <- t(CLOSE(t(dat)))
    }

    if(method=="CLR"){
        dat <- apply(dat,2,function(x) x/exp(mean(log(x))) )
        dat <- log2(dat)
        return(dat)
    }else if(method=="mdCLR"){
        dat <- apply(dat,2,function(x) x/exp(median(log(x))) )
        dat <- log2(dat)
        return(dat)
    }else if(method=="ILR"){
        dat <- as.matrix(t(coda.base::coordinates(t(dat))))
        dat <- dat[,colnames(rawdat)]
        return(dat)
    }else if(method=="manual"){
        if(is.null(gene)){
            stop("No gene provided")
        }else{
            S_Genes <- which(rownames(dat)%in%gene)
            dat <- apply(dat,2,function(x) x/exp(mean(log(x[S_Genes]))) )
            dat <- log2(dat)
            return(dat)
        }
    }else if(method=="IQLR"){
        tmp1 <- CreateSeuratObject(counts=rawdat,min.cells=0,min.features=0)
        tmp1 <- NormalizeData(tmp1,verbose=FALSE)
        tmp1 <- FindVariableFeatures(tmp1,selection.method="vst",nfeatures=nrow(tmp1),verbose=FALSE)

        ##IQLR
        rankedGenes <- data.frame(gene=VariableFeatures(tmp1),rank=1:length(VariableFeatures(tmp1)))
        rankedGenes$quartile <- cut2(rankedGenes$rank,g=4)
        selectedGenes <- rankedGenes[rankedGenes$quartile%in%levels(rankedGenes$quartile)[2:3],]$gene

        dat <- apply(dat,2,function(x) x/exp(mean(log(x[selectedGenes]))) )
        dat <- log2(dat)
        return(dat)
    }else if(method=="LVHA"){
        tmp1 <- CreateSeuratObject(counts=rawdat,min.cells=0,min.features=0)
        tmp1 <- NormalizeData(tmp1,verbose=FALSE)
        tmp1 <- FindVariableFeatures(tmp1,selection.method="vst",nfeatures=nrow(tmp1),verbose=FALSE)

        ##LVHA
        tmpGenes <- tmp1@assays$RNA@meta.features
        tmpGenes$quartile_m <- cut2(tmpGenes$vst.mean,g=4)
        selectedGenes_m <- rownames(tmpGenes[tmpGenes$quartile_m%in%levels(tmpGenes$quartile_m)[4],])
        tmpGenes$quartile_v <- cut2(tmpGenes$vst.variance.standardized,g=4)
        selectedGenes_v <- rownames(tmpGenes[tmpGenes$quartile_v%in%levels(tmpGenes$quartile_v)[1],])
        selectedGenes_lvha <- selectedGenes_m[selectedGenes_m%in%selectedGenes_v]

        dat <- apply(dat,2,function(x) x/exp(mean(log(x[selectedGenes_lvha]))) )
        dat <- log2(dat)
        return(dat)
    }else if(method=="groupIQLR"){
        if(is.null(group)){
            stop("No group information")
        }else{
            selectedGenes1 <- getGroupIQLRGenes(rawdat,group)
            dat <- apply(dat,2,function(x) x/exp(mean(log(x[selectedGenes1]))) )
            dat <- log2(dat)
            return(dat)
        }
    }else if(method=="groupLVHA"){
        if(is.null(group)){
            stop("No group information")
        }else{
            selectedGenes2 <- getGroupLVHAGenes(rawdat,group)
            dat <- apply(dat,2,function(x) x/exp(mean(log(x[selectedGenes2]))) )
            dat <- log2(dat)
            return(dat)
        }
    }else{
        stop("No method available")
    }
}

#' Run various CoDA LR transformations on scRNA-seq with prior log-normalization
#'
#' This function transforms the prior-log-normalized scRNA-seq matrix to different CoDA log-ratio-based matrix as an approximation.
#' @param dat The raw or log-normalized count scRNA-seq matrix (gene(row) x cell(column)).
#' @param method The CoDA log-ratio method used. Options included: CLR(default); IQLR; LVHA; mdCLR; manual(should provide gene list as denominator); groupIQLR(should provide group information); groupLVHA(should provide group information).
#' @param gene The gene list used as denominator. Should be used with method="manual".
#' @param group The group information used for groupIQLR or groupLVHA. Should be used with method="groupIQLR" or "groupLVHA".
#' @param log_normalized Whether the input data has been log-normalized. If not, the raw data will be log2-normalized followed by LR transformations. Default is True.
#' @return The output log-ratio matrix.
#' @examples
#' CLR_dat <- runCoDAhd_LogNorm(raw_dat,method="CLR",log_normalized=FALSE);
#'
#' norm_dat <- t(CLOSE(t(raw_dat))*10000);
#' norm_dat <- log2(norm_dat+1);
#' IQLR_dat <- runCoDAhd_LogNorm(norm_dat,method="IQLR",log_normalized=TRUE);
#'
#' norm_dat <- CreateSeuratObject(counts=raw_dat,min.cells=0,min.features=0);
#' norm_dat <- NormalizeData(norm_dat,verbose=FALSE); #natural log
#' norm_dat <- as.matrix(norm_dat@assays$RNA@data);
#' groupIQLR_dat <- runCoDAhd_LogNorm(norm_dat,method="groupIQLR",group=metadata$celltype,log_normalized=TRUE)
#'
#' HKGLR_dat <- runCoDAhd_LogNorm(raw_dat,method="manual",gene=c("GAPDH","UBC"),log_normalized=FALSE);
#' @export
runCoDAhd_LogNorm <- function(dat,method="CLR",gene=NULL,group=NULL,log_normalized=TRUE){

    rawdat <- dat[,]
    dat <- as.matrix(dat)

    if(!log_normalized){
        dat <- t(CLOSE(t(dat))*10000)
        dat <- log2(dat+1)
    }

    if(method=="CLR"){
        dat <- apply(dat,2,function(x) x-mean(x))
        return(dat)
    }else if(method=="mdCLR"){
        dat <- apply(dat,2,function(x) x-median(x))
        return(dat)
    }else if(method=="manual"){
        if(is.null(gene)){
            stop("No gene provided")
        }else{
            S_Genes <- which(rownames(dat)%in%gene)
            dat <- apply(dat,2,function(x) x-mean(x[S_Genes]))
            return(dat)
        }
    }else if(method=="IQLR"){
        tmp1 <- CreateSeuratObject(counts=rawdat,min.cells=0,min.features=0)
        tmp1 <- NormalizeData(tmp1,verbose=FALSE)
        tmp1 <- FindVariableFeatures(tmp1,selection.method="vst",nfeatures=nrow(tmp1),verbose=FALSE)

        ##IQLR
        rankedGenes <- data.frame(gene=VariableFeatures(tmp1),rank=1:length(VariableFeatures(tmp1)))
        rankedGenes$quartile <- cut2(rankedGenes$rank,g=4)
        selectedGenes <- rankedGenes[rankedGenes$quartile%in%levels(rankedGenes$quartile)[2:3],]$gene

        dat <- apply(dat,2,function(x) x-mean(x[selectedGenes]))
        return(dat)
    }else if(method=="LVHA"){
        tmp1 <- CreateSeuratObject(counts=rawdat,min.cells=0,min.features=0)
        tmp1 <- NormalizeData(tmp1,verbose=FALSE)
        tmp1 <- FindVariableFeatures(tmp1,selection.method="vst",nfeatures=nrow(tmp1),verbose=FALSE)

        ##LVHA
        tmpGenes <- tmp1@assays$RNA@meta.features
        tmpGenes$quartile_m <- cut2(tmpGenes$vst.mean,g=4)
        selectedGenes_m <- rownames(tmpGenes[tmpGenes$quartile_m%in%levels(tmpGenes$quartile_m)[4],])
        tmpGenes$quartile_v <- cut2(tmpGenes$vst.variance.standardized,g=4)
        selectedGenes_v <- rownames(tmpGenes[tmpGenes$quartile_v%in%levels(tmpGenes$quartile_v)[1],])
        selectedGenes_lvha <- selectedGenes_m[selectedGenes_m%in%selectedGenes_v]

        dat <- apply(dat,2,function(x) x-mean(x[selectedGenes_lvha]))
        return(dat)
    }else if(method=="groupIQLR"){
        if(is.null(group)){
            stop("No group information")
        }else{
            selectedGenes1 <- getGroupIQLRGenes(rawdat,group)
            dat <- apply(dat,2,function(x) x-mean(x[selectedGenes1]))
            return(dat)
        }
    }else if(method=="groupLVHA"){
        if(is.null(group)){
            stop("No group information")
        }else{
            selectedGenes2 <- getGroupLVHAGenes(rawdat,group)
            dat <- apply(dat,2,function(x) x-mean(x[selectedGenes2]))
            return(dat)
        }
    }else{
        stop("No method available")
    }
}

#' Select groupIQLR genes as denominators for groupIQLR transformations
#'
#' This function selects group-based IQLR(inter-quartile log-ratio) genes as denominators for groupIQLR transformations.
#' @param dat The raw count scRNA-seq matrix (gene(row) x cell(column)).
#' @param group The group information used for groupIQLR.
#' @return The selected genes.
#' @examples
#' selected_genes <- getGroupIQLRGenes(raw_dat,group=metadata$celltype);
#' @export
getGroupIQLRGenes <- function(x,group){
    uniq_group <- unique(group)

    selected1 <- which(group==uniq_group[1])
    tmp <- CreateSeuratObject(counts=x[,selected1],meta.data=data.frame(groups=group[selected1],row.names=colnames(x[,selected1])),min.cells=0,min.features=0)
    tmp <- NormalizeData(tmp,verbose=FALSE)
    tmp <- FindVariableFeatures(tmp,selection.method="vst",nfeatures=nrow(tmp),verbose=FALSE)

    ##IQLR
    rankedGenes <- data.frame(gene=VariableFeatures(tmp),rank=1:length(VariableFeatures(tmp)))
    rankedGenes$quartile <- cut2(rankedGenes$rank,g=4)
    AllselectedGenes <- rankedGenes[rankedGenes$quartile%in%levels(rankedGenes$quartile)[2:3],]$gene

    for(g in uniq_group[-1]){
        selectedx <- which(group==g)
        tmp <- CreateSeuratObject(counts=x[,selectedx],meta.data=data.frame(groups=group[selectedx],row.names=colnames(x[,selectedx])),min.cells=0,min.features=0)
        tmp <- NormalizeData(tmp,verbose=FALSE)
        tmp <- FindVariableFeatures(tmp,selection.method="vst",nfeatures=nrow(tmp),verbose=FALSE)

        rankedGenes <- data.frame(gene=VariableFeatures(tmp),rank=1:length(VariableFeatures(tmp)))
        rankedGenes$quartile <- cut2(rankedGenes$rank,g=4)
        selectedGenes <- rankedGenes[rankedGenes$quartile%in%levels(rankedGenes$quartile)[2:3],]$gene
        AllselectedGenes <- AllselectedGenes[AllselectedGenes%in%selectedGenes]
    }
    return(AllselectedGenes)
}

#' Select groupLVHA genes as denominators for groupLVHA transformations
#'
#' This function selects group-based LVHA(low-variance-high-abundance) genes as denominators for LVHA transformations.
#' @param dat The raw count scRNA-seq matrix (gene(row) x cell(column)).
#' @param group The group information used for groupLVHA.
#' @return The selected genes.
#' @examples
#' selected_genes <- getGroupLVHAGenes(raw_dat,group=metadata$celltype);
#' @export
getGroupLVHAGenes <- function(x,group){
    uniq_group <- unique(group)

    selected1 <- which(group==uniq_group[1])
    tmp <- CreateSeuratObject(counts=x[,selected1],meta.data=data.frame(groups=group[selected1],row.names=colnames(x[,selected1])),min.cells=0,min.features=0)
    tmp <- NormalizeData(tmp,verbose=FALSE)
    tmp <- FindVariableFeatures(tmp,selection.method="vst",nfeatures=nrow(tmp),verbose=FALSE)

    ##LVHA
    tmpGenes <- tmp@assays$RNA@meta.features
    tmpGenes$quartile_m <- cut2(tmpGenes$vst.mean,g=4)
    selectedGenes_m <- rownames(tmpGenes[tmpGenes$quartile_m%in%levels(tmpGenes$quartile_m)[4],])
    tmpGenes$quartile_v <- cut2(tmpGenes$vst.variance.standardized,g=4)
    selectedGenes_v <- rownames(tmpGenes[tmpGenes$quartile_v%in%levels(tmpGenes$quartile_v)[1],])
    AllselectedGenes <- selectedGenes_m[selectedGenes_m%in%selectedGenes_v]

    for(g in uniq_group[-1]){
        selectedx <- which(group==g)
        tmp <- CreateSeuratObject(counts=x[,selectedx],meta.data=data.frame(groups=group[selectedx],row.names=colnames(x[,selectedx])),min.cells=0,min.features=0)
        tmp <- NormalizeData(tmp,verbose=FALSE)
        tmp <- FindVariableFeatures(tmp,selection.method="vst",nfeatures=nrow(tmp),verbose=FALSE)

        tmpGenes <- tmp@assays$RNA@meta.features
        tmpGenes$quartile_m <- cut2(tmpGenes$vst.mean,g=4)
        selectedGenes_m <- rownames(tmpGenes[tmpGenes$quartile_m%in%levels(tmpGenes$quartile_m)[4],])
        tmpGenes$quartile_v <- cut2(tmpGenes$vst.variance.standardized,g=4)
        selectedGenes_v <- rownames(tmpGenes[tmpGenes$quartile_v%in%levels(tmpGenes$quartile_v)[1],])
        selectedGenes <- selectedGenes_m[selectedGenes_m%in%selectedGenes_v]
        AllselectedGenes <- AllselectedGenes[AllselectedGenes%in%selectedGenes]
    }
    return(AllselectedGenes)
}

