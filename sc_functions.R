library(Seurat)
library(ggplot2)
library(ggthemes)

parallelize <- function(n.workers) {
  library(future)
  library(Seurat)
  library(doParallel)
  plan("multiprocess", workers = n.workers)
  options(future.globals.maxSize= Inf)
}

load.and.aggregate <- function(data.loc, run.ids, group.name, normalize=TRUE, n.delim = ">", n.field=2) {
    parallelize(10)
    d10x.data <- sapply(run.ids, function(i){
        print(file.path(data.loc,i,"outs/filtered_feature_bc_matrix/"))
        d10x <- Read10X(file.path(data.loc,i,"outs/filtered_feature_bc_matrix/"))
        colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep=">")
        d10x <- d10x[grep('^mt-', rownames(d10x), invert=TRUE),]
        d10x
    })

    print('loading complete')

    print('begin aggregating and splitting...')

    experiment.data <- do.call("cbind", d10x.data)


    experiment.aggregate <- CreateSeuratObject(
        experiment.data,
        min.cells = 10,
        names.delim= n.delim,
        names.field= n.field)

    experiment.aggregate <- SplitObject(experiment.aggregate, split.by='ident')

    experiment.aggregate
}
normalize.aggregated.object <- function(experiment.aggregate) {
    library(Seurat)
    library(future)
    plan("multiprocess", workers = 6)
    options(future.globals.maxSize= Inf)
    print('splitting complete')

    for (i in 1:length(x = experiment.aggregate)) {
        print(i)
        experiment.aggregate[[i]] <- NormalizeData(object = experiment.aggregate[[i]], verbose = TRUE)
        print('normalization complete')
        experiment.aggregate[[i]] <- FindVariableFeatures(object = experiment.aggregate[[i]],
            selection.method = "vst", nfeatures = 2000, verbose = TRUE)
        print('variable feature selection complete')
    }

    experiment.aggregate
}
reset.ids <- function(exp.object, ident1.list, ident2.label, ident2.newnames) {
  orig.object.idents <- Idents(exp.object)
  temp <- exp.object
  cell.list <- c()
  for (i in 1:length(ident1.list)) {
    cell.list <- WhichCells(temp, idents=ident1.list[[i]])
    temp <- SetIdent(temp, cells=cell.list, value=ident2.newnames[i])
    temp <- AddMetaData(temp, metadata=Idents(temp), col.name=ident2.label)
  }
  Idents(temp) <- orig.object.idents
  temp
}

reset.ids.cell.ids <- function(exp.object, cellvector.list, ident2.label, ident2.newnames) {
  orig.object.idents <- Idents(exp.object)
  temp <- exp.object
  cell.list <- c()
  for (i in 1:length(cellvector.list)) {
    cell.list <- WhichCells(temp, cells=cellvector.list[[i]])
    temp <- SetIdent(temp, cells=cell.list, value=ident2.newnames[i])
    temp <- AddMetaData(temp, metadata=Idents(temp), col.name=ident2.label)
  }
  Idents(temp) <- orig.object.idents
  temp
}

integrate.exp.list <- function(experiment.aggregated, nns=30,knns=50, knn.score=30, features.to.keep = NULL, future.workers=6, ref.index=NULL) {
    library(Seurat)
    library(future)
    plan("multiprocess", workers = future.workers)
    options(future.globals.maxSize= Inf)


    print('Beginning normalization...')
    experiment.list.anchors <- FindIntegrationAnchors(object.list = experiment.aggregated, dims = 1:nns,k.score=knn.score, k.filter=knns, verbose=T, reference=ref.index)

    print('anchors found')
    print(1)
    experiment.list.integrated <- IntegrateData(anchorset = experiment.list.anchors, dims = 1:nns, features.to.integrate= features.to.keep, verbose=T)
    print(2)
    print('integration complete')
    DefaultAssay(object = experiment.list.integrated) <- "integrated"
    print('script ')

    experiment.list.integrated
}

transfer.data <- function(ref.integrated, query.integrated) {
    library(Seurat)
    ref.anchors <- FindTransferAnchors(normalization.method='LogNormalize', reference = ref.integrated, query = query.integrated,
        dims = 1:30, reduction='pcaproject', project.query=TRUE, l2.norm=TRUE)
    projected.ids <- TransferData(l2.norm=TRUE, anchorset = ref.anchors, refdata = Idents(ref.integrated),
        dims = 1:30,  weight.reduction='pcaproject',k.weight=100)
    query.integrated <- AddMetaData(object = query.integrated, metadata = projected.ids)

    query.integrated
}

find.pcs.and.clustids <- function(exp.data, pc.dims=30, scale.data=TRUE, feature.list =NULL) {
    library(future)
    library(doParallel)
    plan("multiprocess", workers = 8)
    options(future.globals.maxSize= 4194304000)

    if (scale.data==TRUE) {
        exp.data <- ScaleData(object = exp.data, verbose = TRUE)
    }
    exp.data <- RunPCA(object = exp.data, npcs = pc.dims, verbose = TRUE, features=feature.list)
    exp.data <- RunUMAP(object = exp.data, reduction = "pca",
    dims = 1:pc.dims)
    exp.data <- FindNeighbors(exp.data, dims = 1:pc.dims)
    exp.data <- FindClusters(exp.data, resolution = 0.5)

    exp.data
}

euc.dist <- function(data, point) {
    apply(data, 1, function (row) sqrt(sum((point - row) ^ 2)))
}

getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

project.pca.space <- function(ctl.object, exp.object) {
    pc.matrix <- Loadings(ctl.object[['pca']])
    feature.names <- rownames(pc.matrix)
    ctl.embeddings <- Embeddings(ctl.object[['pca']])
    exp.raw.matrix <- GetAssayData(exp.object, slot='scale.data')
    print('beginning matrix multiplication')
    print(dim(t(exp.raw.matrix)))
    print(dim(pc.matrix))
    exp.projected <- t(exp.raw.matrix) %*% pc.matrix
    exp.object[['pca']] <- CreateDimReducObject(embeddings = exp.projected, key = "pca_", assay = DefaultAssay(exp.object))
    exp.object
}




# make ml density map
kernel.project <- function(kernel.ggplot, scat.project, ident.legend=TRUE, kernel.legend=FALSE) {
    p1 <- kernel.ggplot
    p2 <- scat.project
    if(length(p1$data[[1]]) > length(p2$data[[1]])) {
        zeros <- matrix(0.0,ncol=length(p2$data), nrow=(length(p1$data[[1]])-length(p2$data[[1]])))
        colnames(zeros) <- names(p2$data)
        p2$data <- rbind(p2$data, zeros)
    } else {
        zeros <- matrix(NA,ncol=length(p1$data), nrow=(length(p2$data[[1]])-length(p1$data[[1]])))
        colnames(zeros) <- names(p1$data)
        p1$data <- rbind(p1$data, zeros)
    }

    p2$data$UMAP_1_ctl <- p1$data$UMAP_1
    p2$data$UMAP_2_ctl <- p1$data$UMAP_2

    x.lim <- range(p1$data$UMAP_1, na.rm=T)
    y.lim <- range(p1$data$UMAP_2, na.rm=T)

    if(ident.legend) {
        ident.legend <- p2$guides$colour
    }



    x.lim <- x.lim + x.lim*.5
    y.lim <- y.lim + y.lim*.5


    p2 <- p2 + stat_density_2d(aes(x=UMAP_1_ctl, y=UMAP_2_ctl, fill = stat(level)),
    geom = "polygon") + lims(x = x.lim,y = y.lim) + guides(colour=p2$guides$colour)

    if(!(kernel.legend)) {
        p2 <- p2 + guides(fill=FALSE)
    }

    p2$layers <- rev(p2$layers)
    p2
}

plot.many.markers <- function(s.object, markers, file.folder, file.prefix, range=length(markers)) {
    for(i in 1:range) {
        temp.plot <- FeaturePlot(s.object, reduction='umap', features=markers[[i]], split.by='stage')
        ggsave(paste(file.folder,file.prefix, '_cluster_',i-1, '_featureplot.pdf', sep=''), plot=temp.plot)
    }
}


kevinify <- function(gg.object) {
    gg.object <- gg.object + xlab('Umap 1') + ylab('Umap 2') +
    theme_tufte() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_gradient(low='grey50', high='grey15') + scale_colour_fivethirtyeight()
}


find.top.markers <- function(exp.object, id.list=levels(exp.object$seurat_clusters), n.markers=6, future.workers=10) {
  parallelize(10)
  top.markers <- c()
  top.few.markers <- c()
  for (i in id.list) {
      id.name <- paste('cluster',i,'.markers', sep="")
      top.markers[[id.name]] <- FindMarkers(exp.object, ident.1 = i, only.pos=TRUE)
      top.few.markers[[id.name]] <- rownames(top.markers[[id.name]])[1:n.markers]
      print(i)
  }
  list(top.markers,top.few.markers)
}
find.stage.diff <- function(exp.object, stage1.list, stage2.list) {
  markerlist <- c()
  for (i in 1:length(levels(exp.object$seurat_clusters))) {
    try( {
      id.name <- paste('cluster',i-1,'.markers', sep="")
      cluster.subset <- subset(exp.object, idents=c(i-1))
      Idents(cluster.subset) <- cluster.subset$stage
      temp.markers <- FindMarkers(cluster.subset, ident.1=stage1.list, ident.2=stage2.list, assay='RNA', max.cells.per.ident=50)
      markerlist[[id.name]] <- temp.markers
      print(temp.markers)
    })

  }
  markerlist
}

write.markers.tocsv <- function(marker.object, filename) {
  markers.out <- c()
  for (i in 1:length(marker.object)) {
    temp <- marker.object[[i]]
    temp <- cbind(rownames(temp), temp)
    rownames(temp) <- NULL
    temp <- rbind(paste('cluster_',i-1,'_markers', sep=''), temp)
    if(i==1) {
      markers.out <- temp[1:50,]
    }
    else {
      markers.out <- cbind(markers.out, temp[1:50,])
    }
  }
  write.csv(markers.out, filename)
}
