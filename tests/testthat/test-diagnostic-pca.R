.pca_fixture <- function(no.depth = FALSE, missing.sample = FALSE) {
  path <- tempfile(fileext=".gds"); vcf <- tempfile(fileext=".vcf")
  x <- rbind(c(0,0,1,1,2,2),c(0,1,0,2,1,2),c(2,1,2,0,1,0),
             c(NA,0,1,1,2,NA),rep(0,6),rep(NA,6),rep(1,6))
  if (missing.sample) x[,1] <- NA_real_
  writeLines(c("##fileformat=VCFv4.2","##contig=<ID=1,length=10000>",
    "##contig=<ID=2,length=10000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",paste0("s",1:6)),collapse="\t"),
    vapply(seq_len(nrow(x)),function(j) {
      gt <- c("0/0","0/1","1/1")[x[j,]+1]; gt[is.na(gt)] <- "./."
      paste(c(if(j<=3) "1" else "2",j*100,paste0("v",j),"A","G",".","PASS",".",
              if(no.depth) "GT" else "GT:DP",
              if(no.depth) gt else paste0(gt,":",10+1:6)),collapse="\t")
    },character(1))),vcf)
  SeqArray::seqVCF2GDS(vcf,path,verbose=FALSE)
  list(gds=SeqArray::seqOpen(path),files=c(path,vcf),x=t(x))
}

test_that("block PCA matches direct decomposition and preserves the selection", {
  f <- .pca_fixture(); on.exit({SeqArray::seqClose(f$gds);unlink(f$files)})
  meta <- data.frame(INDIVIDUALS=paste0("s",6:1),STRATA=rep(c("A","B"),each=3),
                     PLATE=c("x",NA,"y","y","x","x"))
  before <- SeqArray::seqGetFilter(f$gds)
  result <- diagnostic_pca(f$gds,meta,n.components=3,chunk.size=2,write.files=FALSE,verbose=FALSE)
  expect_identical(SeqArray::seqGetFilter(f$gds),before)
  expect_true(result$active.selection.restored)
  expect_identical(result$samples$INDIVIDUALS,paste0("s",1:6))
  expect_equal(result$markers$STATUS,c(rep("USED",4),"NO_VARIATION","ALL_MISSING","NO_VARIATION"))
  g <- f$x[,1:4]; p <- colMeans(g,na.rm=TRUE)/2
  z <- .diagnostic_pca_standardize(g,p); direct <- svd(z)
  expect_equal(result$components$VARIANCE_PROPORTION,direct$d[1:3]^2/sum(direct$d^2),tolerance=1e-8)
  scores <- as.matrix(result$samples[paste0("PC",1:3)])
  loadings <- as.matrix(result$markers[1:4,paste0("PC",1:3)])
  expect_equal(unname(scores %*% t(loadings)),z,tolerance=1e-8)
  expect_equal(result$samples$PCA_DEPTH,11:16)
  expect_equal(result$samples$PCA_MISSING,rowMeans(is.na(f$x)))
  expect_identical(result$depth.status,"AVAILABLE")
  expect_no_warning(lapply(result$plots,ggplot2::ggplot_build))
  # A subset remains a subset: all dimensions and metadata follow active order.
  SeqArray::seqSetFilter(f$gds,sample.id=c("s2","s3","s4","s5"),variant.id=1:4,verbose=FALSE)
  before <- SeqArray::seqGetFilter(f$gds)
  subset <- diagnostic_pca(f$gds,meta,write.files=FALSE,verbose=FALSE)
  expect_identical(SeqArray::seqGetFilter(f$gds),before)
  expect_identical(subset$samples$INDIVIDUALS,c("s2","s3","s4","s5"))
  expect_error(diagnostic_pca(f$gds,meta[1:2,],write.files=FALSE,verbose=FALSE),"cover every")
  expect_identical(SeqArray::seqGetFilter(f$gds),before)
  expect_error(diagnostic_pca(f$gds,max.cells=2,write.files=FALSE,verbose=FALSE),"exceeds")
  expect_identical(SeqArray::seqGetFilter(f$gds),before)
})

test_that("missing depth and samples without usable calls are explicit", {
  f <- .pca_fixture(no.depth=TRUE)
  on.exit({SeqArray::seqClose(f$gds);unlink(f$files)})
  result <- diagnostic_pca(f$gds,n.components=1,write.files=FALSE,verbose=FALSE)
  expect_true(all(is.na(result$samples$PCA_DEPTH)))
  expect_false("scores_1_PCA_DEPTH" %in% names(result$plots))
  expect_no_warning(lapply(result$plots,ggplot2::ggplot_build))
  expect_equal(nrow(result$components),1L)
  f2 <- .pca_fixture(missing.sample=TRUE)
  on.exit({SeqArray::seqClose(f2$gds);unlink(f2$files)},add=TRUE)
  before <- SeqArray::seqGetFilter(f2$gds)
  expect_error(diagnostic_pca(f2$gds,write.files=FALSE,verbose=FALSE),"No usable genotype calls for: s1")
  expect_identical(SeqArray::seqGetFilter(f2$gds),before)
})

test_that("unmapped marker diagnostics do not imply chromosome structure", {
  f <- .pca_fixture(); on.exit({SeqArray::seqClose(f$gds);unlink(f$files)})
  result <- diagnostic_pca(f$gds,reference.genome=FALSE,write.files=FALSE,verbose=FALSE)
  expect_equal(unique(as.character(result$plots$loadings$data$CHROM)),"Unmapped markers")
  expect_equal(unique(result$markers$CHROM),c("1","2"))
  expect_equal(result$plots$loadings$data$LOCATION,result$plots$loadings$data$ACTIVE_ORDER)
  expect_match(result$plots$loadings$labels$subtitle,"does not imply")
  # Entirely absent chromosome labels must not suppress the loading plot.
  markers <- result$markers; markers$CHROM <- NA_character_
  plots <- .diagnostic_pca_plots(result$samples,markers,result$components,"PCA_MISSING",FALSE)
  expect_true("loadings" %in% names(plots))
  expect_no_warning(ggplot2::ggplot_build(plots$loadings))
})

test_that("PCA files and figures are written without modifying source data", {
  f <- .pca_fixture(); on.exit({SeqArray::seqClose(f$gds);unlink(f$files)})
  folder <- tempfile(); dir.create(folder); on.exit(unlink(folder,recursive=TRUE),add=TRUE)
  hash <- tools::md5sum(f$files[1])
  result <- diagnostic_pca(f$gds,strata=data.frame(INDIVIDUALS=paste0("s",1:6)),
    n.components=2,path.folder=folder,verbose=FALSE,reference.genome=TRUE,
    sensitivity=list(reviewed=list(variant.id=1:3)))
  expect_true(file.exists(file.path(result$path.folder,"pca_samples.tsv")))
  expect_true(file.exists(file.path(result$path.folder,"loadings.png")))
  expect_identical(tools::md5sum(f$files[1]),hash)
  expect_true(file.exists(file.path(result$path.folder,"pca_sensitivity_summary.tsv")))
  expect_true(file.exists(file.path(result$path.folder,"pca_sensitivity_samples.tsv")))
  expect_true(file.exists(file.path(result$path.folder,"sensitivity.png")))
  expect_identical(readRDS(file.path(result$path.folder,"pca_sensitivity_selections.rds")),
    result$sensitivity$selections)
})

test_that("explicit sensitivity scenarios preserve GDS and report comparisons", {
  f <- .pca_fixture(); on.exit({SeqArray::seqClose(f$gds);unlink(f$files)})
  before <- SeqArray::seqGetFilter(f$gds)
  result <- diagnostic_pca(f$gds,n.components=2,write.files=FALSE,verbose=FALSE,
    sensitivity=list(same=list(variant.id=1:7),
      samples=list(sample.id=paste0("s",2:6)),markers=list(variant.id=1:3)))
  expect_identical(SeqArray::seqGetFilter(f$gds),before)
  expect_equal(result$sensitivity$summary$SIMILARITY[1],1,tolerance=1e-10)
  expect_equal(result$sensitivity$summary$RMS_DISPLACEMENT[1],0,tolerance=1e-10)
  expect_equal(result$sensitivity$summary$N_SAMPLES,c(6,5,6))
  expect_equal(result$sensitivity$summary$N_MARKERS,c(7,7,3))
  expect_no_warning(ggplot2::ggplot_build(result$plots$sensitivity))
  expect_error(diagnostic_pca(f$gds,write.files=FALSE,verbose=FALSE,
    sensitivity=list(no_variation=list(variant.id=c(5,7)))),"No usable genotype")
  expect_identical(SeqArray::seqGetFilter(f$gds),before)
  SeqArray::seqSetFilter(f$gds,variant.id=1:4,verbose=FALSE)
  before <- SeqArray::seqGetFilter(f$gds)
  expect_error(diagnostic_pca(f$gds,write.files=FALSE,verbose=FALSE,
    sensitivity=list(bad=list(variant.id=1:5))),"unique active")
  expect_identical(SeqArray::seqGetFilter(f$gds),before)
})

test_that("alignment ignores flips rotations scale and sample order, not changes", {
  a <- data.frame(INDIVIDUALS=letters[1:5],PC1=c(-2,-1,0,1,2),PC2=c(1,-1,0,-1,1))
  b <- a; b[c("PC1","PC2")] <- 3*as.matrix(a[c("PC1","PC2")]) %*% matrix(c(0,1,-1,0),2)
  b <- b[5:1,]
  cmp <- function(x) .diagnostic_pca_compare(a,x,c("PC1","PC2"),c("PC1","PC2"),"test")
  expect_equal(cmp(b)$summary$SIMILARITY,1,tolerance=1e-10)
  expect_equal(cmp(b)$samples$DISPLACEMENT,rep(0,5),tolerance=1e-10)
  b$PC1[1] <- b$PC1[1]+10
  expect_lt(cmp(b)$summary$SIMILARITY,1)
  expect_gt(cmp(b)$summary$MAX_DISPLACEMENT,0)
  expect_equal(cmp(b[1:2,])$summary$STATUS,"INSUFFICIENT_SHARED_SAMPLES")
  b$PC2 <- 0
  expect_equal(cmp(b)$summary$STATUS,"RANK_DEFICIENT")
})
