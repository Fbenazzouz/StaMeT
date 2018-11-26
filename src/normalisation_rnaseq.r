# Script normalisation données de comptage RNA-Seq

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
#loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
sink(stdout(), type = "message")


pack_dispo <- row.names(installed.packages())
if(! "optparse" %in% pack_dispo) install.packages("optparse", repos="https://cloud.r-project.org/")

library(optparse)

option_list = list(
	make_option(c("-count", "--count_file"), type="character", default=NULL,
              help="Path to text file containing raw count data. See documentation for file structure details. [default: %default]"),
	make_option(c("-des", "--design"), type="character", default=NULL,
              help="Path to text file containing design information. See documentation for file structure details. [default: %default]"),
	make_option(c("-rseq_n", "--rnaseq_norm"), type="character", default="DESeq2",
              help="Normalisation method for count data, possible values: DESeq2, edgeR, VOOM. [default: %default]")
)

arg_parser = OptionParser(option_list=option_list)
les_args = parse_args(arg_parser)

counts <- read.table(les_args$count_file, header=TRUE, sep="\t", as.is=TRUE, row.names=1)
if(!is.null(les_args$design)) condition <- read.table(les_args$design, header=TRUE, as.is=TRUE) else condition <- matrix(rep(1, ncol(counts)), ncol=1)

if(!require("edgeR", quietly=TRUE, character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    library("edgeR", quietly=TRUE, character.only=TRUE)
}

if(!require("DESeq2", quietly=TRUE, character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    library("DESeq2", quietly=TRUE, character.only=TRUE)
}

library(MASS)


Normalization <- function(counts, condition, Norm=c("DESeq2", "edgeR", "VOOM")){

	switch(Norm,
	      "DESeq2"={colData <- data.frame (condition,row.names=colnames(counts) )
			dds <- DESeqDataSetFromMatrix(counts, colData, design=if(any(condition!=condition[1])) ~condition else ~1)
					cds=estimateSizeFactors(dds)
					#sizeFactors(cds)
					Data_Norm <- counts(cds, normalized=TRUE)
					Norm_DESeq2_log2 <- log2(Data_Norm+1)
					output <- Norm_DESeq2_log2},
          "edgeR"={cds <- DGEList(counts, group=condition )
				   cds <- calcNormFactors(cds, method="TMM")
				   lcpm <- cpm(cds, log=TRUE)
				   output <- lcpm},
		  "VOOM"={design=model.matrix(~0+condition)
				  voom_trans<-voom(counts, design,span = 0.5, plot = FALSE, save.plot = FALSE )
				  voom_matrix <-  voom_trans$E
				  output <- voom_matrix},
		  {message("Unrecognised normalisation method: ", Norm, ". Normalisation method must be one of 'DESeq2', 'edgeR' or 'VOOM'") ; return()})
	return(output)
}

Norm_count <- Normalization(counts, condition, Norm=les_args$rnaseq_norm)

write.table(data.frame(Gene=row.names(Norm_count), Norm_count), file="RNAseq_normalization.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
