# Script simulateur RNA-Seq

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
#loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
sink(stdout(), type = "message")


pack_dispo <- row.names(installed.packages())
if(! "optparse" %in% pack_dispo) install.packages("optparse", repos="https://cloud.r-project.org/")

library(optparse)

option_list = list(
    make_option(c("-gn", "--gene_number"), type="numeric", default=10000,
              help="Total number of genes that are simulated [default: %default]"),
	make_option(c("-sn1", "--samples_n1"), type="numeric", default=75,
              help="Number of samples with phenotype 1 [default: %default]"),
	make_option(c("-sn2", "--samples_n2"), type="numeric", default=75,
              help="Number of samples with phenotype 2 [default: %default]"),
	make_option(c("-diff", "--diff_genes_ratio"), type="numeric", default=0.1,
              help="Proportion of differentially expressed genes (genes related to phenotype) within all genes [default: %default]"),
	make_option(c("-up", "--up_ratio"), type="numeric", default=0.5,
              help="Proportion of up-regulated genes within differentially-expressed genes [default: %default]"),
	make_option(c("-fc", "--fc_file"), type="character", default=NULL,
              help="Path to text file with Fold-Change values for the simulated genes. See documentation for file structure details. [default: %default]"),
	make_option(c("-rseq_n", "--rnaseq_norm"), type="character", default="DESeq2",
              help="Normalisation method for count data, possible values: DESeq2, edgeR, VOOM. [default: %default]"),
    make_option(c("-s", "--seed"), type="numeric", default=NULL,
              help="Seed value: can be set to make the simulation reproducible. [default: %default]")
)

arg_parser = OptionParser(option_list=option_list)
les_args = parse_args(arg_parser)

# création du vecteur des foldchange :
# - si aucun n'est fourni, récupération des fold-change par défaut, soit tels quels, soit avec échantillonnage pour avoir les bons nombres de DE et de up
# - si le fichier est fourni, import du fichier
if (is.null(les_args$fc_file)) {
    FC=NULL	    
    }
} else {
    FC <- read.table(les_args$fc_file, as.is=TRUE, header=TRUE, sep="\t")[, 1, drop=TRUE]
}



if(!require("edgeR", quietly=TRUE, character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("edgeR")
    library("edgeR", quietly=TRUE, character.only=TRUE)
}

if(!require("MASS", quietly=TRUE, character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("MASS")
    library("MASS", quietly=TRUE, character.only=TRUE)
}
#library(MASS)


counts.simulation <- function(nGenes, n1, n2, pi0, up, fc, seed=NULL){

    ## au cas où on aurait autre chose que des nombres entiers pour nGenes, n1 et n2, on les arrondit pour avoir des entiers
    nGenes <- round(nGenes) ; n1 <- round(n1) ; n2 <- round(n2)

    ## Vérification que tous les paramètres numériques le sont bien
    check_num <- sapply(list(nGenes, n1, n2, pi0, up), is.numeric)
    if(any(!check_num)) {message("All parameters must be numeric: ", paste(sapply(which(!check_num), switch, "gene number", "samples with phenotype 1",  "samples with phenotype 2", "DE proportion", "up DE proportion"), collapse=", "), " was(were) not."); return()}

    ## Vérification que tous les paramètres numériques sont positifs
    check_pos <- c(nGenes, n1, n2, pi0, up)>=0
    if(any(!check_pos)) {message("Parameter(s) ", paste(sapply(which(!check_pos), switch, "gene number", "samples with phenotype 1",  "samples with phenotype 2", "DE proportion", "up DE proportion"), collapse=", "), " should be positive."); return()}

    ## Vérification que l'on a au moins un gène demandé et au moins un patient
    if(any(c(nGenes, n1+n2)<=0)) {message("You need at least one gene and one patient to compute the simulation"); return()}

    ## Vérification que les 2 paramètres de proportions sont bien compris entre 0 et 1
    if(pi0<0 | pi0>1) {message("Proportion of DE genes should be a value between 0 ans 1."); return()}
    if(up<0 | up>1) {message("Proportion of up-regulated genes should be a value between 0 ans 1."); return()}


	## si un seed a été entré, on le fixe
    if(!is.null(seed) & is.numeric(seed)) set.seed(seed)

	## moyenne et dispersion
	mu <- runif(nGenes, 500, 3000); disp=rlnorm(nGenes, -1.13, 1)
	replace <- TRUE
	## Nombre des vrais positifs
	TP <- round(nGenes * pi0)
	## Nombre des Faux positifs
	FP <- nGenes - TP
	## types des TP (up & down)
	TP_up <- round(TP * up)
	TP_down <- TP - TP_up

	## vecteur des DE: 1 pour les gènes Up, -1 pour les gènes down et 0 pour les gènes non DE
	DE <- c(rep(1, TP_up), rep(-1, TP_down), rep(0, FP))

    ## LFC, environ moitié positif, moitié négatif (pour les gènes différentiellement exprimés (DE))
	delta <- rep(0, nGenes)
		if (!is.null(fc)){
	# Verification qu'il n'y a pas de NA et que l'on a des valeurs numériques uniquement
	# sinon, on sort de la fonction avec message d'erreur
	if(any(is.na(fc))) {message("NA are not allowed in fc_file.") ; return()}
	if(!is.numeric(fc)) {fc <- as.numeric(fc) ; if(any(is.na(fc))) {message("Fold change values should be numeric"); return()}}
	if((length(fc)==TP) & (TP_up==0 | all(fc[1:TP_up]>1)) & all(fc[(TP_up+1):TP]<1)){
     	lfc <- log(fc)
		delta[DE != 0] <- lfc[DE != 0]
	}else{
		fc=exp(c(sort(rnorm(TP_up,2*log(2),0.3*log(2)),decreasing = TRUE)[1:(TP_up/2)],sort(rnorm(TP_up,log(2),0.5*log(2)),decreasing = TRUE)[1:(TP_up/2)],
		sort(rnorm(TP_down,-2*log(2),0.3*log(2)),decreasing = FALSE)[1:(TP_down/2)],sort(rnorm(TP_down,-log(2),0.5*log(2)),decreasing = FALSE)[1:(TP_down/2)]))
			lfc <- log(fc)
		delta[DE != 0] <- lfc[DE != 0]
	}}
	if(is.null(fc)){
		fc=exp(c(sort(rnorm(TP_up,2*log(2),0.3*log(2)),decreasing = TRUE)[1:(TP_up/2)],sort(rnorm(TP_up,log(2),0.5*log(2)),decreasing = TRUE)[1:(TP_up/2)],
		sort(rnorm(TP_down,-2*log(2),0.3*log(2)),decreasing = FALSE)[1:(TP_down/2)],sort(rnorm(TP_down,-log(2),0.5*log(2)),decreasing = FALSE)[1:(TP_down/2)]))
			lfc <- log(fc)
		delta[DE != 0] <- lfc[DE != 0]
		}
	

	 # initialisation
	 h <- rep(TRUE, nGenes)
	 counts <- matrix(0, nrow = nGenes, ncol = n1+n2)
	 selected_genes <- true_means <- true_disps <- rep(0, nGenes)
	 left_genes <- 1:length(mu)
	 lambda <- phi <- matrix(0, nrow=nGenes, ncol=n1+n2)

	while(any(h)){
		temp <- sample.int(length(left_genes), sum(h), replace)
		temp <- temp[order(temp)]
		selected_genes[h] <- left_genes[temp]
		if (!replace){
		  left_genes <- left_genes[-temp]
		}

		true_means[h] <- mu[selected_genes[h]]
		true_disps[h] <- disp[selected_genes[h]]

		lambda[h, ] <- matrix(true_means[h], ncol = 1) %*%
					   matrix(rep(1, n1+n2), nrow = 1) *
					   cbind(matrix(rep(1, sum(h) * n1), ncol = n1, nrow=sum(h)),
							 matrix(rep(exp(delta[h]), n2), ncol = n2, nrow=sum(h)))

		## moyenne des comptages
		phi[h, ] <- matrix(rep(true_disps[h],n1+n2 ), ncol = n1+n2)
		## dispersion des comptages
		counts[h, ] <- rnegbin(sum(h) * (n1+n2), lambda[h, ], 1/phi[h, ])
		h <- (rowSums(cpm(counts) > 2) < 3)
	}

	rownames(counts) <- grep("\\d$", c(paste("Gene.up", seq_len(TP_up), sep="_"), paste("Gene.down", seq_len(TP_down), sep="_"), paste("Gene" , seq_len(FP), sep="_")), value=TRUE)
	delta <- delta/log(2)

	colnames(counts) <- grep("\\d$", c(paste("cond1", seq_len(n1), sep="_"), paste("cond2", seq_len(n2), sep="_")), value=TRUE)
	list(counts = counts, DE=DE)
}



RNAseq_counts <- counts.simulation(nGenes=les_args$gene_number, n1=les_args$samples_n1, n2=les_args$samples_n2, pi0=les_args$diff_genes_ratio, up=les_args$up_ratio,
                                   fc=les_args$fc_file, seed=les_args$seed)
counts=RNAseq_counts$counts

# NORM
if(!require("DESeq2", quietly=TRUE, character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
    library("DESeq2", quietly=TRUE, character.only=TRUE)
}


Normalization <- function(counts, n1, n2, Norm=c("DESeq2", "edgeR", "VOOM")){

    ## Si on a qu'une condition (n1 ou n2 nul), on met un design: ~1 sinon on le crée suivant les conditions
    if(n1 & n2){
        exist.cond <- TRUE
    	condition <- factor(rep(paste0("Cond", 1:2), c(n1, n2)))
	} else {
	    exist.cond <- FALSE
	    condition <- factor(rep(1, n1+n2))
	}

	# normalisation
	switch(Norm,
	      "DESeq2"={colData <- data.frame(condition, row.names=colnames(counts))
					dds <- DESeqDataSetFromMatrix(counts, colData, design=if(exist.cond) ~condition else ~1)
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

Norm_count <- Normalization(counts, n1=les_args$samples_n1, n2=les_args$samples_n2, Norm=les_args$rnaseq_norm)

write.table(data.frame(Gene=row.names(Norm_count), Norm_count), file="RNAseq_simulation.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
