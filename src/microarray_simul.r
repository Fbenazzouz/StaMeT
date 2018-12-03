# Script simulateur microarray

# setup R error handling to go to stderr
options( show.error.messages=FALSE, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, FALSE ) } )
#loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
sink(stdout(), type = "message")

pack_dispo <- row.names(installed.packages())
if(! "optparse" %in% pack_dispo) install.packages("optparse", repos="https://cloud.r-project.org/")

library(optparse)

option_list = list(
    make_option(c("-gn", "--gene_number"), type="numeric", default=100, 
              help="Total number of genes that are simulated [default: %default]"),
	make_option(c("-sn1", "--samples_n1"), type="numeric", default=7, 
              help="Number of samples with phenotype 1 [default: %default]"),
	make_option(c("-sn2", "--samples_n2"), type="numeric", default=7, 
              help="Number of samples with phenotype 2 [default: %default]"),
	make_option(c("-diff", "--diff_genes_ratio"), type="numeric", default=0.1, 
              help="Proportion of differentially expressed genes (genes related to phenotype) within all genes, between 0 and 1 [default: %default]"),
	make_option(c("-up", "--up_ratio"), type="numeric", default=0.5, 
              help="Proportion of up-regulated (phenotype 2 / phenotype 1) genes within differentially-expressed genes, between 0 and 1 [default: %default]"),
	make_option(c("-m1", "--m1"), type="numeric", default=1.4, 
              help="average difference between global mean and phenotype mean for highly differentially expressed genes [default: %default]"),			  
	make_option(c("-m2", "--m2"), type="numeric", default=0.8, 
              help="average difference between global mean and phenotype mean for weakly differentially expressed genes [default: %default]"),
    make_option(c("-s", "--seed"), type="numeric", default=NULL, 
              help="Seed value: can be set to make the simulation reproducible. Non-numeric values are considered NULL. [default: %default]")		
)

arg_parser = OptionParser(option_list=option_list)
les_args = parse_args(arg_parser)

Simulation.microarray=function(nGenes, n1, n2, pi0, up, muminde1, muminde2, ratio=FALSE, seed=NULL){
    
    ## au cas où on aurait autre chose que des nombres entiers pour nGenes, n1 et n2, on les arrondit pour avoir des entiers
    nGenes <- round(nGenes) ; n1 <- round(n1) ; n2 <- round(n2)
    
    ## Vérification que tous les paramètres numériques le sont bien
    check_num <- sapply(list(nGenes, n1, n2, pi0, up, muminde1, muminde2), is.numeric)
    if(any(!check_num)) {message("All parameters must be numeric: ", paste(sapply(which(!check_num), switch, "gene number", "samples with phenotype 1",  "samples with phenotype 2", "DE proportion", "up DE proportion", "high mean difference", "low mean difference"), collapse=", "), " was(were) not."); return()}
    
    ## Vérification que tous les paramètres numériques sont positifs
    check_pos <- c(nGenes, n1, n2, pi0, up, muminde1, muminde2)>=0
    if(any(!check_pos)) {message("Parameter(s) ", paste(sapply(which(!check_pos), switch, "gene number", "samples with phenotype 1",  "samples with phenotype 2", "DE proportion", "up DE proportion", "high mean difference", "low mean difference"), collapse=", "), " should be positive."); return()}
    
    ## Vérification que l'on a au moins un gène demandé et au moins un patient
    if(any(c(nGenes, n1+n2)<=0)) {message("You need at least one gene and one patient to compute the simulation"); return()}
    
    ## Vérification que les 2 paramètres de proportions sont bien compris entre 0 et 1
    if(pi0<0 | pi0>1) {message("Proportion of DE genes should be a value between 0 ans 1."); return()}
    if(up<0 | up>1) {message("Proportion of up-regulated genes should be a value between 0 ans 1."); return()}
    
    ## Initialisation des paramètres
	shape2 = 4; lb = 4; ub = 14; lambda1 = 0.13; sdde = 0.5; sdn = 0.4; N=n1+n2+1; lambda2=4

    ## si un seed a été entré et qu'il est bien numérique, on le fixe
    if(!is.null(seed) & is.numeric(seed)) set.seed(seed)

	## Nombre de vrais positifs																			
	TP <- round(nGenes * pi0)		
	## Nombre de Faux positifs 																			
	FP <- nGenes - TP 																			
	## types des TP (up & down)																		
	TP_up <- round(TP * up)
	TP_down <- TP - TP_up 
	## TP avec des LFC élevés en valeur absolue
	TP_up_up=round((TP_up/2))
	TP_down_down=round((TP_down/2))
		
	## simulation de n valeurs selon la loi bêta, une par gène, d'où découleront ensuite les valeurs pour les différents patients.
    x <- rbeta(nGenes, 2, shape2) 
	## mise à l'échelle (de lb à (lb+ub)): 
	## On transforme ces valeurs pour correspondre aux données réelles variant entre la borne inférieure  lb=4 et la borne supérieure ub=14. x2=lb+ub*x.
	## Ces valeurs représentent le niveau d’expression moyen pour les n gènes.
    x2 <- lb + (ub * x) 
	## initialisation de la matrice avec la bonne taille (ou 1 de trop pour l'éch ref si on est pas en ratio)
    xdat <- matrix(c(rep(0, nGenes * N)), ncol = N) 
	## "id" des gènes, initialisées à 0 pour tous => pour récupérer info up/down/non de ensuite
	xid <- matrix(c( rep(1, TP_up), rep(-1, TP_down), rep(0, FP)), ncol = 1)		
	LFC <- matrix(c(rep(0, nGenes)), ncol = 1) 
	## on fait une boucle sur les gènes
	for (i in 1:nGenes) {
		alpha <- lambda1 * exp(-lambda1 * x2[i]) 
		## x2[i] contient la valeur "échantillonnée" pour le gène i; on génère aléatoirement N valeurs depuis la loi uniforme (cf algo)
		## Pour chaque valeur de x2, on génère N valeurs uniformément réparties sur un intervalle centré
		xi_val <- runif(N, min = (1 - alpha) * x2[i], max = (1 + alpha) * x2[i]) 
		## Si le gène n'est pas DE
		if (i > TP) { 
            xdat[i, ] <- xi_val 
			LFC[i,]=0
		} else {	## si le gène DE
			## on garde les (n1+1) premières valeurs, qui reflètent les valeurs dans des conditions "normales"/ première condition
			xi1 <- xi_val[1:(n1 + 1)] 
			## on tire au sort la diff de moyenne, selon muminde et lambda2
			if(( i<(TP_up_up+1)) | ((i>TP_up) & (i<TP_up+TP_down_down+1))) {
				mude <- muminde1 + rexp(1, lambda2) 
				## gène up ou down, up pour le if, down pour le else
				if (i <= TP_up_up) { 
				    fc <- rnorm(n2, mean = mude, sd = sdde)
				    LFC[i, ] <- mean(fc)
				    ## on ajoute la partie "up", distribuée selon loi N autour de mude avec toujours la même dispersion sdde
                    xi2 <- xi_val[(n1 + 2):N] + fc 
                    xid[i] <- 1 # le gène est up
				} else {
					fc <- rnorm(n2, mean = mude, sd = sdde)
					LFC[i, ] <- mean(fc)
					## idem que pour les up mais on retranche pour les down
					xi2 <- xi_val[(n1 + 2):N] - fc 
					# gène est down
					xid[i] <- -1 
				}
		
				xdat[i, ] <- c(xi1, xi2)
			} else {
			    mude <- muminde2 + rexp(1, lambda2) 
				## gène up ou down, up pour le if, down pour le else
				if ((i >TP_up_up)&(i<=TP_up)) { 
				fc <- rnorm(n2, mean = mude, sd = sdde)
				LFC[i, ] <- mean(fc)
				## on ajoute la partie "up", distribuée selon loi N autour de mude avec toujours la même dispersion sdde
                xi2 <- xi_val[(n1 + 2):N] + fc 
                xid[i] <- 1 # le gène est up
				} else {
					fc<-rnorm(n2, mean = mude, sd = sdde)
					LFC[i, ] <- mean(fc)
					## idem que pour les up mais on retranche pour les down
					xi2 <- xi_val[(n1 + 2):N] - fc 
					# gène est down
					xid[i] <- -1 
				}
    		}
	        xdat[i, ] <- c(xi1, xi2)
		}
	}
				
	LFC=LFC*xid

	## annotation des colonnes et des lignes
    colnames(xdat) <- grep("\\d$", c("V1", paste("cond1", seq_len(n1), sep="_"), paste("cond2", seq_len(n2), sep="_")), value=TRUE)
	rownames(xdat) <- grep("\\d$", c(paste("Gene.up", seq_len(TP_up), sep="_"), paste("Gene.down", seq_len(TP_down), sep="_"), paste("Gene" , seq_len(FP), sep="_")), value=TRUE)
	## sd de la 1e colonne = "individu" ref pour les ratios
    xsd <- sd(xdat[, 1]) 
	## si sd est strictement positif, on utilise cette valeur pour recréer une matrice de la même taille que xdat, de moyenne 0 et sd=xsd (pour faire le bruit additif)
    if (sdn > 0) { 
        ndat <- matrix(c(rnorm(nGenes * N, mean = 0, sd = sdn)), ncol = N)
		## On ajoute du bruit additif       
	   xdat <- xdat + ndat 
    }
	
	## si on est en simulation de ratio, on retranche la ref (correspond à un ratio car on simule données en log)
    if (ratio) { 
        xdata <- xdat[, 2:N, drop=FALSE] - xdat[, 1]
    } else { # sinon on supprime la 1e colonne
        xdata <- xdat[, 2:N, drop=FALSE]
    }
	
	list(xdata=xdata, xid=xid, xsd=xsd, delta=LFC)
}

Array_data <- Simulation.microarray(nGenes=les_args$gene_number, n1=les_args$samples_n1, n2=les_args$samples_n2, pi0=les_args$diff_genes_ratio, up=les_args$up_ratio, muminde1=les_args$m1, muminde2=les_args$m2, seed=les_args$seed)

write.table(data.frame(Gene=row.names(Array_data$xdata), Array_data$xdata), file="MicroArray_simulation.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
