
# Script simulateur microarray

library(optparse)

option_list = list(
    make_option(c("-gn", "--gene_number"), type="numeric", default=10000, 
              help="Total number of genes that are simulated [default: %default]"),
	make_option(c("-sn1", "--samples_n1"), type="numeric", default=75, 
              help="Number of samples with phenotype 1 [default: %default]"),
	make_option(c("-sn2", "--samples_n2"), type="numeric", default=75, 
              help="Number of samples with phenotype 2 [default: %default]"),
	make_option(c("-up", "--up_ratio"), type="numeric", default=0.5, 
              help="Proportion of up-regulated genes within differentially-expressed genes [default: %default]"),
	make_option(c("-diff", "--diff_genes_ratio"), type="numeric", default=0.1, 
              help="Proportion of differentially expressed genes (genes related to phenotype) within all genes [default: %default]"),
	make_option(c("-m1", "--m1"), type="numeric", default=1.4, 
              help="average difference between global mean and phenotype mean for highly differentially expressed genes [default: %default]"),			  
	make_option(c("-m2", "--m2"), type="numeric", default=0.8, 
              help="average difference between global mean and phenotype mean for weakly differentially expressed genes [default: %default]"),		
)

arg_parser = OptionParser(option_list=option_list)
les_args = parse_args(arg_parser)


Simulation.microarray=function(nGenes, n1, n2, pi0, up, muminde1, muminde2, ratio=FALSE){
## Initialisation des param�tres
	  shape2 = 4; lb = 4; ub = 14; lambda1 = 0.13; sdde = 0.5;sdn = 0.4;N=n1+n2+1;lambda2=4; 

	  ## Nombre des Faux positifs 																			
	  FP <- round(nGenes * pi0) 																			
	  ## Nombre des vrais positifs																			
	  TP <- nGenes - FP 																					
	  ## types des TP ( up & down)																		
	  TP_up <- round(TP * up)
	  TP_down <- TP - TP_up 
	  ## TP avec des LFC �lev�s en valeur absolue
	  TP_up_up=round((TP_up/2))
	  TP_down_down=round((TP_down/2))
		## simulation de n valeurs selon la loi b�ta, une par g�ne, d'o� d�couleront ensuite les valeurs pour les diff�rents patients.
        x <- rbeta(nGenes, 2, shape2) 
		## mise � l'�chelle (de lb � (lb+ub)): 
		## On transforme ces valeurs pour correspondre aux donn�es r�elles variant entre la borne inf�rieure  lb=4 et la borne sup�rieure ub=14. x2=lb+ub*x.
		## Ces valeurs repr�sentent le niveau d�expression moyen pour les n g�nes.
        x2 <- lb + (ub * x) 
		## initialisation de la matrice avec la bonne taille (ou 1 de trop pour l'�ch ref si on est pas en ratio)
    	xdat <- matrix(c(rep(0, nGenes * N)), ncol = N) 
		## "id" des g�nes, initialis�es � 0 pour tous => pour r�cup�rer info up/down/non de ensuite
		xid <- matrix(c( rep(1, TP_up), rep(-1, TP_down),rep(0, FP)), ncol = 1)		
		LFC <- matrix(c(rep(0, nGenes)), ncol = 1) 
		## on fait une boucle sur les g�nes
		for (i in 1:nGenes) {
			alpha <- lambda1 * exp(-lambda1 * x2[i]) 
			## x2[i] contient la valeur "�chantillonn�e" pour le g�ne i; on g�n�re al�atoirement N valeurs depuis la loi uniforme (cf algo)
			## Pour chaque valeur de x2, on g�n�re N valeurs uniform�ment r�parties sur un intervalle centr�
			xi_val <- runif(N, min = (1 - alpha) * x2[i], max = (1 + alpha) * x2[i]) 
			## Si le g�ne n'est pas DE
			if (i > TP) { 
            xdat[i, ] <- xi_val 
			LFC[i,]=0
			} else {	## si le g�ne DE
				## on garde les (n1+1) premi�res valeurs, qui refl�tent les valeurs dans des conditions "normales"/ premi�re condition
				xi1 <- xi_val[1:(n1 + 1)] 
				## on tire au sort la diff de moyenne, selon muminde et lambda2
				 if(( i<(TP_up_up+1)) | ((i>TP_up)&(i<TP_up+TP_down_down+1))) {

				mude <- muminde1 + rexp(1, lambda2) 
				## on tire au sort si g�ne up ou down, up pour le if, down pour le else
				if (i <= TP_up_up) { 
				fc<-rnorm(n2, mean = mude, sd = sdde)
				LFC[i,]<-mean(fc)
				## on ajoute la partie "up", distribu�e selon loi N autour de mude avec toujours la m�me dispersion sdde
                xi2 <- xi_val[(n1 + 2):N] + fc 
                xid[i] <- 1 # le g�ne est up
				} else {
					fc<-rnorm(n2, mean = mude, sd = sdde)
					LFC[i,]<-mean(fc)
					## idem que pour les up mais on retranche pour les down
					xi2 <- xi_val[(n1 + 2):N] - fc 
					# g�ne est down
					xid[i] <- -1 
				}
		
				 xdat[i, ] <- c(xi1, xi2)

			} else {
			
			mude <- muminde2 + rexp(1, lambda2) 
				## on tire au sort si g�ne up ou down, up pour le if, down pour le else
				if ((i >TP_up_up)&(i<=TP_up)) { 
				fc<-rnorm(n2, mean = mude, sd = sdde)
				LFC[i,]<-mean(fc)
				## on ajoute la partie "up", distribu�e selon loi N autour de mude avec toujours la m�me dispersion sdde
                xi2 <- xi_val[(n1 + 2):N] + fc 
                xid[i] <- 1 # le g�ne est up
				} else {
					fc<-rnorm(n2, mean = mude, sd = sdde)
					LFC[i,]<-mean(fc)
					## idem que pour les up mais on retranche pour les down
					xi2 <- xi_val[(n1 + 2):N] - fc 
					# g�ne est down
					xid[i] <- -1 
				}}
	            xdat[i, ] <- c(xi1, xi2)
				}
				}
				
	LFC=LFC*xid

	## annotation des colonnes et des lignes
    colnames(xdat) <- c("V1", paste("cond1", 1:n1, sep="_"), paste("cond2", 1:n2, sep="_")) 
	rownames(xdat) <- c(paste("Gene.up", 1:TP_up), paste("Gene.down", 1:TP_down), paste("Gene" , 1:FP))
	## sd de la 1e colonne = "individu" ref pour les ratios
    xsd <- sd(xdat[, 1]) 
	## si sd est strictement positif, on utilise cette valeur pour recr�er une matrice de la m�me taille que xdat, de moyenne 0 et sd=xsd (pour faire le bruit additif)
    if (sdn > 0) { 
        ndat <- matrix(c(rnorm(nGenes * N, mean = 0, sd = sdn)), ncol = N)
		## On ajoute du bruit additif       
	   xdat <- xdat + ndat 
    }
	
	## si on est en simulation de ratio, on retranche la ref (correspond � un ratio car on simule donn�es en log)
    if (ratio) { 
        xdata <- xdat[, 2:N] - xdat[, 1]
    } else { # sinon on supprime la 1e colonne
        xdata <- xdat[, 2:N]
    }
	
	list(xdata=xdata, xid=xid, xsd=xsd, delta=LFC)
}

Array_data <- Simulation.microarray(nGenes=les_args$gene_number, n1=les_args$samples_n1, n2=les_args$samples_n2, pi0=1-les_args$diff_genes_ratio, up=les_args$up_ratio, muminde1=les_args$m1, muminde2=les_args$m2)

write.table(data.frame(Gene=row.names(Array_data$xdata), Array_data$xdata), file="MicroArray_simulation.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)













































