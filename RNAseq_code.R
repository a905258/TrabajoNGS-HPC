
#SUPER IMPORTANTE ANTES DE EMPEZAR A CORRER TODO: HACER UNA SESION DE IJOB PARA QUE LAS LIBRERIAS SE QUEDEN CARGADAS
#EN LOS NODOS DE COMPUTO PORQUE SI NO ESTAN GUARDADAS EN OTRO LADO Y NO FUNCIONA, ASI QUE HAY QUE PONER:
#cd /scratch/arubio/bin
#./ijob
#LO QUE VIENE DESPUES ESTA BIEN PARA ASEGURARSE QUE ESTAM TODOS LOS PROGRAMAS DE NGS
#cd /scratch/arubio/bin
#. moduleloadNGS
#module load Java
library("clusterProfiler")
library("org.Hs.eg.db")
library("biomaRt")
library("sleuth")
library(readr)
library(ggplot2)
library(gplots)

samples <- read_csv("/home/osboxes/Desktop/NGSFolder/FinalWork_NGS_HPC/FinalWork_NGS_HPC/samples_finalwork.csv",col_names = FALSE) #matriz ADI todas las samples

colnames(samples) <- c("ID", "library_strategy", "base_pair_fragment", "Bioproject", "Biosample", "Accession", "X7", "library_selection", "library_source", "library_name", "end_date", "X12", "X13", "Run", "SRS", "SRA_study", "sample_name", "read_length", "center_name", "data_type", "sra", "ncbi", "library_instrument", "library_layout", "organism", "seq_type", "Published")

samples$sample_ID <- sapply(strsplit(samples$sample_name, "\\.|_"), function(X) X[1]) #selecciona la primera parte de antes de la "_"
samples$sample_type <- sapply(strsplit(samples$sample_name, "\\.|_"), function(X) X[2]) # selecciona la segunda parte de antes de la "_"

model_matrix <- model.matrix(~ sample_ID + sample_type, data = samples)


#loading data RNAseq isoform expression calculated with kallisto
#en kallisto output tenemos: samples de personas y cada una tiene un archivo abundacne.tsv en el cual aparecen los transcritos que tiene ese sample
#y cuanta cantidad de cada trnascrito hay: es una tabla tipo:
# target_id	        length	   eff_length	 est_counts	tpm
#ENST00000456328.2	1657	       1475.84	   0	       0
#en donde length es el tamaño del transcrito y est_counts te dice cuantas copias de ese trnascrito hay
#cada sample tiene todos los transcritos que se han encontrado? 

base_dir <- "/home/osboxes/Desktop/NGSFolder/FinalWork_NGS_HPC/FinalWork_NGS_HPC/samples" #base_dir = directory where the output of kallisto is stored
sample_id <- dir(base_dir) #sample_id = the names of the folders.

#if in the base_dir directory there are more files, we have to remove them from this variable. For example, if there is an R file:
#  sample_id<-sample_id[-grep(".R",sample_id,fixed=T)]

dirsample<-(paste0(base_dir,"/",sample_id[1]))
dirtoload <- paste0(dirsample,"/","abundance.tsv")
RNASeq <- read.delim(dirtoload,sep = "\t", colClasses = c(NA,"NULL","NULL","NULL",NA)) 
#te quedas con el target_id y los est_counts 
for (n in 2:length(sample_id)){
  dirsample<-(paste0(base_dir,"/",sample_id[n]))
  dirtoload <- paste0(dirsample,"/","abundance.tsv")
  RNASeq[,n+1] <- read.delim(dirtoload,sep = '\t', colClasses = c('NULL','NULL','NULL','NULL',NA))
}

rownames(RNASeq) <- sapply(strsplit(as.character(RNASeq[,1]),"\\|"),function(X) return(X[1]))
RNASeq<-RNASeq[,-1]
RNASeq2<-as.matrix(RNASeq)
colnames(RNASeq2)<-sample_id 

#RNASeq2 es una matriz que muestra en la primera columna los trnascrits ID y en la 
#primera fila los nombres de todos los samples, y te dice los est_counts de cada trnascrito en cada sample 
#para todos los samples

#AHORA SLEUTH: PARA VER DIFFERENTIAL EXPRESSION ENTRE LA MUESTRA  NORMAL Y TUMORAL DE CADA PACIENTE, 
#ASI PODEMOS VER QUE TRANSCRITOS SE EXPRESAN EN CANTIDADES DIFERENTES EN LA MUESTRA DEL TUMOR Y LA NORMAL
#PARA CADA PACIENTE

#en Atlas no esta la library de rhdf5 que se necesita para cargar sleuth, 
#entonces hemos descargado rhdf5 de https://bioconductor.org/packages/release/bioc/html/rhdf5.html

#cuando ya tenemos la library(rhdf5) ya podemos descargar sleuth con el comando -> devtools::install_github("pachterlab/sleuth")


library("sleuth")

samples_red <- samples[which(samples$Run %in% sample_id),] #buscar enla matriz de samples de adi cuales ids coinciden con los samples que se descargo laura (tenemos 32 y en ADI hay 60)
samples_red <- samples_red[match(sample_id, samples_red$Run),] #ordenar para que en el samples_red este ordenado de acuerdo a nuestros sample_id (orden ascendente)


#With the next line we obtain the complete path:
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir,id)) #añadir el path individual de cada sample

# 2.2- Build a data.frame with the condition and the kallisto directories for each sample. 
s2c <- data.frame(sample = sample_id,
                  sample_ID = samples_red$sample_ID,
                  sample_type = samples_red$sample_type) #data frame con ID, nombres de muestras:B28..., y luego cancer/normal
#check that it is ordered (en realidad el nuestro ya esta ordenado)
if (!identical(as.vector(s2c$sample),names(kal_dirs))){
  iix <- match(s2c$sample,names(kal_dirs))
  kal_dirs<-kal_dirs[iix]
}
s2c <- dplyr::select(s2c, sample = sample, sample_ID = sample_ID, sample_type = sample_type)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

#Assign the condition value to the intercept
s2c$sample_type <- relevel(as.factor(s2c$sample_type), "Normal")#convertimos a factor ambas columnas y cogemos como referencia a los samples normales
s2c_ordered <- s2c[with(s2c, order(sample_type)),]

# Load the kallisto processed data into the object
s1 <- sleuth_prep(s2c_ordered, extra_bootstrap_summary=T) 
# see ?sleuth_prep

# Estimate parameters for the sleuth response error measurement (full) model
so <- sleuth_fit(s1, ~ sample_ID + sample_type, fit_name =  'full') #matriz de diseño entera

# Estimate parameters for the sleuth reduced model
so <- sleuth_fit(so, ~ sample_ID, fit_name = 'reduced') #matriz de diseño reducida que pongo solo lo que NO quiero medir (en este caso el sample_ID)

# Perform differential analysis (testing)
so <- sleuth_lrt(so,'reduced', 'full')

# In order to see the design matrix
mod <- so$fits$full$design_matrix
mod
# 3. Read the Sleuth object ##################################################
#tpm not normalized
RNASeq_tpm_raw <- sleuth_to_matrix(so,"obs_raw","tpm") #los de nuestra matriz inicial

#counts not normalized
RNASeq_counts_raw <- sleuth_to_matrix(so,"obs_raw","est_counts")

#tpm normalized
RNASeq_tpm_norm <- sleuth_to_matrix(so,"obs_norm","tpm")

#counts normalized
RNASeq_counts_norm <- sleuth_to_matrix(so,"obs_norm","est_counts")

# Statistical Analysis result
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.001) #hemos bajado el qvalue de 0.05 a 0.001 porque con 0.05 habia 45000 trnascritos significativos, el valor de 0.001 lo hemos ido ajustando
head(sleuth_significant, 20)

sig_transcripts <- sapply(strsplit(sleuth_significant$target_id,"\\."),function(X) return(X[1]))
#el sig_transcripts es para quedarme solo con la parte antes del punto 
#del nombre del trnascritoo: tipo de ENS562343.1 quedarme con ENS562343

#SE SUPONE QUE AHROA DESPUES DEL SLEUTH TE HAS QUEDADO CON LOS TRNASCRITOS DIFERENCIALMENTE EXPRESADOS ENTRE TUMOR Y NORMAL

#NOW DOWNLOAD LIBRARY BIOMATRT:
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("biomaRt")

library("biomaRt") 

# Annotation

# Get all the transcripts, genes, go_id and domains

mart <- biomaRt::useMart(host = 'may2021.archive.ensembl.org',
                         biomart = "ensembl",
                         dataset =   "hsapiens_gene_ensembl") #sacar un dataframe de la web de ensemble de la version 2021 # hsapiens for human 
Atributos <- listAttributes(mart) #ver que caracteristicas tiene la matriz
#To obtain Affy HG U133-PLUS-2 probesets referred to YES1, MAPT and TP53 genes

#to use getBm necesitamos usar una funcion que se llama as_overscope, esta funcion esta en l libreria rlang pero solo
#en la version 0.2.2 y nosotros en atlbas tenemos la version de rlang (1.0.2), asi que hemos eliminado esa version de rlang (remove) y 
#hemos descargado la version 0.2.2 para poder usar funciones
GeneTransciptLocation <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","chromosome_name","transcript_start","transcript_end"),mart=mart)

#quedarme las caracteristocas data frame con la ubicación de los transcritos

#OJO: puede pasar que r te diga que una funcion no hace lo que tu quieres: por ejemplo la funcion
#amplify. la cosa es que puede haber en varios paquetes la misma funcion entonces tienes que especificar
#de que paquete es la funcion amplify en tu caso. por ejemplo biomatRT::amplify
#asi especificas que quieres la funcion amplify del paquete biomatRT

iij <- which(GeneTransciptLocation$ensembl_transcript_id %in% sig_transcripts)
sig_genes <- GeneTransciptLocation$ensembl_gene_id[iij]

universe <- unique(GeneTransciptLocation$ensembl_gene_id)


## GENE ONTOLOGY #PARA VER QUE FUNCIONES BIOLOGICAS SE ENCUENTRAN EN COMUN ENTRE LOS GENES QUE NOS HAN SALIDO DIFERENCIALMENTE EXPRESADOS
#ENTRE LAS MUESTRAS DE TUMOR Y NORMAL, ASI PODEMOS VER QUE FUNCION ESTA AFECTADA EN LOS PACIENTES

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#packageurl <- "https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#remove.packages("clusterProfiler")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")


library("clusterProfiler")
library("org.Hs.eg.db")


enrichGO <- clusterProfiler::enrichGO(gene = sig_genes, OrgDb = "org.Hs.eg.db",
                                      keyType ='ENSEMBL',  ont ='BP', pvalueCutoff = 0.01,
                                      pAdjustMethod = 'fdr', universe = universe,
                                      qvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 600,
                                      readable ='FALSE', pool = 'FALSE')	

database_enrichgo <- data.frame(enrichGO)


# We use simplify to remove redundant GO terms

enrichGO@result <- enrichGO@result[which(enrichGO@result$p.adjust < 0.05),]

enrichGO_simplified <- simplify(enrichGO, cutoff=0.7, by="p.adjust", select_fun=min) #te elemina los terminos que se parezcan igual o mas de un 70%

database_enrichGO_simplified <- data.frame(enrichGO_simplified)


# Dotplots

dotplot_enrichGO_simplified <- enrichplot::dotplot(enrichGO_simplified, showCategory=20) + ggtitle("Dotplot")
dotplot_enrichGO_simplified #gene ratio: cuanto % de los genes tienen esa funcion, y luego el count es el numero de genes que tienen esa funcion, y gene rtio es el count/genes universo


# Enrichment Maps

emapplot_enrichGO_simplified <- emapplot(enrichGO_simplified) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Enrichment Map")
emapplot_enrichGO_simplified 


# Other plots

cnetplot(enrichGO_simplified)


###Ver si los genes que nos salen coinciden con el paper

#BiocManager::install("clusterProfiler")
library(EnsDb.Hsapiens.v79)

sig_genes_symbol_ens <- ensembldb::select(EnsDb.Hsapiens.v79, keys = sig_genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
sig_genes_symbol <- sig_genes_symbol_ens[,1]

genes_paper_RNA <- c("STAG2","ESPL1","FGFR3","TACC3") #genes que aparecen en el abstract. STAG2 y ESPL1 están relacionados con la función de "sister chromatid cohesion and segregation"
genes_conocidos <- c("TP53", "HRAS", "FGFR3", "PIK3CA", "RB1", "KARS", "TSC1") 
chrom_rem_genes <- c("UTX", "ARID1A", "MLL-MLL3", "CREBBP-EP300", "NCOR1", "CHD6")

#PAra el variant calling debe haber mutaciones en los genes: UTX, MLL-MLL3, CREBBP-EP300, NCOR1,ARID1A, CHD6 (estos genes son los que se han descubierto)
#-Otros genes conocidos que estan mutados en bladder cancer: TP53, HRAS, FGFR3, PIK3CA, RB1, KARS, TSC1

which(genes_paper_RNA %in% sig_genes_symbol) #STAG2, ESPL1 y TACC3
which(genes_conocidos %in% sig_genes_symbol) #KRAS y HRAS
which(chrom_rem_genes %in% sig_genes_symbol) #NCOR1

## Hace falta comparar con la lista del choromosome segregation

genes_go <- strsplit(database_enrichgo_simplified[""]$geneID[j],"\\/")[[1]]


############ PCA ###########

rownames(RNASeq2) <- sapply(strsplit(rownames(RNASeq2), "\\."), function(X) X[1])

RNASeq2_filtered <- RNASeq2[sig_transcripts,]

col.order <- s2c_ordered$sample
RNASeq2_filt_ord <- RNASeq2_filtered[, col.order]

PCs <- prcomp(t(log2(1+RNASeq2_filt_ord)), center = T, retx = T)

data_PCA <- data.frame(samples = rownames(PCs$x),
                       type = s2c_ordered$sample_type,
                       PC1 = PCs$x[,1],
                       PC2 = PCs$x[,2])


library(ggplot2)

Fig_PCA <- ggplot(data_PCA) + geom_point(aes(x=PC1, y=PC2, color = type)) +
  xlab('PC1') + ylab ('PC2') + 
  theme_classic() + theme(legend.title = element_blank())
Fig_PCA


###HEAT MAP

num_cancer <- length(which(s2c_ordered$sample_type=="Cancer"))
num_normal <- length(which(s2c_ordered$sample_type=="Normal"))



heatmap_2 <- heatmap.2(x = log2(1+RNASeq2_filt_ord[1:50,]),
                       col="bluered",
                       trace="none",
                       ColSideColors = c(rep("orange", num_normal), rep("dark green", num_cancer)),
                       main="Normal vs Cancer",
                       margins = c(4,10)
)
legend(0.9,1.2, legend = c("Normal", "Cancer"),col = c("orange", "dark green"), lty= 1, lwd = 10, cex = 0.7)

##BOXPLOTS
sig_transcripts_sample <- sig_transcripts[1:10]
samples <- col.order
RNASeq2_filt_ord_sample <- log2(1+RNASeq2_filt_ord[sig_transcripts_sample,])

data_boxplot <- data.frame(transcripts = rep(sig_transcripts_sample, length(samples)), sample_id = rep(samples, each=length(sig_transcripts_sample)), sample_type = gl(2,length(sig_transcripts_sample)*num_normal,labels = c('Normal','Tumor')), values = as.vector(RNASeq2_filt_ord_sample))
data_boxplot$transcripts = as.factor(data_boxplot$transcripts)

Fig_boxplot <- ggplot(data = data_boxplot, aes(x = transcripts, y = values, fill = sample_type)) +
  geom_boxplot(width=0.5,position=position_dodge(width=0.5),outlier.shape = NA) +
  coord_cartesian(ylim =c(0,10))+ theme_classic() +
  stat_boxplot(geom="errorbar",position = position_dodge(width=0.5),width=0.2)+
  theme(plot.title = element_text(hjust = 0.5,margin=margin(0,0,10,0), size = 11),
        axis.text.x = element_text(hjust=0.9,vjust=0.9,angle=45),
        axis.text.y = element_text(margin=margin(0,0,0,10)),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()) +
  ggtitle('!0 most significant transcripts ')
Fig_boxplot
