#library("dplyr")

setwd("/media/leon/Novo volume/RABICO/GB-500-EstatisticaBioinformaticaR/trabProva/teste/dataset/")

sample.01="sample_N_Rep1.csv"
sample.02="sample_N_Rep2.csv"
sample.03="sample_DS_Rep1.csv"
sample.04="Sample_DS_Rep2.csv"

#1 ABRIR OS ARQUIVOS “csv” EM R ATRAVÉS DE UMA MATRIZ EXPERIMENTAL
#GEMEO 1
data01.01 = read.table(sample.01, header = T, sep = "\t", stringsAsFactors = F)
data02.01 = read.table(sample.02, header = T, sep = "\t", stringsAsFactors = F)

#GEMEO 2
data01.02 = read.table(sample.03, header = T, sep = "\t", stringsAsFactors = F)
data02.02 = read.table(sample.04, header = T, sep = "\t", stringsAsFactors = F)

#2 CRIAR UMA FUNÇÃO QUE FILTRE OS DADOS PELO NUMERO DE READS DE CADA VARIANTE
#para filtrar contigs
#Ncontig = 6
tCount = 1000
#para filtrar colunas
desiredCol = c("contig", "position", "variantID", "refAllele", "altAllele", "refCount", "altCount", "totalCount")

#criando subset
#GEMEO 1
sub01.01 <- subset(data01.01, totalCount > tCount ,select = c(desiredCol))
sub02.01 <- subset(data02.01, totalCount > tCount ,select = c(desiredCol))

#GEMEO 2
sub01.02 <- subset(data01.02, totalCount > tCount ,select = c(desiredCol))
sub02.02 <- subset(data02.02, totalCount > tCount ,select = c(desiredCol))

#3 CRIAR UMA FUNÇÃO DE NORMALIZAÇÃO DOS DADOS

#normalizando
#GEMEO 1
sub01.01$REF_RATIO <- sub01.01$refCount/sub01.01$totalCount
sub02.01$REF_RATIO <- sub02.01$refCount/sub02.01$totalCount

#GEMEO 2
sub01.02$REF_RATIO <- sub01.02$refCount/sub01.02$totalCount
sub02.02$REF_RATIO <- sub02.02$refCount/sub02.02$totalCount

#4  JUNTAR OS DOIS ARQUIVOS DE CADA PACIENTE UTILIZANDO AS COLUNAS COMUNS

#juntar dados
colCom01.01 = c("contig", "position", "variantID", "refAllele", "altAllele")

#GEMEO 1                 
merge1 = merge(x=sub01.01, y=sub02.01, by=colCom01.01, all=F)

#GEMEO 2
merge2 = merge(x=sub01.02, y=sub02.02, by=colCom01.01, all=F)

#renomeando colunas
#GEMEO 1
colnames(merge1)[9] <- "REF_RATIO_Rep1"
colnames(merge1)[13] <- "REF_RATIO_Rep2"

#GEMEO 2
colnames(merge2)[9] <- "REF_RATIO_Rep1"
colnames(merge2)[13] <- "REF_RATIO_Rep2"

#media entre replicas
#GEMEO 1
merge1$mean <- rowMeans(merge1[c('REF_RATIO_Rep1', 'REF_RATIO_Rep2')], na.rm=TRUE)

#GEMEO 2
merge2$mean <- rowMeans(merge2[c('REF_RATIO_Rep1', 'REF_RATIO_Rep2')], na.rm=TRUE)

#renomeando coluna da media
#GEMEO 1
colnames(merge1)[14] <- "REF_RATIO_mean"

#GEMEO 2
colnames(merge2)[14] <- "REF_RATIO_mean"

#5 CRIAR UMA FUNÇÃO DE CLASSIFICAÇÃO DAS VARIANTES DENTRO DO MESMO PACIENTES
#GEMEO 1

dim(merge1)[1]
merge1$REF_Class=""
for(i in 1:dim(merge1)[1]){
  if(merge1$REF_RATIO_mean[i] == 0 || merge1$REF_RATIO_mean[i] == 1 ){
    merge1$REF_Class[i] = "Strictly_monoallelic"
  }else{
    if(merge1$REF_RATIO_mean[i] > 0 &  merge1$REF_RATIO_mean[i] <= 0.5 || merge1$REF_RATIO_mean[i] >= 0.85 &  merge1$REF_RATIO_mean[i] < 1 ){
      merge1$REF_Class[i] = "Consistent_with_monoallelic"
    }else{
      if(merge1$REF_RATIO_mean[i] > 0.15 &  merge1$REF_RATIO_mean[i] < 0.5 || merge1$REF_RATIO_mean[i] < 0.85 &  merge1$REF_RATIO_mean[i] > 0.5 ){
        merge1$REF_Class[i] = "Consistent_with_biallelic"
      }
    }
  }
}

#GEMEO 2
dim(merge2)[1]
merge2$REF_Class=""
for(i in 1:dim(merge2)[1]){
  if(merge2$REF_RATIO_mean[i] == 0 || merge2$REF_RATIO_mean[i] == 1 ){
    merge2$REF_Class[i] = "Strictly_monoallelic"
  }else{
    if(merge2$REF_RATIO_mean[i] > 0 &  merge2$REF_RATIO_mean[i] <= 0.5 || merge2$REF_RATIO_mean[i] >= 0.85 &  merge2$REF_RATIO_mean[i] < 1 ){
      merge2$REF_Class[i] = "Consistent_with_monoallelic"
    }else{
      if(merge2$REF_RATIO_mean[i] > 0.15 &  merge2$REF_RATIO_mean[i] < 0.5 || merge2$REF_RATIO_mean[i] < 0.85 &  merge2$REF_RATIO_mean[i] > 0.5 ){
        merge2$REF_Class[i] = "Consistent_with_biallelic"
      }
    }
  }
}


#6 JUNTAR AS TABELAS DOS DOIS PACIENTES DIFERENTES
variant.merge_1_2 = merge(x=merge1, y=merge2, by=colCom01.01, all=F, suffixes = c("_N","_DS"))

#7 ENCONTRAR OS PONTOS DE MÁXIMA VARIAÇÃO ENTRE OS PACIENTES
variant.merge_1_2$N_DS_Dif=""

dim(variant.merge_1_2)[1]

for(i in 1:dim(variant.merge_1_2)[1]){
  if((variant.merge_1_2$REF_RATIO_mean_N[i] - variant.merge_1_2$REF_RATIO_mean_DS[i]) < 0 ){
    variant.merge_1_2$N_DS_Dif[i] = "DS significativo"
  }else{
    if((variant.merge_1_2$REF_RATIO_mean_N[i] - variant.merge_1_2$REF_RATIO_mean_DS[i]) > 0 ){
      variant.merge_1_2$N_DS_Dif[i] = "N significativo"
    }else{
      if((variant.merge_1_2$REF_RATIO_mean_N[i] - variant.merge_1_2$REF_RATIO_mean_DS[i]) == 0 ){
        variant.merge_1_2$N_DS_Dif[i] = "Genes Identicos"
      }
  }
  
  }
}
#write.table(variant.merge_1_2, "Leon_Result_Genes_N_DS.csv", sep = ";", row.names = F)


#para filtrar colunas
desiredCol = c("REF_RATIO_Rep1_N", "REF_RATIO_Rep2_N", "REF_RATIO_Rep1_DS", "REF_RATIO_Rep2_DS")

#criando subset
#GEMEO 1
subX <- subset(variant.merge_1_2,select = c(desiredCol))

#write.table(subX, "REF_RATIO.csv", sep = ";", row.names = F)

t.test(subX$REF_RATIO_Rep1_N, subX$REF_RATIO_Rep1_DS)

TukeyHSD()
