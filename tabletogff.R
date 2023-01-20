#!/usr/bin/env Rscript

################################################
## Resume                                     ##
################################################

# Crear un gff a partir de una tabla de anotaciones

################################################
## LOAD LIBRARIES                             ##
################################################

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
#library(openxlsx, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)

### Not in o in ----

`%notin%` <- Negate(`%in%`)

################################################
## ANALYSIS     ################################
################################################

#mod_data<- openxlsx::read.xlsx("./data/total_data.xlsx", 1, startRow = 2)
#metadata<- read.table("/srv/www/NOBACKUP/JudiSeq/JudiSeq-Genes-Scaffolds-relation.tsv", header = F); colnames(metadata)<- c("Contig", "Scaffold.name")
#data_conjunto<- merge(mod_data, metadata, by = c("Scaffold.name"))
mod_data<- read.csv2("./data/JudiSeq-A25_v_JoinedAnnotations.tsv", header = T, sep = "\t")
names(mod_data)[names(mod_data) == 'X.ID'] <- 'ID'
contigs<- c("JudiSeq-A25_scaf_1", "JudiSeq-A25_scaf_2", "JudiSeq-A25_scaf_3")
mod_data$Scaffold.name[mod_data$Scaffold.name %in% contigs[1]]<- "JudiSeq-A25_scaf_1|200333bp|contig_3368:F:-126:contig_1760:F,contig_2176"
mod_data$Scaffold.name[mod_data$Scaffold.name %in% contigs[2]]<- "JudiSeq-A25_scaf_2|1710977bp|contig_1165:F:1448:contig_3339:R:2:contig_295:F:2930:contig_2332:F,contig_1864"
mod_data$Scaffold.name[mod_data$Scaffold.name %in% contigs[3]]<- "JudiSeq-A25_scaf_3|1934252bp|contig_735:F:2148:contig_109:F,contig_3194:R:791:contig_1726:R:153:contig_3106:R:98:contig_3476:F"
metadata<- read.table("./data/metadata.tsv", sep = "\t", header = T)
data_conjunto<- merge(mod_data[,!names(mod_data) %in% c("X")], metadata, by = c("Scaffold.name"))

################################################
## FUNCTIONS     ###############################
################################################

create_ind<- function(dataframe) {
    ind<- data.frame(
        variables = colnames(dataframe),
        n = seq(1:length(colnames(dataframe)))
    )
    return(ind)
}

list_with_values<- function(dataframe) {
    lista<- list()
    for (i in 1:length(dataframe)){
        columnas<- colnames(dataframe)[i]
        values<- dataframe[1,i]
        lista[[i]]<- paste0(columnas,"=", values)
        }
    return(lista)
}

get_atributtes<- function(dataframe) {
    final_list<- list()
    final_vector<- 0
    for (j in 1:nrow(dataframe)) {
    filter_df<- dataframe[j, c(2,9:38)]
    colnames(filter_df)<- c("gene_id", "sequence", tolower(colnames(filter_df[,3:31])))
    #filter_df$pfam<- gsub(";", ",", filter_df$pfam)
    #filter_df$interpro<- gsub(";", ",", filter_df$interpro)
    columnas_cambio <- colnames(filter_df[,3:31])
    filter_df[columnas_cambio] <- lapply(filter_df[columnas_cambio], gsub, pattern = ";", replacement = ",")
    filter_df$gene_name<- str_split(as.character(filter_df[,1]), ".p", simplify = T)[,1]
    filter_df<- filter_df[,c(1,32,2:31)]
    filter_lista<- list_with_values(filter_df)
    final_list[[j]]<- filter_lista
    final_vector[j]<- paste(unlist(final_list[[j]]), collapse = ";")
    }
    return(final_vector)
}

create_df<- function(dataframe) {
    df_gff<- data.frame(
        seqid = dataframe[,1],
        source = rep("AG", nrow(dataframe)),
        type = str_split(as.character(dataframe[,2]), "\\.", simplify = T)[,2],
        #start = dataframe[,3] + (dataframe[,6]-1),
        start = dataframe[,3],
        #end = dataframe[,4] + (dataframe[,7]-1),
        end = dataframe[,4],
        score = rep(".", nrow(dataframe)),
        #orf_length = dataframe[,5],
        #start_aac = dataframe[,6],
        #end_aac = dataframe[,7],
        strand = dataframe[,8],
        phase = rep(0, nrow(dataframe)),
        attributes = get_atributtes(dataframe)
    )
    return(df_gff)
}

create_gff<- function(dataframe1, dataframe2, file) {
    lista<- list()
    partes<- unique(dataframe2$Scaffold.name)
    for (i in 1:length(unique(dataframe2$Scaffold.name))){
        df_filter<- subset(dataframe2, Scaffold.name == partes[i])
        nombres<- partes[i]
        minvalue<- 1
        maxvalue<- df_filter[,2]
        lista[[i]]<- paste0("##sequence-region", " ", nombres, " ", minvalue, " ", maxvalue, "\n")
    }
    encabezado1<- paste0("##gff-version 3")
    cat(encabezado1, "\n", file = file, append = T, sep = "")
    cat(unlist(lista), file = file, append = T, sep = "")
    write.table(dataframe1, file = file, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
}

################################################
## USE     #####################################
################################################

df_final<- create_df(data_conjunto)
data_prueba<- head(df_final, 300)
create_gff(data_prueba, df_fasta, "output_prueba.gff")
create_gff(df_final, df_fasta, "judiseq.gff")
#write.table(unique(data_prueba$seqid), file = "index_fasta_prueba.txt", col.names = F, row.names = F, quote = F, sep = "\t")
