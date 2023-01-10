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
mod_data<- read.csv2("./data/prueba_data.csv", header = T, sep = ",")
df_fasta<- read.table("./data/length_fasta.tsv", header = T); colnames(df_fasta)<- c("Scaffold.name", "length_contig")
data_conjunto<- merge(mod_data, df_fasta, by = c("Scaffold.name"))

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
    filter_df$pfam<- gsub(";", ",", filter_df$pfam)
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
        start = dataframe[,3],
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

create_gff<- function(dataframe, file) {
    lista<- list()
    for (i in 1:length(unique(dataframe$seqid))){
        partes<- unique(dataframe$seqid)
        df_filter<- subset(dataframe, seqid == partes[i])
        nombres<- partes[i]
        minvalue<- min(df_filter[,4])
        maxvalue<- max(df_filter[,5])
        lista[[i]]<- paste0("##sequence-region", " ", nombres, " ", minvalue, " ", maxvalue, "\n")
    }
    encabezado1<- paste0("##gff-version 3")
    cat(encabezado1, "\n", file = file, append = T, sep = "")
    cat(unlist(lista), file = file, append = T, sep = "")
    write.table(dataframe, file = file, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
}

################################################
## USE     #####################################
################################################

df_final<- create_df(data_conjunto)
View(head(df_final))
data_prueba<- head(df_final, 300)
create_gff(data_prueba, "output_prueba.gff")
create_gff(df_final, "judiseq.gff")
#write.table(unique(data_prueba$seqid), file = "index_fasta_prueba.txt", col.names = F, row.names = F, quote = F, sep = "\t")
