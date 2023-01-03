#!/usr/bin/env Rscript

################################################
## Resume                                     ##
################################################

# Revisar los resultados de la tabla de virus IMG/VR

################################################
## LOAD LIBRARIES                             ##
################################################

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(openxlsx, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)

### Not in o in ----

`%notin%` <- Negate(`%in%`)

################################################
## ANALYSIS     ################################
################################################

mod_data<- openxlsx::read.xlsx("./data/mod_data.xlsx")

################################################
## FUNCTIONS     ###############################
################################################

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
    for (j in 1:nrow(dataframe)) {
    filter_df<- dataframe[j, c(1,2,9:10)]
    colnames(filter_df)<- c("gene_id", "gene_name", "sequence", "pfam")
    filter_df$pfam<- gsub(";", ",", filter_df$pfam)
    filter_lista<- list_with_values(filter_df)
    final_list[[j]]<- filter_lista
    final_vector[j]<- paste(unlist(final_list[[j]]), collapse = ";")
    }
    return(final_vector)
}

create_df<- function(dataframe) {
    df_gff<- data.frame(
        seqid = str_split(as.character(dataframe[,1]), "\\.", simplify = T)[,3],
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

write_gff <- function(dataframe, file){
    encabezado1<- paste0("##gff-version 3")
    encabezado2<- paste0("##sequence-region p1 1 17139")
    cat(encabezado1, "\n", encabezado2, "\n", file = file, sep = "")
    write.table(dataframe, file, append = T, col.names = F, row.names = F, quote = F, sep = "\t")
}

################################################
## USE     #####################################
################################################

df_final<- create_df(mod_data)
write_gff(df_final, "prueba.gff")

################################################
## Pruebas #####################################
################################################
#
#library(ape, quietly = TRUE, warn.conflicts = FALSE)
#
#prueba_gff<- read.gff("./data/genomic.gff", na.strings = c(".", "?"), GFF3 = TRUE); str(prueba_gff)
#prueba_gff_propio<- read.gff("./prueba.gff", na.strings = c(".", "?"), GFF3 = TRUE); str(prueba_gff_propio)
#
#prueba_gff$score
#prueba_gff_propio$score


