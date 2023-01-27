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
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(optparse, quietly = TRUE, warn.conflicts = FALSE)

### Not in o in ----

`%notin%` <- Negate(`%in%`)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list = list(
    make_option(c("-i", "--input_table"), action="store", default=NA, type='character',
              help="Path to input table with annotations"),
    make_option(c("-m", "--path_to_metadata"), action="store", default=NA, type='character',
              help="Path to input table with metadata")
              )

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$input_table)){
  print_help(opt_parser)
  stop("Please provide annotation table.", call.=FALSE)
}
if (is.null(opt$input_table)) {
  print_help(opt_parser)
  stop("Please provide a metadata table consists in: scaffold name, new id and scaffold length", call.=FALSE)
}

################################################
## ANALYSIS     ################################
################################################

mod_data<- read.csv2(opt$input_table, header = T, sep = "\t")
mod_data[mod_data == ""]<- "-"
names(mod_data)[names(mod_data) == 'X.ID'] <- 'ID'
contigs<- c("JudiSeq-A25_scaf_1", "JudiSeq-A25_scaf_2", "JudiSeq-A25_scaf_3")
mod_data$Scaffold.name[mod_data$Scaffold.name %in% contigs[1]]<- "JudiSeq-A25_scaf_1|200333bp|contig_3368:F:-126:contig_1760:F,contig_2176"
mod_data$Scaffold.name[mod_data$Scaffold.name %in% contigs[2]]<- "JudiSeq-A25_scaf_2|1710977bp|contig_1165:F:1448:contig_3339:R:2:contig_295:F:2930:contig_2332:F,contig_1864"
mod_data$Scaffold.name[mod_data$Scaffold.name %in% contigs[3]]<- "JudiSeq-A25_scaf_3|1934252bp|contig_735:F:2148:contig_109:F,contig_3194:R:791:contig_1726:R:153:contig_3106:R:98:contig_3476:F"
metadata<- read.table(opt$path_to_metadata, sep = "\t", header = T)
metadata$final_id<- paste0("judiseq_", metadata$new_id)
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

get_attributes<- function(dataframe) {
    final_list<- list()
    final_vector<- 0
    for (j in 1:nrow(dataframe)) {
        filter_df<- dataframe[j, c("ID", "Scaffold.length", "Pfam", "InterPro", "GO.Terms", "GENENAME", "DESCRIPTION")]
        filter_df$gene_name<- str_split(as.character(filter_df$ID), ".p", simplify = T)[,1]
        filter_df<- filter_df[ ,c("ID", "gene_name", "Scaffold.length", "Pfam", "InterPro", "GO.Terms", "GENENAME", "DESCRIPTION")]
        colnames(filter_df)<- c("gene name", "gene id", "scaffold length", tolower(colnames(filter_df[,4:ncol(filter_df)])))
        filter_df[,4:ncol(filter_df)] <- lapply(filter_df[4:ncol(filter_df)], gsub, pattern = ";", replacement = ",")
        filter_lista<- list_with_values(filter_df)
        final_list[[j]]<- filter_lista
        final_vector[j]<- paste(unlist(final_list[[j]]), collapse = ";")
    }
    return(final_vector)
}

create_df<- function(dataframe) {
    df_gff<- data.frame(
        seqid = dataframe[,41],
        source = rep(".", nrow(dataframe)),
        type = str_split(as.character(dataframe[,2]), "\\.", simplify = T)[,2],
        #start = dataframe[,3] + (dataframe[,6]-1),
        start = dataframe[,3],
        #end = dataframe[,3] + (dataframe[,7]-1),
        end = dataframe[,4],
        score = rep(".", nrow(dataframe)),
        strand = dataframe[,8],
        phase = rep(0, nrow(dataframe)),
        attributes = get_attributes(dataframe)
    )
    return(df_gff)
}

create_gff<- function(dataframe1, dataframe2, file) {
    lista<- list()
    partes<- unique(dataframe2$final_id)
    for (i in 1:length(unique(dataframe2$final_id))){
        df_filter<- subset(dataframe2, final_id == partes[i])
        nombres<- partes[i]
        minvalue<- 1
        maxvalue<- df_filter[,3]
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
create_gff(data_prueba, metadata, "output_prueba.gff")
create_gff(df_final, metadata, "judiseq.gff")
