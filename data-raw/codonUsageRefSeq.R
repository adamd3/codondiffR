library(taxonomizr)
library(stringr)
library(dplyr)

## Prepare SQL database of NCBI taxon names and nodes, using taxonomizr
## functionality
prepareDatabase(sqlFile = "nameNode.sql",)

## The RefSeq codon usage table comes from the Codon Usage Table Database
## available at: hive.biochemistry.gwu.edu/review/codon

RefSeqCU <- read.table(
    "../inst/extdata/o569942-refseq_species.tsv",
    sep = "\t",
    header = TRUE,
    comment.char = "?",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
)


table(RefSeqCU$Organelle)
# genomic mitochondrion       plastid
#  130044           289            66

## Keep only genomic entries
RefSeqGenomic <- subset(RefSeqCU, Organelle == "genomic")

nrow(RefSeqCU)
##[1] 130399
nrow(RefSeqGenomic)
##[1] 130044



##------------------------------------------------------------------------------
## Annotate CUT database entries with taxonomic ranks
##------------------------------------------------------------------------------
taxa_domains_RefSeq <- getTaxonomy(
    RefSeqGenomic$Taxid, "../nameNode.sql",
    desiredTaxa = c("superkingdom")
)
taxa_kingdoms_RefSeq <- getTaxonomy(
    RefSeqGenomic$Taxid, "../nameNode.sql",
    desiredTaxa = c("kingdom")
)
taxa_phyla_RefSeq <- getTaxonomy(
    RefSeqGenomic$Taxid, "../nameNode.sql",
    desiredTaxa = c("phylum")
)


table(taxa_domains_RefSeq, useNA = "always")
# Archaea  Bacteria Eukaryota   Viruses      <NA>
#     721    120860       882      7441       140


prop.table(table(taxa_domains_RefSeq, useNA = "always"))
#     Archaea    Bacteria   Eukaryota     Viruses        <NA>
# 0.005544277 0.929377749 0.006782320 0.057219095 0.001076559


table(taxa_kingdoms_RefSeq, useNA = "always")
# Fungi       Metazoa Viridiplantae          <NA>
#   277           412           102        129253


prop.table(table(taxa_kingdoms_RefSeq, useNA = "always"))
#        Fungi       Metazoa Viridiplantae          <NA>
# 0.0021300483  0.0031681585  0.0007843499  0.9939174433


table(taxa_phyla_RefSeq, useNA = "always")
# Abditibacteriota               Acidobacteria
#                1                          39
#   Actinobacteria                    Annelida
#            13624                           1
#      Apicomplexa                   Aquificae
#               33                          28
#  Armatimonadetes                  Arthropoda
#                8                         126
#       Ascomycota             Bacillariophyta
#              208                           2
#    Bacteroidetes                Balneolaeota
#             2445                           9
#    Basidiomycota                 Brachiopoda
#               53                           1
#      Caldiserica             Calditrichaeota
#                1                           2
# candidate division NC10     Candidatus Atribacteria
#                1                           2
# Candidatus Cloacimonetes     Candidatus Korarchaeota
#                1                           1
# Candidatus Kryptonia  Candidatus Melainabacteria
#               13                           1
# Candidatus Micrarchaeota Candidatus Saccharibacteria
#                1                           2
# Candidatus Tectomicrobia                  Chlamydiae
#                2                         372
#         Chlorobi                 Chloroflexi
#               19                          75
#      Chlorophyta                    Chordata
#               11                         252
#   Chrysiogenetes             Chytridiomycota
#                2                           2
#         Cnidaria        Coprothermobacterota
#                6                           4
#    Crenarchaeota               Cyanobacteria
#              131                         500
#  Deferribacteres         Deinococcus-Thermus
#                7                          90
#      Dictyoglomi               Echinodermata
#                2                           2
#    Elusimicrobia               Euryarchaeota
#                5                         567
#    Fibrobacteres                  Firmicutes
#               43                       38399
#     Fusobacteria            Gemmatimonadetes
#              200                           3
#     Hemichordata             Ignavibacteriae
#                1                           2
# Kiritimatiellaeota               Lentisphaerae
#                1                           2
#    Microsporidia                    Mollusca
#               10                           8
#     Mucoromycota                    Nematoda
#                4                           8
#      Nitrospinae                 Nitrospirae
#                1                          27
#         Placozoa              Planctomycetes
#                1                          52
#  Platyhelminthes                    Porifera
#                4                           1
#       Priapulida              Proteobacteria
#                1                       63328
#  Rhodothermaeota                Spirochaetes
#                5                         810
#     Streptophyta               Synergistetes
#               91                          28
#      Tenericutes              Thaumarchaeota
#              485                          21
# Thermodesulfobacteria                 Thermotogae
#               11                          81
#  Verrucomicrobia                        <NA>
#              104                        7660


RefSeqGenomic$Domain <- taxa_domains_RefSeq
RefSeqGenomic$Kingdom <- taxa_kingdoms_RefSeq
RefSeqGenomic$Phylum <- taxa_phyla_RefSeq

## remove viruses
nrow(RefSeqGenomic)
# 130044

RefSeqGenomic <- subset(RefSeqGenomic, !(Domain == "Viruses"))

nrow(RefSeqGenomic)
# 122463



##------------------------------------------------------------------------------
## Get relative frequency of each codon
##------------------------------------------------------------------------------
codon_cols <- names(RefSeqGenomic)[13:76]

RefSeqCodonSort <- RefSeqGenomic[13:76] %>%
    mutate(sum_freq_codpos = rowSums(.[codon_cols])) %>%
    mutate_at(codon_cols, funs(./sum_freq_codpos)) %>%
    dplyr::select(-c(sum_freq_codpos)) %>%
    dplyr::select(sort(current_vars())) ## sort columns alphabetically


## Select only required cols
RefSeqCUTD <- cbind(RefSeqGenomic[c(3,4,8,77,78,79)], RefSeqCodonSort)


# ## Remove all missing data
# RefSeqCUTD <- na.omit(RefSeqCUTD)


table(RefSeqCUTD$Domain)
# Archaea  Bacteria Eukaryota
#     721    120860       882

prop.table(table(RefSeqCUTD$Domain))
#     Archaea    Bacteria   Eukaryota
# 0.005887493 0.986910332 0.007202175


## Count numbers/proportions for a given length threshold:
table(subset(RefSeqCUTD, X..Codons >= 600)$Domain, useNA="always")
# Archaea  Bacteria Eukaryota      <NA>
#     721    120860       882         0

prop.table(table(subset(RefSeqCUTD, X..Codons >= 600)$Domain))
#     Archaea    Bacteria   Eukaryota
# 0.005887493 0.986910332 0.007202175



table(subset(RefSeqCUTD, X..Codons >= 600)$Phylum, useNA="always")
# Abditibacteriota               Acidobacteria
#                1                          39
#   Actinobacteria                    Annelida
#            13624                           1
#      Apicomplexa                   Aquificae
#               33                          28
#  Armatimonadetes                  Arthropoda
#                8                         126
#       Ascomycota             Bacillariophyta
#              208                           2
#    Bacteroidetes                Balneolaeota
#             2445                           9
#    Basidiomycota                 Brachiopoda
#               53                           1
#      Caldiserica             Calditrichaeota
#                1                           2
# candidate division NC10     Candidatus Atribacteria
#                1                           2
# Candidatus Cloacimonetes     Candidatus Korarchaeota
#                1                           1
# Candidatus Kryptonia  Candidatus Melainabacteria
#               13                           1
# Candidatus Micrarchaeota Candidatus Saccharibacteria
#                1                           2
# Candidatus Tectomicrobia                  Chlamydiae
#                2                         372
#         Chlorobi                 Chloroflexi
#               19                          75
#      Chlorophyta                    Chordata
#               11                         252
#   Chrysiogenetes             Chytridiomycota
#                2                           2
#         Cnidaria        Coprothermobacterota
#                6                           4
#    Crenarchaeota               Cyanobacteria
#              131                         500
#  Deferribacteres         Deinococcus-Thermus
#                7                          90
#      Dictyoglomi               Echinodermata
#                2                           2
#    Elusimicrobia               Euryarchaeota
#                5                         567
#    Fibrobacteres                  Firmicutes
#               43                       38399
#     Fusobacteria            Gemmatimonadetes
#              200                           3
#     Hemichordata             Ignavibacteriae
#                1                           2
# Kiritimatiellaeota               Lentisphaerae
#                1                           2
#    Microsporidia                    Mollusca
#               10                           8
#     Mucoromycota                    Nematoda
#                4                           8
#      Nitrospinae                 Nitrospirae
#                1                          27
#         Placozoa              Planctomycetes
#                1                          52
#  Platyhelminthes                    Porifera
#                4                           1
#       Priapulida              Proteobacteria
#                1                       63328
#  Rhodothermaeota                Spirochaetes
#                5                         810
#     Streptophyta               Synergistetes
#               91                          28
#      Tenericutes              Thaumarchaeota
#              485                          21
# Thermodesulfobacteria                 Thermotogae
#               11                          81
#  Verrucomicrobia                        <NA>
#              104                          79

##------------------------------------------------------------------------------
## Get the Standard Genetic Code
##------------------------------------------------------------------------------
stdgc <- data.frame(
    AA = as.character(GENETIC_CODE),
    codon = names(GENETIC_CODE),
    stringsAsFactors = FALSE
) %>% arrange(codon)

rownames(stdgc) <- stdgc$codon
stdgc$codon <- NULL



##------------------------------------------------------------------------------
## Normalise the reference database values
##------------------------------------------------------------------------------
## Normalise such that the codon frequencies per amino acid sum to 1
RefSeqAnn <- cbind(stdgc, t(RefSeqCUTD[7:ncol(RefSeqCUTD)]))
RefSeqAnn$AA <- as.factor(RefSeqAnn$AA)
gblist <- split(RefSeqAnn, RefSeqAnn$AA)
gblist <- lapply(gblist, function(x) {
    sweep(
        x[2:(ncol((gblist)[[1]]))],
        2,
        colSums(x[2:ncol((gblist)[[1]])]),
        FUN = "/"
    )
})
refseqnorm <- do.call("rbind", gblist)

## For the M and W amino acids, there is only a single codon; therefore, must
## re-insert the codons to rownames
for (uniqaa in c("M", "W")) {
    rownames(refseqnorm)[grepl(uniqaa, rownames(refseqnorm))] <- paste0(
        rownames(refseqnorm)[grepl(uniqaa, rownames(refseqnorm))],
        ".",
        rownames(stdgc)[grepl(uniqaa, stdgc$AA)]
    )
}
rownames(refseqnorm) <- str_split_fixed(rownames(refseqnorm), ".", 3)[,3]
refseqnorm <- t(refseqnorm[order(row.names(refseqnorm)), ])

((length(refseqnorm[is.nan(refseqnorm)]))/(length(refseqnorm)))*100
# [1] 0

refseqnorm <- cbind(RefSeqCUTD[1:6], refseqnorm)

## Sequences with no annotated Kingdom (Kingdom == NA) are mostly prokaryotes:
table(subset(refseqnorm, is.na(Kingdom))$Domain)
# Archaea  Bacteria Eukaryota
#     721    120860        91

## In contrast, fewer than 100 sequences have no annotated Phylum:
table(subset(refseqnorm, is.na(Phylum))$Domain)
# Bacteria Eukaryota
#       23        56


##------------------------------------------------------------------------------
## Save objects for use in package
##------------------------------------------------------------------------------
# saveRDS(refseqnorm, file = "refseq_database.rds")
# saveRDS(stdgc, file = "refseq_database.rds")

save.image(file = "refseq_database.RData")

gbnorm <- refseqnorm
usethis::use_data(gbnorm, stdgc, internal = TRUE, overwrite = TRUE)
