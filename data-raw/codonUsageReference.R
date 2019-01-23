#' @import taxonomizr
#' @import stringr
#' @import dplyr

## Prepare SQL database of NCBI taxon names and nodes, using taxonomizr
## functionality
prepareDatabase(sqlFile = "../inst/extdata/nameNode.sql",)

## There are two sets of data available from the Codon Usage Table Database:
## those for NCBI RefSeq and NCBI Genbank entries, respectively.
## See: hive.biochemistry.gwu.edu/review/codon
## The latter dataset contains more sequences and is used here.

RefSeqCU <- read.table(
    "../inst/extdata/o569942-refseq_species.tsv",
    sep = "\t",
    header = TRUE,
    comment.char = "?",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
)

GenbankCU <- read.table(
    "../inst/extdata/o569942-genbank_species.tsv",
    sep = "\t",
    header = TRUE,
    comment.char = "?",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
)


table(RefSeqCU$Organelle)
## genomic mitochondrion       plastid
##  130044           289            66


table(GenbankCU$Organelle)
#         genomic mitochondrion       plastid
# 1        404993        468875        130606

## Keep only genomic entries
RefSeqGenomic <- subset(RefSeqCU, Organelle == "genomic")
GenbankGenomic <- subset(GenbankCU, Organelle == "genomic")

nrow(RefSeqCU)
##[1] 130399
nrow(RefSeqGenomic)
##[1] 130044

nrow(GenbankCU)
##[1] 1004475
nrow(GenbankGenomic)
##[1] 404993



##------------------------------------------------------------------------------
## Annotate CUT database entries with taxonomic ranks
##------------------------------------------------------------------------------
taxa_domains_RefSeq <- getTaxonomy(
    RefSeqGenomic$Taxid, "../inst/extdata/nameNode.sql",
    desiredTaxa = c("superkingdom")
)
taxa_phyla_RefSeq <- getTaxonomy(
    RefSeqGenomic$Taxid, "../inst/extdata/nameNode.sql",
    desiredTaxa = c("phylum")
)

taxa_domains_Genbank <- getTaxonomy(
    GenbankGenomic$Taxid, "../inst/extdata/nameNode.sql",
    desiredTaxa = c("superkingdom")
)
taxa_phyla_Genbank <- getTaxonomy(
    GenbankGenomic$Taxid, "../inst/extdata/nameNode.sql",
    desiredTaxa = c("phylum")
)


table(taxa_domains_RefSeq)
# Archaea  Bacteria Eukaryota   Viruses
#     721    120860       882      7441

prop.table(table(taxa_domains_RefSeq))
#     Archaea    Bacteria   Eukaryota     Viruses
# 0.005550252 0.930379357 0.006789629 0.057280761

table(taxa_phyla_RefSeq)
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
#  Verrucomicrobia
#              104


table(taxa_domains_Genbank)
# Archaea  Bacteria Eukaryota   Viruses
#    1327     71061    170191    161797

prop.table(table(taxa_domains_Genbank))
#     Archaea    Bacteria   Eukaryota     Viruses
# 0.003281599 0.175730014 0.420873148 0.400115239

table(taxa_phyla_Genbank)
# Acanthocephala               Acidobacteria
#              6                          22
# Actinobacteria                    Annelida
#           9927                        1164
#    Apicomplexa                   Aquificae
#            908                         125
# Armatimonadetes                  Arthropoda
#              4                       58549
#     Ascomycota             Bacillariophyta
#          22389                         238
#  Bacteroidetes               Basidiomycota
#           1332                        7064
# Blastocladiomycota                 Brachiopoda
#             13                          19
#        Bryozoa                 Caldiserica
#             28                           1
# Calditrichaeota     candidate division NC10
#              1                           2
# candidate division WWE3    Candidatus Cloacimonetes
#              1                           1
# Candidatus Gracilibacteria     Candidatus Korarchaeota
#              1                           1
# Candidatus Melainabacteria    Candidatus Micrarchaeota
#              1                           1
# Candidatus Omnitrophica    Candidatus Parvarchaeota
#              1                           1
# Candidatus Saccharibacteria    Candidatus Tectomicrobia
#              3                           2
#   Chaetognatha                  Chlamydiae
#              9                         195
#       Chlorobi                 Chloroflexi
#             34                          87
#    Chlorophyta                    Chordata
#            604                       37273
#     Chromerida              Chrysiogenetes
#              4                           2
# Chytridiomycota                    Cnidaria
#             89                         969
#   Colponemidia        Coprothermobacterota
#              3                           1
#  Crenarchaeota                Cryptomycota
#            140                           2
#     Ctenophora               Cyanobacteria
#             42                        6224
#    Cycliophora             Deferribacteres
#              3                           4
# Deinococcus-Thermus                 Dictyoglomi
#            104                           3
#  Echinodermata               Elusimicrobia
#            451                           6
#     Entoprocta                   Euglenida
#              6                          69
#  Euryarchaeota               Fibrobacteres
#           1146                           3
#     Firmicutes                Fusobacteria
#           5552                          65
#   Gastrotricha            Gemmatimonadetes
#              8                           3
# Gnathostomulida               Haplosporidia
#             16                          13
#   Hemichordata             Ignavibacteriae
#             18                           3
#    Kinorhyncha          Kiritimatiellaeota
#              7                           2
#  Lentisphaerae                  Loricifera
#              1                           1
#  Microsporidia                    Mollusca
#            100                        5090
#   Mucoromycota               Nanoarchaeota
#            344                           2
#       Nematoda                Nematomorpha
#            578                           3
#       Nemertea                 Nitrospinae
#            157                           1
#    Nitrospirae               Olpidiomycota
#             22                           2
#    Onychophora                    Placozoa
#             35                          15
# Planctomycetes             Platyhelminthes
#            118                         390
#       Porifera                  Priapulida
#            247                           6
# Proteobacteria                   Rhombozoa
#          43641                           6
#       Rotifera                Spirochaetes
#            102                        1179
#   Streptophyta               Synergistetes
#          29362                           8
#     Tardigrada                 Tenericutes
#             35                        1041
# Thaumarchaeota       Thermodesulfobacteria
#             32                           9
#    Thermotogae             Verrucomicrobia
#             59                          50
# Xenacoelomorpha               Zoopagomycota
#             52                         184


RefSeqGenomic$Domain <- taxa_domains_RefSeq
RefSeqGenomic$Phylum <- taxa_phyla_RefSeq

GenbankGenomic$Domain <- taxa_domains_Genbank
GenbankGenomic$Phylum <- taxa_phyla_Genbank

##------------------------------------------------------------------------------
## Get relative frequency of each codon
##------------------------------------------------------------------------------
codon_cols <- names(RefSeqGenomic)[13:76]

RefSeqGenomic[13:76] <- RefSeqGenomic[13:76] %>%
    mutate(sum_freq_codpos = rowSums(.[codon_cols])) %>%
    mutate_at(codon_cols, funs(./sum_freq_codpos)) %>%
    select(-c(sum_freq_codpos)) %>%
    select(sort(current_vars())) ## sort codons alphabetically

codon_cols <- names(GenbankGenomic)[13:76]

GenbankCodonSort <- GenbankGenomic[13:76] %>%
    mutate(sum_freq_codpos = rowSums(.[codon_cols])) %>%
    mutate_at(codon_cols, funs(./sum_freq_codpos)) %>%
    select(-c(sum_freq_codpos)) %>%
    select(sort(current_vars())) ## sort codons alphabetically


## Select only required cols
RefSeqCUTD <- cbind(RefSeqGenomic[c(3,4,8,77,78)], RefSeqCodonSort)
GenbankCUTD <- cbind(GenbankGenomic[c(3,4,8,77,78)], GenbankCodonSort)

## Remove all missing data
RefSeqCUTD <- na.omit(RefSeqCUTD)
GenbankCUTD <- na.omit(GenbankCUTD)

table(GenbankCUTD$Domain)
# Archaea  Bacteria Eukaryota
#    1323     69841    166673

prop.table(table(GenbankCUTD$Domain))
#     Archaea    Bacteria   Eukaryota
# 0.005562633 0.293650694 0.700786673



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
GenbankAnn <- cbind(stdgc, t(GenbankCUTD[6:69]))
GenbankAnn$AA <- as.factor(GenbankAnn$AA)
gblist <- split(GenbankAnn, GenbankAnn$AA)
gblist <- lapply(gblist, function(x) {
    sweep(
        x[2:(ncol((gblist)[[1]]))],
        2,
        colSums(x[2:ncol((gblist)[[1]])]),
        FUN = "/"
    )
})
gbnorm <- do.call("rbind", gblist)

## For the M and W amino acids, there is only a single codon; therefore, must
## re-insert the codons to rownames
for (uniqaa in c("M", "W")) {
    rownames(gbnorm)[grepl(uniqaa, rownames(gbnorm))] <- paste0(
        rownames(gbnorm)[grepl(uniqaa, rownames(gbnorm))],
        ".",
        rownames(stdgc)[grepl(uniqaa, stdgc$AA)]
    )
}
rownames(gbnorm) <- str_split_fixed(rownames(gbnorm), ".", 3)[,3]
gbnorm <- t(gbnorm[order(row.names(gbnorm)), ])

((length(gbnorm[is.nan(gbnorm)]))/(length(gbnorm)))*100
# [1] 2.308731
## 2.3% of values are NaN: indicates that the amino acid in question is not used

## e.g. to check sequences with no Asparagine (N) amino acid [AAC/AAT codon]
head(subset(GenbankCUTD, AAC == 0 & AAT == 0))
#       Taxid                           Species X..Codons superkingdom
# 60  1198031           Pasteurella sp. medi 11       300     Bacteria
# 109 1391678            Brachyspira sp. KL-194       420     Bacteria
# 146  343139             Comamonas sp. R-28235       144     Bacteria
# 148  537985 Staphylococcus sp. 020902-022-273       123     Bacteria
# 154  907839            Legionella sp. ST19309       208     Bacteria
# 155   45495               Methylophaga marina        89     Bacteria
#             phylum        AAA AAC         AAG AAT ACA         ACC         ACG
# 60  Proteobacteria 0.00000000   0 0.026666667   0   0 0.003333333 0.043333333
# 109   Spirochaetes 0.00000000   0 0.030952381   0   0 0.000000000 0.000000000
# 146 Proteobacteria 0.09722222   0 0.000000000   0   0 0.020833333 0.006944444
# 148     Firmicutes 0.00000000   0 0.000000000   0   0 0.000000000 0.081300813
# 154 Proteobacteria 0.01442308   0 0.004807692   0   0 0.004807692 0.004807692
# 155 Proteobacteria 0.05617978   0 0.000000000   0   0 0.000000000 0.022471910

gbnorm[is.nan(gbnorm)] <- NA
gbnorm <- cbind(GenbankCUTD[1:5], gbnorm)



##------------------------------------------------------------------------------
## Save objects for use in package
##------------------------------------------------------------------------------
usethis::use_data(GenbankCUTD, gbnorm, stdgc, internal = TRUE, overwrite = TRUE)
