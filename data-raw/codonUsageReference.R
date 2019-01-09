#' @import taxonomizr
#' @import dplyr

## prepare SQL database of NCBI taxon names and nodes
prepareDatabase(sqlFile = "../inst/extdata/nameNode.sql",)


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

RefSeqSubset <- subset(RefSeqCU, Organelle == "genomic")
GenbankSubset <- subset(GenbankCU, Organelle == "genomic")

nrow(RefSeqCU)
##[1] 130399
nrow(RefSeqSubset)
##[1] 130044

nrow(GenbankCU)
##[1] 1004475
nrow(GenbankSubset)
##[1] 404993


## get taxonomic ranks of taxids
taxa_domains_RefSeq <- getTaxonomy(
    RefSeqSubset$Taxid, "../inst/extdata/nameNode.sql",
    desiredTaxa = c("superkingdom")
)
taxa_domains_Genbank <- getTaxonomy(
    GenbankSubset$Taxid, "../inst/extdata/nameNode.sql",
    desiredTaxa = c("superkingdom")
)

taxa_phyla_RefSeq <- getTaxonomy(
    RefSeqSubset$Taxid, "../inst/extdata/nameNode.sql",
    desiredTaxa = c("phylum")
)
taxa_phyla_Genbank <- getTaxonomy(
    GenbankSubset$Taxid, "../inst/extdata/nameNode.sql",
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


RefSeqSubset$Domain <- taxa_domains_RefSeq
RefSeqSubset$Phylum <- taxa_phyla_RefSeq

GenbankSubset$Domain <- taxa_domains_Genbank
GenbankSubset$Phylum <- taxa_phyla_Genbank


## get proportion of each codon
# [13] "TTT"               "TTC"               "TTA"
# [16] "TTG"               "CTT"               "CTC"
# [19] "CTA"               "CTG"               "ATT"
# [22] "ATC"               "ATA"               "ATG"
# [25] "GTT"               "GTC"               "GTA"
# [28] "GTG"               "TAT"               "TAC"
# [31] "TAA"               "TAG"               "CAT"
# [34] "CAC"               "CAA"               "CAG"
# [37] "AAT"               "AAC"               "AAA"
# [40] "AAG"               "GAT"               "GAC"
# [43] "GAA"               "GAG"               "TCT"
# [46] "TCC"               "TCA"               "TCG"
# [49] "CCT"               "CCC"               "CCA"
# [52] "CCG"               "ACT"               "ACC"
# [55] "ACA"               "ACG"               "GCT"
# [58] "GCC"               "GCA"               "GCG"
# [61] "TGT"               "TGC"               "TGA"
# [64] "TGG"               "CGT"               "CGC"
# [67] "CGA"               "CGG"               "AGT"
# [70] "AGC"               "AGA"               "AGG"
# [73] "GGT"               "GGC"               "GGA"
# [76] "GGG"

names_to_transform <- names(RefSeqSubset)[13:76]


## get percentage of each codon
colProps <- RefSeqSubset[13:76] %>%
  mutate(sum_freq_codpos = rowSums(.[names_to_transform])) %>%
  mutate_at(names_to_transform, funs((./sum_freq_codpos)*100))


RefSeqSubset %>%
  group_by(TTT, TTC) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))

RefSeqSubset[13:76]

prop.table(RefSeqSubset[13:76])

apply(RefSeqSubset[13:76], 1, function(x) x/rowSums(RefSeqSubset[13:76]))

apply(names(RefSeqSubset)[13:76], 1, function(x)) {
  RefSeqSubset[paste0(x, "_pct")] <<- RefSeqSubset[x] / sum(RefSeqSubset[13:76])
}

for (row in 1:nrow(RefSeqSubset)) {
  RefSeqSubset[paste0(col, "_pct")] = fruits[col] / sum(fruits[col])
}

usethis::use_data(RefSeqSubset, GenbankSubset, internal = TRUE)

## TO DO: put data in R/sysdata.rda -- you can use usethis::use_data() to create
##  this file with the argument internal = TRUE:
## (see: http://r-pkgs.had.co.nz/data.html; section "Internal data")

## see also https://coolbutuseless.github.io/2018/12/10/r-packages-internal-and-external-data/

## NB: devtools::use_data() is DEPRECATED!
