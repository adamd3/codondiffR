#' @import taxonomizr
#' @import stringr
#' @import dplyr

## Prepare SQL database of NCBI taxon names and nodes, using taxonomizr
## functionality
prepareDatabase(sqlFile = "nameNode.sql",)

## There are two sets of data available from the Codon Usage Table Database:
## those for NCBI RefSeq and NCBI Genbank entries, respectively.
## The latter dataset contains more sequences and is used here.
## (See: hive.biochemistry.gwu.edu/review/codon)

GenbankCU <- read.table(
    "../o569942-genbank_species.tsv",
    sep = "\t",
    header = TRUE,
    comment.char = "?",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
)


table(GenbankCU$Organelle)
#         genomic mitochondrion       plastid
# 1        404993        468875        130606

## Keep only genomic entries
GenbankGenomic <- subset(GenbankCU, Organelle == "genomic")

nrow(GenbankCU)
##[1] 1004475
nrow(GenbankGenomic)
##[1] 404993



##------------------------------------------------------------------------------
## Annotate CUT database entries with taxonomic ranks
##------------------------------------------------------------------------------
taxa_domains_Genbank <- getTaxonomy(
    GenbankGenomic$Taxid, "../nameNode.sql",
    desiredTaxa = c("superkingdom")
)
taxa_kingdoms_Genbank <- getTaxonomy(
    GenbankGenomic$Taxid, "../nameNode.sql",
    desiredTaxa = c("kingdom")
)
taxa_phyla_Genbank <- getTaxonomy(
    GenbankGenomic$Taxid, "../nameNode.sql",
    desiredTaxa = c("phylum")
)


table(taxa_domains_Genbank, useNA = "always")
# Archaea  Bacteria Eukaryota   Viruses      <NA>
#    1327     71061    170191    161797       617

prop.table(table(taxa_domains_Genbank, useNA = "always"))
#     Archaea    Bacteria   Eukaryota     Viruses        <NA>
# 0.003276600 0.175462292 0.420231955 0.399505670 0.001523483

table(taxa_kingdoms_Genbank, useNA = "always")
# Fungi       Metazoa Viridiplantae          <NA>
# 30241        105286         29966        239500

prop.table(table(taxa_kingdoms_Genbank, useNA = "always"))
#      Fungi       Metazoa Viridiplantae          <NA>
# 0.07467043    0.25996993    0.07399140    0.59136825

table(taxa_phyla_Genbank, useNA = "always")
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
#           <NA>
#         167156


GenbankGenomic$Domain <- taxa_domains_Genbank
GenbankGenomic$Kingdom <- taxa_kingdoms_Genbank
GenbankGenomic$Phylum <- taxa_phyla_Genbank

## remove viruses
nrow(GenbankGenomic)
# 404993

GenbankGenomic <- subset(GenbankGenomic, !(Domain == "Viruses"))

nrow(GenbankGenomic)
# 242579



##------------------------------------------------------------------------------
## Get relative frequency of each codon
##------------------------------------------------------------------------------
codon_cols <- names(GenbankGenomic)[13:76]

GenbankCodonSort <- GenbankGenomic[13:76] %>%
    mutate(sum_freq_codpos = rowSums(.[codon_cols])) %>%
    mutate_at(codon_cols, funs(./sum_freq_codpos)) %>%
    select(-c(sum_freq_codpos)) %>%
    select(sort(current_vars())) ## sort columns alphabetically


## Select only required cols
GenbankCUTD <- cbind(GenbankGenomic[c(3,4,8,77,78,79)], GenbankCodonSort)


# ## Remove all missing data
# GenbankCUTD <- na.omit(GenbankCUTD)


table(GenbankCUTD$Domain)
#   Archaea  Bacteria Eukaryota
#      1327     71061    170191

prop.table(table(GenbankCUTD$Domain))
#     Archaea    Bacteria   Eukaryota
# 0.005470383 0.292939620 0.701589997

## Count numbers/proportions for a given length threshold:
table(subset(GenbankCUTD, X..Codons >= 600)$Domain)
# Archaea  Bacteria Eukaryota
#    1050     27402     72372

prop.table(table(subset(GenbankCUTD, X..Codons >= 600)$Domain))
#    Archaea   Bacteria  Eukaryota
# 0.01041419 0.27178053 0.71780528



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
GenbankAnn <- cbind(stdgc, t(GenbankCUTD[7:ncol(GenbankCUTD)]))
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
# [1] 5.424251
## 5.4% of values are NaN: indicates that the amino acid in question is not used

## e.g. to check sequences with no Asparagine (N) amino acid [AAC/AAT codon]
head(subset(GenbankCUTD, AAC == 0 & AAT == 0))
# Taxid                        Species X..Codons superkingdom kingdom
# 357  1321993                Vibrio sp. 3127        26     Bacteria    <NA>
# 609  1835519        Frankia sp. AcoItVII.11       190     Bacteria    <NA>
# 655   362757              Vibrio sp. CCH-64        41     Bacteria    <NA>
# 1031  362570        Klebsiella sp. SO-Y1-40        63     Bacteria    <NA>
# 1443 1221370    Sediminicola sp. YIK-SED-23       118     Bacteria    <NA>
# 1592  212171 Actinomycetales bacterium GP-7        69     Bacteria    <NA>
#       phylum         AAA AAC        AAG AAT         ACA        ACC
# 357  Proteobacteria 0.038461538   0 0.00000000   0 0.000000000 0.00000000
# 609  Actinobacteria 0.005263158   0 0.00000000   0 0.000000000 0.03157895
# 655  Proteobacteria 0.000000000   0 0.00000000   0 0.000000000 0.00000000
# 1031 Proteobacteria 0.000000000   0 0.00000000   0 0.079365079 0.01587302
# 1443  Bacteroidetes 0.033898305   0 0.01694915   0 0.008474576 0.01694915
# 1592 Actinobacteria 0.000000000   0 0.00000000   0 0.000000000 0.01449275


gbnorm[is.nan(gbnorm)] <- NA
gbnorm <- cbind(GenbankCUTD[1:6], gbnorm)



# ##------------------------------------------------------------------------------
# ## Handle missing data (NAs)
# ##------------------------------------------------------------------------------
# ## Certain taxa, even with large "X..Codons" counts, only utilise a small number
# ## of codons -- for example, Taxid 1497514 (Colletotrichum sp. IVS-2014c)
# ## has numerous Genbank entries of highly similar short nucleotide sequences
# ## under different accession numbers (KJ490360, KJ490361, KJ490362, etc )
# gbnorm$NAcount <- rowSums(is.na(gbnorm[7:ncol(gbnorm)]))
#
# nrow(gbnorm)
# # 242579
#
# nrow(subset(gbnorm, NAcount < 20))
# # 239610
#
# gbnorm <- subset(gbnorm, NAcount < 20)
#
# gbnorm$NAcount <- NULL



##------------------------------------------------------------------------------
## Save objects for use in package
##------------------------------------------------------------------------------
usethis::use_data(gbnorm, stdgc, internal = TRUE, overwrite = TRUE)
