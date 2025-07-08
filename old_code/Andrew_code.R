
Trypsin_path = "Data/all_search.tsv"


columns_to_read <-c('Genes','Run','Lib.PG.Q.Value','RT','Precursor.Id','Stripped.Sequence','Precursor.Mz',
                    'Precursor.Charge','Precursor.Quantity','Ms1.Area','Protein.Group','Translated.Q.Value','Channel.Q.Value')

Bulk_trypsin <- data.table::fread(Trypsin_path,select = columns_to_read)
# make unique precursors
Bulk_trypsin$seqcharge <- paste0(Bulk_trypsin$Stripped.Sequence,Bulk_trypsin$Precursor.Charge)


# Pulse was only with Lysine, so separate R and K peptides
Bulk_trypsin$type = str_sub(Bulk_trypsin$Stripped.Sequence, -1)
Bulk_R <- Bulk_trypsin %>% filter(type != 'K')
Bulk_K <- Bulk_trypsin %>% filter(type == 'K')


# Translation filtering for Heavy/Light peptides
Bulk_K <- Bulk_K %>% filter(Translated.Q.Value < .01)


# Filter for Heavy and light
Bulk_K$lab = str_sub(Bulk_K$Precursor.Id,-3, -3)
Bulk_K_L <- Bulk_K %>% filter(lab == 'L')
Bulk_K_H <- Bulk_K %>% filter(lab == 'H')



# Make Heavy and Light peptide X sample matricies 
Bulk_K_H <- dcast(Bulk_K_H,seqcharge + Protein.Group ~ Run, value.var = 'Ms1.Area')
Bulk_K_L <- dcast(Bulk_K_L,seqcharge + Protein.Group ~ Run, value.var = 'Ms1.Area')

sect = intersect(Bulk_K_H$seqcharge, Bulk_K_L$seqcharge)

Bulk_K_H <- Bulk_K_H %>% filter(seqcharge %in% sect)
Bulk_K_H <- Bulk_K_H %>% distinct(seqcharge,.keep_all = T)
Bulk_K_L <- Bulk_K_L %>% filter(seqcharge %in% sect)
Bulk_K_L <- Bulk_K_L %>% distinct(seqcharge,.keep_all = T)

rownames(Bulk_K_L) <- Bulk_K_L$seqcharge
rownames(Bulk_K_H) <- Bulk_K_H$seqcharge

# save peptide to protein mapping
PGsL <- cbind(Bulk_K_L$seqcharge,Bulk_K_L$Protein.Group)


Bulk_K_L$Protein.Group <- NULL
Bulk_K_L$seqcharge <- NULL
Bulk_K_H$seqcharge <- NULL
Bulk_K_H$Protein.Group <- NULL




## Compute abundance matrices for modeling turnover

# Total K peptide amount
Bulk_all_Kpep <- as.matrix(Bulk_K_H) + as.matrix(Bulk_K_L)
rownames(Bulk_all_Kpep) <- PGsL[,1]


# Take top 3 most abundant peptides for each protein to reduce noise
#Top3_abundant <- FindTop3(Bulk_all_Kpep,PGsL)





# Save fract Light over total and filtering out peptides with more heavy than
# should be possible with availible AA fraction
Bulk_half_life = log(Bulk_K_L/(Bulk_K_L+Bulk_K_H) + 1)

### You need to then divide this by the time point^^^
