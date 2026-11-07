###############################################################################
## OneK1K_nodes.R  --  candidate genes: KEGG pathway hsa04662 "B cell receptor signaling pathway" (Homo sapiens)
##
## Source: KEGG PATHWAY hsa04662 (https://www.kegg.jp/entry/hsa04662), gene list retrieved with the KEGG REST API
## (https://rest.kegg.jp/link/hsa/hsa04662 and /get/hsa04662) on the date in KEGG_ACCESSED: 91 entries, of which 82 are named
## genes (9 unnamed LOC entries and the P3R3URF-PIK3R3 readthrough are omitted). Genes that share a KEGG orthology
## group (KO) form one node; the node name is the KO description (or the gene symbol for a singleton).
## Sourced by OneK1K_PreProcessing.R and OneK1K_ModelFitting.R. A node with several genes
## contributes every gene that has an instrument; figures and tables label genes by symbol.
##
## Reference: Kanehisa, M. & Goto, S. (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30.
###############################################################################

KEGG_PATHWAY <- "hsa04662"
KEGG_ACCESSED <- "2026-10-26"

nodes <- list(
  "AKT" = c("AKT3", "AKT1", "AKT2"),   # K04456
  "LILR (inhibitory/activating receptors)" = c("LILRB2", "LILRB1", "LILRB5", "LILRB4", "LILRA1", "LILRB3", "LILRA3", "LILRA2", "LILRA4", "LILRA5", "LILRA6"),   # K06512
  "VAV" = c("VAV3", "VAV1", "VAV2"),   # K05730
  "MALT1" = c("MALT1"),   # K07369
  "CHUK" = c("CHUK"),   # K04467
  "PIK3AP1" = c("PIK3AP1"),   # K12230
  "CR2" = c("CR2"),   # K04012
  "FCGR2B" = c("FCGR2B"),   # K12560
  "FOS" = c("FOS"),   # K04379
  "RASGRP3" = c("RASGRP3"),   # K12362
  "DAPP1" = c("DAPP1"),   # K12229
  "GRB2" = c("GRB2"),   # K04364
  "GSK3B" = c("GSK3B"),   # K03083
  "BLNK" = c("BLNK"),   # K07371
  "HRAS" = c("HRAS"),   # K02833
  "IKBKB" = c("IKBKB"),   # K07209
  "INPP5D" = c("INPP5D"),   # K03084
  "INPPL1" = c("INPPL1"),   # K15909
  "JUN" = c("JUN"),   # K04448
  "KRAS" = c("KRAS"),   # K07827
  "LYN" = c("LYN"),   # K05854
  "NFATC1" = c("NFATC1"),   # K04446
  "NFATC2" = c("NFATC2"),   # K17332
  "NFATC3" = c("NFATC3"),   # K17333
  "NFKB1" = c("NFKB1"),   # K02580
  "NFKBIA" = c("NFKBIA"),   # K04734
  "NFKBIB" = c("NFKBIB"),   # K02581
  "NFKBIE" = c("NFKBIE"),   # K05872
  "NRAS" = c("NRAS"),   # K07828
  "PI3K p110" = c("PIK3CA", "PIK3CB", "PIK3CD"),   # K00922
  "PI3K p85" = c("PIK3R1", "PIK3R2", "PIK3R3"),   # K02649
  "PLCG2" = c("PLCG2"),   # K05859
  "BANK1" = c("BANK1"),   # K26239
  "Calcineurin A" = c("PPP3CA", "PPP3CB", "PPP3CC"),   # K04348
  "Calcineurin B" = c("PPP3R1", "PPP3R2"),   # K06268
  "PRKCB" = c("PRKCB"),   # K19662
  "ERK1/2" = c("MAPK1", "MAPK3"),   # K04371
  "MEK1" = c("MAP2K1"),   # K04368
  "MEK2" = c("MAP2K2"),   # K04369
  "PTPN6" = c("PTPN6"),   # K05697
  "RAC1" = c("RAC1"),   # K04392
  "RAC2" = c("RAC2"),   # K07860
  "RAC3" = c("RAC3"),   # K07861
  "RAF1" = c("RAF1"),   # K04366
  "RELA" = c("RELA"),   # K04735
  "BLK" = c("BLK"),   # K08890
  "SOS" = c("SOS1", "SOS2"),   # K03099
  "SYK" = c("SYK"),   # K05855
  "BTK" = c("BTK"),   # K07370
  "CARD11" = c("CARD11"),   # K07367
  "IKBKG" = c("IKBKG"),   # K07210
  "IFITM1" = c("IFITM1"),   # K19831
  "BCL10" = c("BCL10"),   # K07368
  "CD19" = c("CD19"),   # K06465
  "CD22" = c("CD22"),   # K06467
  "CD72" = c("CD72"),   # K06504
  "CD79A" = c("CD79A"),   # K06506
  "CD79B" = c("CD79B"),   # K06507
  "CD81" = c("CD81")   # K06508
)

## Entrez IDs and KO groups (for the record)
kegg_genes <- data.frame(
  entrez = c(10000, 10288, 10451, 10859, 10892, 10990, 11006, 11024, 11025, 11026, 11027, 1147, 118788, 1380, 207, 208, 2213, 2353, 23547, 25780, 27071, 2885, 2932, 29760, 3265, 353514, 3551, 3635, 3636, 3725, 3845, 4067, 4772, 4773, 4775, 4790, 4792, 4793, 4794, 4893, 5290, 5291, 5293, 5295, 5296, 5336, 55024, 5530, 5532, 5533, 5534, 5535, 5579, 5594, 5595, 5604, 5605, 5777, 5879, 5880, 5881, 5894, 5970, 640, 6654, 6655, 6850, 695, 7409, 7410, 79168, 84433, 8503, 8517, 8519, 8915, 930, 933, 971, 973, 974, 975),
  symbol = c("AKT3", "LILRB2", "VAV3", "LILRB1", "MALT1", "LILRB5", "LILRB4", "LILRA1", "LILRB3", "LILRA3", "LILRA2", "CHUK", "PIK3AP1", "CR2", "AKT1", "AKT2", "FCGR2B", "FOS", "LILRA4", "RASGRP3", "DAPP1", "GRB2", "GSK3B", "BLNK", "HRAS", "LILRA5", "IKBKB", "INPP5D", "INPPL1", "JUN", "KRAS", "LYN", "NFATC1", "NFATC2", "NFATC3", "NFKB1", "NFKBIA", "NFKBIB", "NFKBIE", "NRAS", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2", "PLCG2", "BANK1", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "PRKCB", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "PTPN6", "RAC1", "RAC2", "RAC3", "RAF1", "RELA", "BLK", "SOS1", "SOS2", "SYK", "BTK", "VAV1", "VAV2", "LILRA6", "CARD11", "PIK3R3", "IKBKG", "IFITM1", "BCL10", "CD19", "CD22", "CD72", "CD79A", "CD79B", "CD81"),
  KO     = c("K04456", "K06512", "K05730", "K06512", "K07369", "K06512", "K06512", "K06512", "K06512", "K06512", "K06512", "K04467", "K12230", "K04012", "K04456", "K04456", "K12560", "K04379", "K06512", "K12362", "K12229", "K04364", "K03083", "K07371", "K02833", "K06512", "K07209", "K03084", "K15909", "K04448", "K07827", "K05854", "K04446", "K17332", "K17333", "K02580", "K04734", "K02581", "K05872", "K07828", "K00922", "K00922", "K00922", "K02649", "K02649", "K05859", "K26239", "K04348", "K04348", "K04348", "K06268", "K06268", "K19662", "K04371", "K04371", "K04368", "K04369", "K05697", "K04392", "K07860", "K07861", "K04366", "K04735", "K08890", "K03099", "K03099", "K05855", "K07370", "K05730", "K05730", "K06512", "K07367", "K02649", "K07210", "K19831", "K07368", "K06465", "K06467", "K06504", "K06506", "K06507", "K06508"),
  stringsAsFactors = FALSE)
stopifnot(identical(sort(unlist(nodes, use.names = FALSE)), sort(kegg_genes$symbol)), !anyDuplicated(kegg_genes$symbol))
