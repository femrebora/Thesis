# PID of current job: 3818662
mSet<-InitDataObjects("conc", "msetqea", FALSE)
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "colu", "cont");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "(-)-Gallocatechin-3-O-gallate");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "(2R)-2-Hydroxy-3-butenyl glucosinolate");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "(9Z, 12Z)-Octadecadienoate");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "1,2-Disinapoyl diglucoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "102FTA");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "10-Deacetyl-7-xylosylpaclitaxel");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "13-Methylberberine");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "1-O-p-Coumaroylglycerol");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "2',4,5'-Trihydroxychalcone");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "2-[(4,6-diamino-1,3,5-triazin-2-yl)thio]-Ethanesulfonic acid");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "4-Methylsulfinyl-3-butenyl glucosinolate");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "62FTA");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "73FTA");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Acetylendicarboxylate");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Acetylendicarboxylate", "Acetylenedicarboxylate");
mSet<-PerformDetailMatch(mSet, "Apigenin 7-O-glucuronide");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Apiopaeonoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Baccatin I");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Bestatin");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Brevifoliol");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Catechin-(4alpha->6)-gallocatechin");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Chloratranol");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Cyanidin 3-(6''-acetyl)galactoside");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Cyanidin 3-(6''-acetyl)galactoside", "Cyanidin 3-(6''-acetyl-galactoside)");
mSet<-PerformDetailMatch(mSet, "Cyanidin 3-(6''-p-coumarylsophoroside)-5-glucoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Cyclo(Pro-Phe)");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Dimethyl lithospermate");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "discorhabdin G");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Galacturonate 1-phosphate");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Galacturonate 1-phosphate", "1-Phospho-alpha-D-galacturonate");
mSet<-PerformDetailMatch(mSet, "Gangaleoidin");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Geyerline");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "GlcNAcThrNAc");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Gluconapoleiferin");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Glutathione (oxidised form)");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Glycyrrhizate");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Hedysarimcoumestan A");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Herbacetin 8-O-xylose-3-O-glucose");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Kutkoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Microcolin B");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Naringenin-7-O-neohesperidoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "nPFOA");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Octopine");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Pinoresinol-4-O-glucoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Posthumulone");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Quercetin-3,4'-O-di-beta-glucopyranoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Quercetin-3-O-beta-galactoside");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Quercetin-3-O-beta-galactoside", "Quercetin 3-galactoside");
mSet<-PerformDetailMatch(mSet, "Quercetin-3,4'-O-di-beta-glucopyranoside");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Scopularide D");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Taxinine A");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Taxinine M");
mSet<-GetCandidateList(mSet);
mSet<-PerformDetailMatch(mSet, "Zearalanone");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Zearalanone", "(S,E)-Zearalenone");
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "CrNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_1_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_1_", "png", 72, width=NA)
mSet<-Normalization(mSet, "MedianNorm", "SrNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_2_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_2_", "png", 72, width=NA)
mSet<-Normalization(mSet, "MedianNorm", "CrNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_3_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_3_", "png", 72, width=NA)
mSet<-Normalization(mSet, "MedianNorm", "CrNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_4_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_4_", "png", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_0_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_0_", "png", 72, width=NA)
mSet<-PlotQEA.MetSet(mSet, "Glycine and Serine Metabolism", "png", 72, width=NA)
mSet<-PlotQEA.MetSet(mSet, "Glycine and Serine Metabolism", "png", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_1_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_1_", "png", 72, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "RaMP_pathway", 2);
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_2_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_2_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)
