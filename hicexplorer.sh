hicSumMatrices --matrices AS-937495.h5 AS-937497.h5 AS-894137.h5 --outFileName E14.h5
hicSumMatrices --matrices AS-937499.h5 AS-937501.h5 AS-894139.h5 --outFileName E18.h5
hicSumMatrices --matrices AS-937523.h5 AS-937503.h5 AS-894116.h5 --outFileName P0.h5
hicSumMatrices --matrices AS-937505.h5 AS-937525.h5 AS-894118.h5 --outFileName P60.h5
hicSumMatrices --matrices AS-937527.h5 AS-937529.h5  --outFileName HNF4A.h5
hicSumMatrices --matrices AS-937531.h5 AS-937533.h5 --outFileName CEBPA.h5

hicMergeMatrixBins --matrix E14.h5 --numBins 2 --outFileName E14.20kb.h5
hicMergeMatrixBins --matrix E18.h5 --numBins 2 --outFileName E18.20kb.h5
hicMergeMatrixBins --matrix P0.h5 --numBins 2 --outFileName P0.20kb.h5
hicMergeMatrixBins --matrix P60.h5 --numBins 2 --outFileName P60.20kb.h5
hicMergeMatrixBins --matrix CEBPA.h5 --numBins 2 --outFileName CEBP.20kb.h5
hicMergeMatrixBins --matrix HNF4A.h5 --numBins 2 --outFileName HNF4A.20kb.h5

hicNormalize -m E14.h5  E18.h5  P0.h5 P60.h5 --normalize smallest -o E14.10kb.norm.h5 E18.10kb.norm.h5 P0.10kb.norm.h5 P60.10kb.norm.h5
hicNormalize -m P60.h5  CEBP.h5  --normalize smallest -o P60v2.10kb.norm.h5 CEBP.10kb.norm.h5
hicNormalize -m P60.h5  HNF4A.h5   --normalize smallest -o P60v3.10kb.norm.h5 HNF4A.10kb.norm.h5

hicConvertFormat -m E14.10kb.norm.h5 -o E14.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m E14.10kb.cool -o E14loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8

hicConvertFormat -m E18.10kb.norm.h5 -o E18.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m E18.10kb.cool -o E18loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8

hicConvertFormat -m P0.10kb.norm.h5 -o P0.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m P0.10kb.cool -o P0loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8

hicConvertFormat -m P60.10kb.norm.h5 -o P60.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m P60.10kb.cool -o P60loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8

# Diagnostic plot
hicCorrectMatrix diagnostic_plot -m E14.norm.20kb.0_1.h5 -o diagnosticPlot_E14.20kb.pdf --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
# Correct Matrix (-1.5 - 2 for 250kb) (-1.5 - 3 for 10kb)
hicCorrectMatrix correct -m E14.norm.20kb.0_1.h5 --filterThreshold -1.5 2 -o E14.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
hicCorrectMatrix correct -m E18.norm.20kb.0_1.h5 --filterThreshold -1.5 2 -o E18.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
hicCorrectMatrix correct -m P0.norm.20kb.0_1.h5 --filterThreshold -1.5 2 -o P0.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
hicCorrectMatrix correct -m P60.norm.20kb.0_1.h5 --filterThreshold -1.5 2 -o P60.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
hicCorrectMatrix correct -m CEBP.norm.20kb.0_1.h5 --filterThreshold -1.5 2 -o CEBP.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
hicCorrectMatrix correct -m HNF4A.norm.20kb.0_1.h5 --filterThreshold -1.5 2 -o HNF4A.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY

hicAdjustMatrix --matrix E14.norm.20kb.0_1.cor.h5 -r centromeres_mm10.bed -o E14.norm.cor.adj.20kb.h5 -a mask
hicAdjustMatrix --matrix E18.norm.20kb.0_1.cor.h5 -r centromeres_mm10.bed -o E18.norm.cor.adj.20kb.h5 -a mask
hicAdjustMatrix --matrix P0.norm.20kb.0_1.cor.h5 -r centromeres_mm10.bed -o P0.norm.cor.adj.20kb.h5 -a mask
hicAdjustMatrix --matrix P60.norm.20kb.0_1.cor.h5 -r centromeres_mm10.bed -o P60.norm.cor.adj.20kb.h5 -a mask
hicAdjustMatrix --matrix P60v2.norm.20kb.0_1.cor.h5 -r centromeres_mm10.bed -o P60v2.norm.cor.adj.20kb.h5 -a mask
hicAdjustMatrix --matrix CEBP.norm.20kb.0_1.cor.h5 -r centromeres_mm10.bed -o CEBP.norm.cor.adj.20kb.h5 -a mask
hicAdjustMatrix --matrix HNF4A.norm.20kb.0_1.cor.h5 -r centromeres_mm10.bed -o HNF4A.norm.cor.adj.20kb.h5 -a mask

hicFindTADs -m E14.norm.20kb.0_1.cor.h5 --outPrefix E14.20kbv2 --numberOfProcessors 8 --correctForMultipleTesting fdr --threshold 0.001
hicFindTADs -m E18.norm.20kb.0_1.cor.h5 --outPrefix E18.20kbv2 --numberOfProcessors 8 --correctForMultipleTesting fdr --threshold 0.001
hicFindTADs -m P0.norm.20kb.0_1.cor.h5 --outPrefix P0.20kbv2 --numberOfProcessors 8 --correctForMultipleTesting fdr --threshold 0.001
hicFindTADs -m P60.norm.20kb.0_1.cor.h5 --outPrefix P60.20kbv2 --numberOfProcessors 8 --correctForMultipleTesting fdr --threshold 0.001
hicFindTADs -m CEBP.norm.20kb.0_1.cor.h5 --outPrefix CEBP.20kbv2 --numberOfProcessors 8 --correctForMultipleTesting fdr --threshold 0.001
hicFindTADs -m HNF4A.norm.20kb.0_1.cor.h5 --outPrefix HNF4A.20kbv2 --numberOfProcessors 8 --correctForMultipleTesting fdr --threshold 0.001

hicCompareMatrices --matrices E18.norm.20kb.0_1.cor.h5 E14.norm.20kb.0_1.cor.h5 --outFileName E14_E18.h5 --operation log2ratio
hicCompareMatrices --matrices P0.norm.20kb.0_1.cor.h5 E18.norm.20kb.0_1.cor.h5 --outFileName E18_P0.h5 --operation log2ratio
hicCompareMatrices --matrices P60.norm.20kb.0_1.cor.h5 P0.norm.20kb.0_1.cor.h5 --outFileName P0_P60.h5 --operation log2ratio
hicCompareMatrices --matrices P60.norm.20kb.0_1.cor.h5 E18.norm.20kb.0_1.cor.h5 --outFileName E18_P60.h5 --operation log2ratio
hicCompareMatrices --matrices CEBP.norm.20kb.0_1.cor.h5 P60.norm.20kb.0_1.cor.h5 --outFileName P60_CEBP.h5 --operation log2ratio
hicCompareMatrices --matrices HNF4A.norm.20kb.0_1.cor.h5 P60.norm.20kb.0_1.cor.h5 --outFileName P60_HNF4A.h5 --operation log2ratio

hicPlotMatrix -m E14_E18.h5 --clearMaskedBins --region chr8:105000000-115000000 -o E14_E18_log2ratiocontactm_chr8.pdf --colorMap RdBu_r  --title "E14-E18 log2ratio" --vMin -4 --vMax 4
hicPlotMatrix -m E14_E18.h5 --clearMaskedBins --region chr11:96000000-106000000 -o E14_E18_log2ratiocontactm_chr11.pdf --colorMap RdBu_r  --title "E14-E18 log2ratio" --vMin -4 --vMax 4
hicPlotMatrix -m E18_P0.h5 --clearMaskedBins --region chr8:105000000-115000000 -o E18_P0_log2ratiocontactm_chr8.pdf --colorMap RdBu_r  --title "E18-P0 log2ratio" --vMin -4 --vMax 4
hicPlotMatrix -m E18_P0.h5 --clearMaskedBins --region chr11:96000000-106000000 -o E18_P0_log2ratiocontactm_chr11.pdf --colorMap RdBu_r  --title "E18-P0 log2ratio" --vMin -4 --vMax 4
hicPlotMatrix -m E18_P60.h5 --clearMaskedBins --region chr8:105000000-115000000 -o E18_P60_log2ratiocontactm_chr8.pdf --colorMap RdBu_r  --title "E18-P60 log2ratio" --vMin -4 --vMax 4
hicPlotMatrix -m E18_P60.h5 --clearMaskedBins --region chr11:96000000-106000000 -o E18_P60_log2ratiocontactm_chr11.pdf --colorMap RdBu_r  --title "E18-P60 log2ratio" --vMin -4 --vMax 4
hicPlotMatrix -m P0_P60.h5 --clearMaskedBins --region chr8:105000000-115000000 -o P0_P60_log2ratiocontactm_chr8.pdf --colorMap RdBu_r  --title "P0-P60 log2ratio" --vMin -4 --vMax 4
hicPlotMatrix -m P0_P60.h5 --clearMaskedBins --region chr11:96000000-106000000 -o P0_P60_log2ratiocontactm_chr11.pdf --colorMap RdBu_r  --title "P0-P60 log2ratio" --vMin -4 --vMax 4

hicConvertFormat -m E14.20kb.h5 -o E14.20kb.cool --inputFormat h5 --outputFormat cool
hicConvertFormat -m E14.h5 -o E14.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m E14.10kb.cool -o E14loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E14.10kb.cool -o E14loops2mb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E14.20kb.cool -o E14loops1mb20kb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E14.20kb.cool -o E14loops2mb20kb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E14.20kb.cool -o E14.20kbdefault.loop -t 8

hicConvertFormat -m E18.20kb.h5 -o E18.20kb.cool --inputFormat h5 --outputFormat cool
hicConvertFormat -m E18.h5 -o E18.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m E18.10kb.cool -o E18loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E18.10kb.cool -o E18loops2mb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E18.20kb.cool -o E18loops1mb20kb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E18.20kb.cool -o E18loops2mb20kb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m E18.20kb.cool -o E18.20kbdefault.loop -t 8

hicConvertFormat -m P0.20kb.h5 -o P0.20kb.cool --inputFormat h5 --outputFormat cool
hicConvertFormat -m P0.h5 -o P0.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m P0.10kb.cool -o P0loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P0.10kb.cool -o P0loops2mb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P0.20kb.cool -o P0loops1mb20kb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P0.20kb.cool -o P0loops2mb20kb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P0.20kb.cool -o P0.20kbdefault.loop -t 8

hicConvertFormat -m P60.20kb.h5 -o P60.20kb.cool --inputFormat h5 --outputFormat cool
hicConvertFormat -m P60.h5 -o P60.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m P60.10kb.cool -o P60loops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P60.10kb.cool -o P60loops2mb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P60.20kb.cool -o P60loops1mb20kb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P60.20kb.cool -o P60loops2mb20kb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m P60.20kb.cool -o P60.20kbdefault.loop -t 8

hicConvertFormat -m CEBP.20kb.h5 -o CEBP.20kb.cool --inputFormat h5 --outputFormat cool
hicConvertFormat -m CEBPA.h5 -o CEBP.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m CEBP.10kb.cool -o CEBPloops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m CEBP.10kb.cool -o CEBPloops2mb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m CEBP.20kb.cool -o CEBPloops1mb20kb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m CEBP.20kb.cool -o CEBPloops2mb20kb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m CEBP.20kb.cool -o CEBP.20kbdefault.loop -t 8

hicConvertFormat -m HNF4A.20kb.h5 -o HNF4A.20kb.cool --inputFormat h5 --outputFormat cool
hicConvertFormat -m HNF4A.h5 -o HNF4A.10kb.cool --inputFormat h5 --outputFormat cool
hicDetectLoops -m HNF4A.10kb.cool -o HNF4Aloops1mb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m HNF4A.10kb.cool -o HNF4Aloops2mb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m HNF4A.20kb.cool -o HNF4Aloops1mb20kb.bedgraph --maxLoopDistance 1000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m HNF4A.20kb.cool -o HNF4Aloops2mb20kb.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 -t 8
hicDetectLoops -m HNF4A.20kb.cool -o HNF4A.20kbdefault.loop -t 8

hicPCA -m E14.norm.20kb.0_1.cor.h5 -o E14_pca1chr8.bw E14_pca2chr8.bw E14_pca3chr8.bw E14_pca4chr8.bw E14_pca5chr8.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr8
hicPlotMatrix -m E14.norm.20kb.0_1.cor.h5 --outFileName E14_20kb_chr8.pdf --region chr8:107000000-114500000 --bigwig E14_pca1chr8.bw --log --colorMap afmhot_r --title "E14 20kb" --loops E14loops1mb.bedgraph
hicPCA -m E14.norm.20kb.0_1.cor.h5 -o E14_pca1chr11.bw E14_pca2chr11.bw E14_pca3chr11.bw E14_pca4chr11.bw E14_pca5chr11.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr11
hicPlotMatrix -m E14.norm.20kb.0_1.cor.h5 --outFileName E14_20kb_chr11.pdf --region chr11:97000000-102000000 --bigwig E14_pca1chr11.bw --log --colorMap afmhot_r --title "E14 20kb" --loops E14loops1mb.bedgraph

hicPCA -m E18.norm.20kb.0_1.cor.h5 -o E18_pca1chr8.bw E18_pca2chr8.bw E18_pca3chr8.bw E18_pca4chr8.bw E18_pca5chr8.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr8
hicPlotMatrix -m E18.norm.20kb.0_1.cor.h5 --outFileName E18_20kb_chr8.pdf --region chr8:107000000-114500000 --bigwig E18_pca1chr8.bw --log --colorMap afmhot_r --title "E18 20kb" --loops E18loops1mb.bedgraph
hicPCA -m E18.norm.20kb.0_1.cor.h5 -o E18_pca1chr11.bw E18_pca2chr11.bw E18_pca3chr11.bw E18_pca4chr11.bw E18_pca5chr11.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr11
hicPlotMatrix -m E18.norm.cor.adj.20kb.h5 --outFileName E18_20kb_chr11.pdf --region chr11:97000000-102000000 --bigwig E18_pca1chr11.bw --log --colorMap afmhot_r --title "E18 20kb" --loops E18loops1mb.bedgraph

hicPCA -m P0.norm.20kb.0_1.cor.h5 -o P0_pca1chr8.bw P0_pca2chr8.bw P0_pca3chr8.bw P0_pca4chr8.bw P0_pca5chr8.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr8
hicPlotMatrix -m P0.norm.20kb.0_1.cor.h5 --outFileName P0_20kb_chr8.pdf --region chr8:107000000-114500000 --bigwig P0_pca1chr8.bw --log --colorMap afmhot_r --title "P0 20kb" --loops P0loops1mb.bedgraph
hicPCA -m P0.norm.20kb.0_1.cor.h5 -o P0_pca1chr11.bw P0_pca2chr11.bw P0_pca3chr11.bw P0_pca4chr11.bw P0_pca5chr11.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr11
hicPlotMatrix -m P0.norm.cor.adj.20kb.h5 --outFileName P0_20kb_chr11.pdf --region chr11:97000000-102000000 --bigwig P0_pca1chr11.bw --log --colorMap afmhot_r --title "P0 20kb" --loops P0loops1mb.bedgraph


hicPCA -m P60.norm.20kb.0_1.cor.h5 -o P60_pca1chr8.bw P60_pca2chr8.bw P60_pca3chr8.bw P60_pca4chr8.bw P60_pca5chr8.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr8
hicPlotMatrix -m P60.norm.20kb.0_1.cor.h5 --outFileName P60_20kb_chr8.pdf --region chr8:107000000-114500000 --bigwig P60_pca1chr8.bw --log --colorMap afmhot_r --title "P60 20kb" --loops P60loops1mb.bedgraph
hicPCA -m P60.norm.20kb.0_1.cor.h5 -o P60_pca1chr11.bw P60_pca2chr11.bw P60_pca3chr11.bw P60_pca4chr11.bw P60_pca5chr11.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr11
hicPlotMatrix -m P60.norm.cor.adj.20kb.h5 --outFileName P60_20kb_chr11.pdf --region chr11:97000000-102000000 --bigwig P60_pca1chr11.bw --log --colorMap afmhot_r --title "P60 20kb" --loops P60loops1mb.bedgraph

hicPCA -m CEBP.norm.20kb.0_1.cor.h5 -o CEBP_pca1chr8.bw CEBP_pca2chr8.bw CEBP_pca3chr8.bw CEBP_pca4chr8.bw CEBP_pca5chr8.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr8
hicPlotMatrix -m CEBP.norm.20kb.0_1.cor.h5 --outFileName CEBP_20kb_chr8.pdf --region chr8:107000000-114500000 --bigwig CEBP_pca1chr8.bw --log --colorMap afmhot_r --title "CEBP 20kb" --loops CEBPloops1mb.bedgraph
hicPCA -m CEBP.norm.cor.adj.20kb.h5 -o CEBP_pca1chr11.bw CEBP_pca2chr11.bw CEBP_pca3chr11.bw CEBP_pca4chr11.bw CEBP_pca5chr11.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr11
hicPlotMatrix -m CEBP.norm.20kb.0_1.cor.h5 --outFileName CEBP_20kb_chr11.pdf --region chr11:97000000-102000000 --bigwig CEBP_pca1chr11.bw --log --colorMap afmhot_r --title "CEBP 20kb" --loops CEBPloops1mb.bedgraph


hicPCA -m HNF4A.norm.20kb.0_1.cor.h5 -o HNF4A_pca1chr8.bw HNF4A_pca2chr8.bw HNF4A_pca3chr8.bw HNF4A_pca4chr8.bw HNF4A_pca5chr8.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr8
hicPlotMatrix -m HNF4A.norm.cor.adj.20kb.h5 --outFileName HNF4A_20kb_chr8.pdf --region chr8:107000000-114500000 --bigwig HNF4A_pca1chr8.bw --log --colorMap afmhot_r --title "HNF4A 20kb" --loops HNF4Aloops1mb.bedgraph
hicPCA -m HNF4A.norm.20kb.0_1.cor.h5 -o HNF4A_pca1chr11.bw HNF4A_pca2chr11.bw HNF4A_pca3chr11.bw HNF4A_pca4chr11.bw HNF4A_pca5chr11.bw --format bigwig --whichEigenvectors 1 2 3 4 5 --chromosomes chr11
hicPlotMatrix -m HNF4A.norm.20kb.0_1.cor.h5 --outFileName HNF4A_20kb_chr11.pdf --region chr11:97000000-102000000 --bigwig HNF4A_pca1chr11.bw --log --colorMap afmhot_r --title "HNF4A 20kb" --loops HNF4Aloops1mb.bedgraph


hicPlotTADs --tracks track_TADS_P0testv1.ini.txt --region chr8:107000000-114500000 --outFileName TADS_E14v2chr8v2.pdf
hicPlotTADs --tracks track_TADS_P0testv1.ini.txt --region chr11:97000000-102000000 --outFileName TADS_E14v2chr1v2.pdf

hicPlotTADs --tracks track_TADS_P0testv2.ini.txt --region chr8:107000000-114500000 --outFileName TADS_E18v2chr8.pdf
hicPlotTADs --tracks track_TADS_P0testv2.ini.txt --region chr11:97000000-102000000 --outFileName TADS_E18v2chr11.pdf

hicPlotTADs --tracks track_TADS_P0testv3.ini.txt --region chr8:107000000-114500000 --outFileName TADS_P0v2chr8.pdf
hicPlotTADs --tracks track_TADS_P0testv3.ini.txt --region chr11:97000000-102000000 --outFileName TADS_P0v2chr11.pdf

hicPlotTADs --tracks track_TADS_P0testv4.ini.txt --region chr8:107000000-114500000 --outFileName TADS_P60v2chr8.pdf
hicPlotTADs --tracks track_TADS_P0testv4.ini.txt --region chr11:97000000-102000000 --outFileName TADS_P60v2chr11.pdf

hicPlotTADs --tracks track_TADS_P0testv5.ini.txt --region chr8:107000000-114500000 --outFileName TADS_CEBPv2chr8.pdf
hicPlotTADs --tracks track_TADS_P0testv5.ini.txt --region chr11:97000000-102000000 --outFileName TADS_CEBPv2chr11.pdf

hicPlotTADs --tracks track_TADS_P0testv6.ini.txt --region chr8:107000000-114500000 --outFileName TADS_HNF4Av2chr8.pdf
hicPlotTADs --tracks track_TADS_P0testv6.ini.txt --region chr11:97000000-102000000 --outFileName TADS_HNF4Av2chr11.pdf

hicPlotTADs --tracks track_log2ratio_Stages.ini --region chr8:107000000-114500000 --outFileName TADS_stageslog2ratioschr8.pdf
hicPlotTADs --tracks track_log2ratio_Stages.ini --region chr11:97000000-102000000 --outFileName TADS_stageslog2ratioschr11.pdf
hicPlotTADs --tracks track_insulation_Heatmap_Stages.ini --region chr8:107000000-114500000 --outFileName TADS_stageslog2ratioschr8Heatmap.pdf
hicPlotTADs --tracks track_insulation_Heatmap_Stages.ini --region chr11:97000000-102000000 --outFileName TADS_stageslog2ratioschr11Heatmap.pdf

hicPlotTADs --tracks track_log2ratio_Stages.ini --region chr8:107000000-114500000 --outFileName TADS_stageslog2ratioschr8.pdf
hicPlotTADs --tracks track_log2ratio_Stages.ini --region chr11:97000000-102000000 --outFileName TADS_stageslog2ratioschr11.pdf
hicPlotTADs --tracks track_insulation_Heatmap_Stages.iniversion3.txt --region chr8:107000000-114500000 --outFileName TADS_stageslog2ratioschr8HeatmapE14vsE18test2.pdf
hicPlotTADs --tracks track_insulation_Heatmap_Stages.iniversion3.txt --region chr11:97000000-102000000 --outFileName TADS_stageslog2ratioschr11HeatmapE14vsE18test2.pdf

track_insulation_Heatmap_Stages

hicPlotTADs --tracks track_log2ratio_Stagesv2.ini --region chr8:107000000-114500000 --outFileName TADS_stageslog2ratioschr8v2.pdf
hicPlotTADs --tracks track_log2ratio_Stagesv2.ini --region chr11:97000000-102000000 --outFileName TADS_stageslog2ratioschr11v2.pdf
hicPlotTADs --tracks track_insulation_Heatmap_Stages.iniv2.txt --region chr8:107000000-114500000 --outFileName TADS_stageslog2ratioschr8Heatmapv2.pdf
hicPlotTADs --tracks track_insulation_Heatmap_Stages.iniv2.txt --region chr11:97000000-102000000 --outFileName TADS_stageslog2ratioschr11Heatmapv2.pdf

####################################dcHiC

hicAdjustMatrix -m E14.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chrX chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 --action keep -o E14.norm.20kb.0_1.cor_chr.h5
hicConvertFormat -m E14.norm.20kb.0_1.cor_chr.h5 --inputFormat h5 --outputFormat cool -o E14.norm.20kb.0_1.cor_chr.cool
python /home/orsalia/Downloads/dcHiC-master/utility/preprocess.py -input cool -file  E14.norm.20kb.0_1.cor_chr.cool -genomeFile mm10.chrom.sizes2.txt -res 20000 -prefix E14

hicAdjustMatrix -m E18.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chrX chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 --action keep -o E18.norm.20kb.0_1.cor_chr.h5
hicConvertFormat -m E18.norm.20kb.0_1.cor_chr.h5 --inputFormat h5 --outputFormat cool -o E18.norm.20kb.0_1.cor_chr.cool
python /home/orsalia/Downloads/dcHiC-master/utility/preprocess.py -input cool -file  E18.norm.20kb.0_1.cor_chr.cool -genomeFile mm10.chrom.sizes2.txt -res 20000 -prefix E18

hicAdjustMatrix -m P0.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chrX chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 --action keep -o P0.norm.20kb.0_1.cor_chr.h5
hicConvertFormat -m P0.norm.20kb.0_1.cor_chr.h5 --inputFormat h5 --outputFormat cool -o P0.norm.20kb.0_1.cor_chr.cool
python /home/orsalia/Downloads/dcHiC-master/utility/preprocess.py -input cool -file  P0.norm.20kb.0_1.cor_chr.cool -genomeFile mm10.chrom.sizes2.txt -res 20000 -prefix P0

hicAdjustMatrix -m P60.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chrX chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 --action keep -o P60.norm.20kb.0_1.cor_chr.h5
hicConvertFormat -m P60.norm.20kb.0_1.cor_chr.h5 --inputFormat h5 --outputFormat cool -o P60.norm.20kb.0_1.cor_chr.cool
python /home/orsalia/Downloads/dcHiC-master/utility/preprocess.py -input cool -file  P60.norm.20kb.0_1.cor_chr.cool -genomeFile mm10.chrom.sizes2.txt -res 20000 -prefix P60

hicAdjustMatrix -m CEBP.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chrX chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 --action keep -o CEBP.norm.20kb.0_1.cor_chr.h5
hicConvertFormat -m CEBP.norm.20kb.0_1.cor_chr.h5 --inputFormat h5 --outputFormat cool -o CEBP.norm.20kb.0_1.cor_chr.cool
python /home/orsalia/Downloads/dcHiC-master/utility/preprocess.py -input cool -file  CEBP.norm.20kb.0_1.cor_chr.cool -genomeFile mm10.chrom.sizes2.txt -res 20000 -prefix CEBP

hicAdjustMatrix -m HNF4A.norm.20kb.0_1.cor.h5 --chromosomes chr1 chr2 chr3 chrX chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 --action keep -o HNF4A.norm.20kb.0_1.cor_chr.h5
hicConvertFormat -m HNF4A.norm.20kb.0_1.cor_chr.h5 --inputFormat h5 --outputFormat cool -o HNF4A.norm.20kb.0_1.cor_chr.cool
python /home/orsalia/Downloads/dcHiC-master/utility/preprocess.py -input cool -file  HNF4A.norm.20kb.0_1.cor_chr.cool -genomeFile mm10.chrom.sizes2.txt -res 20000 -prefix HNF4A

###############################################################################################################dcHiC
cd /home/orsalia/Desktop/HiC_PHiC/P60/HiC/again/again/plots/downsample/matrices/E14vsE18/
conda activate dchic
Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E14_E18.txt --pcatype cis --dirovwt T

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E14_E18.txt --pcatype select --dirovwt T --genome mm10

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E14_E18.txt --pcatype analyze --dirovwt T --diffdir E18_vs_E14

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E14_E18.txt --pcatype subcomp --dirovwt T --diffdir E18_vs_E14

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E14_E18.txt --pcatype viz --diffdir E18_vs_E14 --genome mm10

cd /home/orsalia/Desktop/HiC_PHiC/P60/HiC/again/again/plots/downsample/matrices/E18vsP0/
conda activate dchic
Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P0.txt --pcatype cis --dirovwt T

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P0.txt --pcatype select --dirovwt T --genome mm10

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P0.txt --pcatype analyze --dirovwt T --diffdir E18_vs_P0

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P0.txt --pcatype subcomp --dirovwt T --diffdir E18_vs_P0

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P0.txt --pcatype viz --diffdir E18_vs_P0 --genome mm10


cd /home/orsalia/Desktop/HiC_PHiC/P60/HiC/again/again/plots/downsample/matrices/P0vsP60/
conda activate dchic
Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P0_P60.txt --pcatype cis --dirovwt T

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P0_P60.txt --pcatype select --dirovwt T --genome mm10

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P0_P60.txt --pcatype analyze --dirovwt T --diffdir P0_vs_P60

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P0_P60.txt --pcatype subcomp --dirovwt T --diffdir P0_vs_P60

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P0_P60.txt --pcatype viz --diffdir P0_vs_P60 --genome mm10

cd /home/orsalia/Desktop/HiC_PHiC/P60/HiC/again/again/plots/downsample/matrices/E18vsP60/
conda activate dchic
Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P60.txt --pcatype cis --dirovwt T

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P60.txt --pcatype select --dirovwt T --genome mm10

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P60.txt --pcatype analyze --dirovwt T --diffdir E18_vs_P60

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P60.txt --pcatype subcomp --dirovwt T --diffdir E18_vs_P60

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.E18_P60.txt --pcatype viz --diffdir E18_vs_P60 --genome mm10


cd /home/orsalia/Desktop/HiC_PHiC/P60/HiC/again/again/plots/downsample/matrices/P60vsCEBP
conda activate dchic
Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P60_CEBP.txt --pcatype cis --dirovwt T

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P60_CEBP.txt --pcatype select --dirovwt T --genome mm10

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P60_CEBP.txt --pcatype analyze --dirovwt T --diffdir P60_vs_CEBP

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P60_CEBP.txt --pcatype subcomp --dirovwt T --diffdir P60_vs_CEBP

Rscript /home/orsalia/Downloads/dcHiC-master/dchicf.r --file input.P60_CEBP.txt --pcatype viz --diffdir P60_vs_CEBP --genome mm10
