##agregate plots around CTCF binding sites from chip-seq data macs II

hicAggregateContacts --matrix P60v2.norm.20kb.0_1.cor.h5 --BED P60commonHNF4AcommonCEBP_CTCF.bed --outFileName P60_com3.pdf --numberOfBins 50 --operationType mean --transform obs/exp --mode intra-chr --dpi 300 --keep_outlier --range 20000:15000000 --vMin 0.8 --vMax 1.6 

hicAggregateContacts --matrix HNF4A.norm.20kb.0_1.cor.h5 --BED P60commonHNF4AcommonCEBP_CTCF.bed --outFileName HNF4A_com3.pdf --numberOfBins 50 --operationType mean --transform obs/exp --mode intra-chr --dpi 300 --keep_outlier --range 20000:15000000 --vMin 0.8 --vMax 1.6 

hicAggregateContacts --matrix CEBP.norm.20kb.0_1.cor.h5 --BED P60commonHNF4AcommonCEBP_CTCF.bed --outFileName CEBP_com3.pdf --numberOfBins 50 --operationType mean --transform obs/exp --mode intra-chr --dpi 300 --keep_outlier --range 20000:15000000 --vMin 0.8 --vMax 1.6 




hicAggregateContacts --matrix P60v2_COMMON.h5 --BED P60commonHNF4AcommonCEBP_CTCF.bed --outFileName P60_com4.pdf --numberOfBins 50 --operationType mean --transform obs/exp --mode intra-chr --dpi 300 --keep_outlier --range 20000:15000000 --vMin 0.8 --vMax 1.6 

hicAggregateContacts --matrix HNF4A_COMMON.h5 --BED P60commonHNF4AcommonCEBP_CTCF.bed --outFileName HNF4A_com4.pdf --numberOfBins 50 --operationType mean --transform obs/exp --mode intra-chr --dpi 300 --keep_outlier --range 20000:15000000 --vMin 0.8 --vMax 1.6 

hicAggregateContacts --matrix CEBP_COMMON.h5 --BED P60commonHNF4AcommonCEBP_CTCF.bed --outFileName CEBPA_com4.pdf --numberOfBins 50 --operationType mean --transform obs/exp --mode intra-chr --dpi 300 --keep_outlier --range 20000:15000000 --vMin 0.8 --vMax 1.6 

