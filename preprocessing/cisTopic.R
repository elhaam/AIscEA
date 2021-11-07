library(Matrix)
suppressWarnings(library(cisTopic, lib.loc ="~/.local/bin"))

atac_data <- read.csv(file = 'ATACseq.csv', row.names = 1, header= TRUE, sep = ',')


atac_mat <- data.matrix(atac_data)
RowName <- rownames(atac_data)
ColName <- colnames(atac_data)
Indrow = row(atac_data)[which(!atac_data == 0)]
Indcol = col(atac_data)[which(!atac_data == 0)]
Val = atac_mat[ which(!atac_mat == 0)]
atac_SpaMat <-sparseMatrix(
  i = Indrow, 
  j = Indcol, 
  x = Val,
  dims = c(nrow(atac_data), ncol(atac_data)), 
  dimnames = list(RowName, ColName)
)

cisTopicObject <- createcisTopicObject(atac_SpaMat, project.name='Mus')
rm(atac_SpaMat)
rm(atac_mat)

cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(2:15, 20, 25, 35, 40, 45, 50), seed=123, nCores=17, addModels=FALSE)

cisTopicObject <- selectModel(cisTopicObject, type='derivative')

#cisTopicObject <- runUmap(cisTopicObject, target='cell')
#plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', dim=2)
#ATAC_Topics <- cisTopicObject@selected.model$document_expects

topics <- cisTopicObject@selected.model$document_expects
write.table(topics,'atac_topics.tsv', sep = '\t', row.names = FALSE, quote=FALSE)

pred.matrix <- predictiveDistribution(cisTopicObject)
write.table(pred.matrix,'pred_atac.csv', sep = '\t', row.names = TRUE, quote=FALSE)
