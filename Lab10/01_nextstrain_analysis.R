library(dplyr)
library(data.table)
library(reshape2)

df <- fread("nextstrain_ncov_metadata_old_lab.tsv") %>% data.frame()

# Get some univariate statistics
table(df$Clade)
table(df$Sex)
table(df$Country)

# See how many of each clade is in each country
summary_df <- df %>% group_by(Clade, Country) %>% summarize(count = n()) 

# Determine the country that has the most cases per clade
df %>% group_by(Clade, Country) %>% summarize(count = n()) %>%
  ungroup() %>% group_by(Clade) %>% top_n(n = 1, wt = count)

# Do a test of association
count_reshape <- reshape2::dcast(summary_df, Country ~ Clade, value.var = "count", fill = 0)
count_mat <- data.matrix(count_reshape[,c(-1)])
rownames(count_mat) <- count_reshape$Country
count_mat[1:5,1:5]

# Perform association test
chisq.test(count_mat)
dim(count_mat)
