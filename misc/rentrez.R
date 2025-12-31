library(rentrez)


hox_paper <- entrez_search(
  db = "pubmed",
  term = "10.1016/j.neuron.2024.01.013[doi]"
)

hox_paper$ids
hox_data <- entrez_link(db = "all", id = hox_paper$ids, dbfrom = "pubmed")
hox_data$links$pubmed_sra

hox_sra <- entrez_fetch(
  db = "sra",
  id = hox_data$links$pubmed_sra,
  rettype = "runinfo",
  retmode = "csv"
)


library(GEOquery)
library(rentrez)

gse_id <- "GSE226602"

# Retrieve the metadata
gse_metadata <- getGEO(gse_id)

# Print the metadata
print(gse_metadata)

get_sra_runs_for_geo


entrez_search(
  db
)


geo_search <- entrez_search(db = "gds", term = gse_id)
geo_search$ids
geo_summary <- entrez_summary(db = "gds", id = geo_search$ids)
geo_data <- entrez_fetch(
  db = "gds",
  id = geo_search$ids,
  rettype = "tsv"
)


cat(geo_data)


install.packages("rentrez")
library(rentrez)


gse_id <- "GSE226602"
geo_search <- entrez_search(db = "gds", term = gse_id)

geo_summary <- entrez_summary(db = "sra", id = geo_search$ids)


geo_data <- entrez_fetch(db = "gds", id = geo_search$ids, rettype = "xml")

# Load XML package to parse the data
library(XML)
geo_xml <- xmlParse(geo_data)

# Extract Study Title and BioProject ID
study_title <- xpathSApply(geo_xml, "//Title", xmlValue)
bioproject_id <- xpathSApply(geo_xml, "//Bioproject//Accession", xmlValue)
library(rentrez)

# 获取项目信息
project_id <- "PRJNA616380"
search_res <- entrez_search(db = "sra", term = project_id)
summary_res <- entrez_summary(db = "sra", id = search_res$ids)

# 获取相关的 Run 信息
sra_xml <- entrez_fetch(
  db = "sra",
  id = search_res$ids,
  rettype = "xml",
  parsed = TRUE
)

# 打印摘要
summary_res
# 查看 summary_res 的结构
str(summary_res)

# 提取某些关键字段
project_title <- summary_res$title
study_accession <- summary_res$study_acc
submission_date <- summary_res$submission_date
organism <- summary_res$organism
experiments <- summary_res$experiment

# 打印提取的信息
print(project_title)
print(study_accession)
print(submission_date)
print(organism)
print(experiments)
