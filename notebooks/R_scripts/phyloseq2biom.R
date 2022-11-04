##

source("src/_load_packages.R")
base::load("data/Phyloseq_Objects/Phyloseq_all.RData")
source(paste0(fs::path_home(), "/git/microfiltR/microfiltR_source_code.R"))


for (ref_DB in names(refDB_phyloseq)) {
  for (level in names(refDB_phyloseq[[ref_DB]])) {
    e3f.f <- write.dataset(ps = dat.obj,
                           filePATH = fp,
                           filePREFIX = "e3.f",
                           writeFASTA=F, rename = F)

  }
}

dat.obj <- refDB_phyloseq[["UHGG"]][["Species"]]


e3f.f <- write.dataset(ps = dat.obj,
                       filePATH = "data/biom",
                       filePREFIX = "testpref_",
                       writeFASTA=F, rename = F)
