#Solution copied and adapted from following website:
#http://biostar.stackexchange.com/questions/15016/get-gene-symbol-synonyms-aliases-with-biomart

#INSTALL
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")

#To cite this package in a publication, start R and enter:
citation("org.Hs.eg.db")

library(org.Hs.eg.db)

# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()

# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
#give a list with aliases
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)


# load genenames you want to find aliases for from file
queryGeneNames <- as.character(read.table("89GenesWithPossibleAliases.txt",header=T)$Gene)



#Get all aliases (if any) for all genenames in a list
allAliases <- list()
for (i in 1:length(queryGeneNames))
{
  #from aliasSymbol obtained by SQL query find for genename at hand all aliases
  #this assumes that upper case and lower case doesn't change meaning
  #can pick up extra aliases due to inconsistent casing
  aliases <-aliasSymbol[which(tolower(aliasSymbol[,2]) == tolower(queryGeneNames[i])),5] 
  if(length(aliases)==1){print(paste("Warning: No alias found for gene name",queryGeneNames[i]))}
  allAliases[[i]] <- aliases
}
allAliases

#Now of course we want to know for each genename in the list which gene
#it is in h18list
h18List <- read.table("GeneList_hg18.txt",header=F)

genePairs <- matrix(,ncol=4,nrow=length(queryGeneNames))
colnames(genePairs) <- c("queryGeneId","queryGeneName","refGeneId","refGeneName")

for (i in 1:length(queryGeneNames))
{
  genePairs[i,1:2] <- c(i,queryGeneNames[i])
  #all known names for the current gene
  aliases <- c(queryGeneNames[i],allAliases[[i]])
  #or alternatively only look for alias matches not the original query name
  #aliases <- allAliases[[i]]
  
  for (j in 1:length(aliases))
  {
    matchid <- which(h18List[,4]==aliases[j])
    if(length(matchid)>0)  genePairs[i,3:4] <- c(matchid,h18List[matchid,4])
  }
  
}

genePairs

