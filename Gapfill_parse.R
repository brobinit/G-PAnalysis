# Gapfill Parsing
# Blaire Robinson
# Created 1 May 2016; Updated 1 May 2016

#load data
frame <- read.table("gapfilled_roles_ED252.tsv", sep="\t", header=TRUE, quote="")

#Split media lists and create new data frame
s<-strsplit(as.character(frame$media),";")
medias <- data.frame(media=unlist(s),
                     reaction_IDs=rep(frame$reaction_id, sapply(s, FUN=length)), 
                     roles=rep(frame$roles, sapply(s, FUN=length)),
                     equation=rep(frame$equation, sapply(s, FUN=length)))

# Conbining like media names to get concatenated list of media~reaction pairs
# wide_format <- unstack(medias, reaction_IDs~media) #simple method for two columns
# wide_format
# wide_format1 <- data.frame(
#   Media=rep(names(wide_format), lapply(wide_format, length)), 
#   Reactions = unlist(wide_format))
# rownames(wide_format1) <-NULL
# attach(wide_format1)
# write.csv(wide_format1, file="Medias_ED144.csv")

# rearrange data and write to .csv file
library(plyr)
medias1<-arrange(medias, media,reaction_IDs,roles,equation)
attach(medias1)
write.csv(medias1, file="Medias_ED252.csv")

