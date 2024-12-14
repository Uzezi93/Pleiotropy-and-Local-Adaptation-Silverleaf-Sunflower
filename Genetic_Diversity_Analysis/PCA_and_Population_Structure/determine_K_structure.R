#read in the data

data <-list.files("~/admixture/", pattern = ".log", full.names = T)

#look at data to make sure it only has the log files

data

#use lapply to read in all our log files at once
bigData <-lapply(1:7, FUN = function(i) readLines(data[i]))


# find the line that starts with "best like=" or just "b"
library(stringr)


#this will pull out the line that starts with "b" from each file and return it as a list
foundset<-sapply(1:7, FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])

foundset

#we're getting there!

#now we need to pull out the first number in the string, we'll do this with the function sub

as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )


#now lets store it in a dataframe
#make a dataframe with an index 1:7, this corresponds to our K values
logs <-data.frame(K = rep(1:7, each=3))

#add to it our likelihood values

logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

#and now we can calculate our delta K and probability

tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

#----------------------------------------------------------------------------------

par(mfrow=c(1,1))
par(mfrow=c(3,1)) # Adjust the margin to reduce space

# Adjusting the margins to be even smaller
par(mar = c(3, 4, 0, 2) + 0.2)

# Increase y-axis label size using par()
par(cex.lab = 1.5)  # Increase y-axis label size to 1.5 times the default


#par(mfrow=c(3,1))

# Plot admixture proportions
#read in the data for k=2
q<-read.table("~/admixture/str_k_2.qopt")

#read in the population information file located in the bam.filelist
pop<-read.table("~/admixture/samples.txt")

q <- data.frame(cbind(q,pop))
names(q) <- c("K1","K2","individuals")

q <- q %>%
  filter(!individuals %in% c("ARG1820.Aligned.sortedByCoord.out.bam", "ARG1805.Aligned.sortedByCoord.out.bam", "ARG1834.Aligned.sortedByCoord.out.bam", "Ames695.Aligned.sortedByCoord.out.bam", "Ames449.Aligned.sortedByCoord.out.bam", "btm7B-14.Aligned.sortedByCoord.out.bam", "btm5-1.Aligned.sortedByCoord.out.bam", "arg11B-11.Aligned.sortedByCoord.out.bam"))

# pop <- q %>% dplyr::select(individuals) %>% separate(individuals, into = c("pop", NA))

pop <- c("Coastal", "North", "North", "Coastal", "North", "Coastal", "Coastal", "North", "Coastal", "Coastal", "North", "North", "Coastal", "Coastal", "Coastal", "North", "North", "North", "North")

q <- cbind(q, pop)

#order it by population
ord<-order(q$pop)

# Custom colors
custom_colors <- c("darkgrey", "brown")

barplot(t(q)[,ord],
        col=custom_colors,
        names.arg = NULL,
        las=2,
        space=0,
        border="black",
        yaxt="n",
        # xlab="Individuals",
        ylab="K=2")


#-------------------------------------------------------------------------------
# Plot admixture proportions
#read in the data for k=3
q <-read.table("~/admixture/str_k_3.qopt")

#read in the population information file located in the bam.filelist
pop<-read.table("~/admixture/samples.txt")

q <- data.frame(cbind(q,pop))
names(q) <- c("K1","K2","K3","individuals")

q <- q %>%
  filter(!individuals %in% c("ARG1820.Aligned.sortedByCoord.out.bam", "ARG1805.Aligned.sortedByCoord.out.bam", "ARG1834.Aligned.sortedByCoord.out.bam", "Ames695.Aligned.sortedByCoord.out.bam", "Ames449.Aligned.sortedByCoord.out.bam", "btm7B-14.Aligned.sortedByCoord.out.bam", "btm5-1.Aligned.sortedByCoord.out.bam", "arg11B-11.Aligned.sortedByCoord.out.bam"))

# pop <- q %>% dplyr::select(individuals) %>% separate(individuals, into = c("pop", NA))

pop <- c("Coastal", "North", "North", "Coastal", "North", "Coastal", "Coastal", "North", "Coastal", "Coastal", "North", "North", "Coastal", "Coastal", "Coastal", "North", "North", "North", "North")


q <- cbind(q, pop)

#order it by population
ord<-order(q$pop)

custom_colors2 <- c("darkgrey", "brown", "aquamarine4")

barplot(t(q)[,ord],
        col=custom_colors2,
        names.arg = NULL,
        las=2,
        space=0,
        border="black",
        yaxt="n",
        # xlab="Individuals",
        ylab="K=3")

#-------------------------------------------------------------------------------
# Plot admixture proportions
#read in the data for k=4
q <-read.table("~/admixture/str_k_4.qopt")

#read in the population information file located in the bam.filelist
pop<-read.table("~/admixture/samples.txt")

q <- data.frame(cbind(q,pop))
names(q) <- c("K1","K2","K3","K4","individuals")

q <- q %>%
  filter(!individuals %in% c("ARG1820.Aligned.sortedByCoord.out.bam", "ARG1805.Aligned.sortedByCoord.out.bam", "ARG1834.Aligned.sortedByCoord.out.bam", "Ames695.Aligned.sortedByCoord.out.bam", "Ames449.Aligned.sortedByCoord.out.bam", "btm7B-14.Aligned.sortedByCoord.out.bam", "btm5-1.Aligned.sortedByCoord.out.bam", "arg11B-11.Aligned.sortedByCoord.out.bam"))

# pop <- q %>% dplyr::select(individuals) %>% separate(individuals, into = c("pop", NA))

pop <- c("Coastal", "North", "North", "Coastal", "North", "Coastal", "Coastal", "North", "Coastal", "Coastal", "North", "North", "Coastal", "Coastal", "Coastal", "North", "North", "North", "North")


q <- cbind(q, pop)

#order it by population
ord<-order(q$pop)

custom_colors3 <- c("darkgrey", "brown", "aquamarine4", "gold")

# Adjusting the margins to be even smaller
par(mar = c(4, 4, 0, 2) + 0.2)

barplot(t(q)[,ord],
        col=custom_colors3,
        names=q$pop[ord],
        las=2,
        space=0,
        border="black",
        yaxt="n",
        # xlab="Individuals",
        ylab="K=4")

#-------------------------------------------------------------------------------

