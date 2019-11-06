#Script for calculating continental syntopy using BBS 50-stop data
	
rm(list=ls())   # Remove all variables from active memory

#Merge the BBS 50-Stop data files into a single data.frame	
setwd("~/BBS/BBS data 2018/50-StopData")
x <- list.files()
data.50 <- read.csv(x[1], header=TRUE, colClasses="numeric", na.strings="*")
for(i in 2:length(x)) {
    newdata <- read.csv(x[i], header=TRUE, colClasses="numeric", na.strings="*")
    data.50 <- rbind(data.50, newdata)
}

#Open the species pair list and create list of unique AOU numbers
    #absence of a species in the BBS list was assumed to mean that the species has not been recorded in BBS surveys; species pairs including those species are tagged with In.BBS = 0

wc.pairs <- read.csv("~/BBS/Species_pairs.csv", header=TRUE)
wc.pairs <- subset(wc.pairs, wc.pairs$In.BBS==1)
	#(removes species wc.pairs that are not in BBS Species.List (i.e., one or both species missing))
species <- append(wc.pairs$AOU.1, wc.pairs$AOU.2)
species <- unique(species)

#Create a subset of the 50-stop data that includes only the species in the species list
data.50 <- subset(data.50, data.50$AOU %in% species)

#Further restrict the data to runs with 1 observer
data.50 <- subset(data.50, data.50$RPID!=203)
  #(RPID 203 is the only 2-observer code used in the entire 50-stop dataset)

#Make a list of routes in the 50-stop data subset
routes <- unique(data.50[,c(2:4)])

#Cycle through the species pairs calculating syntopy and add syntopy measures to the species pairs data frame
for(j in 1:length(wc.pairs[,1])) {
	
#Make a list of routes where both species were observed
for(i in 1:length(routes[,1])) {
    i.state <- routes$statenum[i]
    i.route <- routes$Route[i]
    
    #make a file containing only runs of route i
        route.i <- subset(data.50, data.50$statenum==i.state & data.50$Route==i.route)
    #(although there is also a country column, the state numbers are unique, i.e., not repeated in US and Canada)
    
    #generate a list of species observed on the route (already reestricted to the species on the wc taxa list)
    sp.on.route <- unique(route.i$AOU)
    
    # If both species were observed on this route, tag the route TRUE otherwise FALSE
    if(wc.pairs$AOU.1[j] %in% sp.on.route & wc.pairs$AOU.2[j] %in% sp.on.route) {
        routes$both.sp.present[i] <- TRUE
        
    } else {
        routes$both.sp.present[i] <- FALSE
    }
} #end of i loop

routes.j <- subset(routes, routes$both.sp.present==TRUE) # the subset of routes where both species were observed

if (length(routes.j[,1])==0) {  #if there are no routes on which both species were found, reconstitute routes.j with a row of NAs and zeros for this species pair and skip the k loop

routes.j <-routes[1,]
routes.j$countrynum <- rep(NA, length(routes.j[,1]))
routes.j$statenum <- rep(NA, length(routes.j[,1]))
routes.j$Route <- rep(NA, length(routes.j[,1]))
routes.j$pair.n <- rep(wc.pairs$pair.n[j], length(routes.j[,1]))
routes.j$sp1.english <- rep(wc.pairs$sp.1.english[j], length(routes.j[,1]))
routes.j$sp2.english <- rep(wc.pairs$sp.2.english[j], length(routes.j[,1])) 
routes.j$sp1.AOU <- rep(wc.pairs$AOU.1[j], length(routes.j[,1])) 
routes.j$sp2.AOU <- rep(wc.pairs$AOU.2[j], length(routes.j[,1]))  
routes.j$total.stops.1 <- rep(NA, length(routes.j[,1]))                       # Total number of stops (runs*50) for route k
routes.j$total.stops.2 <- rep(NA, length(routes.j[,1]))                       # Total number of stops on runs of route k where >=1 member of sp. pair was seen
routes.j$sp1.stops <- rep(NA, length(routes.j[,1]))
routes.j$sp2.stops <- rep(NA, length(routes.j[,1]))
routes.j$shared.stops <- rep(NA, length(routes.j[,1]))
routes.j$sp1.indiv <- rep(NA, length(routes.j[,1]))
routes.j$sp2.indiv <- rep(NA, length(routes.j[,1]))
routes.j$sp1.indiv.avg <- rep(NA, length(routes.j[,1]))
routes.j$sp2.indiv.avg <- rep(NA, length(routes.j[,1]))
routes.j$sp1.indiv.var <- rep(NA, length(routes.j[,1]))
routes.j$sp2.indiv.var <- rep(NA, length(routes.j[,1]))
  
print(routes.j) #for testing
}
  else {
# Go through routes on the list and calculate syntopy measures for species pair j
print(routes.j)  #for testing
routes.j$pair.n <- rep(wc.pairs$pair.n[j], length(routes.j[,1]))
routes.j$sp1.english <- rep(wc.pairs$sp.1.english[j], length(routes.j[,1]))
routes.j$sp2.english <- rep(wc.pairs$sp.2.english[j], length(routes.j[,1])) 
routes.j$sp1.AOU <- rep(wc.pairs$AOU.1[j], length(routes.j[,1])) 
routes.j$sp2.AOU <- rep(wc.pairs$AOU.2[j], length(routes.j[,1]))  
routes.j$total.stops.1 <- rep(NA, length(routes.j[,1]))                       # Total number of stops (runs*50) for route k
routes.j$total.stops.2 <- rep(NA, length(routes.j[,1]))                       # Total number of stops on runs of route k where >=1 member of sp. pair was seen
routes.j$sp1.stops <- rep(NA, length(routes.j[,1]))
routes.j$sp2.stops <- rep(NA, length(routes.j[,1]))
routes.j$shared.stops <- rep(0, length(routes.j[,1]))
routes.j$sp1.indiv <- rep(NA, length(routes.j[,1]))
routes.j$sp2.indiv <- rep(NA, length(routes.j[,1]))
routes.j$sp1.indiv.avg <- rep(NA, length(routes.j[,1]))
routes.j$sp2.indiv.avg <- rep(NA, length(routes.j[,1]))
routes.j$sp1.indiv.var <- rep(NA, length(routes.j[,1]))
routes.j$sp2.indiv.var <- rep(NA, length(routes.j[,1]))


for(k in 1:length(routes.j[,1])) {
	  # Create a subset of the 50-stop data for the current route
    route.sub <- subset(data.50, data.50$statenum==routes.j$statenum[k] & data.50$Route==routes.j$Route[k])
       
    # Create a variable that is unique for each combination of year and RPID
    route.sub$yearRPID <- paste(route.sub$year, route.sub$RPID, sep=":")

    # Subset data by species, returning only the rows for a species in question.
    sp1.sp <- subset(route.sub, route.sub$AOU==wc.pairs$AOU.1[j])
    sp2.sp <- subset(route.sub, route.sub$AOU==wc.pairs$AOU.2[j])
    
    #Calculate the possible number of stops where the species could have been observed
    routes.j$total.stops.1[k] <- length(unique(route.sub$yearRPID)) * 50 						# Total number of stops (runs*50) for route k
    routes.j$total.stops.2[k] <- length(unique(c(sp1.sp$yearRPID, sp2.sp$yearRPID))) * 50   # Total number of stops on runs of route k where >=1 member of sp. pair was seen
    	#(using yearRPID instead of year allows for multiple runs per year)
    	
   # Create datasets containing only the 50-stop data for each species on route k
    sp1.counts <- sp1.sp[,8:57]
    sp2.counts <- sp2.sp[,8:57]
    
   # Calculate the number of stops on route k where each species was observed
    routes.j$sp1.stops[k] <- sum(sp1.counts > 0)
    routes.j$sp2.stops[k] <- sum(sp2.counts > 0)
    
    # Calculate the number of stops on route k where both species were observed
    routes.j$shared.stops[k] <- 0
    # For each unique run, add up the total number of shared stops between the two species
    for(l in 1:length(sp1.sp[,1])) {                                          
        yR <- sp1.sp$yearRPID[l]
        if(yR %in% sp2.sp$yearRPID) {
            routes.j$shared.stops[k] <- routes.j$shared.stops[k] +
                sum((sp2.sp[sp2.sp$yearRPID == yR, 8:57] * sp1.sp[l,8:57]) != 0)
                #(multiplies the vectors of 50 stop counts of the sp1 and sp2 species and counts the number of products > 0)
        }
    } #end of l loop
    
  # Calculate the total number of individuals of each species observed on route k
    routes.j$sp1.indiv[k] <- sum(sp1.counts)
    routes.j$sp2.indiv[k] <- sum(sp2.counts)
    
    # Calculate the mean of the nonzero counts of each species
    routes.j$sp1.indiv.avg[k] <- mean(sp1.counts[sp1.counts > 0])
    routes.j$sp2.indiv.avg[k] <- mean(sp2.counts[sp2.counts > 0])
    
    # Calculate the variance of the nonzero counts of each species
    routes.j$sp1.indiv.var[k] <- var(sp1.counts[sp1.counts > 0])
    routes.j$sp2.indiv.var[k] <- var(sp2.counts[sp2.counts > 0])   
    
} #end of k loop
}

    if(j==1) {
        write.table(routes.j, file="~/BBS/Continental.syntopy.counts.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
        } else {
        write.table(routes.j, file="~/BBS/Continental.syntopy.counts.csv", append=TRUE, quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)
        }
        print(paste("Saved after ", j, " wc.pairs on ", date(), sep=""))
   
} #end of j loop

############################################################################################################################################################
    
#Calculate indices of syntopy for each species pair across all routes where both species were observed

rm(list=ls())   # Remove all variables from active memory
counts <- read.csv("~/BBS/Continental.syntopy.counts.csv", header=TRUE)
  
  # Generate a list of unique species pairs and include english, latin and AOU names
  u.pairs <- unique(counts[,c(4:9)])
  names(u.pairs)[1] <- "on.same.route"  #changing the name of both.sp.present
  
  # Calculate indices of syntopy for each species pair
  for(i in 1:length(u.pairs[,1])) {
    # Subset the data by species pair
    sub <- subset(counts, pair.n == u.pairs$pair.n[i])
  
  # Add variables to the species pair list
  u.pairs$num.routes[i] <- length(sub[,1])                             # Number of routes in the subset created above
  if (u.pairs$on.same.route[i]==FALSE) {u.pairs$num.routes[i] <- 0}   #If the species were not found on the same route, reset number of routes to 0
  u.pairs$tot.stops[i] <- sum(sub$total.stops.1)                       # Total number of stops on all routes in the subset created above
  u.pairs$tot.sp1[i] <- sum(sub$sp1.stops)                         # Total number of stops at which sp1 species was observed
  u.pairs$tot.sp2[i] <- sum(sub$sp2.stops)                         # Total number of stops at which sp2 species was observed
  u.pairs$tot.both[i] <- sum(sub$shared.stops)                         # Total number of stops at which both species were observed
  u.pairs$sp1.indiv[i] <- sum(sub$sp1.indiv)                         # Total number of individuals of sp 1 observed
  u.pairs$sp2.indiv[i] <- sum(sub$sp2.indiv)                         # Total number of individuals of sp 2 observed
  u.pairs$sp1.indiv.avg[i] <- sum(sub$sp1.indiv.avg)                         # Average number individuals of sp 1 observed per stop
  u.pairs$sp2.indiv.avg[i] <- sum(sub$sp2.indiv.avg)                         # Total number of individuals of sp 2 observed
  
  if(u.pairs$num.routes[i]==0) {  #if there are no routes on which both species were observed, set syntopy to NA
    u.pairs$syntopy1w[i] <-NA
    u.pairs$syntopy2w[i] <-NA
    u.pairs$syntopy1u[i] <- NA
    u.pairs$syntopy2u[i] <- NA
    u.pairs$syntopyOw[i] <- NA
    u.pairs$syntopyOu[i] <- NA
  }
  
  else {
  # Proportion metric of syntopy, weighted (i.e. calculated from totals, which means that routes with more runs are weighted more than routes with fewer runs)
  # This measure of syntopy is based on the proportion of stops that the less frequent species occupied, that the more common species also occupied

  if(u.pairs$tot.sp1[i] <= u.pairs$tot.sp2[i]) {
    u.pairs$syntopy1w[i] <- u.pairs$tot.both[i] / u.pairs$tot.sp1[i]   # If sp1 was seen at fewer stops than sp2
  } else {
    u.pairs$syntopy1w[i] <- u.pairs$tot.both[i] / u.pairs$tot.sp2[i]  # If sp2 was seen at fewer stops than sp1
  }
  
  # Observed vs. expected metric of syntopy, weighted (i.e. calculated from totals, which means that routes with more runs are weighted more than routes with fewer runs)
  u.pairs$syntopy2w[i] <- u.pairs$tot.both[i] / ((u.pairs$tot.sp1[i]/u.pairs$tot.stops[i]) * u.pairs$tot.sp2[i])
    #(the denominator is npq where q = proportion of stops with sp 1, q = proportion of stops with sp 2 and n = total number of stops)
  
  #Overall syntopy metric (weighted)
    u.pairs$syntopyOw[i] <-u.pairs$tot.both[i] / (u.pairs$tot.sp1[i] + u.pairs$tot.sp2[i] - u.pairs$tot.both[i])
    #denominator is the number of stops where one or both species were recorded, pooled across all runs of all routes
  
  # Proportion metric of syntopy, unweighted (i.e. calculated for each route separately and then averaged across routes)
  for(j in 1:length(sub[,1])) {
    # Proportion metric of syntopy
    if(sub$sp1.stops[j] <= sub$sp2.stops[j]) {
      sub$syntopy1u[j] <- sub$shared.stops[j] / sub$sp1.stops[j]
    } else {
      sub$syntopy1u[j] <- sub$shared.stops[j] / sub$sp2.stops[j]
    }
    
    # Observed vs. expected metric of syntopy, unweighted (i.e. calculated for each route separately and then averaged across routes)
    sub$syntopy2u[j] <- sub$shared.stops[j] / (sub$sp1.stops[j] * sub$sp2.stops[j] / sub$total.stops.1[j])
    
    #Overall syntopy (unweighted)
    sub$syntopyOu[j] <- sub$shared.stops[j] / (sub$sp1.stops[j] + sub$sp2.stops[j] - sub$shared.stops[j])
  }
  # Average across routes to get unweighted average syntopy estimates
  u.pairs$syntopy1u[i] <- mean(sub$syntopy1u)
  u.pairs$syntopy2u[i] <- mean(sub$syntopy2u)
  u.pairs$syntopyOu[i] <- mean(sub$syntopyOu)
  
  if(i %% 10 == 0) {
    print(paste("Finished ", i, " pairs on ", date(), sep=""))            # Progress report
  }
  }
  }
write.table(u.pairs, file="~/BBS/Continental_syntopy.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

