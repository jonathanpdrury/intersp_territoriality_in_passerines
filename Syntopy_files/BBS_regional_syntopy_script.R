# Step 1: Get list of routes within 250km of study location that were run in the time period specified
# Step 2: Narrow these routes to the runs that have 50-stop data
# Step 3: Identify routes where both IT species were observed and create a list of species pairs involving one or both IT species
# Step 4: For each species pair in each route, get counts of total, focal, other, shared stops
# Step 5: For each species pair at each location, calculate indices of syntopy
# Step 6: Average the syntopy estimates available for each species pair of interest

# --------------------------------------------------
# Step 1: Get list of routes within 250km of study location that were run in the time period specified

install.packages("fossil")

rm(list=ls())
library(fossil)
weather <- read.csv("~/BBS/BBS data 2018/weather.csv", header=TRUE)
routes <- read.csv("~/BBS/BBS data 2018/routes.csv", header=TRUE)
locations <- read.csv("~/BBS/locations.csv", header=TRUE)
#(locations.csv file is a list of the locations where interspecific territoriality was reported; it also lists the IT species pairs)

#Make a list of all routes that were run from 1982 to 2017, including their lat long coordinates

year.lower<-1982
year.upper<-2017


    # Define variables
    run.lower <- rep(NA, length(routes$Route))
    run.upper <- rep(NA, length(routes$Route))
    run.tot.yr <- rep(0, length(routes$Route))
    
    # For each route in the "routes" data.frame...
    for(i in 1:length(routes$Route))
    {
        state.id <- routes$statenum[i]                   # Get state ID number for the route
        route.id <- routes$Route[i]                      # Get route ID number for the route
        for(j in year.lower:year.upper)                  # For each year in the period of interest...
        {
            if(length(weather$Route[weather$statenum==state.id & weather$Route==route.id & weather$Year==j & weather$RunType==1]) > 0) # Only proceed if there's data for the given route in the given year

            {
                
                # Also, increment run.tot.yr by 1
                run.tot.yr[i] <- run.tot.yr[i]+1
                
                # If we haven't set run.lower yet (i.e. if it's still at its default NA)
                # Then set run.lower for route i to year j
                if(is.na(run.lower[i]))
                {
                    run.lower[i] <- j

                }
                if(is.na(run.upper[i])) run.upper[i] <- j
                if(run.upper[i] < j) run.upper[i] <- j
            }
        }
    }
       
    rts <- data.frame(cbind(routes$statenum, routes$Route, routes$Longi, routes$Lati, run.lower, run.upper, run.tot.yr))
    names(rts)[1] <- "statenum"
    names(rts)[2] <- "Route"
    names(rts)[3] <- "Long"
    names(rts)[4] <- "Lat"
    
    #remove routes that were not run during the specified time period
	rts<-rts[rts$run.tot.yr>0,]
	
# Make a list of routes within 250km of each location

for(j in 1:length(locations$FocalSp)) {                                         

	#add distance to current location from each of the routes that were run in the year range specified above
    rts$route.dist <- NA
    for(k in 1:length(rts[,1])) {
    	rts$route.dist[k] <- deg.dist(-1*locations$Long[j], locations$Lat[j], rts$Long[k], rts$Lat[k])
    }
    
    #remove routes that are farther than 250 km from the study location
    routes.250 <- rts[rts$route.dist < 250,]
    
    locations$n.routes.250[j]<-length(routes.250[,1])
    
    if(length(routes.250[,1])>0) { #skip locations with no routes within 250 km
        
    # Add variables for species pair, study, location, focal species to the routes.250 dataset
    	#(routes.250 contains the routes < 250 for the current location)
    	
    for(k in 1:length(routes.250[,1])) {
        routes.250$SpPair[k] <- locations$SpPair[j]
        routes.250$StudyID[k] <- locations$StudyID[j]
        routes.250$LocationID[k] <- locations$LocationID[j]
        routes.250$FocalSp[k] <- locations$FocalSp[j]
        routes.250$ITSp[k] <- locations$ITSp[j]
    }
    
    # Status update
    print(paste(j, " locations completed.", sep=""))
    
    if(j==1) {
        write.table(routes.250, file="~/BBS/Routes250.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
    } else {
        write.table(routes.250, file="~/BBS/Routes250.csv", append=TRUE, quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)
    }
}
}
#save the locations file with the new routes.250 column
write.table(locations, file="~/BBS/locations.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)



# --------------------------------------------------
# Step 2: Make a subset of the BBS 50-stop data that includes runs for the routes in Routes250.csv  
	
# Read in the routes.250 dataset generated in step 1:
routes.250 <- read.csv("~/BBS/Routes250.csv", header=TRUE)

setwd("~/BBS/BBS data 2018/50-StopData")
x <- list.files()

# Read in the 50-stop BBS dataset (comes in several separate .csv files)
data.50 <- read.csv(x[1], header=TRUE, colClasses="numeric", na.strings="*")
for(i in 2:length(x)) {
    newdata <- read.csv(x[i], header=TRUE, colClasses="numeric", na.strings="*")
    data.50 <- rbind(data.50, newdata)
}

unique.runs.50 <- unique(data.50[,3:6])

unique.routes.250 <- unique(routes.250[,c(1:2,9:13)])
	#(this is not really a unique list of routes but a unique list of "statenum Route SpPair StudyID LocationID FocalSp ITSp")
	
# Generate a data.frame containing a list of all runs of routes within 250km of study locations
# For each route in unique.routes.250...
for(b in 1:length(unique.routes.250[,1])) {
    # If this is the first route...
    if(b==1) {
        # Create a subset of unique.runs.50 restricted to the current route, omitting 2-observer counts
        runs.250 <- unique.runs.50[unique.runs.50$statenumstatenum==unique.routes.250$statenum[b]
            & unique.runs.50$Route==unique.routes.250$Route[b]
            & unique.runs.50$RPID!=203,]
        # Add columns representing SpPair, StudyID, etc.
        for(c in 1:5) {
            runs.250[,c+4] <- rep(unique.routes.250[b,c+2], length(runs.250[,1]))
        }
    # Otherwise, if this is not the first route...
    } else {
        # Do the same thing, but append the data set to the end of the existing runs.250 data set
        tempRuns <- unique.runs.50[unique.runs.50$statenum==unique.routes.250$statenum[b]
            & unique.runs.50$Route==unique.routes.250$Route[b]
            & unique.runs.50$RPID!=203,]
        for(c in 1:5) {
            tempRuns[,c+4] <- rep(unique.routes.250[b,c+2], length(tempRuns[,1]))
        }
        runs.250 <- rbind(runs.250, tempRuns)
    }
}
# Name the new columns
names(runs.250)[5:9] <- c("SpPair", "StudyID", "LocationID", "FocalSp", "ITSp")     
 # Create a dummy matching variable by concatenating several variables
runs.250$matching.factor <- paste(runs.250[,1],runs.250[,2],runs.250[,3],runs.250[,4],sep=":")
          
# Generate a list of unique "matching factors" from the previous line
unique.matching.factor <- unique(runs.250$matching.factor)         
                                      
data.50$matching.factor <- paste(data.50[,3],data.50[,4],data.50[,5],data.50[,6],sep=":")

# Subset the 50-stop dataset so it only includes routes for which we will
# require data at a later stage (only those within 250km of a study location)
data.50.final <- subset(data.50, data.50$matching.factor %in% unique.matching.factor)

write.table(data.50.final, file="~/BBS/Data_50-stop.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(runs.250, file="~/BBS/Runs_250.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

# --------------------------------------------------
# Step 3: Identify routes where both IT species were observed and create a list of species pairs involving one or both IT species

species.list <- read.csv("~/BBS/SpeciesList.csv", header=TRUE)
data.50.final <- read.csv("~/BBS/Data_50-stop.csv", header=TRUE)
routes.250 <- read.csv("~/BBS/Routes250.csv", header=TRUE)

# Generate a list of unique routes within 250km of study locations
unique.routes.250 <- unique(routes.250[,c(1:2,9:13)])
	#(this is not really a unique list of routes but a unique list of "statenum Route SpPair StudyID LocationID FocalSp ITSp")

#First make a list 

# Define variables that will be used to build a species-pair-list dataset
sp.pair <- numeric()
study.ID <- numeric()
location.ID <- numeric()
state <- numeric()
route <- numeric()
IT.sp.1 <- numeric()                              # Sp 1 of IT pair
IT.sp.2 <- numeric()                              # Sp 2 of IT pair
IT.species <- numeric()                           # IT species of focal pair
other.species <- numeric()                        # other species of focal pair

# For each route within 250km of a study location...
#for(i in 3:100) {                                                               # Shortened loop for testing
for(i in 1:length(unique.routes.250[,1])) {
    i.state <- unique.routes.250$statenum[i]
    i.route <- unique.routes.250$Route[i]
    
    # Subset the big 50-stop dataset by route
    route.sub <- subset(data.50.final, data.50.final$statenum==i.state
        & data.50.final$Route==i.route)
    
    # Generate a list of species observed in the route and included on species.list
    sp.by.route <- unique(route.sub$AOU)
    sp.by.route <- subset(sp.by.route, sp.by.route %in% species.list$AOU)
    
#For the routes where both IT species were recorded, construct a list of species pairs involving the IT species and all other species of interest (as defined by specieslist) recorded on the route
      
    if(unique.routes.250$FocalSp[i] %in% sp.by.route & unique.routes.250$ITSp[i] %in% sp.by.route) {
        unique.routes.250$both.sp.present[i] <- TRUE  #G tags route as TRUE if both IT species were seen on the route
        
        num.sp <- length(sp.by.route)                                           # Number of species observed on a route
        
 # Generate a list of routes for each non-IT species pair where one of the species is one of the IT species
        
        sp.pair <- append(sp.pair, rep(unique.routes.250$SpPair[i], 2*(num.sp-2)))
        study.ID <- append(study.ID, rep(unique.routes.250$StudyID[i], 2*(num.sp-2)))
        location.ID <- append(location.ID, rep(unique.routes.250$LocationID[i], 2*(num.sp-2)))
        state <- append(state, rep(unique.routes.250$statenum[i], 2*(num.sp-2)))
        route <- append(route, rep(unique.routes.250$Route[i], 2*(num.sp-2)))
        IT.sp.1 <- append(IT.sp.1, rep(unique.routes.250$FocalSp[i], 2*(num.sp-2)))
        IT.sp.2 <- append(IT.sp.2, rep(unique.routes.250$ITSp[i], 2*(num.sp-2)))
        IT.species <- append(IT.species, c(rep(unique.routes.250$FocalSp[i], (num.sp-2)), rep(unique.routes.250$ITSp[i], (num.sp-2))))
        other.species <- append(other.species, rep(sp.by.route[sp.by.route!=unique.routes.250$FocalSp[i] &
            sp.by.route!=unique.routes.250$ITSp[i]], 2))
    } else {
        unique.routes.250$both.sp.present[i] <- FALSE #tags route as FALSE if both IT species were not seen on the route, and leaves this route out of sp.pair.list

    }
}
sp.pair.list <- data.frame(sp.pair, study.ID, location.ID, state, route, IT.sp.1, IT.sp.2, IT.species, other.species)

 # Generate a list of routes for each IT species pair 
 
# Define variables that will be used to build a species-pair-list dataset
sp.pair <- numeric()
study.ID <- numeric()
location.ID <- numeric()
state <- numeric()
route <- numeric()
IT.sp.1 <- numeric()
IT.sp.2 <- numeric()
IT.species <- numeric()
other.species <- numeric()

# For each route within 250km of a study location...
# for(i in 1:100) {                                                               # Shortened loop for testing
for(i in 1:length(unique.routes.250[,1])) {
    i.state <- unique.routes.250$statenum[i]
    i.route <- unique.routes.250$Route[i]
    
    # Subset the big 50-stop dataset by route
    route.sub <- subset(data.50.final, data.50.final$statenum==i.state
        & data.50.final$Route==i.route)
        
    # Generate a list of unique species observed in the route
    sp.by.route <- unique(route.sub$AOU)
    sp.by.route <- subset(sp.by.route, sp.by.route %in% species.list$AOU)
       
    # IF both of the focal species were observed on this route...
    if(unique.routes.250$FocalSp[i] %in% sp.by.route & unique.routes.250$ITSp[i] %in% sp.by.route) {
        
        # Generate a list of two species pairs. In each pair, both members of the pair are
        # the focal (e.g., IT) species, but the order differs
        
        sp.pair <- append(sp.pair, rep(unique.routes.250$SpPair[i], 2))
        study.ID <- append(study.ID, rep(unique.routes.250$StudyID[i], 2))
        location.ID <- append(location.ID, rep(unique.routes.250$LocationID[i], 2))
        state <- append(state, rep(unique.routes.250$statenum[i], 2))
        route <- append(route, rep(unique.routes.250$Route[i], 2))
        IT.sp.1 <- append(IT.sp.1, rep(unique.routes.250$FocalSp[i], 2))
        IT.sp.2 <- append(IT.sp.2, rep(unique.routes.250$ITSp[i], 2))
        IT.species <- append(IT.species, c(unique.routes.250$FocalSp[i], unique.routes.250$ITSp[i]))
        other.species <- append(other.species, c(unique.routes.250$ITSp[i], unique.routes.250$FocalSp[i]))
    }
    
    if(i %% 100 == 0) {
        print(paste("Finished ", i, " pairs on ", date(), sep=""))                   # Progress report
    }
}
sp.pair.list2 <- data.frame(sp.pair, study.ID, location.ID, state, route, IT.sp.1, IT.sp.2, IT.species, other.species)
# ____________________________________________________________

write.table(sp.pair.list, file="~/BBS/50-stop-pairs.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(sp.pair.list2, file="~/BBS/50-stop-pairs-focal.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(unique.routes.250, file="~/BBS/50-stop-routes.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

# --------------------------------------------------


# Step 4: For each species pair in each route, get counts of total, focal, other, shared stops

rm(list=ls())
setwd("~/BBS/")
unique.routes.250 <- read.csv("50-stop-routes.csv", header=TRUE)
runs.250 <- read.csv("Runs_250.csv", header=TRUE)
pairs <- read.csv("50-stop-pairs.csv", header=TRUE)
pairs.focal <- read.csv("50-stop-pairs-focal.csv", header=TRUE)
data.50 <- read.csv("Data_50-stop.csv", header=TRUE)

pairs$stateroute <- paste(pairs$state, pairs$route, sep=":")
pairs$total.stops.1 <- rep(NA, length(pairs[,1]))
pairs$total.stops.2 <- rep(NA, length(pairs[,1]))
pairs$focal.stops <- rep(NA, length(pairs[,1]))
pairs$other.stops <- rep(NA, length(pairs[,1]))
pairs$shared.stops <- rep(0, length(pairs[,1]))
pairs$focal.indiv <- rep(NA, length(pairs[,1]))
pairs$other.indiv <- rep(NA, length(pairs[,1]))
pairs$focal.indiv.avg <- rep(NA, length(pairs[,1]))
pairs$other.indiv.avg <- rep(NA, length(pairs[,1]))
pairs$focal.indiv.var <- rep(NA, length(pairs[,1]))
pairs$other.indiv.var <- rep(NA, length(pairs[,1]))
cur.route <- "start"
print(paste("Start time: ", date(), sep=""))

# For each species pair on each route
for(i in 1:length(pairs[,1])) {
    # Check if cur.route corresponds to the appropriate route under consideration
    if(pairs$stateroute[i] != cur.route) {  #this avoids having to subset data.50 every time through the loop
        cur.route <- pairs$stateroute[i]    #only re-subsets the data when the route number has changed
        # Subset 50-stop data based on current route
        route.sub <- subset(data.50, data.50$statenum==pairs$state[i]
            & data.50$Route==pairs$route[i])
        # Create a dummy variable that is unique for each combination of year and RPID
        route.sub$yearRPID <- paste(route.sub$year, route.sub$RPID, sep=":")
    }
    
    # Subset data for one route by species, returning only the rows for a species in question.
    focal.sp <- subset(route.sub, route.sub$AOU==pairs$IT.species[i]) #runs of this route on which focal sp was recorded
    other.sp <- subset(route.sub, route.sub$AOU==pairs$other.species[i]) #runs of this route on which other sp was recorded
    
    # Calculate the number of stops at which a species could have been observed
    pairs$total.stops.1[i] <- length(unique(route.sub$yearRPID)) * 50 # runs of this route times the number of stops per run
    pairs$total.stops.2[i] <- length(unique(c(focal.sp$yearRPID, other.sp$yearRPID))) * 50 # runs of this route on which either species was recorded times the number of stops per run

    # Create a dataset containing only the 50-stop data for one species.
    focal.counts <- focal.sp[,8:57]
    other.counts <- other.sp[,8:57]
    
    # Calculate the number of stops where each species was observed
    pairs$focal.stops[i] <- sum(focal.counts > 0)
    pairs$other.stops[i] <- sum(other.counts > 0)
    
    # Calculate the number of stops where both species were observed 
    pairs$shared.stops[i] <- 0
    # For each unique run, add up the total number of shared stops between the two species
    for(j in 1:length(focal.sp[,1])) {                                          
    	        yR <- focal.sp$yearRPID[j]
        if(yR %in% other.sp$yearRPID) {
            pairs$shared.stops[i] <- pairs$shared.stops[i] +
                sum((other.sp[other.sp$yearRPID == yR, 8:57] * focal.sp[j,8:57]) != 0)
        }
    }
    
    # Calculate the total number of each species observed on the current route
    pairs$focal.indiv[i] <- sum(focal.counts)
    pairs$other.indiv[i] <- sum(other.counts)
    
    # Calculate the mean of the nonzero counts of each species
    pairs$focal.indiv.avg[i] <- mean(focal.counts[focal.counts > 0])
    pairs$other.indiv.avg[i] <- mean(other.counts[other.counts > 0])
    
    # Calculate the variance of the nonzero counts of each species
    pairs$focal.indiv.var[i] <- var(focal.counts[focal.counts > 0])
    pairs$other.indiv.var[i] <- var(other.counts[other.counts > 0])
    
    # Output and diagnostics
    if(i %% 100 == 0) {
        print(paste("Finished ", i, " pairs on ", date(), sep=""))
    }
    if(i %% 1000 == 0) {
        write.table(pairs, file="~/BBS/50-stop-pairs2.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
        print(paste("Saved after ", i, " pairs on ", date(), sep=""))
    }  
}
write.table(pairs, file="~/BBS/50-stop-pairs2.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


# ____________________________________________________________
# Now, do the same thing for the species pairs involving the two IT species

pairs.focal$stateroute <- paste(pairs.focal$state, pairs.focal$route, sep=":")
pairs.focal$total.stops.1 <- rep(NA, length(pairs.focal[,1]))                   # Total number of runs for route in question * 50
pairs.focal$total.stops.2 <- rep(NA, length(pairs.focal[,1]))                   # Total number of runs for route in question where >=1 member of sp. pair was seen
pairs.focal$focal.stops <- rep(NA, length(pairs.focal[,1]))                     # Number of stops at which focal sp. observed
pairs.focal$other.stops <- rep(NA, length(pairs.focal[,1]))                     # Number of stops at which other sp. observed
pairs.focal$shared.stops <- rep(0, length(pairs.focal[,1]))                     # Number of stops at which both sp. observed
pairs.focal$focal.indiv <- rep(NA, length(pairs.focal[,1]))                     # Number of individuals of the focal sp. observed
pairs.focal$other.indiv <- rep(NA, length(pairs.focal[,1]))                     # Number of individuals of the other sp. observed
pairs.focal$focal.indiv.avg <- rep(NA, length(pairs.focal[,1]))                 # Mean # of individuals of focal sp. at stops where it's observed
pairs.focal$other.indiv.avg <- rep(NA, length(pairs.focal[,1]))                 # Mean # of individuals of other sp. at stops where it's observed
pairs.focal$focal.indiv.var <- rep(NA, length(pairs.focal[,1]))                 # Variance of individuals of focal sp. at stops where it's observed
pairs.focal$other.indiv.var <- rep(NA, length(pairs.focal[,1]))                 # Variance of individuals of other sp. at stops where it's observed
cur.route <- "start"
print(paste("Start time: ", date(), sep=""))
for(i in 1:length(pairs.focal[,1])) {                                                            
    if(pairs.focal$stateroute[i] != cur.route) {
        cur.route <- pairs.focal$stateroute[i]
        route.sub <- subset(data.50, data.50$statenum==pairs.focal$state[i]
            & data.50$Route==pairs.focal$route[i])
        route.sub$yearRPID <- paste(route.sub$year, route.sub$RPID, sep=":")
    }
    focal.sp <- subset(route.sub, route.sub$AOU==pairs.focal$IT.species[i])
    other.sp <- subset(route.sub, route.sub$AOU==pairs.focal$other.species[i])
    pairs.focal$total.stops.1[i] <- length(unique(route.sub$yearRPID)) * 50
    pairs.focal$total.stops.2[i] <- length(unique(c(focal.sp$yearRPID, other.sp$yearRPID))) * 50
    focal.counts <- focal.sp[,8:57]
    other.counts <- other.sp[,8:57]
    pairs.focal$focal.stops[i] <- sum(focal.counts > 0)
    pairs.focal$other.stops[i] <- sum(other.counts > 0)
    pairs.focal$shared.stops[i] <- 0
    for(j in 1:length(focal.sp[,1])) {
        yR <- focal.sp$yearRPID[j]
        if(yR %in% other.sp$yearRPID) {
            pairs.focal$shared.stops[i] <- pairs.focal$shared.stops[i] +
                sum((other.sp[other.sp$yearRPID == yR, 8:57] * focal.sp[j,8:57]) != 0)
        }
    }
    pairs.focal$focal.indiv[i] <- sum(focal.counts)
    pairs.focal$other.indiv[i] <- sum(other.counts)
    pairs.focal$focal.indiv.avg[i] <- mean(focal.counts[focal.counts > 0])
    pairs.focal$other.indiv.avg[i] <- mean(other.counts[other.counts > 0])
    pairs.focal$focal.indiv.var[i] <- var(focal.counts[focal.counts > 0])
    pairs.focal$other.indiv.var[i] <- var(other.counts[other.counts > 0])
    if(i %% 100 == 0) {
        print(paste("Finished ", i, " pairs on ", date(), sep=""))
    }
    if(i %% 200 == 0) {
        write.table(pairs.focal, file="~/BBS/50-stop-pairs-focal2.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
        print(paste("Saved after ", i, " pairs on ", date(), sep=""))
    }  
}
write.table(pairs.focal, file="~/BBS/50-stop-pairs-focal2.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

# --------------------------------------------------
# Step 5: Calculate indices of syntopy for each species pair at each location

pairs <- read.csv("~/BBS/50-stop-pairs2.csv", header=TRUE)
pairs.focal <- read.csv("~/BBS/50-stop-pairs-focal2.csv", header=TRUE)

pairs <- rbind(pairs, pairs.focal)

# Generate a list of unique species pairs at each location
u.pair.loc <- unique(pairs[,c(1:3,8:9)])

# Calculate indices of syntopy for each species pair at each location
for(i in 1:length(u.pair.loc[,1])) {
# for(i in 1:50) {
    # Subset the pairs data set by study, location, and species pair
    sub <- subset(pairs, sp.pair==u.pair.loc$sp.pair[i] &
        study.ID==u.pair.loc$study.ID[i] & location.ID==u.pair.loc$location.ID[i] &
        IT.species==u.pair.loc$IT.species[i] & other.species==u.pair.loc$other.species[i])
    
    # Filling in some variables
    u.pair.loc$num.routes[i] <- length(sub[,1])                             # Number of routes in the subset created above
    u.pair.loc$tot.stops[i] <- sum(sub$total.stops.1)                       # Total number of stops on all routes in the subset created above
    u.pair.loc$tot.focal[i] <- sum(sub$focal.stops)                         # Total number of stops at which focal species was observed
    u.pair.loc$tot.other[i] <- sum(sub$other.stops)                         # Total number of stops at which other species was observed
    u.pair.loc$tot.both[i] <- sum(sub$shared.stops)                         # Total number of stops at which both species were observed
    
    # If the focal species was seen at fewer stops than the other species...
    # Proportion metric of syntopy
    # This measure of syntopy is based on the proportion of stops that the less frequent species occupied, that the more common species also occupied
    if(u.pair.loc$tot.focal[i] <= u.pair.loc$tot.other[i]) {
        u.pair.loc$syntopy1w[i] <- u.pair.loc$tot.both[i] / u.pair.loc$tot.focal[i]
    } else {
        u.pair.loc$syntopy1w[i] <- u.pair.loc$tot.both[i] / u.pair.loc$tot.other[i]
    }
    
    # Observed vs. expected (weighted, i.e. calculated by run)
    u.pair.loc$syntopy2w[i] <- u.pair.loc$tot.both[i] / (u.pair.loc$tot.focal[i]*(u.pair.loc$tot.other[i] / u.pair.loc$tot.stops[i]))
        #G the denominator is npq where q = proportion of stops with sp 1, q = proportion of stops with sp 2 and n = total number of stops
        #(this order of calculation avoids integer overflow)
        
    #Overall syntopy metric (weighted)
    u.pair.loc$syntopyOw[i] <-u.pair.loc$tot.both[i] / (u.pair.loc$tot.focal[i] + u.pair.loc$tot.other[i] - u.pair.loc$tot.both[i])
    #denominator is the number of stops where one or both species were recorded, pooled across all runs of all routes
    
    # For each route in subset "sub"...
    	#G Some routes were run more than once; this way of doing the calculation gives equal weight to all routes but still uses all of the runs per route
    for(j in 1:length(sub[,1])) {
        # If the focal species was seen at fewer stops than the other species...
        # Proportion metric of syntopy
        if(sub$focal.stops[j] <= sub$other.stops[j]) {
            sub$syntopy1u[j] <- sub$shared.stops[j] / sub$focal.stops[j]
        } else {
            sub$syntopy1u[j] <- sub$shared.stops[j] / sub$other.stops[j]
        }
        
        # Observed vs. expected (unweighted, i.e. calculated by route)
        sub$syntopy2u[j] <- sub$shared.stops[j] / (sub$focal.stops[j] * sub$other.stops[j] / sub$total.stops.1[j])
        
        #Overall syntopy (unweighted)
        sub$syntopyOu[j] <- sub$shared.stops[j] / (sub$focal.stops[j] + sub$other.stops[j] - sub$shared.stops[j])

    }
    # Average across routes to get unweighted average syntopy estimates
    u.pair.loc$syntopy1u[i] <- mean(sub$syntopy1u)
    u.pair.loc$syntopy2u[i] <- mean(sub$syntopy2u)
    u.pair.loc$syntopyOu[i] <- mean(sub$syntopyOu)
    
    if(i %% 100 == 0) {
        print(paste("Finished ", i, " pairs on ", date(), sep=""))            # Progress report
    }
}
write.table(u.pair.loc, file="~/BBS/Syntopy-pair-location.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


#Step 6: Average the syntopy estimates available for each species pair of interest

#(the species pairs of interest could include all or a subset of the species pairs in Syntopy-pair-location.csv)

rm(list=ls())

data <- read.csv("~/BBS/Syntopy-pair-location.csv", header=TRUE)
pairs <-read.csv("~/BBS/Species_pairs.csv", header=TRUE)

#add unique species pair and location by species pair columns to syntopy data
for(i in 1:length(data[,1])){
if(data$IT.species[i]<data$other.species[i]) {data$spsp[i]<-as.character(paste(data$IT.species[i],data$other.species[i],sep=""))}
else {data$spsp[i]<-as.character(paste(data$other.species[i],data$IT.species[i],sep=""))}
}
data$loc.spsp<-as.character(paste(data$location.ID,data$spsp,sep="."))
data$IT.spsp<-as.character(paste(data$sp.pair,data$spsp,sep="."))
data$IT.loc.spsp<-as.character(paste(data$sp.pair,data$loc.spsp,sep="."))

#add unique species pair column to species pairs list
for(i in 1:length(pairs[,1])){
if(pairs$AOU.1[i]<pairs$AOU.2[i]) {pairs$spsp[i]<-as.character(paste(pairs$AOU.1[i],pairs$AOU.2[i],sep=""))}
else {pairs$spsp[i]<-as.character(paste(pairs$AOU.2[i],pairs$AOU.1[i],sep=""))}
}
#subset data to species pairs on list
data <- subset(data, data$spsp %in% pairs$spsp)

#remove duplicates
data$dup<-0
data <- data[order(data$IT.loc.spsp),]
for(i in 2:length(data[,1])){
	if(data$IT.loc.spsp[i]==data$IT.loc.spsp[i-1]) {data$dup[i]<-1}
	}
data<-subset(data,data$dup!=1)  

#First average within species pairs over locations within IT pairs
#make a list of IT.spsp to store interim sums and averages
int<- unique(data[,c(1,17,19)])  
for (i in 1:length(int[,1])) {
	sub<-subset(data,data$IT.spsp==int$IT.spsp[i])
	int$locations[i] <-length(sub$sp.pair)
	int$num.routes[i] <-sum(sub$num.routes)
	int$tot.stops[i] <-sum(sub$tot.stops)
	int$tot.focal[i]<-sum(sub$tot.focal)
	int$tot.other[i]<-sum(sub$tot.other)
	int$tot.both[i]<-sum(sub$tot.both)
	int$syntopy1w[i]<-mean(sub$syntopy1w)
	int$syntopy2w[i]<-mean(sub$syntopy2w)
	int$syntopy1u[i]<-mean(sub$syntopy1u)
	int$syntopy2u[i]<-mean(sub$syntopy2u)}

#Next average within species pairs over IT pairs

#make a list of spsp to store final sums and averages
final<- data.frame(unique(int$spsp))
names(final)[1] <- "spsp"

for (i in 1:length(final[,1])) {
	sub<-subset(int,int$spsp==final$spsp[i])
	final$num.ITpairs[i]<-length(sub$sp.pair)
	final$num.locations[i] <-sum(sub$locations)
	final$num.routes[i] <-sum(sub$num.routes)
	final$tot.stops[i] <-sum(sub$tot.stops)
	final$tot.focal[i]<-sum(sub$tot.focal)
	final$tot.other[i]<-sum(sub$tot.other)
	final$tot.both[i]<-sum(sub$tot.both)
	final$syntopy1w[i]<-mean(sub$syntopy1w)
	final$syntopy2w[i]<-mean(sub$syntopy2w)
	final$syntopy1u[i]<-mean(sub$syntopy1u)
    final$syntopy2u[i]<-mean(sub$syntopy2u)}		

write.table(final, file="~/BBS/Syntopy.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

