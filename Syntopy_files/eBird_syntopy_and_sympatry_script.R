#This script generates syntopy and sympatry estimates using data from eBird dataset, which can be requested from ebird.org

#load packages
library(sp)
library(rgeos)
library(maptools)
library(ggmap)
library(rgdal)
library(raster)
library(geosphere)
library(tools)

#open species pair list
taxa <- read.csv("Species_pairs.csv", header=TRUE)

#Create output variables if starting over with new species pair list
if(length(taxa$ebird.min.syntopy)==0) {
  taxa$ebird.min.syntopy<-NA
  taxa$ebird.sp1.syntopy<-NA
  taxa$ebird.sp2.syntopy<-NA
  taxa$ebird.mean.syntopy<-NA
  taxa$ebird.overall.syntopy<-NA
  taxa$ebird.sp1.sympatry<-NA
  taxa$ebird.sp2.sympatry<-NA
  taxa$ebird.min.sympatry<-NA
  taxa$ebird.mean.sympatry<-NA
  taxa$ebird.overall.sympatry<-NA
  
  	#Rename version of species names with spaces
	taxa$sp.1.sci_spc <- taxa$sp.1.sci
	taxa$sp.2.sci_spc <- taxa$sp.2.sci
	
	#Change spaces to underscore in species names to match eBird 
	taxa$sp.1.sci <- gsub(" ", "_", taxa$sp.1.sci)
	taxa$sp.2.sci <- gsub(" ", "_", taxa$sp.2.sci)
  }

for(i in 1:length(taxa$pair.n)){

    #Find eBird files	
    ebirdpath1U<-list.files(path="US_Canada", pattern=as.character(taxa$sp.1.sci[i]), full.names='TRUE')
    ebirdpath2U<-list.files(path="US_Canada", pattern=as.character(taxa$sp.2.sci[i]), full.names='TRUE')
    ebirdpath1M<-list.files(path="Mexico", pattern=as.character(taxa$sp.1.sci[i]), full.names='TRUE')
    ebirdpath2M<-list.files(path="Mexico", pattern=as.character(taxa$sp.2.sci[i]), full.names='TRUE')
    dummyfile<-"eBird_file_used_by_syntopy_script.txt"
    
    #Read dummy files in case some ebird files are missing or empty (allows for merging later)
    ebird1U <-read.delim(dummyfile, header=FALSE, sep="\t")
    ebird2U <-read.delim(dummyfile, header=FALSE, sep="\t")
    ebird1M <-read.delim(dummyfile, header=FALSE, sep="\t")
    ebird2M <-read.delim(dummyfile, header=FALSE, sep="\t")
    
    #Try to read eBird files
    try(ebird1U <- read.delim(ebirdpath1U, header=FALSE, sep="\t"),silent=TRUE)
    try(ebird2U <- read.delim(ebirdpath2U, header=FALSE, sep="\t"),silent=TRUE)
    try(ebird1M <- read.delim(ebirdpath1M, header=FALSE, sep="\t"),silent=TRUE)
    try(ebird2M <- read.delim(ebirdpath2M, header=FALSE, sep="\t"),silent=TRUE)

	#Remove uneeded columns
	ebird1U<-subset(ebird1U, select = c(V1,V5:V8))
	ebird2U<-subset(ebird2U, select = c(V1,V5:V8))
	ebird1M<-subset(ebird1M, select = c(V1,V5:V8))
	ebird2M<-subset(ebird2M, select = c(V1,V5:V8))

    #Remove US entries from Mexico files (e.g., Mexico, NY)
    ebird1M<-subset(ebird1M,V5=="MX")
    ebird2M<-subset(ebird2M,V5=="MX")
    
    #Merge US and Mexico files
    ebird1<-rbind(ebird1U,ebird1M)
    ebird2<-rbind(ebird2U,ebird2M)
    
    #Remove uneeded files from memory
    rm(ebird1U,ebird1M,ebird2U,ebird2M)
       
    #Clean up data
      #The eBird was subsetted using a pattern search, so the single species files can contain rows for other species if the focal species was mentioned in notes
      #This also removes unidentified species names like "Baeolophus bicolor/atricristatus" but retains hybrids iff the focal species' name comes before " x", as in "Baeolophus bicolor x atricristatus”
      #Currently, in the subsequent code, all "Baeolophus bicolor x“ are grouped with "Baeolophus bicolor” and, likewise, any "Baeolophus atricristatus x” are grouped with "Baeolophus atricristatus”.  
      #A potential problem is that if one person recorded an individual hybrid as "Baeolophus bicolor x atricristatus” and another person recorded it as "Baeolophus atricristatus x bicolor” in the same year,
      #this would be counted as an instance of the two species being found in syntopy (assuming the locations are <402 m apart).
    
    ebird1$sp <- as.character(ebird1$V1)
    ebird1<-subset(ebird1,sp==taxa$sp.1.sci_spc[i] | sp %in% grep(paste(taxa$sp.1.sci_spc[i],"x",sep=" "), ebird1$sp, value=TRUE))
    ebird2$sp <- as.character(ebird2$V1)
    ebird2<-subset(ebird2,sp==taxa$sp.2.sci_spc[i] | sp %in% grep(paste(taxa$sp.2.sci_spc[i],"x",sep=" "), ebird2$sp, value=TRUE))
    
    #Filter by peak breeding months
    first<-taxa$first_breeding_month[i] #the first month that substantially overlaps breeding season of both species and peak breeding of at least one of the two species, according to BNA
    last<-taxa$last_breeding_month[i] #the last month that substantially overlaps breeding season of both species and peak breeding of at least one of the two species, according to BNA
    
    ebird1$date <- as.character(ebird1$V8) 
    ebird1$date <- as.Date(ebird1$date, "%Y-%m-%d")
    ebird1$month<-format(ebird1$date, format="%m")
    ebird1$month<- as.numeric(ebird1$month) 
    ebird1<-subset(ebird1,(month>=first & month <=last))
    
    ebird2$date <- as.character(ebird2$V8) 
    ebird2$date <- as.Date(ebird2$date, "%Y-%m-%d")
    ebird2$month<-format(ebird2$date, format="%m")
    ebird2$month<- as.numeric(ebird2$month) 
    ebird2<-subset(ebird2,(month>=first & month <=last))

    
    if(length(ebird1)>0 & length(ebird2)>0) {#skips species pairs if one or both species have no eBird records 
    	 
    #Make list of unique years in ebird files
    ebird1$year<-format(ebird1$date, format="%Y")
    ebird2$year<-format(ebird2$date, format="%Y")
    sp1years<-unique(ebird1$year)
    sp2years<-unique(ebird2$year)
    years<-c(sp1years,sp2years)
    years<-unique(years)
    years<-as.numeric(years)
    years<-subset(years, years>1950)  #for ignoring older eBird records 
    years<-as.character(years)

    #Cycle through years; break p1 and p2 into subsets by year; remove locations within <50 m from both p1 and p2; count numnber of shared locations between p1 and p2
    #calculate syntopy metrics for each year; after finishing year loop, calculate mean syntopy across years
    
    #Create vectors
    sp1.syntopy<-rep(NA, length(years))
    sp2.syntopy<-rep(NA, length(years))
    mean.syntopy<-rep(NA, length(years))
    overall.syntopy<-rep(NA, length(years))
    sp1.total.count<-rep(NA, length(years))
    sp2.total.count<-rep(NA, length(years))
    sp1.sympatry.count<-rep(NA, length(years))
    sp2.sympatry.count<-rep(NA, length(years))    
    sp1.syntopy.count<-rep(NA, length(years))
    sp2.syntopy.count<-rep(NA, length(years))
    min.syntopy<-rep(NA, length(years))
    sp1.name<-rep(NA, length(years))
    sp2.name<-rep(NA, length(years))
    min.sympatry<-rep(NA, length(years))
    sp1.sympatry<-rep(NA, length(years))
    sp2.sympatry<-rep(NA, length(years))
    mean.sympatry<-rep(NA, length(years))
    overall.sympatry<-rep(NA, length(years))
    
    print("entering year loop")
    print(length(years))
    
#####################YEAR LOOP#############################    
    for(k in 1:length(years)) {
      p1<-subset(ebird1,year==years[k])
      p2<-subset(ebird2,year==years[k])
      
      #Get species names
    	sp1.name[k]<-as.character(taxa$sp.1.sci[i])
    	sp2.name[k]<-as.character(taxa$sp.2.sci[i])
    
    if(length(p1[,1])>0 & length(p2[,1])>0) {  #skip the rest of the year loop if either species has no records in year k
    
    #Round off coordinates to 3 decimal places to make locations within about 70 m of each other identical
    p1$V6<-round(p1$V6, digits=3)
    p1$V7<-round(p1$V7, digits=3)
    p2$V6<-round(p2$V6, digits=3)
    p2$V7<-round(p2$V7, digits=3)
    
	#Remove records with identical coordinates and drop uneeded columnns
    p1<-unique(p1[,c(3:4)])
    p2<-unique(p2[,c(3:4)])      
    p1$keep<-1
    p2$keep<-1
    
    #Convert into SpatialPointsDataFrames (required by distGeo)
    coordinates(p1) <- ~V7+V6
    coordinates(p2) <- ~V7+V6
    
    #Apply the datum (required by distGeo)
    crs(p1)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    crs(p2)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

    #Store the total counts of non-redundant records of each species in this year k
    N1<-length(p1)
    N2<-length(p2)
      
    #make a subset of records that are within within 24.5 miles of where the other species was recorded in the same year 
      #(24.5 miles is the length of a BBS route)
      n1<-0
      n2<-0
      if(length(p1)>0 & length(p2)>0) {
        for(j in 1:length(p1)) {
          if (min(distGeo(p1[j,],p2))>39428.93) {p1$keep[j]<-0}}
            
        for(j in 1:length(p2)) {
          if (min(distGeo(p2[j,],p1))>39428.93) {p2$keep[j]<-0}}
          
          p1<-subset(p1,keep==1)
          p2<-subset(p2,keep==1)
             #Store the counts of sympatric records
      			n1<-length(p1)
      			n2<-length(p2)
      }  
      
    #find the number of shared locations, within 1/4 mile, where both species were recorded, from the perspective of each species
      #(1/4 mile is the BBS stop search radius)
      #these represent the number of records where a species was found in syntopy
    
    sp1.shared<-0
    sp2.shared<-0
      
    if(n1>0 & n2>0) { #if either species has no records in sympatry then they can't be syntopic
    for(j in 1:n1) {
    if (min(distGeo(p1[j,],p2))<402.34) {sp1.shared<-sp1.shared+1}
    }

    for(j in 1:n2) {
          if (min(distGeo(p2[j,],p1))<402.34) {sp2.shared<-sp2.shared+1}
    }
      all.shared<-sp1.shared+sp2.shared
      }
    
    #calculate syntopy metrics
    if(n1+n2>10){ #(requires there to be at least 10 observations of the species in sympatry in a given year, which should screen out some bad records)
    	if(n1<n2){
        min.syntopy[k]<-sp1.shared/n1 
      } else {
        min.syntopy[k]<-sp2.shared/n2
      }
    sp1.syntopy[k]<-sp1.shared/n1 #proportion of species 1 records in syntopy
    sp2.syntopy[k]<-sp2.shared/n2 #proportion of species 2 records in syntopy
    mean.syntopy[k]<-(sp1.syntopy[k]+sp2.syntopy[k])/2  #mean of the species specific syntopy measures
    overall.syntopy[k]<-all.shared/(n1+n2) #proportion of all sympatric records in syntopy
    }
    
       #calculate eBird sympatry metrics
    if(N1>0 & N2>0){ #leave sympatry variables NA if one or both species was not recorded in year k, which leaves year k out of the sympatry calculation
      if(N1<N2){
        min.sympatry[k]<-n1/N1 
      } else {
        min.sympatry[k]<-n2/N2
      }
      sp1.sympatry[k]<-n1/N1
      sp2.sympatry[k]<-n2/N2
      mean.sympatry[k]<-(n1/N1+n2/N2)/2
      overall.sympatry[k]<-(n1+n2)/(N1+N2)
    }
    
    #Store counts for each year in vectors
    sp1.total.count[k]<-N1
    sp2.total.count[k]<-N2
    sp1.sympatry.count[k]<-n1
    sp2.sympatry.count[k]<-n2
    sp1.syntopy.count[k] <-sp1.shared
    sp2.syntopy.count[k] <-sp2.shared
   
print(k)  
} #closes first if in k loop
} #closes year loop

  #combine year vectors into a dataframe and save it, so that the script above doesn't have to be re-run to calculate different metrics
  
  year.data <- data.frame(sp1.name, sp2.name, years, sp1.total.count, sp2.total.count, sp1.sympatry.count, sp2.sympatry.count, sp1.syntopy.count, sp2.syntopy.count, min.syntopy, sp1.syntopy, sp2.syntopy, mean.syntopy, overall.syntopy, min.sympatry, sp1.sympatry, sp2.sympatry, mean.sympatry, overall.sympatry)
  
  if(i == 1){
  write.table(year.data, file="eBird_syntopy_counts.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)	
  } else {  
  write.table(year.data, file="eBird_syntopy_counts.csv", append=TRUE, quote=FALSE, sep=",", row.names=FALSE, col.names=FALSE)
  }
    
  taxa$ebird.min.syntopy[i]<-mean(min.syntopy, na.rm=TRUE)
  taxa$ebird.sp1.syntopy[i]<-mean(sp1.syntopy, na.rm=TRUE)
  taxa$ebird.sp2.syntopy[i]<-mean(sp2.syntopy, na.rm=TRUE)
  taxa$ebird.mean.syntopy[i]<-mean(mean.syntopy, na.rm=TRUE)
  taxa$ebird.overall.syntopy[i]<-mean(overall.syntopy, na.rm=TRUE)
  taxa$ebird.sp1.sympatry[i]<-mean(sp1.sympatry, na.rm=TRUE)
  taxa$ebird.sp2.sympatry[i]<-mean(sp2.sympatry, na.rm=TRUE)
  taxa$ebird.min.sympatry[i]<-mean(min.sympatry, na.rm=TRUE)
  taxa$ebird.mean.sympatry[i]<-mean(mean.sympatry, na.rm=TRUE)
  taxa$ebird.overall.sympatry[i]<-mean(overall.sympatry, na.rm=TRUE)
  
  print("Species pair")
  print(i)
  print(taxa$ebird.overall.syntopy[i])

  
} #closes if for skipping species pairs with no overlap of eBird records in region of overlap between breeding season shapefiles
#} #closes if for skipping species not recognized by BirdLife 
  
#save the data after finishing each species pair
write.table(taxa, file="eBird_syntopy.csv", append=FALSE, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)	
} #closes species pair loop
