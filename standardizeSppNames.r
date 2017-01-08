# Function to standardize species names across regions
# Originally from 5_select_species.r


standardizeSppNames <- function(spp){

	# see spptaxonomy_2015-02-09_plusManual.csv for a useful conversion table
	Spptax<-read.csv("data/spptaxonomy_plusManual.csv") #note: new column in CSV file 
	spptax<-apply(Spptax,2,tolower)
	sppu <- unique(spp)
	sppl<-tolower(sppu)
	datspp<-unique(sppl)

	spptax<-as.data.frame(spptax)
	sum(datspp %in% spptax$taxon)

	#Match when possible and assign new genus species, otherwise keep old name
	#Find taxa with matches in the taxonomy table
	#In some cases, the "name" field of taxonomy has spelling errors (missing "n"s). In this case, genus+species should be used. But in some cases, "genus species" are not given if full taxonomy is missing. The spptaxonomy.csv file now has a "newname" column which merges these optimally.
	matches<-datspp %in% spptax$taxon + datspp %in% spptax$name + datspp %in% spptax$common + datspp %in% spptax$genusspecies  #826 matches/4937

	#For those species datspp[matches>0], replace current name with spptax$newname
	newnames=character(length(datspp))
	for (i in 1:length(datspp)){
		if(matches[i]==0){
			newnames[i]<-datspp[i]
		} else if(datspp[i] %in% spptax$taxon){
			ind<-match(datspp[i],spptax$taxon)
			newnames[i]<-as.character(spptax$newname[ind])
		} else if(datspp[i] %in% spptax$name){
			ind<-match(datspp[i],spptax$name)
			newnames[i]<-as.character(spptax$newname[ind])
		} else if(datspp[i] %in% spptax$genusspecies){
			ind<-match(datspp[i],spptax$genusspecies)
			newnames[i]<-as.character(spptax$newname[ind])
		} else if(datspp[i] %in% spptax$common){ # is this needed? (Malin 5/4/2015)
			ind<-match(datspp[i],spptax$common)
			newnames[i]<-as.character(spptax$newname[ind])
		}
	}

	length(unique(datspp))
	length(unique(newnames)) # down to 4831 unique spp, from 4937

	oldnewnames<-cbind(sppl=datspp,sppnew=newnames)
	oldnewnames2 <- merge(data.frame(spp=sppu, sppl=sppl), oldnewnames)
	return(oldnewnames2[,c('spp', 'sppnew')])	#Can merge these new names with old names in dat file
}