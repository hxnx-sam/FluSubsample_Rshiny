# Rcode_utils is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.
# FluSubsample_Rshiny is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Rcode_utils.  
# If not, see <http://www.gnu.org/licenses/>.

# functions to classify bird names to latin species names and countries to continents
# S. J. Lycett
# 14 Dec 2011 - extract from from matrixTranslate.R
# 3 April 2012 - added great cormorant to pel, and Bulgaria to Europe
# 23 May 2012 - added more birds
# 24 May 2012 - added more countries
# 13 July 2012 - changed egret from ans to cic
# 5 Dec 2012 - more countries
# 25 Mar 2013 - more birds and Iceland in Europe
# 5 April 2013 - more birds, and corrected pel and cic orders (egret is a pel, not cic)
# 9 April 2013 - added countryClass - change place to country
# 25 April 2013 - getWildDomestic
# 26 April 2013 - problem with grep and owl
# 1  May   2013 - domestic_green-winged_teal => Domestic
# 7  May   2013 - Wenzhou & Zhangzhou in China
# 8  May   2013 - Cuckoo & abbreviations (mostly for HKU)
# 9  May   2013 - CU = chukar
# 5  Aug   2013 - african stonechat
# 2  Dec   2013 - finch (pas)
# 12 Feb   2014 - dove (col)
# 11 Nov   2014 - from NS allele B work
# 17 Nov   2014 - from NS allele B work
# 22 Nov   2014 - H5N8 work, include Kurdistan in MiddleEast (parts are in Turkey & Iraq)
# 27 Nov   2014 - un geoscheme
# 15 Jan   2015 - crane (Gruiformes)
# 20 Feb   2015 - anas_ = anseriformes (e.g. Anas_americana, Anas_discors)
# 25 Mar   2015 - breeder_chicken = domestic
# 12 June  2015 - Chinese francolin = galliformes
# 23 June  2015 - black-legged kittywake = cha
# 17 July  2015 - added option on unGeoScheme to return the full list
# 20 July  2015 - horned puffin = cha
# 13 Oct   2015 - Bahrain, Oman in middle east, Paraguay in South America
# 15 Oct   2015 - WallisandFutuna, Yucatan
# 13 Nov   2015 - Anthropoides virgo this is a crane (Gru)
# 18 Nov   2015 - getLongShort
# 18 Nov   2015 - added some more places to countryClass (note this is not exhaustive, just what req for seqs)
# 25 Nov   2015 - Snowy Owl (stg)
# 10 Feb   2016 - more countries and places (H5N8 work using GISAID)
# 19 Feb   2016 - added more long range migrants for N5
# 20 May   2016 - updated path
# 6  Sept  2017 - updated for owl owl, anseriformes and phalacrocorax (H5N8 work)
# 18 Feb 2018 - common porchard = common pochard
# 19 Feb 2018 - Anser cygnoides = ans (Swan goose), bulbul, thrush = pas, condor = acc

######################################################################
# function to classify birds in H5N1 - to bird order

# Bird Orders

# Acc	Accipitriformes	http://en.wikipedia.org/wiki/Accipitriformes	Eagles, Hawks (falcons are separate)
# Ans	Anseriformes	http://tolweb.org/Anseriformes	Ducks, Geese, Swans
# Cha	Charadriiformes	http://tolweb.org/Charadriiformes	Shorebirds
# Cic	Ciconiiformes	http://tolweb.org/Ciconiiformes/26332	Storks, Herons
# Cor Coraciiformes	http://en.wikipedia.org/wiki/Coraciiformes - like pas.
# Col Columbiformes	http://en.wikipedia.org/wiki/Columbiformes, doves and pigeons
# Cuc Cuculiformes	http://en.wikipedia.org/wiki/Cuculiformes	Cuckoos
# Gal	Galliformes	http://tolweb.org/Galliformes	Fowl, Quail
# Gru	Gruiformes	http://en.wikipedia.org/wiki/Gruiformes	Coot, Moorhen
# Fal	Falconiformes	http://tolweb.org/Falconiformes/26379	Falcons
# Pas	Passeriformes	http://tolweb.org/Passeriformes	Perching birds (e.g Sparrows)
# Pel Pelecaniforms	http://en.wikipedia.org/wiki/Pelecaniformes Pelican
# Pic Piciformes		http://en.wikipedia.org/wiki/Piciformes - like pas, includes woodpeckers
# Pro	Procellariiformes	http://en.wikipedia.org/wiki/Procellariiformes	Shearwater
# Psi Psittaciformes	http://en.wikipedia.org/wiki/Parrot	Parrot
# Rhe	Rheiformes	http://en.wikipedia.org/wiki/Rhea_(bird)	Rhea (looks like an emu)
# Str	Struthioniformes	http://en.wikipedia.org/wiki/Ostrich	Ostrich, Emu
# Sph	Sphenisciformes	http://en.wikipedia.org/wiki/Penguin	Penguin
# Tin Tinamiformes	http://en.wikipedia.org/wiki/Red-winged_Tinamou Tinamou (Brazil)

# birdNames is a vector

# updated 25 Mar 2013
# updated 26 April 2013
# updated 9  May   2013
# updated 12 Feb   2014
# updated 11 Nov   2014
# updated 22 Nov   2014 - Ch = chicken, Larus michahellis (yellow-legged gull) = cha
# updated 20 Feb	 2015 - "anas_" = anseriformes (e.g. Anas_americana, Anas_discors)
# updated 12 June 2015 - Chinese francolin = galliformes
# updated 24 June 2015 - black-legged kittywake = cha
# updated 20 July 2015 - horned puffin = cha
# updated 13 Nov  2015 - Anthropoides virgo = gru
# updated 25 Nov  2015 - Strigiformes = stg, Snowy Owl
# updated 29 Nov  2016 - cockatoo = psi
# updated 9 Dec 2016 - 	pte = Pteroclidiformes, Syrrhaptes_paradoxus = sand grouse
# updated 12 june 2017 - stg owls
# updated 14 july 2017 - apo - swiftlet
# updated 6 Sept 2017 - owl owl = stg, anseriformes = ans, Phalacrocorax = cormorant = pel
# updated 19 Feb 2018
birdClass	<- function( birdNames ) {

	birdNames	<- tolower(birdNames)

	# nasty coding of ck for chicken messes up duck !!
	inds		<- which((birdNames=="ck") | (birdNames=="chick") | (birdNames=="ch") )
	if (length(inds) > 0) {
		birdNames[inds] <- "chicken"
	}	

	inds		<- which(birdNames=="gf")
	if (length(inds) > 0) {
		birdNames[inds] <- "guineafowl"
	}

	inds		<- which(birdNames=="ph")
	if (length(inds) > 0) {
		birdNames[inds] <- "pheasant"
	}

	inds		<- which(birdNames=="qa")
	if (length(inds) > 0) {
		birdNames[inds] <- "quail"
	}

	inds		<- which( (birdNames=="mllard") | (birdNames=="msllard") )
	if (length(inds) > 0) {
		birdNames[inds] <- "mallard"
	}

	# crowned gets synoned with crow - so remove all crowneds

	#inds		<- grep("crowned", birdNames)
	#if (length(inds) > 0) {
	#	for (i in 1:length(inds)) {
	#		birdNames[inds[i]] = sub("crowned", "", birdNames[inds[i]])
	#		}
	#}

	# quicker ! not that it matters terribly much, but still
	birdNames <- gsub("crowned", "", birdNames)

	# owl might get synoned with other things (e.g. fowl) so do it double
	inds		<- which(birdNames=="owl")
	if (length(inds) > 0) {
		birdNames[inds] <- "owl owl"
	}

	birdOrders	<- array("-", length(birdNames))

	apo   <- c("swiftlet")

	acc	<- c("eagle","hawk","buzzard","bussard","hooded vulture","condor")

	ans	<- c("dk","duck","mallard","teal","goose","gs","swan","cygnus","wigeon",
			"munia","grebe","goldeneye","pochard","goosander",
			"waterfowl","pintail","gadwall","garganey","platyrhynchos","crecca",
			"plathyrhynchos","shoveler","aquatic bird","aquaticbird",
			"widgeon","bufflehead","eider","scoter","redhead","water fowl","watercock","anas acuta",
			"greater scaup","mergus albellus","smew",
			"canvasback","hooded merganser","lesser scaup","northern shoverl","brant",
			"anas angustirostris","anas querquedula","ruddy sheldrake","scaup",
			"anas_","anas ","anser fabalis","shoveller","tadorna","aythya fuligula","dendrocygna viduata","mergus albellus",
			"anseriformes","muscovy",
			"common porchard","anser cygnoides")

	cic	<- c("stork","open-billed ibis","ibis") #19 feb 2018 - why is condor in here ? #c("stork","condor","open-billed ibis","ibis")

	cor	<- c("rollers")

	col	<- c("peaceful dove","rock dove", "dove")

	cha	<- c("shorebird","gull","tern","turnstone","dunlin","sandpiper","sanderling",
			"red knot","redknot","knot","plover","red-necked stint","redneckedstint","curlew",
			"rufous-necked stint","red-necked_stint","rufous-necked_stint","black-legged kittiwake","common murre","murre",
			"arenaria interpres","arenaria_interpres","thick-billed murre","arenaria-interpres",
			"eurasian woodcock","lapwing","oystercatcher","common snipe","snipe",
			"guillemot","larus argentatus","razorbill","larus michahellis","yellow-legged gull",
			"larus","bar-tailed godwit","black-tailed godwit","black-legged kittywake","horned puffin","puffin",
			"common redshank","black-necked stilt")

	cuc	<- c("cuckoo")

	gal	<- c("chicken","pheasant","quail","turkey","sck","guinea","guineafowl","guinea_fowl",
			"chukar","bustard","peacock","peahen","peafowl","partridge",
			"poultry","chukkar","chukka","silky fowl","bantam","rooster","numida meleagris",
			"chinese francolin","pavo cristatus")

	gru	<- c("coot","moorhen","hooded crane","black-neck crane","crane","anthropoides virgo","gallinula chloropus")

	fal	<- c("falcon","kestrel","peregrine","harrier")

	pas	<- c("sparrow","crow","magpie","pigeon","blackbird","myna",
			"starling","wild bird","wildbird","shrike","robin","softbill",
			"japanese white eye","japanesewhiteeye","common iora","commoniora",
			"fairy bluebird","fairybluebird","rook",
			"babbler","black bulbul","golden mountain thrush",
			"japanese white-eye","japanese_white-eye","silver-eared mesia","brambling","chinese hwamei",
			"african stonechat","finch","raven","barn swallow","barnswallow","swallow",
			"alder flycatcher","black-headed weaver","scaly thrush","campylorhamphus pucherani","copsychus saularis",
			"bulbul","thrush")

	psi	<- c("parrot","conure","parakeet","psittacine","macaw","yellow-headed amazon","budgerigar","cockatoo")

	pro	<- c("shearwater")	
	
	rhe	<- c("rhea")

	str	<- c("ostrich","emu")

	sph	<- c("penguin")

	stg	<- c("snowy owl","snowy_owl","snowyowl","long-eared owl","ural owl","eurasian eagel owl","owl owl")

	tin	<- c("red-winged tinamou", "redwingedtinamou")

	pel	<- c("pelican","great cormorant","greatcormorant", "cormorant", "phalacrocorax",
			"ardea cinerea", "grey heron", "egret", "heron","ardea cinerea")

	pic	<- c("great barbet")

	pte	<- c("syrrhaptes paradoxus","sandgrouse")

	orders <- list(	acc=acc,ans=ans,apo=apo,cic=cic,cor=cor,col=col,
				cha=cha,cuc=cuc,gru=gru,gal=gal,fal=fal,
				pas=pas,pel=pel,pro=pro,pic=pic,pte=pte,rhe=rhe,str=str,sph=sph,stg=stg,psi=psi,tin=tin)

	for (i in 1:length(orders)) {
		oo	<- unlist(orders[i])
		for (j in 1:length(oo)) {
			if (oo[j] != "owl") {
				inds	<- grep(oo[j], birdNames)
			} else {
				inds  <- which(birdNames==oo[j])
			}
			if (length(inds) > 0) {
				birdOrders[inds] <- attributes(orders[i])$names
			}
		}
	}	

	inds	<- grep("feces", birdNames)
	if (length(inds) > 0) {
		birdOrders[inds] <- "-"
	}

	return( birdOrders )

}

# 2 Sept 2013 - added turkey (not nec. previously), and muscovy_duck as Domestic
# see also http://www.ncbi.nlm.nih.gov/pubmed/20521681 (Qinghai lake and teals)
# During 2007-08 we marked wild ducks at Poyang Lake with satellite transmitters to examine the location and timing of spring migration 
# and identify any spatiotemporal relationship with HPAI H5N1 outbreaks. 
# Species included the 
# Eurasian wigeon (Anas penelope), 
# northern pintail (Anas acuta), 
# common teal (Anas crecca), 
# falcated teal (Anas falcata), 
# Baikal teal (Anas formosa), 
# mallard (Anas platyrhynchos), 
# garganey (Anas querquedula), and 
# Chinese spotbill (Anas poecilohyncha). 
# These wild ducks (excluding the resident mallard and Chinese spotbill ducks) ..

# updates 17 nov 2014
# updates 22 nov 2014 - Ch = chicken
# updates 25 mar 2015 - breeder_chicken => domestic
# updated 6  Sept 2017 - mulard ducks are domestic
# updated 19 Feb 2018 - Anser cygnoides are probably domestic (Swan goose)
getWildDomestic <- function( bird ) {
	bird2			<- tolower(bird)
	bird2			<- gsub(" ","_",bird2)			# 2 Sept 2013

	wildDomestic	<- array("Wild", length(bird2))
	inds			<- which( 	(bird2=="ck") | 
					(bird2=="ch") |
					(bird2=="chicken") | 
					(bird2=="guinea_fowl") | 
					(bird2=="guineafowl") |
					(bird2=="sck") |
					(bird2=="silky_chicken") |
					(bird2=="silkie_chicken") |
					(bird2=="silkie") |
					(bird2=="silky_fowl") |
					(bird2=="fowl") |
					(bird2=="poultry") |
					(bird2=="quail") |
					(bird2=="gs") |
					(bird2=="goose") |
					(bird2=="dk") |
					(bird2=="duck") |
					(bird2=="pheasant") |
					(bird2=="ph") |
					(bird2=="partridge") |
					(bird2=="Gf") |
					(bird2=="turkey") |
					(bird2=="village_chicken") |
					(bird2=="muscovy_duck") | 
					(bird2=="muscovy") |
					(bird2=="bantam") |
					(bird2=="breeder_duck") |
				      (bird2=="broiler_duck") |
					(bird2=="rooster") |
					(bird2=="breeder_chicken") |
					(bird2=="mulard") |
					(bird2=="mulard_duck") |
					(bird2=="korean_native_chicken") |
					(bird2=="bronze_turkey") |
					(bird2=="black_chicken") |
					(bird2=="anser_cygnoides"))

	wildDomestic[inds] <- "Domestic"
	inds			 <- which(bird2=="mammal")
	wildDomestic[inds] <- "Domestic"
	inds			 <- which(bird2=="human")
	wildDomestic[inds] <- "Domestic"
	inds			 <- which(bird2=="environment")
	wildDomestic[inds] <- "Domestic"
	inds			 <- grep("domestic",bird2)
	wildDomestic[inds] <- "Domestic"

	inds			 <- grep("wild",bird2)
	wildDomestic[inds] <- "Wild"

	return (wildDomestic) 
}

unabbreviatedBird <- function( bird ) {
	bird2		<- toupper(bird)
	abbrev	<- c("CK","CU","DK","GS","PG","PH","SCK","WDK","WWF")
	unabbrev	<- c("chicken","chukar","duck","goose","pigeon","pheasant","chicken","wild duck","wild water fowl")
	minds		<- match(bird2, abbrev)
	jj		<- which(is.finite(minds))
	ubird		<- bird
	ubird[jj]	<- unabbrev[minds[jj]]
	return( ubird )
}

# function to classify birds into long-range or short-range (or sedentary) species
# really only applicable for the wild birds anseriformes, so should run getWildDomestic first
# 18 nov 2015

getLongShort <- function( bird ) {

#Dear Sam, 
#for your information:
#- the "large" distance migratory birds in which HPAI H5N8 virus was detected in Korea/Japan are: baikal teal, bean goose, greater white fronted goose, tundra swan, common teal
#- the "large" distance migratory birds in which HPAI H5N8 virus was detected in Europe are: Eurasian wigeon, common teal
#- the "large" distance migratory birds in which HPAI H5N8/2 virus was detected in the US are: American green-winged teal, northern pintail, American wigeon
#Please let me know if you need additional information.
#Best regards, Rogier

	long <- c(	"baikal teal",
			"bean goose",
			"greater white fronted goose",
			"tundra swan",
			"eurasian wigeon",
			"common teal",
			"american green-winged teal",
			"northern pintail",
			"american wigeon",
			"wigeon",
			"white fronted goose",
			"anas crecca" )

	# 1 feb 2016 - SJL added due to more USA sequences in larger data set (seg6)
	long <- c(long,"blue-winged teal","snow goose","green-winged teal")

	# 19 feb 2016 - SJL added due to more sequences in N5 data set (seg6)
	long <- c(long,"bewicks swan","greylag goose","migratory duck")

	# 21 feb 2016 - SJL widgeon is same as wigeon
	long <- c(long, "widgeon")


	long2		<- gsub(" ","",long)
	long2		<- gsub("-","",long2)

	bird2		<- gsub(" ","",bird)
	bird2		<- gsub("_","",bird2)
	bird2		<- gsub("-","",bird2)
	bird2		<- tolower(bird2)

	minds		<- match(bird2,long2)
	jj		<- which(is.finite(minds))
	all( bird2[jj] == long2[minds[jj]] )
	mig		<- array("Short",length(bird2))
	mig[jj] 	<- "Long"

	return( mig )
}


#####################################################################################################
# function to classify countries to continents

# check this before use - it only has limited values
# updated 16 Nov 2011
# updated 31 Jan 2012
# updated 14 Feb 2012
# updated 17 Feb 2012
# updated 24 May 2012
# updated 9  Jul 2012
# updated 5  Dec 2012
# updated 25 Mar 2013
# updated 17 Nov 2014
# updated 22 Nov 2014
# updated 26 Nov 2014 - scale="main" => 7 continents, scale="sub" => finer scale

continentClass	<- function( countries, scale="main" ) {
	Europe		<- c(	"UK","France","Austria","Belgium","Bolivia","CzechRepublic","Denmark","Estonia",
					"Germany","Greece","Italy","Luxembourg","Netherlands",
					"Norway","Poland","Turkey","Spain","Serbia",
					"UnitedKingdom","Ireland","Hungary","Sweden","Portugal",
					"Slovenia","Croatia","Herzegovina", "Bulgaria",
					"Switzerland","BosniaandHerzegovina","Slovakia","Finland","Iceland",
					"Romania")

	NorthAmerica	<- c("USA","Canada","Mexico")
	CentralAmerica	<- c("ElSalvador","Nicaragua","Haiti","DominicanRepublic",
					"Guatemala","Panama","PuertoRico","Cuba","Barbados","TrinidadandTobago")

	SouthAmerica	<- c("Chile","Colombia","Argentina","Ecuador","Peru","Brazil","Uruguay","Venezuela","Paraguay")

	Africa		<- c("Ethiopia","Nigeria","Mali","Senegal","Djibouti",
					"Africa","Sudan","Egypt","BurkinaFaso", "Niger", "Kenya",
					"IvoryCoast", "CotedIvoire", "Zambia", "SouthAfrica","Ghana","Uganda",
					"Morocco","Tunisia","Togo","Cameroon","Algeria", "Benin","Libya")

	Asia			<- c("India","Malaysia","HongKong","Hong_Kong","China","Singapore","SouthKorea","Thailand","Japan",
					"Cambodia","Taiwan","Myanmar","VietNam","Viet_Nam","Vietnam","Pakistan",
					"Bangladesh","Bhutan", "Indonesia", "Laos","Korea","Malaya",
					"Nepal","Philippines","SriLanka","NorthKorea")

	Oceania		<- c("Australia","NewZealand","New_Zealand","Guam","NewCaledonia","Tonga","Fiji")

	MiddleEast		<- c("Jordan","Israel","Afghanistan","Arabia",
					"Kuwait","GazaStrip", "Iran", "Iraq", "SaudiArabia", "Lebanon",
					"UnitedArabEmirates","Qatar","Azerbaijan","MiddleEast","Kurdistan","PalestinianTerritory","Bahrain","Oman")

	CentralAsia		<- c(	"Kyrgyzstan","Turkmenistan","Russia","USSR",
					"Mongolia","Kazakhstan","Ukraine","Uzbekistan",
					"Georgia","Crimea")

	Antarctica		<- c("Antarctica")

	continents		<- array("-", length(countries))

	if (scale=="main") {
		conts			<- list( 	Europe=Europe, NorthAmerica=NorthAmerica,
						SouthAmerica=c(CentralAmerica,SouthAmerica),
						Africa=Africa, Oceania=Oceania, 
						Asia=c(Asia,CentralAsia,MiddleEast),
						Antarctica=Antarctica )

	} else {
		conts			<- list( Europe=Europe, NorthAmerica=NorthAmerica,
						   SouthAmerica=SouthAmerica, CentralAmeria=CentralAmerica,
						   Africa=Africa, Oceania=Oceania,
						   Asia=Asia, CentralAsia=CentralAsia, MiddleEast=MiddleEast,
						   Antarctica=Antarctica )
	}

	for (i in 1:length(conts)) {
		oo 	<- unlist(conts[i])
		
		for (j in 1:length(oo)) {
			inds 	<- which(countries==oo[j])
			if (length(inds) > 0) {
				continents[inds] <- attributes(conts[i])$names
			}	 
		}
	}
	
	return (continents)
}



countryClass	<- function( places ) {

	USA   <- c( "Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut",
			"Delaware","Florida","Georgia","Hawaii","Idaho","Illinois",
			"Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts",
			"Michigan","Minnesota","Mississippi","Missouri",
			"Montana","Nebraska","Nevada","New Hampshire",
			"New Jersey","New Mexico","New York","North Carolina",
			"North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania",
			"Rhode Island","South Carolina","South Dakota","Tennessee",
			"Texas","Utah","Vermont","Virginia",
			"Washington","West Virginia","Wisconsin","Wyoming","Delaware Bay",
			"Interior_Alaska","Interior Alaska","Southcentral_Alaska","Southcentral Alaska","South Central Alaska")
	Usa2  <- gsub(" ", "_", USA)
	Usa3	<- gsub(" ", "", USA)
	Usa4  <- c("CA","DE","Memphis","MN","NJ","NY","LA")
	USA   <- sort(unique(c(USA,Usa2,Usa3,Usa4)))

	Mexico <- c("Jalisco")

	Canada <- c("Ontario","Quebec","Nova Scotia","New Brunswick",
			"Manitoba","British Columbia","Prince Edward Island","Saskatchewan",
			"Alberta","Newfoundland and Labrador","Newfoundland","Labrador","Edmonton")
	Canada2<- gsub(" ", "_", Canada)
	Canada3<- gsub(" ", "", Canada)
	Canada4<- c("ALB","BC")
	Canada <- sort(unique(c(Canada,Canada2,Canada3,Canada4)))

	Thailand <- c(	"Bangkok","Kanchanaburi","Nakhon Si Thammarat Province","Phra Nakhon Si Ayutthaya","Satun",
				"Amnat Charoen","Khon Kaen","Nan","Phrae","Sing Buri",
				"Ang Thong","Krabi","Narathiwat","Phuket","Sisaket",
				"Bueng Kan","Lampang","Nong Bua Lamphu","Prachinburi","Songkhla",
				"Buriram","Lamphun","Nong Khai","Prachuap Khiri Khan","Sukhothai",
				"Chachoengsao","Loei Province","Nonthaburi","Ranong","Suphan Buri",
				"Chainat","Lopburi Province","Pathum Thani","Ratchaburi","Surat Thani",
				"Chaiyaphum","Mae Hong Son","Pattani","Rayong","Surin",
				"Chanthaburi","Maha Sarakham","Phang Nga","Roi Et","Tak",
				"Chiang Mai","Mukdahan","Phatthalung","Sa Kaeo","Trang",
				"Chiang Rai","Nakhon Nayok","Phayao","Sakon Nakhon","Trat",
				"Chonburi","Nakhon Pathom","Phetchabun","Samut Prakan","Ubon Ratchathani",
				"Chumphon","Nakhon Phanom","Phetchaburi","Samut Sakhon","Udon Thani",
				"Kalasin","Nakhon Ratchasima","Phichit","Samut Songkhram","Uthai Thani",
				"Kamphaeng Phet","Nakhon Sawan","Phitsanulok","Saraburi","Uttaradit",
				"Yala","Yasothon",
				"Phathumthani","Prachinburi_Thailand","NaraThiwat","Nakorn-Patom","Kohn Kaen",
				"Nakhonsawan","Suphanburi" )
	Thailand2 <- gsub(" ", "_", Thailand)
	Thailand3 <- gsub(" ", "", Thailand)
	Thailand  <- sort(unique(c(Thailand,Thailand2,Thailand3)))

	Indonesia <- c("East Java","West Java","East Kalimantan")
	Indonesia2<- gsub(" ","_",Indonesia)
	Indonesia3<- gsub(" ", "", Indonesia)
	Indonesia <- c(Indonesia,Indonesia2,Indonesia3,"Java","Kalimantan")

	China <- c(	"Beijing","Tianjin","Hebei","Shanxi","Mongolia",
			"Liaoning","Jilin","Heilongjiang",
			"Shanghai","Jiangsu","Zhejiang","Anhui","Fujian",
			"Jiangxi","Jiang_Xi","Jiang Xi","Shandong",
			"Henan","Hubei","Hunan","Guangdong","Guangxi",
			"Hainan","Chongqing","Sichuan","Guizhou","Yunnan",
			"Tibet","Xizang","Shaanxi","Gansu",
			"Qinghai","Ningxia","Xinjiang","HongKong","Hong_Kong","HK",
			"Xianggang","Macau","Aomen","Taiwan",
			"Nanchang","Nanjing","Qianzhou","Anyang","Guiyang","Hongze","Huadong","Jiawang",
			"Shantou","Shenzhen","Shijiazhuang","Shuanggou","Taixing","Taizhou","Tongshan",
			"Wuxi","Xianghai","Xiangshui","Xigou","Xuzhou","Yangzhou","Yongcheng","ZhaLong","Dawang","Zibo",
			"jiyuan","Kaifeng","Qixian","Zhuhai","Hangzhou","Wenzhou","Zhangzhou","Guangzhou","Kuming",
			"Dongting","Eastern_China","EasternChina")

	Japan <- c("Kinai","Tokaido","Tosando","Hokurikudo","Sanindo","Sanyodo","Nankaido","Hokkaido","Kumamoto","Miyazaki","Yamagata",
			"Yamaguchi","Osaka","Kobe","Yokohama","Akita","Aomori","Chiba","Miyagi","Shiga","Shimane","Tsukuba","Kagoshima")

	Korea <- c("Geumgang","Nakdonggang","Seongdong","Shihwa","Hadoree")

	Vietnam <- c("Quang_Ninh","QuangNinh","Quanh_Ninh","Ninh_Binh","Nam_Dinh","Mong_Cai","Ha_Tay","Ca_Mau","Hanoi","Soc_Trang","Nha_Trang")
	Vietnam2 <- gsub("_"," ",Vietnam)
	Vietnam3 <- gsub("_","",Vietnam)
	Vietnam  <- sort(unique(c(Vietnam,Vietnam2,Vietnam3)))

	Myanmar <- c("Hmawbi","Thanatpin")

	Kazakhstan <- c("Aktau")

	Australia <- c("Victoria","Western Australia","South Australia","Queensland","New South Wales")
	Australia2<- gsub(" ","_",Australia)
	Australia3<- gsub(" ", "", Australia)
	Australia <- sort(unique(c(Australia,Australia2,Australia3)))

	Russia <- c("Astrakhan","Neimonggu","Primorie","Altai","Gurjev","Siberia","Sakha","Buryatiya","Chabarovsk","Burjatia")

	Germany <- c("Rugen","Berlin","Heinersdorf","Potsdam","Germany-MV","Germany-NI","Stralsund","Bavaria","Germany-RP")

	UnitedKingdom	<- c("England","Scotland")

	countries		<- array("-", length(places))

	conts			<- list( USA=USA, Canada=Canada, Mexico=Mexico, China=China, Japan=Japan, Thailand=Thailand,
					   Indonesia=Indonesia, Korea=Korea, Vietnam=Vietnam, Myanmar=Myanmar,
						Japan=Japan, Kazakhstan=Kazakhstan,
					   	Russia=Russia, Australia=Australia, Germany=Germany, UnitedKingdom=UnitedKingdom )

	for (i in 1:length(conts)) {
		oo 	<- unlist(conts[i])
		
		for (j in 1:length(oo)) {
			inds 	<- which(places==oo[j])
			if (length(inds) > 0) {
				countries[inds] <- attributes(conts[i])$names
			}	 
		}
	}
	
	inds			<- which(countries=="-")
	countries[inds] 	<- places[inds]

	return (countries)

}


isProvinceOfChina	<- function( place, return.places=FALSE ) {

	provinces <- c(	"Beijing","Tianjin","Hebei","Shanxi",
				"Mongolia","Liaoning","Jilin","Heilongjiang",
				"Shanghai","Jiangsu","Zhejiang","Anhui",
				"Fujian","Jiangxi","Shandong","Henan",
				"Hubei","Hunan","Guangdong",
				"Guangxi","Hainan","Chongqing","Sichuan",
				"Guizhou","Yunnan","Tibet","Xizang",
				"Shaanxi","Gansu","Qinghai","Ningxia",
				"Xinjiang","HongKong","Xianggang",
				"Macau","Aomen","Taiwan")
	#provinces 	<- tolower(provinces)

	place2 	<- tolower(place)
	place2	<- gsub("_", "", place2)
	place2	<- gsub(" ", "", place2)
	place2	<- gsub("-", "", place2)
	
	if (length(place2) > 1) {
		minds		<- match(place2, tolower(provinces))
		jj		<- which(is.finite(minds))
		inds		<- array(FALSE, length(place2))
		inds[jj]	<- TRUE

		if (return.places) {
			res 	  <- array("-",length(place))
			res[jj] <- provinces[minds[jj]]
			return ( res )
		} else {
			return( inds )
		}
	} else {
		inds		<- which(tolower(provinces)==place2)
		return ( length(inds) == 1 )
	}
}


# 2 Sept 2013 - updated
# 12 Feb 2014 - updated
# 8 Dec 2014  - updated
# 8 May 2015 - updated
# 4 Feb 2016 - updated
# 25 Mar 2016 - updated
getProvinceOfChina <- function( places ) {
	pc 		<- isProvinceOfChina( places )

	Anhui		<- c("Anhui","Chuzhou","Shitai")
	Guizhou 	<- c("Guiyang","Gui_Yang")
	Zhejiang	<- c("Hangzhou","Huadong","Wenzhou","Shaoxing","Ningbo","Huzhou","Jiaxing")
	HongKong	<- c("HK","Hong_Kong")

	Jiangsu 	<- c(	"Jiawang","Nanjing","Qianzhou",
				"Taixing","Taizhou","Tongshan","Wuxi","Xiangshui","Xuzhou","Hongze",
				"Yangzhou","Shuanggou","Xigou","Xuyi","Suzhou","Changzhou","Sheyang","Yanchen","Yan_chen","Suzhou")		# not sure about Shuanggou or Xigou

	Jiangxi 	<- c("Nanchang","Zhenjiang")
	Henan	  	<- c("Kaifeng","Qixian","Yongcheng","Shangqiu","Xuchang","Anyang")		# not sure about Qixian
	Yunnan  	<- c("Kuming","Kunming")
	Guangdong 	<- c("Shantou","Zhuhai","Shenzhen","Guangzhou","Huizhou","Dongguan")
	Hebei		<- c("Shijiazhuang")
	Hubei		<- c("Chibi","Wuhan")
	Shanghai	<- c("Xianghai")
	Fujian	<- c("Zhangzhou")
	Shandong	<- c("Zibo","Rizhao","Qingdao","Yantai","Jinan","Laiwu","Chiping")
	Shaanxi	<- c("Dawang")
	Guangxi	<- c("SanJiang","Sanjiang","San_Jiang")				# not sure about SanJiang
	Heilongjiang<- c("ZhaLong","Mudanjiang")

	Hunan		<- c("Donting","Dongting","Dongting_Lake","Dongting Lake","DongtingLake",
				"Changsha","Lengshuitan","HN")	# HN is from A/Dk/HN/5806/2003 - checked in paper

	Mongolia	<- c("Neimeng","Neimonggu","Neimenggu")		# not completely sure about these
	Liaoning	<- c("Shenyan","Sheny")					# just guessing Sheny==Shenyan
	Sichuan	<- c("Juxian")
	Gansu		<- c("Sunan")

	provinces	<- list( Anhui=Anhui, Guizhou=Guizhou, Zhejiang=Zhejiang, HongKong=HongKong,
				   Jiangsu=Jiangsu, Jiangxi=Jiangxi, Henan=Henan, Yunnan=Yunnan,
				   Guangdong=Guangdong, Hebei=Hebei, Hubei=Hubei, Shanghai=Shanghai, Fujian=Fujian,
				   Shandong=Shandong, Shaanxi=Shaanxi, 
				   Guangxi=Guangxi, Heilongjiang=Heilongjiang, Hunan=Hunan,
				   Mongolia=Mongolia, Liaoning=Liaoning, Sichuan=Sichuan, Gansu=Gansu )

	provs			 <- array("-", length(places))
	#provs[ which(pc) ] <- places[ which(pc) ]
	provs[ which(pc) ] <-  isProvinceOfChina( places[ which(pc) ], return.places=TRUE )
	for (i in 1:length(provinces)) {
		oo 	<- unlist(provinces[i])
		
		for (j in 1:length(oo)) {
			inds 	<- which(places==oo[j])
			if (length(inds) > 0) {
				provs[inds] <- attributes(provinces[i])$names
			}	 
		}
	}
	
	return( provs )

}

USA_state_abbrev <- function( names, type="toState",
					fname = "usa_state_abbreviations_nov_2014.txt", 
					path = "") {

	#names <- gsub("Interior_Alaska","Alaska",names)
	#names <- gsub("Interior Alaska","Alaska",names)
	#names <- gsub("Southcentral_Alaska","Alaska",names)
	#names <- gsub("Southcentral Alaska","Alaska",names)
	names <- gsub("Barrow","Alaska",names)
	names <- gsub("Delaware Bay","Delaware",names)
	names <- gsub("Delaware_Bay","Delaware",names)

	ainds <- grep("Alaska",names)
	if ( length(ainds) > 0 ) {
		names[ainds] <- "Alaska"
	}

	if (type=="toState") {
		k <- 1
	} else {
		k <- 2
	}

	tbl <- read.table( paste(path,fname,sep=""), header=TRUE, sep="\t")

	minds1 <- match(gsub("[ _]","",names), gsub(" ","",paste(tbl[,1]) ) )
	minds2 <- match(names,tbl[,2])

	states <- array("-",length(names))

	jj1	 <- which(is.finite(minds1))
	if ( length(jj1) > 0 ) {
		states[jj1] <- paste(tbl[minds1[jj1],k])
	}	

	jj2	 <- which(is.finite(minds2))
	if (length(jj2) > 0) {
		states[jj2] <- paste(tbl[minds2[jj2],k])
	}

	return( states )
}

Canada_state_abbrev <- function( names, type="toState",
						fname = "canada_state_abbreviations_nov_2014.txt",
						path = "") {

	names <- gsub("StJohns","Newfoundland and Labrador",names,fixed=TRUE)
	inds	<- which(names=="ALB")
	if (length(inds) > 0) {
		names[inds] <- "Alberta"
	}
	names <- gsub("Nunavet","Nunavut",names)

	if (type=="toState") {
		k <- 1
	} else {
		k <- 2
	}

	tbl <- as.matrix(read.table( paste(path,fname,sep=""), header=TRUE, sep="\t"))

	minds1 <- match(gsub("[ _]","",names), gsub(" ","",paste(tbl[,1]) ) )
	minds2 <- match(names,tbl[,2])

	states <- array("-",length(names))

	jj1	 <- which(is.finite(minds1))
	if ( length(jj1) > 0 ) {
		states[jj1] <- paste(tbl[minds1[jj1],k])
	}	

	jj2	 <- which(is.finite(minds2))
	if (length(jj2) > 0) {
		states[jj2] <- paste(tbl[minds2[jj2],k])
	}

	return( states )
}






# this is like continentClass, but uses UN GeoScheme

# WARNING ONLY FOR EUROPE, ASIA, AFRICA, AMERICA AT THE MOMENT (8 Dec 2014)
# http://en.wikipedia.org/wiki/United_Nations_geoscheme_for_Asia
# http://en.wikipedia.org/wiki/United_Nations_geoscheme_for_Europe
# http://en.wikipedia.org/wiki/United_Nations_geoscheme_for_Africa
# http://en.wikipedia.org/wiki/United_Nations_geoscheme_for_the_Americas
# http://en.wikipedia.org/wiki/United_Nations_geoscheme_for_Oceania
# http://millenniumindicators.un.org/unsd/methods/m49/m49regin.htm
# 17 july 2015 - added option to just return the entire list
# 17 Feb 2018 - problem with locale in R studio with Sao Tome and Principe
unGeoScheme <- function( countries, getFullList=FALSE ) {

	CentralAsia <- c(		"Kazakhstan",
    					"Kyrgyzstan",
    					"Tajikistan",
    					"Turkmenistan",
    					"Uzbekistan" )


	EasternAsia	<- c(		"China","Hong Kong","Macau",
					"Korea","South Korea","North Korea",
					"Japan",
    					"Mongolia",
    					"Taiwan" ) 


	SouthernAsia <- c(	"Afghanistan",
				    	"Bangladesh",
					"Bhutan",
					"India",
					"Iran",
					"Maldives",
					"Nepal",
					"Pakistan",
					"Sri Lanka" )


	SouthEasternAsia <- c(  "Brunei Darussalam","Brunei",
					"Cambodia",
					"Indonesia",
					"Lao","Laos",
					"Malaysia",
					"Myanmar","Burma",
					"Philippines",
					"Singapore",
					"Thailand",
					"Timor-Leste","Timor","East Timor",
					"Viet Nam" )


	WesternAsia		<- c(	"Armenia",
					"Azerbaijan",
    					"Bahrain",
    					"Cyprus",
    					"Georgia","Republic of Georgia",
    					"Iraq",
    					"Israel",
    					"Jordan",
    					"Kuwait",
    					"Lebanon",
    					"Oman",
    					"Qatar",
    					"Saudi Arabia",
    					"State of Palestine","PalestinianTerritory","GazaStrip","Gaza","Gaza Strip","PAR",
    					"Syrian Arab Republic","Syria","MiddleEast",
    					"Turkey",
    					"United Arab Emirates",
    					"Yemen" )


	EasternEurope	<- c(	"Belarus",
					"Bulgaria",
    					"Czech Republic",
    					"Hungary",
    					"Poland",
    					"Republic of Moldova","Moldova",
    					"Romania",
    					"Russian Federation","Russia","USSR",
					"Ukraine",
    					"Slovakia")

	NorthernEurope	<- c( "Aland Islands",
    					"Denmark",
    					"Estonia",
    					"Faroe Islands",
    					"Finland",
    					"Guernsey",
    					"Iceland",
    					"Ireland","Eire","Northern Ireland",
    					"Isle of Man",
    					"Jersey",
    					"Latvia",
    					"Lithuania",
    					"Norway",
    					"Sark",
    					"Svalbard and Jan Mayen Islands","Svalbard",
    					"Sweden",
    					"United Kingdom","UK","GB","Great Britain")


	SouthernEurope	<- c(	"Albania",
    					"Andorra",
    					"Bosnia and Herzegovina",
    					"Croatia",
    					"Cyprus",
    					"Gibraltar",
    					"Greece",
    					"Holy See",
    					"Italy",
    					"Malta",
    					"Montenegro",
    					"Portugal",
    					"San Marino",
    					"Serbia","Kosovo",
    					"Slovenia",
    					"Spain",
    					"The former Yugoslav Republic of Macedonia","Macedonia","Yugoslavia")

	WesternEurope	<- c(	"Austria",
    					"Belgium",
    					"France",
    					"Germany",
    					"Liechtenstein",
    					"Luxembourg",
    					"Monaco",
    					"Netherlands",
    					"Switzerland" )


	EasternAfrica	<- c(	"Burundi",
    					"Comoros",
    					"Djibouti",
    					"Eritrea",
    					"Ethiopia",
					"Kenya",
    					"Madagascar",
    					"Malawi",
    					"Mauritius",
    					"Mayotte",
					"Mozambique",
    					"Reunion",
    					"Rwanda",
    					"Seychelles",
    					"Somalia",
					"South Sudan",
    					"Uganda",
    					"United Republic of Tanzania","Tanzania",
    					"Zambia",
    					"Zimbabwe" )

	CentralAfrica <- c(	"Angola",
    					"Cameroon",
    					"Central African Republic",
    					"Chad",
					"Democratic Republic of the Congo","Congo","DRC",
    					"Equatorial Guinea",
    					"Gabon",
    					"Republic of the Congo",
    					"Sao Tome and Principe","Sao Tome" )


	NorthernAfrica <- c(	"Algeria",
					"Egypt",
    					"Libya",
    					"Morocco",
					"Sudan",
    					"Tunisia",
					"Western Sahara" )

	SouthernAfrica <- c(    "Botswana",
    					"Lesotho",
    					"Namibia",
    					"South Africa",
    					"Swaziland" )

	WesternAfrica <- c(    	"Benin",
    					"Burkina Faso",
    					"Cape Verde",
    					"Cote d'Ivoire","Cote dIvoire","Ivory Coast",
    					"Gambia",
    					"Ghana",
    					"Guinea",
    					"Guinea-Bissau",
    					"Liberia",
    					"Mali",
    					"Mauritania",
    					"Niger",
    					"Nigeria",
    					"Saint Helena",
    					"Senegal",
    					"Sierra Leone",
   		 			"Togo" )



	CaribbeanAmerica <- c(	"Anguilla",
    					"Antigua and Barbuda","Antigua","Barbuda",
    					"Aruba",
    					"The Bahamas","Bahamas",
    					"Barbados",
    					"Bonaire",
					"Sint Eustatius","Sint Eustatius and Saba","Saba",
    					"British Virgin Islands","Virgin Islands",
    					"Cayman Islands",
    					"Cuba",
    					"Curacao",
    					"Dominica",
    					"Dominican Republic",
    					"Grenada",
    					"Guadeloupe",
    					"Haiti",
    					"Jamaica",
    					"Martinique",
    					"Montserrat",
    					"Puerto Rico",
    					"Saint Barthelemy",
    					"Saint Kitts and Nevis","Saint Kitts","Nevis",
    					"Saint Lucia",
    					"Collectivity of Saint Martin","Saint Martin",
    					"Saint Vincent and the Grenadines","Saint Vincent","Grenadines",
    					"Sint Maarten",
    					"Trinidad and Tobago","Trinidad","Tobago",
    					"Turks and Caicos Islands",
    					"United States Virgin Islands" )

	CentralAmerica <- c(	"Belize",
    					"Costa Rica",
    					"El Salvador",
    					"Guatemala",
    					"Honduras",
    					"Mexico","Yucatan",
    					"Nicaragua",
    					"Panama" )

	SouthernAmerica <- c(	"Argentina",
    					"Bolivia",
    					"Brazil",
    					"Chile",
    					"Colombia",
    					"Ecuador",
    					"Falkland Islands","Malvinas",
    					"French Guiana",
    					"Guyana",
    					"Paraguay",
					"Peru",
    					"Suriname",
    					"Uruguay",
    					"Venezuela" )
	# note it says SouthAmerica not SouthernAmerica, not sure which is right

	NorthernAmerica <- c(	"Bermuda",
    					"Canada",
    					"Greenland",
    					"Saint Pierre and Miquelon","Saint Pierre","Miquelon",
    					"United States of America","USA","United States" )


	AustraliaNewZealandOceania <- c(	"Australia",
							"New Zealand",
							"Norfolk Island")


	MelanesiaOceania <- c(	"Fiji",
					"New Caledonia",
    					"Papua New Guinea",
    					"Solomon Islands",
    					"Vanuatu" )

	MicronesiaOceania <- c(	"Guam",
    					"Kiribati",
    					"Marshall Islands",
    					"Micronesia","Federated States of Micronesia",
    					"Nauru",
    					"Northern Mariana Islands",
    					"Palau" )

	PolynesiaOceania <- c(	"American Samoa","Samoa",
    					"Cook Islands",
    					"French Polynesia","Polynesia",
    					"Niue",
    					"Pitcairn",
    					"Tokelau",
    					"Tonga",
    					"Tuvalu",
    					"Wallis and Futuna Islands","WallisandFutuna" )


	Antarctica	<- c("Antarctica")


	continents		<- array("-", length(countries))

	conts			<- list( 	CentralAsia=CentralAsia, EasternAsia=EasternAsia, SouthernAsia=SouthernAsia,
						SouthEasternAsia=SouthEasternAsia, WesternAsia=WesternAsia,
						EasternEurope=EasternEurope, WesternEurope=WesternEurope, 
						SouthernEurope=SouthernEurope, NorthernEurope=NorthernEurope,
						EasternAfrica=EasternAfrica, CentralAfrica=CentralAfrica,
						SouthernAfrica=SouthernAfrica, WesternAfrica=WesternAfrica,
						NorthernAfrica=NorthernAfrica,
						CaribbeanAmerica=CaribbeanAmerica, CentralAmerica=CentralAmerica,
						SouthernAmerica=SouthernAmerica, NorthernAmerica=NorthernAmerica,
						AustraliaNewZealandOceania=AustraliaNewZealandOceania,
						MelanesiaOceania=MelanesiaOceania, MicronesiaOceania=MicronesiaOceania,
						PolynesiaOceania=PolynesiaOceania,
						Antarctica=Antarctica )

  if (getFullList) {
	return( conts )
  } else {

	countries2		<- tolower(countries)
	countries2		<- gsub(" ","",countries2)
	countries2		<- gsub("-","",countries2)

	for (i in 1:length(conts)) {
		oo 	<- unlist(conts[i])
		oo    <- tolower(oo)
		oo	<- gsub(" ","",oo)
		oo	<- gsub("-","",oo)
		
		for (j in 1:length(oo)) {
			inds 	<- which(countries2==oo[j])
			if (length(inds) > 0) {
				continents[inds] <- attributes(conts[i])$names
			}	 
		}
	}
	
	return (continents)

  }


}





