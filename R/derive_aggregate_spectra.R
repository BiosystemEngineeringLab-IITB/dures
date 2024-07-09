#' derive_aggregate_spectra
#'
#' Create aggregate spectra by grouping fragments across multiple spectra for a given feature using a predefined tolerance and calculate the frequencies of fragments
#' @param sps: A spectra object containing top x% TIC spectra concatenated for a given feature.
#' @param mz_tol: mass tolerance (in Da) required when grouping fragments across multiple spectra for a given feature
#' @return a list containing the aggregate spectra as a Spectra object and a dataframe of the aggregate spectra with its mean MZ, Intensity and frequency and the number of fragments before and after inter- and intra-spectra grouping
#' @examples
#' # Example usage of the function
#' derive_aggregate_spectra(sps, 0.05)
#' @export
derive_aggregate_spectra <- function(sps, mz_tol=0.05){
  df=data.frame();
  c=0
  for(i in 1:length(sps)){
    #print(i)
    sps_2=sps[i]
    pks=Spectra::peaksData(sps_2)
    #print(dim(pks[[1]]))
    p=as.data.frame(pks[[1]])
    if(nrow(p) == 0){ #condition to check if dataframe is empty
      c=c+1 #no fragments are there in the scan
      idx = i
    }}
  if(c!=0){
    sps = sps[-idx] #remove that scan which doesn't have fragments
  }
  sps_mod=sps
  for(i in 1:length(sps_mod)){
    sps_2=sps_mod[i]
    pks=Spectra::peaksData(sps_2)
    p=as.data.frame(pks[[1]])
    if(dim(p)[[1]] == 0){
      next
    }
    #add a way to know which spectra (RT range) each row is coming from.
    p$rtime=as.character(round(sps_mod$rtime,1))[i]
    p$fileid=fileID
    df=rbind(df,p)
  }

  #sps = sps[-idx]
  # print("Number of scans with zero fragments")
  # print(c)
  #creating unique id for every scan
  df$unique_id=paste(df$fileid, df$rtime, sep="_")
  df= df[which(df$intensity != 0),] #removing all fragments with intensity 0
  df$mz=round(df$mz,2)
  df$group=rep("", dim(df)[1])
  df$unique_id=factor(df$unique_id)
  #For a given scan
  #If a particular fragment hasn't been assigned to a group yet,
  #look for all fragments within +-0.05 of this fragment's mz
  #value and assign a group to it, do this for all fragments
  #This stage is the grouping stage
  #soring mz values in decreasing order so that the highest mz value is chosen first

  df=df[order(df$mz, decreasing = T),]
  df_merge=data.frame()
  #for(i in 1:length(sps_mod)){ #was giving an error. Number of scans not equal to number of unique ids
  for(i in 1:length(levels(df$unique_id))){
    df_sps=subset(df, df$unique_id == levels(df$unique_id)[i])
    df_sps=df_sps[order(df_sps$mz, decreasing = T),]



    for(j in 1:nrow(df_sps)){
      if(df_sps$group[j] == ""){
        irows=with(df_sps,mz<= df_sps[j,1]+mz_tol & mz >= df_sps[j,1]-mz_tol)
        df_sps$group[irows]=as.character(df_sps$mz[j])
      }}
    df_merge=rbind(df_merge, df_sps)


  }
  #This stage will sum all fragments in a single group
  df_merge$mz_group=paste(df_merge$unique_id, df_merge$group, sep="_")
  #df_merge is the dataframe with the groups formed.
  #Grouping was done based on mz tolerance range within a given scan
  # print(paste("Number of fragments Before intra-scan grouping", dim(df)[1], )
  # print(dim(df))
  before_dim= dim(df)[1]
  #For a given scan, average all mz values and sum all intensity values within 0.05 range of mz
  agg_df <- cbind(aggregate(df_merge$mz, by=list(df_merge$mz_group), FUN=mean), aggregate(df_merge$intensity, by=list(df_merge$mz_group), FUN=sum))
  #delet the third column
  agg_df[,3]=NULL
  colnames(agg_df)=c("Fragment_ID", "MZ", "Intensity")
  # print("After intra-scan grouping")
  # print(dim(agg_df))
  #print("Percentage reduction in data")
  #print(1-dim(agg_df)[1]/dim(df_merge)[1])


  #a=df_table
  #a$mz=as.numeric(as.character(rownames(a)))
  a=agg_df
  a$group=rep("", dim(a)[1])
  for(j in 1:nrow(a)){
    if(a$group[j] == ""){
      irows=with(a,MZ<= a[j,"MZ"]+mz_tol & MZ >= a[j,"MZ"]-mz_tol) #third column in mz
      a$group[irows]=as.character(a$MZ[j])
    }}
  #a$mz_group=paste(a$rt_group, a$group, sep="_")
  agg_df_freq <- cbind(aggregate(a$MZ, by=list(a$group), FUN=mean),
                       aggregate(a$Intensity, by=list(a$group), FUN=mean))
  agg_df_freq[,3]=NULL
  colnames(agg_df_freq)=c("MZ_group", "Mean_MZ", "Mean_Intensity")

  for(i in 1:nrow(agg_df_freq)){
    y=a[which(a$group == agg_df_freq$MZ_group[i]),]
    mat=strsplit(y$Fragment_ID, "_")
    mat=matrix(unlist(mat),ncol=4,byrow=T) #expects input in the format create_consensus_table(sps, "IDA_1") not "IDA-1"
    agg_df_freq$Frequency[i] = length(unique(mat[,3]))

  }
  agg_df_freq$Frequency = agg_df_freq$Frequency / length(sps)

  agg_df_freq$MZ_group=as.numeric(as.character(agg_df_freq$MZ_group))


  #Visualize spectra
  #plotSpectra(sps, main = paste("FileID = ", fileID , "ScanIndex = ", sps$scanIndex, " , ", "RT = ", as.character(round(sps$rtime,2)), sep=""))

  # par(mfrow=c(1,1))
  #
  # pdf(paste(inchikey,"_raw_spectra.pdf", sep=""))
  # #overlay raw spectra
  # plotSpectraOverlay(sps, lwd = 1, main=paste(inchikey, "_Overlaid Raw spectra_", sep=""))
  # dev.off()

  #plotting the reduced spectra
  #first creating the spectra object
  spd <- DataFrame(
    msLevel = c(2L),
    polarity = c(1L),
    name = paste(fileID, range(sps$rtime)[1], range(sps$rtime)[2], sep="-"))

  spd$mz <- list(
    c(agg_df_freq$Mean_MZ))

  spd$intensity <- list(c(agg_df_freq$Mean_Intensity))
  #spd$Frequency = agg_df_freq$Frequency
  sps_mean <- Spectra::Spectra(spd)
  sps_mean$Frequency <- list(c(agg_df_freq$Frequency))


  mzLabel <- function(z) {
    pks <- Spectra::peaksData(z)[[1L]]
    p <- as.data.frame(pks)
    p[,3] = z$Frequency
    colnames(p)[3]="Frequency"
    lbls <- format(p[, "Frequency"], digits = 2)
    lbls[p[, "Frequency"] < 0.5] <- ""
    lbls
  }

  # pdf(paste(inchikey, "_consensus_spectra.pdf", sep=""))
  # plotSpectra(sps_mean, labels = mzLabel,
  #             labelSrt = -30, labelPos = 2, labelOffset = 0.1, main = inchikey)
  # dev.off()
  #
  #print("After inter-scan grouping")
  #print(dim(agg_df_freq)[1])
  after_dim= dim(agg_df_freq)[1]
  #print("Percentage reduction in data")
  #percent_red=c(percent_red, 1-dim(agg_df_freq)[1]/dim(df_merge)[1])
  #print(1-dim(agg_df_freq)[1]/dim(df_merge)[1])
  #agg_df_freq$Identifier=paste("Denoise_", name, sep="_")

  agg = list(Mean = sps_mean, Df = agg_df_freq, before_dim, after_dim)

  return(agg)

  #also return for every precursor mz, the percentage reduction stored in variables percent_red, after_dim, before_dim

}
