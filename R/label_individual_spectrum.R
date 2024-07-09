#' label_individual_spectrum
#'
#' Using the frequencies determined from the call_aggregate(), label individual spectrum
#' @param aggregate_list: list of aggregate spectra objects for all features (output from the call_aggregate())
#' @param folder_path: same path as the one used in preprocess()
#' @param mz_tol: mass tolerance (in Da) required when grouping fragments across multiple spectra for a given feature
#' @return a list containing all spectra for a given feature and the corresponding fragment frequencies
#' @examples
#' # Example usage of the function
#' label_individual_spectrum(l2, folder_path, 0.05)
#' @export

label_individual_spectrum <- function(aggregate_list, folder_path, mz_tol){
  df_sps_mean_agg_df = aggregate_list; freq_df = list()
  for(i in 1:length(aggregate_list)){
    metabolite = names(aggregate_list)[i]
    path = paste(folder_path, "MS2_scans_before_denoising/",metabolite,sep="")
    for(j in 1:length(list.files(path))){
      scan = read.delim(paste(path,"/",list.files(path)[j],sep=""))
      scan_number = list.files(path)[j]
      new_scan = dures::intrascan_grouping(scan, mz_tol)
      if(dim(new_scan)[1] == 0){
        new_scan = scan
      }
      scan_group = df_sps_mean_agg_df[[metabolite]]$Df
      scan_group$Mean_MZ = round(scan_group$Mean_MZ,2)
      #which(between(scan$fragments, scan_group$Mean_MZ - 0.05, scan_group$Mean_MZ + 0.05))
      f = data.frame(); mzd = c(); freqd = c(); d= 1
      #scan = scan[-c(which(scan$intensity == 0)),] no need for this. already doing this in the code before
      new_scan$unique_id=paste(rep("frag",length(new_scan$fragments)), seq(1,length(new_scan$fragments)), sep="_")
      new_scan$fragments=round(new_scan$fragments,2)
      new_scan$group=rep("", dim(new_scan)[1])
      new_scan$unique_id=factor(new_scan$unique_id)
      #For a given scan
      #If a particular fragment hasn't been assigned to a group yet,
      #look for all fragments within +-0.05 of this fragment's mz
      #value and assign a group to it, do this for all fragments
      #This stage is the grouping stage
      #soring mz values in decreasing order so that the highest mz value is chosen first

      new_scan=new_scan[order(new_scan$fragments, decreasing = T),]
      df_merge=data.frame()
      #for(i in 1:length(sps_mod)){ #was giving an error. Number of scans not equal to number of unique ids
      for(i in 1:length(levels(new_scan$unique_id))){
        df_sps=subset(new_scan, new_scan$unique_id == levels(new_scan$unique_id)[i])
        df_sps=df_sps[order(df_sps$fragments, decreasing = T),]

        for(j in 1:nrow(df_sps)){
          if(df_sps$group[j] == ""){
            irows=with(df_sps,fragments<= df_sps[j,1]+mz_tol & fragments >= df_sps[j,1]-mz_tol)
            df_sps$group[irows]=as.character(df_sps$fragments[j])
          }}
        df_merge=rbind(df_merge, df_sps)


      }
      #This stage will sum all fragments in a single group
      df_merge$mz_group=paste(df_merge$unique_id, df_merge$group, sep="_")
      #df_merge is the dataframe with the groups formed.
      #Grouping was done based on mz tolerance range within a given scan
      # print("Before intra-scan grouping")
      # print(dim(new_scan))
      #before_dim=c(before_dim, dim(df)[1])
      #For a given scan, average all mz values and sum all intensity values within 0.05 range of mz
      agg_df <- cbind(aggregate(df_merge$fragments, by=list(df_merge$mz_group), FUN=mean),
                      aggregate(df_merge$intensity, by=list(df_merge$mz_group), FUN=sum))
      #delet the third column
      agg_df[,3]=NULL
      colnames(agg_df)=c("Fragment_ID", "MZ", "Intensity")
      agg_df$Freq = 0
      # print("After intra-scan grouping")
      # print(dim(agg_df)[1])
      for(k in 1:length(agg_df$Fragment_ID)){
        mz = agg_df$MZ[k]
        mzd = c(mzd, mz)
        f = scan_group[which(dplyr::between(scan_group$Mean_MZ, mz - mz_tol, mz + mz_tol) == TRUE),
                       "Frequency"]
        if(length(f) == 0)
        {
          print(mz)
          print("Fragment of this scan couldn't be labelled")
          agg_df$Freq[k] = NA
          n1 = paste(metabolite, scan_number, mz, sep = "_")
          l = c(l, n1)
        }else{
          agg_df$Freq[k] = max(f)
        }
      }
      agg_df[which(agg_df$Freq < 0),"Freq"] = 0
      if(length(which(is.na(agg_df$Freq)) != 0)){
        agg_df = agg_df[-c(which(is.na(agg_df$Freq))),]
      }

      freq_df[[metabolite]][[scan_number]] = agg_df

    }
  }
  return(freq_df)

}
