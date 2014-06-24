# Set up the required packages
require(ggplot2)
require(RANN)


# 
LoadData <- function(data.location="C:/Users/Rasinj/Documents/Work/High-level analysis_tissue vs cellline/151-ml-m_ff/151-ml-m_44/results/Decoding/QT_0.4_0.01_details.csv"){
#   Loads the data set that should be used for analysis
#   The csv needs the following columns:
#     name: The names of the genes 
#     global_X_pos: The X-position of the transcript
#     global_Y_pos: The X-position of the transcript
#
#   Args: 
#     data.location: The total path to the csv-file
#   
#   Returns: 
#     A data frame that can be used for further analysis
  data<-read.csv(data.location)
  return(data)
}     
SubsetData<-function(data.orig,exclude=c("LGR5","ALCAM/CD166","CD44","ALDH1A","MSI1","CTNNB1","OCT4","SOX2","NNNN","RAC1")){
  #   Subsets the data set that is produced by LoadData.
  #   This is for excluding unwanted transcripts from the downstream analysis
  #   
  #   Args: 
  #     data.orig: The original, unsubsetted dataframe, produced by LoadData
  #     exclude: c("VGEF","NNNN") will exclude VGEF and undefined sequences
  #   Returns: 
  #     A subsetted data frame for downstream analysis
  
  indices.exclude<-which(data.orig$name%in% exclude)
  data.include<-data.orig[-indices.exclude,]
}
# Plot the transcripts per name
FacetPlot<-function(data){
#   Produces an overview plot of the different transcripts for a quick overview
#   
#   Args: 
#     data: a dataframe from LoadData or SubsetData
#     
#   Returns: 
#     A transcript-faceted plot of x and y position
  
  overview<-ggplot(data,aes(global_X_pos,global_Y_pos))+
    geom_point(size=0.5,alpha=0.3,aes(color=name))+
    facet_wrap(~name)+
    guides(color=FALSE)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  overview
  #ggsave(filename="overview.png",plot=overview,width=12,height=12)  
} 
GetNameList<-function(data=data.subset){
  #   Gets the list of unique transcript names for the data used in the analysis
  #   This is for downstream analysis by CalculateNearestNeighborNeighbor
  #   
  #   Args: 
  #     data: Any, subsetted of not data
  #     exclude: c("VGEF","NNNN") will exclude VGEF and undefined sequences
  #   Returns: 
  #     A subsetted data frame for downstream analysis
  
  namelist<-sort(unique(data$name))
}

# Loop through the names to check proximity for each transcript
CalculateNearestNeighbor<-function(data,namelist){
  data$closest<-c("")
  for (name in namelist){
    nnqueryids<-which(data$name==name)
    nndata<-data.frame(data[-nnqueryids,])
    nnquery.xy<-data.frame(data[nnqueryids,c(3,4)])
    nndata.xy<-data.frame(data[-nnqueryids,c(3,4)])
    
    # The nn2 calculates the nearest transcripts
    
    nnsd<-nn2(nndata.xy,query=nnquery.xy,k=2)
    closest<-nndata$name[nnsd$nn.idx[,1]]    
    data$closest[nnqueryids]<-as.character(closest)
  }
  data
}

# Plot the proximity matrix - not normalized for transcripts abundance
ProximityHeatmap<-function(data){
  g<-ggplot(data,aes(name,factor(closest)))+geom_point(size=3,alpha=0.2,color="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+ guides(color=FALSE) + #guide_legend(ncol=3))
  xlab(label="Registered transcript")+
  ylab(label="Closest transcript")
  g
  #ggsave(filename="C:/Users/Rasinj/Documents/Work/High-level analysis_tissue vs cellline/distance.png",plot=g,width=8,height=8)
}

# Plot the transcript counts
counts<-ggplot(data,aes(x=name))+geom_histogram()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
#ggsave(filename="C:/Users/Rasinj/Documents/Work/High-level analysis_tissue vs cellline/transcript_counts.png",plot=counts,width=8,height=8)

GetNormalizedPairwiseProximity<-function(data){
  names<-names(table(droplevels(data$name)))
  norm<-as.numeric(table(droplevels(data$name)))
  mat<-array(0,c(length(names),length(names)))
  for (i in seq_along(names)){
    tempdata<-data[data$name==names[i],]
    for (j in seq_along(names)){
      mat[i,j]=length(which(tempdata$closest==names[j]))/(norm[j]+norm[i])
    } 
  }
  df<- data.frame(mat)
  names(df)<-c(names)
  rownames(df)<-c(names)
  df$id<-names
  dfid<-melt(df,id.vars=c("id"))
}

NormalizedPairwiseProximityPlot<-function(normalizedPairwiseProximity,pointsize=4){
  normdistance<-ggplot(normalizedPairwiseProximity,aes(id,variable))+
    geom_point(size=pointsize,aes(color=round(value,2)*10+0.5))+
    scale_color_gradient(low="red",high="green")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
  normdistance
  #ggsave(filename="C:/Users/Rasinj/Documents/Work/High-level analysis_tissue vs cellline/distance_normalized.png",
  #       plot=normdistance,width=10,height=9)
}

# Main setting

data.orig<-LoadData()
data.subset<-SubsetData(data.orig,exclude=c("NNNN"))
FacetPlot(data.subset)
data.namelist<-droplevels(GetNameList(data.subset))
data.subset.nn<-CalculateNearestNeighbor(data.subset,data.namelist)
ProximityHeatmap(data.subset.nn)
normalizedPairwiseProximity<-GetNormalizedPairwiseProximity(data.subset.nn)
k<-names(rev(table(data.orig$name)[order(table(data.orig$name))]))
normalizedPairwiseProximity$id<-with(normalizedPairwiseProximity,
                                     factor(id,levels=k))
normalizedPairwiseProximity$variable<-with(normalizedPairwiseProximity,
                                     factor(variable,levels=k))
NormalizedPairwiseProximityPlot(normalizedPairwiseProximity,pointsize=3)

#overview 
normdistance