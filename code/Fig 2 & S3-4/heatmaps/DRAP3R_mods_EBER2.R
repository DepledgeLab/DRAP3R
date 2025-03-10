library(dplyr)
library(readr)
library(tidyr)
library(Biostrings)
library(tools)
library(ggplot2)
library(gdata)
library(stringr)
Sys.setenv(LAN="en")
setwd("c:/Users/depledgd/Desktop/code/Fig 2 & S3-4/heatmaps/")

##-------------FILE SELECTION---------------------##
all_files <- list.files("./EBER2_data", pattern = ".bed", full.names=T)
filepaths <- all_files[grepl("pseU", basename(all_files), ignore.case =TRUE)]
filepaths
df_list <- list()
##-------------FILE SELECTION---------------------##


##--------TRANSFORM FILES INTO DATA FRAMES--------##
#OPEN FILES
for (i in filepaths) {
  getname <- file_path_sans_ext(basename(i))
  if (file.info(i)$size == 0) {
    warning(paste("The file", getname, "is empty and will be skipped."))
    next  # Skip to the next file if the current file is empty
  }
  #Open file and store in data frame with name of file
  df <- as.data.frame(read.table(i,header = FALSE , sep="\t",stringsAsFactors=FALSE, quote="")) #opens bed file as dataframe
  df$DataBase <- getname #adds column for file name
  colnames(df) <- trimws(colnames(df))
  df <- df %>%
    unite("V10",V10:V18, sep = " ")%>%
    mutate("V12"=as.numeric("6956"), "V13"=as.numeric("7128"))%>%
    select(V1,V2,V3,V6,V10,V12,V13,DataBase)
    dataframe <- paste(getname)
    assign(dataframe,df)
  df_list[[getname]] <- df
}

#MERGE ALL DATA FRAMES  
merged_df <- bind_rows(df_list)
#GET INFO FROM THE DATA FRAMES AS IN DRAP3R_PseudoU.R
merged_df <- merged_df %>% separate(V10, c("coverage","Score","Nmod","Ncanonical", "Nother_mod","Ndelete","Nfail","Ndiff","Nno_call"), " ") 
merged_df <- merged_df %>% mutate(coverage = as.numeric(coverage)) %>% mutate(Nmod = as.numeric(Nmod)) %>% mutate(Ncanonical = as.numeric(Ncanonical))%>% mutate(Nother_mod =as.numeric(Nother_mod))%>% mutate(Ndiff =as.numeric(Ndiff))%>% mutate(Nfail =as.numeric(Nfail))
start=merged_df$V12
end= merged_df$V13
merged_df <- merged_df %>% 
  filter(V2 >= V12 & V3 <= V13)

#FILTERING
#Change to base 1
merged_df <- merged_df %>% mutate(V2=V2+1, V12=V12+1)
#Calculate: Modification Position
negative_strand <- "-"
merged_df <- merged_df %>%
  mutate(Mod_position = ifelse(V6 == negative_strand, V13 - V2+1, V2-V12+1)) %>% mutate(Length=V13-V12+1)%>%
  mutate(Nfail=as.numeric(Nfail))
merged_df <- merged_df%>%rename( "Score"= "mod_freq")
merged_df$highlight <- 1
merged_df$highlight <- ifelse(merged_df$coverage < 50, NA, merged_df$highlight )

#LABELLING
merged_df <- merged_df %>% mutate(DataBase = case_when(
  grepl("CRO-AP5", DataBase, ignore.case = TRUE) & grepl("0.98", DataBase, ignore.case = TRUE) ~ "CRO-AP5 - 0.98",
  grepl("CRO-AP5", DataBase, ignore.case = TRUE) & grepl("noFilt", DataBase, ignore.case = TRUE) ~ "CRO-AP5 - Unfiltered",
  grepl("EBER2-IVT", DataBase, ignore.case = TRUE) & grepl("0.98", DataBase, ignore.case = TRUE) ~ "EBER2-IVT-3 - 0.98",
  grepl("EBER2-IVT", DataBase, ignore.case = TRUE) & grepl("noFilt", DataBase, ignore.case = TRUE) ~ "EBER2-IVT3 - Unfiltered",
  TRUE ~ DataBase # Keep the original value if no match is found
))
count(merged_df,highlight)
##--------TRANSFORM FILES INTO DATA FRAMES--------##




##------------------DATA PREP --------------------##  
matrix <- merged_df %>% select(DataBase,Mod_position,Length, mod_freq, V6, highlight) %>% distinct() 
newdf <- matrix %>% select(DataBase,Length,V6) %>% distinct()
df1<- data.frame()
for (i in 1:nrow(newdf)){
  ids<- c(paste0(rep(newdf$DataBase[i], newdf$Length[i])))
  position <- c(paste0((1:newdf$Length[i])))
  temp<- data.frame(ids,position) #temporary data frame holding the iteration's data
  df1 <- rbind(df1,temp)}

df1$Score <- integer(nrow(df1)) #creates a vector of zeros with the same length as the number of rows in df1.
df1$Mod_freq <- integer(nrow(df1))

for (i in 1:nrow(matrix)){
  name <- matrix$DataBase[i]
  tx <- matrix$V14[i]
  highlight <- matrix$highlight[i]
  score_vector <- integer(nrow(df1)) # Create a temporary score vector
  temp_modfreq <- numeric(nrow(df1))
  for (j in 1:nrow(df1)){
    condition_matched <- (df1$ids[j]== name && df1$position[j]==matrix$Mod_position[i])
    #Update Score vector and mod freq:
    score_vector[j] <- ifelse(condition_matched, 1, 0) #If both conditions are met, 1 is assigned to the corresponding value in 'score_vector'. If not, 0 is assigned.
    temp_modfreq[j] <-ifelse(condition_matched, matrix$mod_freq[i], 0)
  }
  df1$Score = df1$Score | score_vector  # Use bitwise OR to persist 1s in Score column
  df1$Mod_freq = ifelse(score_vector, temp_modfreq, df1$Mod_freq)
}
df1$position <- as.numeric(df1$position)

##Discard modifications 
df1$Highlight <- "MADEUP_POSITION"
for (i in 1:nrow(df1)){
  name <- df1$ids[i]
  position <- df1$position[i]
  for (j in 1:nrow(matrix)){
    name_matrix <- matrix$DataBase[j]
    position_matrix <- matrix$Mod_position[j]
    if (name==name_matrix && position==position_matrix){
      df1$Highlight[i] = matrix$highlight[j]
    }
    
  }}

desired_order <- c( "CRO-AP5 - 0.98",   
                    "EBER2-IVT-3 - 0.98"  ,
                    "CRO-AP5 - Unfiltered",
                    "EBER2-IVT3 - Unfiltered")
##------------------DATA PREP --------------------##  



##---------------------PLOT ----------------------##  
df1$ids <- factor(df1$ids, levels=desired_order)
df1$Mod_freq <- as.numeric(df1$Mod_freq)
df1 <- df1 %>%
  mutate(color = case_when(
    #Highlight == "MADEUP_POSITION" ~ "seashell",          # Color for "Discarded"
    Highlight == "MADEUP_POSITION" ~ "white",
    is.na(Highlight) ~ "#fff",               # Color for "NA values"
    TRUE ~ scales::col_numeric(
      palette = c("#f0f5ff", "darkblue"), domain = c(0, 100)
    )(Mod_freq)                                # Gradient for numeric Mod_freq
  ))

# |  #INFO FOR COLOR PALLETE#   |
# | #min for darkred= #faebeb   |
# | #min for darkblue= #f0f5ff  |


#DISPLAY THE TX IN DIFFERENT ROWS:use the method you prefer 
#df1$position_group <- as.factor(ceiling(as.numeric(df1$position) / as.numeric(merged_df$Length[1]/3)))
df1$position_group <- as.factor(ceiling(as.numeric(df1$position) / 60))

#PLOTTING
pdf("NAME.pdf", width =22.96  , height = 5.5 )
ggplot(df1, aes(x = position, y = ids, colour=color)) +
  geom_tile(aes(fill = color), color = "white", lwd = 2, linetype = 1) +
  scale_fill_identity() +  # Use the precomputed color column directly
  ggtitle("m6A in EBER transcript") +
  xlab("mRNA Position") +
  scale_x_continuous(breaks = seq(0, max(df1$position, na.rm = TRUE), by = 10)) +
  ylab("Data Base") +
  facet_wrap(~ position_group, scales = "free_x", ncol = 1) +
  theme(
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust=0.5, margin=margin(b=20)),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.x = element_text(margin=margin(t=20)),
    axis.title.y =  element_text(margin=margin(r=20)),
    axis.text.y = element_text(margin = margin(r = -10), angle = 0, size = 12),
    axis.text.x = element_text(margin = margin(t = 0), angle = 0, hjust = 0.5, vjust=-0.3, size = 40),
    strip.text = element_blank(),
    legend.position = "right",
  )
dev.off()
##---------------------PLOT ----------------------##  

