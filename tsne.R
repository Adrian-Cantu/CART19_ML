library(tidyverse)
library(Rtsne)
#library(sigmoid)

raw_data <- readRDS('result_full_ALL_CLL_ML.rds')
unlearn_features<- c("specimen","cancer","CellType","Patient","Timepoint","Response_class")
full_data <- raw_data %>% drop_na()
full_data_x <- full_data %>% select(-all_of(c(unlearn_features,'pctSameOrt')))
full_data_y <- full_data %>% select(all_of(unlearn_features))
t_sne <- Rtsne(full_data_x)

df_tmp <- data.frame(
  X=t_sne$Y[,1],
  Y=t_sne$Y[,2]
)

df <- cbind(df_tmp,full_data_y)

fig_c <- ggplot(df, aes(x=X, y=Y, color=cancer)) +
  geom_point() + geom_rug()
ggsave('fig/tsne_c.pdf',plot=fig_c)

fig_res <- ggplot(df, aes(x=X, y=Y, color=Response_class)) +
  geom_point() + geom_rug()
ggsave('fig/tsne_resp.pdf',plot=fig_res)

timepoint_to_days <- function(x) {
  ret <- 0
  if (startsWith(x,'D')) {
    ret <- as.numeric(sub('^D','',x))
  } else if (startsWith(x,'M')) {
    ret <- as.numeric(sub('^M','',x))*30
  } else if (startsWith(x,'Y')) {
    ret <- as.numeric(sub('^Y','',x))*365
  } else if (startsWith(x,'d')) {
    ret <- as.numeric(sub('^d','',x))
  } else if (startsWith(x,'m')) {
    ret <- as.numeric(sub('^m','',x))*30
  } else if (startsWith(x,'y')) {
    ret <- as.numeric(sub('^y','',x))*365
  } 
  return(ret)
}

days_to_range <- function(x) {
  if (x==0){
    return('D0')
  } else if (x<=30){
    return('<1M')
  } else if (x<=180){
    return('1M~6M')
  } else if (x<=366){
    return('6M~1Y')
  } else {
    return('+1Y')
  }
}

timepoint_num <- sapply(full_data_y$Timepoint,timepoint_to_days)
names(timepoint_num) <- NULL
df$timepoint_num <- timepoint_num
fact_order <- c('D0','<1M','1M~6M','6M~1Y','+1Y')
df$timepoint_factor <-factor(sapply(timepoint_num,days_to_range),levels=fact_order)
#df$Timepoint <- timepoint_num/max(timepoint_num)
fig_time <- ggplot(df, aes(x=X, y=Y, color=as.factor(timepoint_num))) +
  geom_point() + geom_rug() +
  scale_color_viridis_d()
ggsave('fig/tsne_time.pdf',plot=fig_time)

fig_time_f <- ggplot(df, aes(x=X, y=Y, color=timepoint_factor)) +
  geom_point() + geom_rug() 
ggsave('fig/tsne_time_f.pdf',plot=fig_time_f)
saveRDS(df,'t_sne.rds')

