library(tidyverse)
library(randomForest)

raw_data <- readRDS('result_full_ALL_CLL_ML.rds')
unlearn_features<- c("specimen","cancer","CellType","Patient","Timepoint","Response_class")

train_data <- raw_data %>% filter(Response_class!='MISS') %>% drop_na()
train_data_x <- train_data %>% select(-all_of(c(unlearn_features,'pctSameOrt')))
train_data_y <- train_data %>% select(all_of(unlearn_features))
train_response <- droplevels(train_data_y$Response_class)

test_data <- raw_data %>% filter(Response_class=='MISS') %>% drop_na()
test_data_x <- test_data %>% select(-all_of(c(unlearn_features,'pctSameOrt')))
test_data_y <- test_data %>% select(all_of(unlearn_features))

rf <- randomForest(x=train_data_x,y=train_response,ntree=600,
             importance = TRUE,proximity = TRUE,classwt=c(1,100))
rf
plot(rf)

top_rf <- rownames(rf$importance[order(rf$importance[,3],decreasing = TRUE),])[1:20]
rf2 <- randomForest(x=train_data_x %>% select(all_of(top_rf))
                      ,y=train_response,ntree=600,
                      importance = TRUE,proximity = TRUE)
rf2
saveRDS(rf2,'rf_resp.rds')
varImpPlot(rf2,type=1)
# predict cancer type ------------


full_data <- raw_data %>% drop_na()
full_data_x <- full_data %>% select(-all_of(c(unlearn_features,'pctSameOrt')))
full_data_y <- full_data %>% select(all_of(unlearn_features))

full_cancer <- as.factor(full_data_y$cancer)

rf_c <- randomForest(x=full_data_x,y=full_cancer,ntree=600,
                  importance = TRUE,proximity = TRUE)
rf_c
plot(rf_c)
varImpPlot(rf_c,type=1)
saveRDS(rf_c,'rf_can.rds')
