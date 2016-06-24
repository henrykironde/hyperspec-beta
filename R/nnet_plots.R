# checking neural network predictions
library(caret)
library(ggrepel)
source('R/fit_neural_net.R')

test_df$scientificname <- as.factor(test_df$scientificname)
test_df$pred <- as.factor(test_df$pred)
levels(test_df$pred) <- levels(test_df$scientificname)
cm <- confusionMatrix(test_df$pred, test_df$scientificname)

# plot commonness vs. sensitivity
commonness <- c(table(test_df$scientificname))
sensitivity <- cm$byClass[, "Sensitivity"]
specificity <- cm$byClass[, "Specificity"]

check_df <- data.frame(Species = names(commonness), 
                       commonness = commonness, 
                       sensitivity = sensitivity, 
                       specificity = specificity)

ggplot(check_df, aes(x = commonness, y = sensitivity)) + 
  geom_point(size = 2) + 
  theme_minimal() + 
  xlab('Commonness') + 
  ylab('Sensitivity')

ggplot(check_df, aes(x = commonness, y = specificity)) + 
  geom_point(size = 2)  + 
  theme_minimal() + 
  xlab('Commonness') + 
  ylab('Specificity')

ggplot(check_df, aes(x = sensitivity, y = specificity)) + 
  geom_point(size = 3, alpha = .3) + 
  coord_equal() + 
  ylim(0, 1) +
  theme_minimal() + 
  xlab('Sensitivity') + 
  ylab('Specificity') + 
  geom_text_repel(aes(label = Species)) +
  ggtitle("Neural network classification results")
