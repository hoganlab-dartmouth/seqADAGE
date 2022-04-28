###############################################################################
# This script compares pathway coverage between ensemble models and individual
# models after removing pathway crosstalk effects.
#
# Usage:
#     Rscript analyze_ensemble_coverage_no_crosstalk.R
#
###############################################################################

pacman::p_load("ggplot2")

########### load in constant

# no_crosstalk_file stores the pathway coverage results of each model after
# crosstalk correction
no_crosstalk_file <- "./models_sigpathways_no_crosstalk/MODEL_PATHWAY_COVERAGE_no_crosstalk.txt"
# summary_file stores the output of this script, it will be overwritten
# in each run.
summary_file = "./analyze_ensemble_coverage_no_crosstalk_output.txt"

########### load in function

analyze_coverage_no_crosstalk <- function(coverage.df, type, summary_file){
  # This function analyzes the pathway coverage result after crosstalk correction
  # for either ensemble ADAGE models built from 100 same individual models
  # or ensemble ADAGE models built from 1000 different individual models.
  #
  # Input:
  # coverage.df: a data.frame that contains pathway coverage results of all
  #              models.
  # type: either "same" or "different", other input will raise error.
  # summary_file: text output of this function (mainly t test results) will be
  #               written to summary_file

  if (type == "same") {
    opp.type = "different"
  } else if (type == "different") {
    opp.type = "same"
  } else {
    stop("model type not recognized!")
  }

  write(paste("comparing ensemble models built from the", type, "100 individual",
               "models with 1000 individual ADAGE models:\n"),
        summary_file, append = TRUE)

  coverage.df.type <- coverage.df[!grepl(opp.type, coverage.df$model), ]

  # make a boxplot of pathway coverage
  ggplot(coverage.df.type, aes(x = model, y = coverage)) +
    geom_boxplot() + labs(x = "Model", y = "Pathway Coverage") + theme_bw()
  ggsave(paste(type, "pathway coverage comparison (no crosstalk).pdf", sep = "_"),
         height = 5, width = 5)

  # t test to compare the number of pathways covered by ADAGE models and the
  # ensemble models
  write("comparing corADAGE with ADAGE:", summary_file, append = TRUE)
  t.result <- t.test(coverage.df.type[
    coverage.df.type$model == paste("corADAGE", type), "coverage"],
    coverage.df.type[coverage.df.type$model == "individual ADAGE", "coverage"],
    alternative = "greater")
  write(capture.output(t.result), summary_file, append = TRUE)
  write("comparing eADAGE with ADAGE:", summary_file, append = TRUE)
  t.result <- t.test(coverage.df.type[
    coverage.df.type$model == paste("eADAGE", type), "coverage"],
    coverage.df.type[coverage.df.type$model == "individual ADAGE", "coverage"],
    alternative = "greater")
  write(capture.output(t.result), summary_file, append = TRUE)
  write("comparing eADAGE with corADAGE:", summary_file, append = TRUE)
  t.result <- t.test(coverage.df.type[
    coverage.df.type$model == paste("eADAGE", type), "coverage"],
    coverage.df.type[coverage.df.type$model == paste("corADAGE", type), "coverage"],
    alternative = "greater")
  write(capture.output(t.result), summary_file, append = TRUE)

}

########### read in and process the pathway coverage result

write("Here are pathway coverage results after crosstalk correction:\n",
      summary_file)
no_crosstalk_result <- read.delim(no_crosstalk_file, header = FALSE,
                                  row.names = 1, stringsAsFactors = F)
colnames(no_crosstalk_result) <- c("Count", "Coverage")

coverage.df.list <- lapply(1:nrow(no_crosstalk_result), function(x) {
  coverage <- as.numeric(unlist(strsplit(no_crosstalk_result[x, "Coverage"], ";")))
  coverage.df <- data.frame(model = rep(rownames(no_crosstalk_result)[x],
                                        length(coverage)),
                            coverage = coverage)
  coverage.df
})

coverage.df <- do.call("rbind", coverage.df.list)
coverage.df$model <- factor(coverage.df$model, levels = c("individual ADAGE",
                                                          "corADAGE same",
                                                          "corADAGE different",
                                                          "eADAGE same",
                                                          "eADAGE different"))

# compare pathway coverage between ensemble models built from same 100
# individual models with 1000 individual models
analyze_coverage_no_crosstalk(coverage.df, type = "same", summary_file)

# compare pathway coverage between ensemble models built from different 100
# individual models with 1000 individual models
analyze_coverage_no_crosstalk(coverage.df, type = "different", summary_file)
