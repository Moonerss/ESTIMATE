#' @name plotPurity
#' @title Plot tumor purity
#' @description Plot tumor purity based on ESTIMATE score.
#' @param scores A data frame object get from \code{estimateScore}.
#' @param samples ector of character strings specifying sample names to be plotted. Defaults to "all_samples", which creates plots for all input samples
#' @param platform character string indicating platform type. Defaults to "affymetrix"
#' @details This function produces scatterplots for each requested sample; it plots tumor purity against ESTIMATE score. At present, only the Affymetrix platform is supported.
#' @return return an ggplot object
#' @import ggplot2 patchwork
#' @export
#' @author Erjie Zhao <2055469819@qq.com>
#' @examples
#'  \dontrun{
#'    file <- system.file("extdata", "sample_input.txt", package="ESTIMATE")
#'    expression <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
#'    merged_expr <- filterCommonGenes(expr = expression, id = "GeneSymbol")
#'    scores <- estimateScore(merged_expr, platform = "affymetrix")
#'    plot_obj <- plotPurity(scores)
#'  }
plotPurity <- function(scores, samples="all_samples",
                       platform=c("affymetrix", "agilent", "illumina")) {
  options(warn = -1)
  ## Check arguments
  platform <- match.arg(platform)

  if (platform != "affymetrix"){
    stop("not implemented")
  }

  ## Begin processing

  ##-------------------------------------------------------------------------
  convert_row_estimate_score_to_tumor_purity <- function(x) {
    stopifnot(is.numeric(x))
    cos(0.6049872018 + 0.0001467884*x)
  }

  ## Read ESTIMATE data file
  estimate.df <- as.data.frame(t(scores))
  samplenames <- rownames(estimate.df)
  Affy.model <- PurityDataAffy
  pred.p <- Affy.model[, 5:7]
  est <- estimate.df[, 3]
  est.new <- estimate.df[, 4]

  ## ESTIMATE based tumor purity in scatterplot with prediction interval
  message("Plotting tumor purity based on ESTIMATE score")

  max.af <- max(Affy.model$ESTIMATEScore)
  min.af <- min(Affy.model$ESTIMATEScore)

  if (samples[1] == "all_samples"){
    Num.S <- nrow(estimate.df)
  } else {
    Num.S <- as.numeric(length(samples))
  }

  ## create plot object
  plot.obj <- list()
  for (i in 1:Num.S) {
    if(samples[1] =="all_samples"){
      samplename <- samplenames[i]
    } else {
      samplename <- samples[i]
    }

    geMin <- est[i] >= min.af
    leMax <- est[i] <= max.af
    withinMinMax <- geMin && leMax

    xlim <- if (!withinMinMax) {
      ## Expands plot boundary
      adjustment <- 500    # Arbitrary
      if (geMin) {
        from <- min.af
        to   <- est[i] + adjustment
      } else {
        from <- est[i] - adjustment
        to   <- max.af
      }
      c(from, to)
    } else {
      NULL
    }

    plot.obj[[i]] <- ggplot(Affy.model) +
      aes(x = Affy.model$ESTIMATEScore, y = Affy.model$tumor.purity) +
      geom_point(size = 1, colour = "lightgrey", fill = "white", shape = 21) +
      annotate("point", x = est[i], y = est.new[i], size = 3, colour = "black") +
      geom_line(mapping = aes(x = Affy.model$ESTIMATEScore, y = Affy.model$fit),
                colour = "darkgrey", linetype = 1, size = 1) +
      geom_line(mapping = aes(x = Affy.model$ESTIMATEScore, y = Affy.model$lwr.p),
                colour = "darkgrey", linetype = 2, size = 1) +
      geom_line(mapping = aes(x = Affy.model$ESTIMATEScore, y = Affy.model$upr.p),
                colour = "darkgrey", linetype = 2, size = 1) +
      geom_vline(xintercept = est[i], colour = "black", linetype = 2, size = 1) +
      geom_hline(yintercept = est.new[i], colour = "black", linetype = 2, size = 1) +
      labs(title = samplename, x = "ESTIMATE score", y = "Tumor purity") +
      ylim(0, 1) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank())

    if (!withinMinMax) {
      plot.obj[[i]] <- plot.obj[[i]] +
        stat_function(fun = convert_row_estimate_score_to_tumor_purity,
                      n = 10000, colour = "grey", linetype = 1, size = 1) +
        xlim(from, to)
    }
  }

  ## merge all ggplot object
  res_obj <- patchwork::wrap_plots(plot.obj)

  return(res_obj)
}

