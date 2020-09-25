# sort level names in order of coefficients
sort.levels.in.order.of.coef <- function(aov_res, factor_name) {

  # fetch level names
  level_name <- aov_res$xlevels[[factor_name]]

  # fetch cofficients and rename
  coefval <- aov_res$coefficients
  coefval <- coefval[grepl(paste0('^', factor_name), names(coefval))]
  names(coefval) <- gsub(factor_name, '', names(coefval))

  # fill all the coefficients
  coefval <- coefval[level_name]
  names(coefval) <- level_name
  coefval[is.na(coefval)] <- 0

  # sort level name
  level_name <- names(coefval)[order(coefval, decreasing = TRUE)]

  # output
  level_name

}

# create initial empty matrix for the algorithm
create.empty.matrix <- function(level_name) {

  # create matrix with NA
  n_level <- length(level_name)
  alphabet_mat <- matrix('-', nrow = n_level, ncol = n_level)
  rownames(alphabet_mat) <- level_name
  colnames(alphabet_mat) <- level_name

  # fill upper triangle with '0
  alphabet_mat[upper.tri(alphabet_mat)] <- '0'

  # output
  alphabet_mat

}

# fill the cells of alphabet_mat
# if the difference of the groups are significant
reflect.significance <- function(alphabet_mat, level0, tukey_res, alpha) {

  # extract result matrix of the target level
  selector1 <- grepl(level0, rownames(tukey_res))
  tukey_res_level0 <- tukey_res[selector1, ]

  # extract level names with significance
  selector2 <- tukey_res_level0[, 'p adj'] < alpha
  signif_level <- rownames(tukey_res_level0)[selector2]
  signif_level <- gsub(level0, '', signif_level)
  signif_level <- gsub('-', '', signif_level)

  # fill matrix
  selector3 <- (alphabet_mat[signif_level, level0] == '-')
  signif_level <- signif_level[selector3]
  alphabet_mat[signif_level, level0] <- '1'

  # output
  alphabet_mat

}


#' Make groups for the results of TukeyHSD() with alphabets
#'
#' @param aov_res An object created with aov().
#' @param factor_name The name of explanatory variable you want to test.
#' @param alpha The significance level in multi-comparison test.
#' @param ... Other parameters for TukeyHSD().
#'
#' @export
#'
TukeyHSD.alphabet <- function(aov_res, factor_name, alpha = 0.05, ...) {

  # factor_name <- 'x1'
  # alpha <- 0.01

  # check class of the object
  is_class_good_1 <- identical(class(aov_res), c('aov', 'lm'))
  is_class_good_2 <- identical(class(factor_name), 'character')
  is_class_good_3 <- identical(class(alpha), 'numeric')
  if (!is_class_good_1) {
    stop('Class of aov_res is not correct')
  } else if (!is_class_good_2) {
    stop('Class of factor_name is not correct')
  } else if (!is_class_good_3) {
    stop('Class of alpha is not correct')
  }

  # extract result matirx of TukeyHSD
  tukey_res <- TukeyHSD(aov_res, which = factor_name, ...)
  tukey_res <- tukey_res[[factor_name]]
  # tukey_res

  # sort estimated coefficients
  level_name <- sort.levels.in.order.of.coef(aov_res, factor_name)

  # create empty matrix
  alphabet_mat <- create.empty.matrix(level_name)

  # fill the cells of alphabet_mat
  # where the differences of the groups are significant
  for (i in 1:length(level_name)) {
    level0 <- level_name[i]
    alphabet_mat <- reflect.significance(
      alphabet_mat, level0, tukey_res, alpha
    )
    # print(alphabet_mat)
  }

  # fill the cells of alphabet_mat with alphabet
  n_alphabet <- 1
  for (i in 1:length(level_name)) {
    # i <- 1
    # i <- 2

    # level name
    level0 <- level_name[i]

    # special procedure for the first loop
    if (i == 1) {
      selector_i1 <- (alphabet_mat[, level0] == '-')
      alphabet_mat[selector_i1, level0] <- 'a'
      next
    }

    # judge if the significance is equal to the left column
    signif1 <- (alphabet_mat[, i - 1] == '1')
    signif2 <- (alphabet_mat[, i] == '1')

    # if significance is different with the left column,
    # change the alphabet
    if (!identical(signif1, signif2)) {
      n_alphabet <- n_alphabet + 1
    }

    # fill matrix with alphabet
    selector_i2 <- (alphabet_mat[, level0] == '-')
    alphabet_mat[selector_i2, level0] <- letters[n_alphabet]

  }

  # create grouping identifiers
  alphabet_val <- apply(alphabet_mat, 1, function(x) {
    x <- x[x != '0']
    x <- x[x != '1']
    x <- unique(x)
    paste(x, collapse = '')
  })

  # output
  alphabet_val <- alphabet_val[order(names(alphabet_val))]
  alphabet_val

}



### test code from here:

# # create simulation data
# create.simulation.data <- function() {
#
#   # factor name
#   name_x1 <- paste0('P', LETTERS[1:5])
#   name_x2 <- paste0('Q', LETTERS[1:3])
#
#   # create data
#   dat <- expand.grid(
#     x1 = name_x1,
#     x2 = name_x2,
#     rep = paste0('rep', 1:3)
#   )
#
#   # set effects
#   eff_x1 <- c(4, 5, 0, 9, 7)
#   eff_x2 <- c(-2, 0, -5)
#   names(eff_x1) <- name_x1
#   names(eff_x2) <- name_x2
#   sd_e <- 2
#
#   # calculate
#   n <- nrow(dat)
#   set.seed(1000)
#   dat$y <- with(
#     dat, eff_x1[x1] + eff_x2[x2] + rnorm(n, 0, sd_e)
#   )
#   dat
#
# }
#
# dat <- create.simulation.data()
# plot(dat$x1, dat$y); points(dat$x1, dat$y)
# plot(dat$x2, dat$y); points(dat$x2, dat$y)
# aov_res <- aov(y ~ x1 + x2, dat)
# TukeyHSD.alphabet(aov_res, 'x1', alpha = 0.05)
# TukeyHSD.alphabet(aov_res, 'x1', alpha = 0.01)
# TukeyHSD.alphabet(aov_res, 'x2', alpha = 0.05)

