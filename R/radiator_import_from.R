# magrittr ---------------------------------------------------------------------
# Forward-pipe operator --------------------------------------------------------
#' @title Forward-pipe operator
#' @description magrittr forward-pipe operator
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

# Exposition pipe-operator ------------------------------------------------------
#' @title Exposition pipe-operator
#' @description magrittr Exposition pipe-operator
#' @name %$%
#' @rdname Exposition_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
NULL

# Compound assignment pipe operator --------------------------------------------
#' @title compound assignment pipe operator
#' @description magrittr compound assignment pipe operator
#' @name %<>%
#' @rdname compound_assignment_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

# dplyr n ----------------------------------------------------------------------
# The number of observations in the current group.
#' @title The number of observations in the current group.
#' @description Check dplyr
#' @name n
#' @rdname n
#' @keywords internal
#' @export
#' @importFrom dplyr n
#' @usage n()
NULL

# importFrom -------------------------------------------------------------------
#' @importFrom pbmcapply pbmclapply
#' @importFrom stringi stri_replace_all_fixed stri_join stri_sub stri_replace_na stri_pad_left
#' @importFrom stats IQR
#' @importFrom utils sessionInfo packageDescription
#' @importFrom parallel detectCores makeCluster clusterExport mclapply parLapply stopCluster
#' @importFrom stringi stri_join stri_replace_all_fixed stri_sub stri_detect_fixed
