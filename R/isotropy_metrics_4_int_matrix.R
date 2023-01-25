
#' Isotropy metrics for a given interaction matrix
#'
#' Given an interaction matrix, this function generates a tibble with S rows,
#' where S is the number of species in the matrix. Each row contains the
#' following information for a given: replicates of its probability of exclusion
#' and the corresponding value of capital Omega for each replicate (i.e., the 
#' proportion of the feasible parameter #' space inside the unit sphere). Those 
#' estimations are computed via a quasi-Monte Carlo method and by 
#' bootstrapping the quasi-Mote Carlo results.
#'
#' @param Parameters: 
#'  \describe {
#'    \item{A_int: SxS interaction matrix, where S is the number of species.}
#'    \item{number_Omega_replicates: a number that specify how many estimations 
#'    of Omega will be calculated by the quasi-Monte Carlo method. By default, 
#'    this parameter is set to 1,000 replicates.}
#'    \item{number_boot_replicates: a number that specify how many bootstrap
#'    estimations of the probability of exclusion and Omega will be calculated
#'    for each species. By default, this parameter is set to 1,000 replicates.}
#'    
#'  }
#' 
#' @section Dependencies: 
#'  \describe {
#'    \item{tidyverse}
#'    \item{foreach}
#'    \item{doParallel}
#'    \item{mvtnorm}
#'    \item{boot}
#'    \item{pracma}
#'    \item{matlib}
#'    \item{zipfR}
#'    \item{vertices_unit_ball (function)}
#'    \item{cone_vertices_director_vertices (function)}
#'    \item{incenter_inradius_isoprob_calculation (function)}
#'    \iten{boot_prob_excl_Omega_raw (function)}
#'  }
#'
#' @return A tibble S rows and 2 x number_boot_replicates columns. Each row 
#' contains the estimated probabilities of exclusion and values of Omega for 
#' each species, when considering the intersection of the feasibility domain 
#' and the unit ball. Specifically, each row contains:
#'     \describe {
#'    \item{Resamples of the estimated probabilities of exclusion for species i 
#'    (prob_excl_rep_x), where x takes values from 1 to number_boot_replicates.}
#'    \item{The corresponding value of capital Omega for species i and 
#'    a given prob_excl_rep_x (Omega_rep_x), where x takes values from 1 to 
#'    number_boot_replicates.}
#'    } 
#'    
#' @examples
#' A_int <- -1*diag(c(3,5,7))
#' boot_prob_excl_Omega_raw(A_int, replicates = 1e3)
#' 
#' @references \url{https://doi.org/10.1016/j.jtbi.2018.04.030}
#'
#' @export


isotropy_metrics_4_int_matrix <- function(A_int, number_Omega_replicates,
                                          number_boot_replicates){
  
  prob_excl_Omega_df <- boot_prob_excl_Omega_raw(A_int, number_Omega_replicates,
                                             number_boot_replicates)
  
  isotropy_agg_info <- tibble(interaction_matrix = 1)
  small_omega_mean <- NA
  small_omega_lowerCI <- NA
  small_omega_upperCI <- NA
  isotropy_index_mean <- NA
  isotropy_index_lowerCI <- NA
  isotropy_index_upperCI <- NA
  
  
  Omega_rep <- prob_excl_Omega_df[1,(1+number_boot_replicates):(2*number_boot_replicates)] %>%
    as.numeric()
  small_omega_rep <- Omega_rep^(1/nrow(A_int))
  isotropy_rep <- NULL
  
  for(i in 1:number_boot_replicates){
    
    probability_vector <- prob_excl_Omega_df[,i+1] %>% pull()
    
    isotropy_rep <-  c(isotropy_rep, isotropy_index(probability_vector))
    
  }
  
  small_omega_rep_CI <- quantile(small_omega_rep, prob=c(.025,.975)) %>% as.numeric()
  isotropy_rep_CI <- quantile(isotropy_rep, prob=c(.025,.975)) %>% as.numeric()
  
  isotropy_agg_info$small_omega_mean <- mean(small_omega_rep, na.rm = T)
  isotropy_agg_info$small_omega_lowerCI <- small_omega_rep_CI[1]
  isotropy_agg_info$small_omega_upperCI <- small_omega_rep_CI[2]
  isotropy_agg_info$isotropy_index_mean <- mean(isotropy_rep, na.rm = T)
  isotropy_agg_info$isotropy_index_lowerCI <- isotropy_rep_CI[1]
  isotropy_agg_info$isotropy_index_upperCI <- isotropy_rep_CI[2]
  
  return(isotropy_agg_info %>% select(-interaction_matrix))
  
}
