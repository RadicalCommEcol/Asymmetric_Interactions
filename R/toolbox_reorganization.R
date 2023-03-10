# R-code of "A guideline to study the feasibility domain of multi-trophic and changing ecological communities" by:
# Chuliang Song, Rudolf P. Rohr, Serguei Saavedra
# published in: Journal of Theoretical Biology

# load necessary packages
library(tidyverse)
library(mvtnorm)
library(mgcv)

Probabilities_all_Alf <- function(alpha,S,sims,delta_t) {
  survival <- rep(0,S)
  # sims <- 250  ### number of simualtions to obtain the probabilities of species persistence. In the text we used 100,000
  y <- 1
  while(y <= sims){
    time_step <- seq(0,300,by=delta_t)
    N0 <- rep(1,S) #initial conditions for species abundances
    N0 <- N0 / sum(N0)
    r <- sphere_sampling_ALF(S) ### sampling random K's in the unit sphere
    parms <- list(r=r, alpha = alpha) ##ODE
    model <- function(t,N,parms){ dN <- N * (parms$r - parms$alpha %*% N); list(dN)}
    sol <- ode(N0,time_step,model,parms)
    
    #plot(sol)
    
    survival_aux <- rep(0,S)
    
    for(z in 1:S){
       survival_aux[z] <- (sol[nrow(sol),z+1]>0.000001)*1 ### check if the species survived in the simulation
    }
      
    if(!all(is.na(survival_aux))){
      survival <- survival + survival_aux
      y <- y + 1
    }
    cat(y,"\n")
    
  }
  survival <- survival / sims ### calculate probabilities
  out <- survival
  return(out)
}

Probabilities_all_Alf_POSITIVE <- function(alpha,S,sims,delta_t) {
  survival <- rep(0,S)
  # sims <- 250  ### number of simualtions to obtain the probabilities of species persistence. In the text we used 100,000
  y <- 1
  while(y <= sims){
    time_step <- seq(0,300,by=delta_t)
    N0 <- rep(1,S) #initial conditions for species abundances
    N0 <- N0 / sum(N0)
    r <- abs(sphere_sampling_ALF(S)) ### sampling random K's in the unit sphere
    parms <- list(r=r, alpha = alpha) ##ODE
    model <- function(t,N,parms){ dN <- N * (parms$r - parms$alpha %*% N); list(dN)}
    sol <- ode(N0,time_step,model,parms)
    
    #plot(sol)
    
    survival_aux <- rep(0,S)
    
    for(z in 1:S){
      survival_aux[z] <- (sol[nrow(sol),z+1]>0.000001)*1 ### check if the species survived in the simulation
    }
    
    if(!all(is.na(survival_aux))){
      survival <- survival + survival_aux
      y <- y + 1
    }
    cat(y,"\n")
    
  }
  survival <- survival / sims ### calculate probabilities
  out <- survival
  return(out)
}

Probabilities_pos <- function(alpha,S) {
  survival <- rep(0,S)
  sims <- 250  ### number of simualtions to obtain the probabilities of species persistence. In the text we used 100,000
  for(y in 1:sims){
    delta_t <- 0.01 #time step
    time_step <- seq(0,300,by=delta_t)
    N0 <- rep(1,S) #initial conditions for species abundances
    N0 <- N0 / sum(N0)
    r <- c_sphere_sampling(S) ### sampling random K's in the positive orthant of the unit sphere
    parms <- list(r=r, alpha = alpha) ##ODE
    model <- function(t,N,parms){ dN <- N * (parms$r - parms$alpha %*% N); list(dN)}
    sol <- ode(N0,time_step,model,parms)
    for(z in 1:S){
      survival[z] <- survival[z] + (sol[nrow(sol),z+1]>0.000001)*1 ### check if the species survived in the simulation
    }
  }
  survival <- survival / sims ### calculate probabilities
  out <- survival
  return(out)
}

# function that computes the NORMALIZED feasibility from an interaction matrix
# inputs: alpha = interaction matrix
# output: out = the normalized feasibility
Omega <- function(alpha) {
  S <- nrow(alpha)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    return(out)
  }
#   if (length(which(diag(alpha) == 0)) == 0) {
#     Sigma <- chol2inv(alpha, size = NCOL(alpha), LINPACK = FALSE)
#     return(omega(S, Sigma))
#   }
#   else {
    f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
    if (f(alpha) == FALSE) {
      return(0)
    }
    else {
      Sigma <- solve(t(alpha) %*% alpha)
      return(omega(S, Sigma))
    }
#   }
# }
}

# Samples m vectors randomly on the positive orthant of a n-dimensional unit sphere
c_sphere_sampling <- function(m) {
  r <- abs(rnorm(m))
  d <- sqrt(sum(r^2))
  r <- r/d
  return(r)
}

# Samples m vectors randomly on a n-dimensional unit sphere
sphere_sampling <- function(m) {
  r <- (rnorm(m))
  r[1:3] <- abs(r[1:3])
  d <- sqrt(sum(r^2))
  r <- r/d
  return(r)
}

# Samples m vectors randomly on a n-dimensional unit sphere
sphere_sampling_ALF <- function(m) {
  r <- (rnorm(m))
  d <- sqrt(sum(r^2))
  r <- r/d
  return(r)
}


# function that normalizes a vector in the L2 norm
# inputs: a = the orignal vector
# output: the normalized vector
normalization <- function(a) {
  a / sqrt(sum(a^2))
}

# function that normalizes the spanning vectors of the feasibility domain in the L2 norm
# inputs: alpha = interaction matrix
# output: Span = the normalized spanning vectors
span_vectors <- function(alpha) {
  Span <- matrix(0, ncol = ncol(alpha), nrow = nrow(alpha))
  for (k in 1:ncol(alpha)) {
    Span[, k] <- -alpha[, k] / sqrt(sum(alpha[, k]^2))
  }
  Span
}

# function that computes all the extreme points that belong to original vertexes
# inputs: A = one interaction matrix, B = another interaction matrix
# output: inside_vertex = all the extreme points that belong to original vertexes
inside_vertex_detection <- function(A, B) {
  SpanA <- span_vectors(A)
  SpanB <- span_vectors(B)
  # to determine whether a vertex of one cone is inside another cone or not.
  inside_detection <- function(Span, vector) {
    lambda <- solve(Span, vector)
    if (sum(lambda >= -1e-10) == length(lambda)) {
      return(1)
    } else {
      return(0)
    }
  }
  inside_vertex <- list()
  l <- 1
  for (i in 1:ncol(B)) {
    auxi <- inside_detection(SpanA, SpanB[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanB[, i]
      l <- l + 1
    }
  }
  for (i in 1:ncol(A)) {
    auxi <- inside_detection(SpanB, SpanA[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanA[, i]
      l <- l + 1
    }
  }
  return(inside_vertex)
}

# function that computes all the extreme points generated that are generated by the intersections of the cones
# inputs: S = one interaction matrix, M = another interaction matrix
# output: intersection_vertex = all the extreme points that are generated by the intersections of the cones
intersection_vertex_detection <- function(S, M) {
  num <- ncol(S)
  combination_S <- combn(1:ncol(S), 2)
  combination_M <- combn(1:ncol(S), (num - 1))
  Span_S <- span_vectors(S)
  Span_M <- span_vectors(M)

  border_M <- list()
  extreme_point_M <- list()
  for (i in 1:ncol(M)) {
    coeff_matrix <- matrix(1, ncol = num, nrow = num)
    for (j in 1:(num - 1))
      coeff_matrix[j, ] <- Span_M[, combination_M[j, i]]
    coeff_vector <- c(rep(0, num - 1), 1)
    border_M[[i]] <- solve(coeff_matrix, coeff_vector)
    extreme_point_M[[i]] <- t(coeff_matrix)[1:(num - 1), 1:(num - 1)]
  }

  inside_face_detection <- function(extreme_point, test_vector) {
    lambda <- solve(extreme_point, test_vector)
    if (sum(lambda >= -1e-10) == length(lambda)) {
      return(1)
    } else {
      return(0)
    }
  }

  l <- 1
  intersection_vertex <- list()
  side <- c()
  for (i in 1:ncol(combination_S)) {
    vertex_1 <- Span_S[, combination_S[1, i]]
    vertex_2 <- Span_S[, combination_S[2, i]]
    for (j in 1:length(border_M)) {
      n1 <- sum(vertex_1 * border_M[[j]])
      n2 <- sum(vertex_2 * border_M[[j]])

      auxi <- n1 * n2
      if (auxi < -1e-10) {
        lambda <- n2 / (n2 - n1)
        possible <- lambda * vertex_1 + (1 - lambda) * vertex_2
        if (det(extreme_point_M[[j]]) != 0) {
          auxi2 <- inside_face_detection(extreme_point_M[[j]], possible[1:(num - 1)])
          if (auxi2 == 1) {
            intersection_vertex[[l]] <- possible
            side[l] <- j
            l <- l + 1
          }
        }
      }
    }
  }

  if (length(intersection_vertex) > 0) {
    for (i in 1:length(intersection_vertex)) {
      intersection_vertex[[i]] <- normalization(intersection_vertex[[i]])
    }
  }

  return(intersection_vertex)
}

# function that computes all the extreme points
# inputs: A = one interaction matrix, B = another interaction matrix
# output: out = all the extreme points that generate the intersection region
vertex_detection <- function(A, B) {
  num <- ncol(A)
  inside_vertex <- inside_vertex_detection(A, B)
  intersection_vertex <- intersection_vertex_detection(A, B)

  # combine the two vertex lists
  if (length(inside_vertex) > 0) {
    vertex <- matrix(unlist(inside_vertex), nrow = num, byrow = FALSE)
  } else {
    vertex <- matrix(0, nrow = num, ncol = 2)
  }
  if (length(intersection_vertex) > 0) {
    vertex <- cbind(vertex, matrix(unlist(intersection_vertex), nrow = num, byrow = FALSE))
  }

  # delete the points that are nonzero due to numerical error
  delete_zeroes <- c()
  for (i in 1:ncol(vertex)) {
    if (near(sum(vertex[, i]^2), 0)) {
      delete_zeroes <- c(delete_zeroes, i)
    }
  }
  if (length(delete_zeroes) > 0) vertex <- vertex[, -delete_zeroes]


  # delete the same ones
  if (length(vertex) > num) {
    for (test in 1:ncol(vertex)) {
      vertex[, test] <- normalization(vertex[, test])
    }
    delete_duplicates <- c()
    for (i in 1:(ncol(vertex) - 1)) {
      for (j in (i + 1):ncol(vertex)) {
        if (sum(near(vertex[, i], vertex[, j])) == nrow(vertex)) {
          delete_duplicates <- c(delete_duplicates, j)
        }
      }
    }
    if (length(delete_duplicates) > 0) vertex <- vertex[, -unique(delete_duplicates)]
  }
  return(vertex)
}

# function that computes all the extreme points
# inputs: p = the list of all extreme points
# output: out = all the extreme points that generate the intersection region
partitionize <- function(p) {
  triangulize <- function(p) {
    set.seed(100)
    inside_detection <- function(Span, vector) {
      lambda <- solve(Span, vector)
      if (sum(lambda >= -1e-10) == length(lambda)) {
        return(1)
      } else {
        return(0)
      }
    }
    border <- function(position) {
      coeff_matrix <- t(p[, position])
      coeff_matrix <- rbind(coeff_matrix, 1)
      coeff_vector <- c(rep(0, NR - 1), 1)
      if (abs(det(coeff_matrix)) > 1e-10) {
        border_M <- solve(coeff_matrix, coeff_vector)
      } else {
        border_M <- rep(0, NR)
      }
      return(border_M)
    }
    side_determine <- function(all_del) {
      # choose a candidate
      candidate <- sample(setdiff(c(1:NC), all_del), NR - 1)
      candidate_face <- border(candidate)

      # check which side other points are
      side <- rep(0, ncol(p))
      for (j in 1:length(side)) {
        side[j] <- sum(p[, j] * candidate_face)
      }
      side <- ifelse(abs(side) < 1e-10, 0, side)
      left_side <- which(side < 0)
      right_side <- which(side > 0)
      left_side <- setdiff(left_side, all_del)
      right_side <- setdiff(right_side, all_del)

      abundant_side <- which(side == 0)
      abundant_side <- setdiff(abundant_side, all_del)
      abundant_side <- setdiff(abundant_side, candidate)

      if (length(left_side) == 1) {
        all_del <- c(all_del, left_side)
        all_del <- unique(c(all_del, abundant_side))
        return(list(all_del, c(candidate, left_side)))
      }
      else if (length(right_side) == 1) {
        all_del <- c(all_del, right_side)
        all_del <- unique(c(all_del, abundant_side))
        return(list(all_del, c(candidate, right_side)))
      }
      else {
        return(c())
      }
    }
    isolation_finder <- function(all_del) {
      try <- side_determine(all_del)
      while (length(try) == 0) {
        try <- side_determine(all_del)
      }
      return(try)
    }
    NR <- nrow(p)
    NC <- ncol(p)
    partition <- matrix(0, nrow = 1, ncol = NR)
    all_del <- c()
    while (length(setdiff(c(1:NC), all_del)) >= NR) {
      result <- isolation_finder(all_del)
      all_del <- result[[1]]
      if (length(result) > 1) partition <- rbind(partition, result[[2]])
    }
    # partition
    partition <- partition[-1, ]
    partition <- matrix(partition, ncol = NR)
    # partition <- partition[!duplicated(partition[,]),]
    return(partition)
  }
  # if the number of vertexes is the same as the dimension
  if (nrow(p) == ncol(p)) {
    matrix(1:(nrow(p)), nrow = 1)
  }
  # else it is a triangulization problem
  else {
    return(triangulize(p))
  }
}

# function that computes the full volume
# inputs: partition = triangulation of the intersection region, vertex = all the extreme points
# output: the normalize feasibility of the intersection region
total_volume <- function(partition, vertex) {
  vol <- c()
  num <- ncol(partition)
  if (length(partition) == num) {
    auxi <- vertex[, partition]
    if (near(det(auxi), 0)) {
      vol <- 0
    } else {
      vol <- Omega(auxi)
    }
  } else {
    for (i in 1:nrow(partition)) {
      auxi <- vertex[, partition[i, ]]
      if (near(det(auxi), 0)) {
        vol[i] <- 0
      } else {
        vol[i] <- Omega(auxi)
      }
    }
  }
  sum(vol^num)^(1 / num)
}

sampling_overlap <- function(A,B){
  mat <- solve(B) %*% A
  Nsample <- 10^5
  abunbance_all <- rmvnorm(n = Nsample, mean = rep(0, nrow(A))) %>% 
    {abs(./sqrt(rowSums(.^2)))}
  get_feasibility <- function(N_A){
    N_B <- mat %*% matrix(N_A,ncol=1) %>% c()
    if_else(sum(N_B >= -1e-10) == length(N_B), 1, 0) 
  }
  percent <- 1:Nsample %>% 
    map_dbl(~get_feasibility(abunbance_all[.x,])) %>% 
    mean()
  Omega(A) * percent^(1/nrow(A))
}


# function that computes the overlap of two feasibility domains
# inputs: A = one interaction matrix, B = another interaction matrix
# output: volume_overlap = the normalize feasibility of the intersection region
Omega_overlap <- function(A, B) {
  num <- nrow(A)

  overlap_vertex <- vertex_detection(A, B)
  if (qr(overlap_vertex)$rank < num) {
    volume_overlap <- 0
  } else {
    partition <- partitionize(overlap_vertex)
    volume_overlap <- total_volume(partition, overlap_vertex)
  }

  volume_overlap <- sampling_overlap(A,B)
  volume_overlap
}

# function that generates a random matrix (as in May, Nature (1972))
# inputs: alpha = interaction matrix; inequality = biological linear ineqaulity
# output: the size of the constrained feasibility domain
Omega_with_inequality_constraints <- function(alpha, inequality){
  beta <- alpha
  for(j in 1:ncol(alpha)){
    beta[,j] <- -alpha[,j] / (-alpha[,j] %*% inequality)
  }
  abs(det(beta))
}

# function that generates a random matrix (as in May, Nature (1972))
# inputs: num = number of species; stren = standard deviation of interaction strength; conne = connectance of the interaction matrix
# output: Inte = the generated random matrix
interaction_matrix_random <- function(num, stren, conne){
  Inte <- rnorm(num*num, mean = 0, sd = stren)
  zeroes <- sample(c(rep.int(1,floor(num*num*conne)),rep.int(0,(num*num-floor(num*num*conne)))))
  Inte[which(zeroes==0)] <- 0
  Inte <- matrix(Inte, ncol = num, nrow = num)
  diag(Inte) <- -1
  return(Inte)
}



# A <- c(-1.00000000,  -0.6381123, -0.7345733, -1.0000000) %>% 
#   matrix(ncol=2)
# get_N <- function(){
#   theta <- runif(1, 0, pi/2)
#   r <- c(cos(theta), sin(theta))
#   # N <- 
#   - A %*% r %>% 
#     t() %>% 
#     as_tibble()
# }
# 1:1000 %>% 
#   map_dfr(~get_N()) %>% 
#   ggplot(aes(V1,V2))+
#   geom_point()

# parameterization_feasible <- function(A, num){
#   lambda <- runif(num, min = 0, max = 1)
#   lambda <- lambda/sum(lambda)
#   -solve(A, lambda)
# }

spanned_vectors <- function(A, num){
  G <- matrix(0, ncol=ncol(A), nrow=nrow(A))
  for(k in 1:num) G[,k] <- -A[,k]/sqrt(sum(A[,k]^2))
  G
}

parameterization_feasible <- function(A, num){
  G <- spanned_vectors(A, num)
  lambda <- runif(num, min = 0, max = 1)
  lambda <- lambda/sum(lambda)
  growth <- G %*% matrix(lambda, ncol=1) %>%
    as.vector()
}
