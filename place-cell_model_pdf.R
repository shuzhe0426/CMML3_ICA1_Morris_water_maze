run_trial <- function(weights0, Wmult, sigma_pc, sigma_ac, PC_x, PC_y,
                      Vdecay, ac_const, beta, etdecay, lrate, discf, noise,
                      platform_x, platform_y, starting_x, starting_y,
                      speed, wall_pun)
{
  # FIXED PARAMETERS OF THE EXPERIMENT
  pool_diameter <- 1.4   # Maze diameter in metres (m)
  platform_radius <- 0.06 # Platform radius
  
  N_pc <- 211  # Population of place cells
  N_ac <- 36   # Population of action cells
  which <- 0
  
  dist <- 0
  wall_zone <- 0
  quadrants <- c(0, 0, 0, 0) # Percentage spent on each quadrant
  
  weights <- weights0
  el_tr <- matrix(rep(0, N_pc * N_ac), nrow = N_pc)
  
  # Initialize trajectories
  track_x <- starting_x
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0
  
  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 +
         (track_y[length(track_y)] - platform_y)^2 > platform_radius^2)
  {
    weights <- weights * (1 - noise) +
      matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult * noise
    
    # Calculate PC activation
    PC_activation <- rep(0, N_pc)
    for (i in 1:N_pc) {
      PC_activation[i] <- exp(-((track_x[length(track_x)] - PC_x[i])^2 +
                                  (track_y[length(track_y)] - PC_y[i])^2) /
                                (2 * sigma_pc^2))
    }
    
    # Calculate AC activation (Q values)
    if (length(track_x) > 1) {
      prevQ <- AC_activation[which]
    }
    
    AC_activation <- PC_activation %*% weights
    
    # Make an action
    ACsel <- AC_activation^beta
    ACsel <- ACsel / sum(ACsel)
    ASrand <- runif(1)
    which <- 1
    ASsum <- ACsel[1]
    
    while (which < N_ac && ASsum < ASrand) {
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }
    
    # Eligibility traces
    el_tr <- el_tr * etdecay
    
    for (j in 1:N_ac) {
      itmp <- min(abs(j - which), N_ac - abs(j - which))
      actgaus <- exp(-(itmp * itmp) / (2 * sigma_ac * sigma_ac))
      el_tr[, j] <- el_tr[, j] + actgaus * AC_activation[j] * t(t(PC_activation))
    }
    
    vel_x <- c(vel_x, (vel_x[length(vel_x)] + ac_const * cos(which / N_ac * 2 * pi)) * Vdecay)
    vel_y <- c(vel_y, (vel_y[length(vel_y)] + ac_const * sin(which / N_ac * 2 * pi)) * Vdecay)
    
    track_x <- c(track_x, track_x[length(track_x)] + vel_x[length(vel_x)])
    track_y <- c(track_y, track_y[length(track_y)] + vel_y[length(vel_y)])
    
    # Check boundaries
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter / 2)^2) {
      ratio <- (track_x[length(track_x)]^2 + track_y[length(track_y)]^2) / ((pool_diameter / 2)^2)
      track_x[length(track_x)] <- track_x[length(track_x)] / sqrt(ratio)
      track_y[length(track_y)] <- track_y[length(track_y)] / sqrt(ratio)
      vel_x[length(vel_x)] <- track_x[length(track_x)] - track_x[length(track_x) - 1]
      vel_y[length(vel_y)] <- track_y[length(track_y)] - track_y[length(track_y) - 1]
    }
    
    if (length(track_x) > 2) {
      if ((track_x[length(track_x)] - platform_x)^2 +
          (track_y[length(track_y)] - platform_y)^2 < platform_radius^2) {
        rew <- 10
      } else if (track_x[length(track_x)]^2 +
                 track_y[length(track_y)]^2 > (0.99 * pool_diameter / 2)^2) {
        rew <- -wall_pun
      } else {
        rew <- 0
      }
      
      currQ <- AC_activation[which]
      tderr <- rew + discf * currQ - prevQ
      weights <- pmax(weights + lrate * tderr * el_tr, 0)
    }
    
    laststep <- sqrt((track_x[length(track_x)] - track_x[length(track_x) - 1])^2 +
                       (track_y[length(track_y)] - track_y[length(track_y) - 1])^2)
    dist <- dist + laststep
    
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > 0.8 * (pool_diameter / 2)^2) {
      wall_zone <- wall_zone + 1
    } else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0) {
      quadrants[1] <- quadrants[1] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0) {
      quadrants[2] <- quadrants[2] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0) {
      quadrants[3] <- quadrants[3] + 1
    } else {
      quadrants[4] <- quadrants[4] + 1
    }
    
    if (length(track_x) > 100) {
      speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2)))
      latency <- (length(track_x) - 1) * speed_ts / speed
      if (latency > 60) {
        break
      }
    }
  }
  
  latency <- length(track_x) - 1
  wall_zone <- wall_zone / latency
  quadrants <- quadrants / latency
  speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2)))
  latency <- latency * speed_ts / speed
  
  return(list(weights, track_x, track_y, vel_x, vel_y,
              dist, wall_zone, quadrants, latency))
}

# =========================
# PARAMETERS
# =========================
N_pc <- 211
N_ac <- 36

plot_trajectories <- 1
plot_cognitive_maps <- 1
pln <- plot_trajectories + plot_cognitive_maps
Nruns <- 10

pool_diameter <- 1.4
platform_radius <- 0.06
sigma_pc <- 0.1
sigma_ac <- 2

etdecay <- 0.83
beta <- 6
alpha <- 0.01
gamma <- 0.85

Vdecay <- 0.82
ac_const <- 0.02
Wnoise <- 0.0004
Wmult <- 0.1
hitwall <- 0.5
speed <- 0.175

Ntrials <- 4
Ndays <- 5

# performance measures: latency, distance, target quadrant, opposite quadrant, wall zone
if (pln > 0.5) {
  PMs <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
} else {
  PMs <- array(rep(0, 5 * Ndays * Ntrials * Nruns), c(5, Ndays, Ntrials, Nruns))
}

# Platform coordinates
platform_x <- cos(-pi / 4) * pool_diameter / 4
platform_y <- sin(-pi / 4) * pool_diameter / 4

# Starting locations
strad <- pool_diameter / 2 * 0.85
starting_xs <- strad * c(cos(pi / 6), cos(pi / 3), cos(7 * pi / 6), cos(4 * pi / 3))
starting_ys <- strad * c(sin(pi / 6), sin(pi / 3), sin(7 * pi / 6), sin(4 * pi / 3))

th <- (0:100) / 50 * pi

# 先关闭旧图形设备，避免冲突
graphics.off()

if (pln > 0.5) {
  
  # ===== 输出到 PDF，避免 RStudio plotting pane 太小 =====
  pdf("trajectories_and_maps.pdf", width = 16, height = 10)
  
  # 缩小边距
  par(mfrow = c(Ndays, Ntrials * pln), mar = c(1.5, 1.5, 2, 1))
  
  # Generate initial weights
  weights <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult
  
  # Generate place cells
  PC_x <- rep(0, N_pc)
  PC_y <- rep(0, N_pc)
  
  for (i in 1:N_pc) {
    PC_x[i] <- (runif(1) - 0.5) * pool_diameter
    PC_y[i] <- (runif(1) - 0.5) * pool_diameter
    while (PC_x[i]^2 + PC_y[i]^2 > (pool_diameter / 2)^2) {
      PC_x[i] <- (runif(1) - 0.5) * pool_diameter
      PC_y[i] <- (runif(1) - 0.5) * pool_diameter
    }
  }
  
  for (day in 1:Ndays) {
    idxs <- sample(4)
    
    for (trial in 1:Ntrials) {
      idx <- idxs[trial]
      starting_x <- starting_xs[idx]
      starting_y <- starting_ys[idx]
      
      modresults <- run_trial(
        weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y,
        Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise,
        platform_x, platform_y, starting_x, starting_y, speed, hitwall
      )
      
      weights <- modresults[[1]]
      track_x <- modresults[[2]]
      track_y <- modresults[[3]]
      
      PMs[1, day, trial] <- modresults[[9]]
      PMs[2, day, trial] <- modresults[[6]]
      PMs[3, day, trial] <- modresults[[8]][4] * 100
      PMs[4, day, trial] <- modresults[[8]][2] * 100
      PMs[5, day, trial] <- modresults[[7]] * 100
      
      if (plot_trajectories == 1) {
        plot(pool_diameter / 2 * cos(th), pool_diameter / 2 * sin(th),
             type = "l",
             xlab = paste("day", day, "trial", trial),
             ylab = "trajectory",
             asp = 1)
        lines(track_x, track_y, type = "l")
        lines(platform_x + platform_radius * cos(th),
              platform_y + platform_radius * sin(th),
              type = "l")
      }
      
      if (plot_cognitive_maps == 1) {
        plot(pool_diameter / 2 * cos(th), pool_diameter / 2 * sin(th),
             type = "l",
             xlab = paste("day", day, "trial", trial),
             ylab = "cognitive map",
             asp = 1)
        
        for (x in (-3:3) * (pool_diameter / 6)) {
          for (y in (-3:3) * (pool_diameter / 6)) {
            if (x^2 + y^2 <= (pool_diameter / 2)^2) {
              x2 <- x
              y2 <- y
              
              for (k in 1:N_ac) {
                PC_activation <- rep(0, N_pc)
                for (i in 1:N_pc) {
                  PC_activation[i] <- exp(-((x - PC_x[i])^2 + (y - PC_y[i])^2) / (2 * sigma_pc^2))
                }
                
                AC_activation <- rep(0, N_ac)
                for (i in 1:N_ac) {
                  for (j in 1:N_pc) {
                    AC_activation[i] <- AC_activation[i] + PC_activation[j] * weights[j, i]
                  }
                }
                
                x2 <- c(x2, x + (AC_activation[k] / 10) * cos(k / N_ac * 2 * pi))
                y2 <- c(y2, y + (AC_activation[k] / 10) * sin(k / N_ac * 2 * pi))
              }
              
              lines(x2, y2, type = "l", col = "blue")
            }
          }
        }
        
        lines(platform_x + platform_radius * cos(th),
              platform_y + platform_radius * sin(th),
              type = "l")
      }
    }
  }
  
  dev.off()
  message("PDF saved as: trajectories_and_maps.pdf")
  
} else {
  
  # run multiple times without plotting
  for (reps in 1:Nruns) {
    weights <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult
    
    PC_x <- rep(0, N_pc)
    PC_y <- rep(0, N_pc)
    
    for (i in 1:N_pc) {
      PC_x[i] <- (runif(1) - 0.5) * pool_diameter
      PC_y[i] <- (runif(1) - 0.5) * pool_diameter
      while (PC_x[i]^2 + PC_y[i]^2 > (pool_diameter / 2)^2) {
        PC_x[i] <- (runif(1) - 0.5) * pool_diameter
        PC_y[i] <- (runif(1) - 0.5) * pool_diameter
      }
    }
    
    for (day in 1:Ndays) {
      idxs <- sample(4)
      
      for (trial in 1:Ntrials) {
        idx <- idxs[trial]
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]
        
        modresults <- run_trial(
          weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y,
          Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise,
          platform_x, platform_y, starting_x, starting_y, speed, hitwall
        )
        
        weights <- modresults[[1]]
        
        PMs[1, day, trial, reps] <- modresults[[9]]
        PMs[2, day, trial, reps] <- modresults[[6]]
        PMs[3, day, trial, reps] <- modresults[[8]][4] * 100
        PMs[4, day, trial, reps] <- modresults[[8]][2] * 100
        PMs[5, day, trial, reps] <- modresults[[7]] * 100
      }
    }
  }
}