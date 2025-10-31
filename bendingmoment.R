# Beam load and bending moment distribution simulation
# www.overfitting.net
# https://www.overfitting.net/2025/10/simulando-cargas-en-una-viga-con-r.html


beam_bending_moments <- function(L, x, m, g = 9.81,
                                 labels = TRUE, png = TRUE, name = "beam.png") {
    require(Cairo)
    
    # Checks
    if (length(x) != length(m)) stop("x and m must have the same length")
    if (any(x < 0 | x > L)) stop("Positions x must be within the span [0, L]")
    neworder=order(x)  # in case the vector is not ordered
    x=x[neworder]
    m=m[neworder]

    fscale <- 3  # scaling factor to plot forces in N
    
    
    ################################
    # 1. CALCULATIONS
    
    # Force calculation
    P <- m * g
    N <- length(x)
    
    # Support reactions (linear effects)
    RA <- sum(P * (L - x)) / L
    RB <- sum(P * x) / L

    # Bending moments at each load
    M <- numeric(N)
    for (i in 1:N) {
        # Forces to the left of point x[i]
        idx_left <- which(x < x[i])
        M[i] <- RA * x[i] - sum(P[idx_left] * (x[i] - x[idx_left]))
    }
    
    
    ################################
    # 2. PLOTTING
    
    if (png) CairoPNG(name, width=1280, height=800)

    xlimits=c(-0.2, L+0.2)
    ylimits=c(-0.05 * max(RA,RB) * fscale, max(M)) * 1.2
    # ylimits=c(-30, 80)  # Example 2: moving loads (animation)
    # ylimits=c(-8, 23)  # Example 3: load distributions
    plot(xlimits, ylimits, type = "n", axes = TRUE,
         xlab = "Position along the beam (m)", ylab = "Bending moment (N·m)",
         main = "Load and bending moment distribution",
         cex.main = 2.5, cex.lab  = 1.8, cex.axis = 1.8)
    
    # Beam
    segments(0, 0, L, 0, lwd = 16, col = rgb(0.5, 0.5, 0.5, 0.5), lend = "butt")
    text(L/2, ylimits[1]/15, paste0("L=", round(L,1), "m beam"), cex = 1.8, pos = 1)
    abline(h = 0, lty = 2)
    
    # Forces
    # Loads
    arrows(x, 0.05 * P * fscale, x, 0, length = 0.3, col = "red", lwd = 2)
    if (labels) text(x, 0.5 * m * fscale, paste0(round(P,1), "N\n(m=", round(m,1), "kg)"),
         col = "red", cex = 1.8, pos = 3)
    
    # Reactions
    arrows(c(0,L), -0.05 * c(RA, RB) * fscale, c(0,L), c(0,0),
           length = 0.3, col = "darkgreen", lwd = 2)
    if (labels) text(c(0, L), -0.05 * c(RA, RB) * fscale,
         c(paste0("RA=", round(RA,1), "N"), paste0("RB=", round(RB,1), "N")),
         col = "darkgreen", cex = 1.8, pos = 1)
    
    # Bending moments
    points(x, M, pch = 19, col = "blue")
    if (labels) text(x, M, paste0(round(M,1), "N·m"), pos = 3, col = "blue", cex = 1.8)
    segments(0, 0, x[1], M[1], col = 'blue', lty = 'dotted')
    for (i in 1:(N-1)) {
        segments(x[i], M[i], x[i+1], M[i+1], col = 'blue', lty = 'dotted')
    }
    segments(x[N], M[N], L, 0, col = 'blue', lty = 'dotted')
    
    if (png) dev.off()
    
    
    return(data.frame(x, m, P, M))
}


# Example 1: check https://www.youtube.com/watch?v=fJOg62E2754
L <- 3.2
x <- cumsum(c(0.6, 0.9, 1.5))
m <- c(40000, 32000, 16000)/9.81
results <- beam_bending_moments(L, x, m, name="beamcheck.png")
print(results)


# Example 2: moving loads (animation)
L <- 2.0
m <- c(9, 4, 6, 12)
NFRAMES=24*4
for (frame in 0:(NFRAMES-1)) {
    name=paste0("beam_", ifelse(frame<10, "0", ""), frame, ".png")
    print(name)
    x <- c((0.3+sin(2*pi*frame/NFRAMES)*0.3/2), (0.6+sin(2*pi*frame/NFRAMES+pi/4)*0.6/2),
           (1.2+sin(2*pi*frame/NFRAMES+pi/2)*0.8/2), (1.8+sin(2*pi*frame/NFRAMES+3*pi/4)*0.2/2))
    results <- beam_bending_moments(L, x, m,
                                    labels=TRUE, png=TRUE, name=name)  
}


# Example 3: load distributions

# Uniformly distributed loads
N=11
L <- 2.0
x <- seq(0, L, length.out=N)
m <- rep(1, N)
results <- beam_bending_moments(L, x, m, labels=FALSE, png=TRUE, name="beamuniform.png")
print(results)

# Optimal distributed loads
N=6
L <- 2.0
x <- seq(0, L, length.out=N)
m <- c(seq(2, 1, length.out=N/2), seq(1, 2, length.out=N/2))
results <- beam_bending_moments(L, x, m, labels=TRUE, png=TRUE, name="beamoptimum.png")
print(results)

# Critically distributed loads
N=6
L <- 2.0
x <- seq(0, L, length.out=N)
m <- c(seq(1, 2, length.out=N/2), seq(2, 1, length.out=N/2))
results <- beam_bending_moments(L, x, m, labels=TRUE, png=TRUE, name="beamcritical.png")
print(results)

# Empty centre
N=6
L <- 2.0
x <- c(0, 0.2, 0.4, 2-0.4, 2-0.2, 2)  # seq(0, L, length.out=N)
m <- c(seq(2, 1, length.out=N/2), seq(1, 2, length.out=N/2))
results <- beam_bending_moments(L, x, m, labels=TRUE, png=TRUE, name="beamprotect.png")
print(results)



