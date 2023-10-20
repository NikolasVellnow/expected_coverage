rm(list=ls())
library(ggplot2)
library(viridis)
library(readxl)

# mac
samples_path <- "~/sciebo/Bioinformatik Bachelor Bielefeld/7_SS_2023/bachelor_project/"

# windows
#samples_path <- "C:/Users/nvellnow/sciebo7/7_SS_2023/bachelor_project/"

#### Model 1 (null model) functions #######

# function to calculate coverage (eq.(2))
coverage_model1 <- function(L, N_m, G){
    c <- (L*N_m)/G
}

log2_c_ratio <- function(N_mA, N_mB, p){
  ratio <- log2(N_mA/N_mB)
}


##### Model 1 example #######

# prepare variables
read_len <- 150
G_RG <- 1.0*10^6
p_value <- 0.3

N_m <- 4.0*10^5

# calculate coverage according to model 1
c <- coverage_model1(read_len, N_m, G_RG)
c


#### Model 2 functions #######

# function to calculate number RG reads that map to RG, eq.(5a))
calc_n_RG_mapped_to_RG_m2 <- function(N_m, p){
    n_RG_mapped_to_RG_m2 <- N_m/(1+p)
}

# function to number of GRC reads that map to RG (eq.(5b))
calc_n_GRC_mapped_to_RG_m2 <- function(N_m, p){
    n_GRC_mapped_to_RG_m2 <- N_m/(1+(1/p))
}

# function to calculate exp.coverage for region where only RG maps to (eq. (7))
coverage_RG_mapped_to_RG_m2 <- function(L, G_RG, N_m, p){
    n_RG_mapped_to_RG_m2 <- calc_n_RG_mapped_to_RG_m2(N_m, p)
    c_RG_mapped_to_RG_m2 <- (L * n_RG_mapped_to_RG_m2) / G_RG
}

# function to calculate exp. coverage for region where RG and GRC map to (eq. (8))
coverage_GRC_mapped_to_RG_m2 <- function(L, G_RG, N_mapped, p){
    n_RG_mapped_to_RG_m2 <- calc_n_RG_mapped_to_RG_m2(N_mapped, p)
    n_GRC_mapped_to_RG_m2 <- calc_n_GRC_mapped_to_RG_m2(N_mapped, p)
    c_GRC_mapped_to_RG_m2 <- (L*((p*n_RG_mapped_to_RG_m2)+n_GRC_mapped_to_RG_m2))/(p*G_RG)
}


# function to calculate exp. log-ratio of sample A (model 2) over sample B (model 1) (eq. (11))
log2_c_ratio_RG_m2 <- function(N_mA, N_mB, p){
    result <- log2((1/(1+p)) * (N_mA/N_mB))
}

# function to calculate eq. (12)
log2_c_ratio_GRC_m2 <- function(N_mA, N_mB, p){
    result <- log2((2/(1+p)) * (N_mA/N_mB))
}

##### Model 2 example #######

# prepare variables
read_len <- 150
G_R <- 1.0*10^6
p_value <- 0.3

N_m <- 4.0*10^5


# calculate number of mapped reads according to model 2 
n_R_m2 <- calc_n_RG_mapped_to_RG_m2(N_m, p_value)
n_R_m2

n_G_m2 <- calc_n_GRC_mapped_to_RG_m2(N_m, p_value)
n_G_m2

# coverage for region where only RG sequences map to according to model 2
c_R_m2  <- coverage_RG_mapped_to_RG_m2(read_len, G_R, N_m, p_value)
c_R_m2


# coverage for region where RG and GRC sequences map to according to model 2
c_G_m2  <- coverage_GRC_mapped_to_RG_m2(read_len, G_R, N_m, p_value)
c_G_m2


#### plot expected coverage of example model2

# range of values for p (proportion of GRC of RG)
p <- seq(0.01,1,0.01)


m2 <- data.frame(p)
m2$c_R_m2 <- coverage_RG_mapped_to_RG_m2(read_len, G_R, N_m, m2$p)
m2$c_G_m2 <- coverage_GRC_mapped_to_RG_m2(L=read_len, G_RG=G_R, N_m=N_m, m2$p)


plot(m2$p, m2$c_R_m2, ylim=c(0,130),
     main="Model 2",
     type="l", col="blue",
     xlab="p (relative size of GRC)", ylab="Expected coverage")
lines(m2$p, m2$c_G, lty=2, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.70,130,legend=c("RG-region","GRC-region"), col=c("blue","red"),
       lty=c(1,2), ncol=1)



#### Model 3 functions #######

# function to calculate exp.coverage for region where only RG maps to (eq. (14))
coverage_RG_mapped_to_RG_m3 <- function(L=L, G_RG=G_RG, N_m=N_m, p=p, g=g){
    c_RG_mapped_to_RG_m3 <- (L / G_RG) * ((N_m) * (1 + p - g * p)/(1 + p))
    c_RG_mapped_to_RG_m3
}

# function to calculate exp.coverage for region where RG and GRC map to (eq. (15)))
coverage_GRC_mapped_to_RG_m3 <- function(L, G_RG, N_m, p, g){
    c_GRC_mapped_to_RG_m3 <- (L / G_RG) * ((2 * g * N_m)/(1+p) + N_m - g * N_m)
    c_GRC_mapped_to_RG_m3
}

# function to calculate exp. log-ratio of sample A (model 3) over sample B (model 1) (eq. (16))
log2_c_ratio_RG_m3 <- function(N_mA, N_mB, p, g){
    ratio <- (N_mA / N_mB) * ((1 + p - g * p)/(1 + p))
    result <- log2(ratio)
    result
}

# function to calculate exp. log-ratio of sample A (model 3) over sample B (model 1) (eq. (17))
log2_c_ratio_GRC_m3 <- function(N_mA, N_mB, p, g){
    ratio <- (N_mA / N_mB) * ((2 * g)/(1+p) + 1 - g)
    result <- log2(ratio)
    result
}

##### Model 3 example #######

# prepare variables
p <- seq(0.01,1,0.01)
g <- seq(0.01,1,0.01)

read_len <- 150
G_R <- 1.0*10^6
p_value <- 0.3
g_value <- 0.5

N_m <- 4.0*10^5


# coverage for region where only RG sequences map to according to model 2
c_R_m3  <- coverage_RG_mapped_to_RG_m3(read_len, G_RG, N_m, p_value, g_value)
c_R_m3


# coverage for region where RG and GRC sequences map to according to model 2
c_G_m3  <- coverage_GRC_mapped_to_RG_m3(read_len, G_RG, N_m=N_m, p=p_value, g=g_value)
c_G_m3



m3 <- expand.grid(p = p, g = g)

m3$c_RG_m3 <- coverage_RG_mapped_to_RG_m3(L=read_len, G_RG=G_R, N_m=N_m, m3$p, m3$g)
m3$c_GRC_m3 <- coverage_GRC_mapped_to_RG_m3(L=read_len, G_RG=G_R, N_m=N_m, m3$p, m3$g)


#### plot expected coverage

plot(m3$g[m3$p=="0.1"], m3$c_RG_m3[m3$p=="0.1"], ylim=c(40,120), xlim=c(0,1),
     main="Model 3",
     type="l",
     lty=1,
     col="blue",
     xlab="g (proportion of PGC)", ylab="Expected coverage")
lines(m3$g[m3$p=="0.1"], m3$c_GRC_m3[m3$p=="0.1"], lty=1, col="red")
lines(m3$g[m3$p=="0.2"], m3$c_RG_m3[m3$p=="0.2"], lty=2, col="blue")
lines(m3$g[m3$p=="0.2"], m3$c_GRC_m3[m3$p=="0.2"], lty=2, col="red")
lines(m3$g[m3$p=="0.3"], m3$c_RG_m3[m3$p=="0.3"], lty=3, col="blue")
lines(m3$g[m3$p=="0.3"], m3$c_GRC_m3[m3$p=="0.3"], lty=3, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.0,120,legend=c("RG-region, p=0.1",
                        "GRC-region, p=0.1",
                        "RG-region, p=0.2",
                        "GRC-region, p=0.2",
                        "RG-region, p=0.3",
                        "GRC-region, p=0.3"),
       col=c("blue",
             "red",
             "blue",
             "red",
             "blue",
             "red"),
       lty=c(1,
             1,
             2,
             2,
             3,
             3), ncol=1)


#### Model 4 functions #######

# function to calculate exp.coverage for region where only RG maps to
coverage_RG_mapped_to_RG_m4 <- function(L=L, G_RG=G_RG, N_m=N_m, p=p){
    c_RG_mapped_to_RG_m4 <- (L / G_RG) * (N_m/(1 + p))
    c_RG_mapped_to_RG_m4
}

# function to calculate exp.coverage for region where RG and GRC map to (eq. (15)))
coverage_GRC_mapped_to_RG_m4 <- function(L, G_RG, N_m, p, d){
    c_GRC_mapped_to_RG_m4 <- (L / G_RG) * (N_m/(1+p)) * (1 + d)
    c_GRC_mapped_to_RG_m4
}

# function to calculate exp. log-ratio of sample A (model 4) over sample B (model 1)
log2_c_ratio_RG_m4 <- function(N_mA, N_mB, p){
    result <- log2((1/(1+p)) * (N_mA/N_mB))
}


# function to calculate exp. log-ratio of sample A (model 3) over sample B (model 1) (eq. (17))
log2_c_ratio_GRC_m4 <- function(N_mA, N_mB, p, d){
    ratio <- (N_mA / N_mB) * ((1+d) / (1+p))
    result <- log2(ratio)
    result
}


###### Model 4 example ######

# prepare variables
p <- seq(0.01,1,0.01)
d <- seq(1,100,1)

read_len <- 150
G_R <- 1.0*10^6
p_value <- 0.3
d_value <- 5

N_m <- 4.0*10^5


# coverage for region where only RG sequences map to according to model 4
c_R_m4  <- coverage_RG_mapped_to_RG_m4(read_len, G_RG, N_m, p_value)
c_R_m4


# coverage for region where RG and GRC sequences map to according to model 2
c_G_m4  <- coverage_GRC_mapped_to_RG_m4(read_len, G_RG, N_m=N_m, p=p_value, d=d_value)
c_G_m4



m4 <- expand.grid(p = p, d = d)

m4$c_RG_m4 <- coverage_RG_mapped_to_RG_m4(L=read_len, G_RG=G_R, N_m=N_m, m4$p)
m4$c_GRC_m4 <- coverage_GRC_mapped_to_RG_m4(L=read_len, G_RG=G_R, N_m=N_m, m4$p, m4$d)


#### plot expected coverage

plot(m4$d[m4$p=="0.1"], m4$c_RG_m4[m4$p=="0.1"],
     ylim=c(40,650),
     xlim=c(1,10),
     main="Model 4",
     type="l",
     lty=1,
     col="blue",
     xlab="d (number of times that GRC maps)", ylab="Expected coverage")
lines(m4$d[m4$p=="0.1"], m4$c_GRC_m4[m4$p=="0.1"], lty=1, col="red")
lines(m4$d[m4$p=="0.2"], m4$c_RG_m4[m4$p=="0.2"], lty=2, col="blue")
lines(m4$d[m4$p=="0.2"], m4$c_GRC_m4[m4$p=="0.2"], lty=2, col="red")
lines(m4$d[m4$p=="0.3"], m4$c_RG_m4[m4$p=="0.3"], lty=3, col="blue")
lines(m4$d[m4$p=="0.3"], m4$c_GRC_m4[m4$p=="0.3"], lty=3, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(1,660,legend=c("RG-region, p=0.1",
                        "GRC-region, p=0.1",
                        "RG-region, p=0.2",
                        "GRC-region, p=0.2",
                        "RG-region, p=0.3",
                        "GRC-region, p=0.3"),
       col=c("blue",
             "red",
             "blue",
             "red",
             "blue",
             "red"),
       lty=c(1,
             1,
             2,
             2,
             3,
             3), ncol=1)




#### Comparison of samples ######

##### Model 2 vs. Model 1 #####

p <- seq(0.01,1,0.01)
read_len <- 150
G_R <- 1.0*10^6
p_value <- 0.3
d_value <- 5

N_mA <- 4.0*10^5
N_mB <- 4.0*10^5


m2 <- data.frame(p)
m2$log2_c_ratio_RG_m2 <- log2_c_ratio_RG_m2(N_mA = N_mA, N_mB = N_mB, m2$p)
m2$log2_c_ratio_GRC_m2 <- log2_c_ratio_GRC_m2(N_mA = N_mA, N_mB = N_mB, m2$p)

plot(m2$p, m2$log2_c_ratio_RG_m2, ylim=c(-2,2),
     main="Model 2 vs. Model 1",
     type="l", col="blue",
     xlab="p (relative size of GRC)", ylab="log2-ratio (Model 2 / Model 1)")
lines(m2$p, m2$log2_c_ratio_GRC_m2, lty=2, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.72,2,legend=c("RG-region","GRC-region"), col=c("blue","red"),
       lty=c(1,2), ncol=1)



##### Model 3 vs. Model 1 #####

p <- seq(0.01,1,0.01)
g <- seq(0.01,1,0.01)
read_len <- 150
G_R <- 1.0*10^6
p_value <- 0.3
d_value <- 5
N_mA <- 4.0*10^5
N_mB <- 4.0*10^5

m3 <- expand.grid(p = p, g = g)

m3$log2_c_ratio_RG_m3 <- log2_c_ratio_RG_m3(N_mA = N_mA, N_mB = N_mB, m3$p, m3$g)
m3$log2_c_ratio_GRC_m3 <- log2_c_ratio_GRC_m3(N_mA = N_mA, N_mB = N_mB, m3$p, m3$g)


# log2-ratio

plot(m3$g[m3$p=="0.1"], m3$log2_c_ratio_RG_m3[m3$p=="0.1"], ylim=c(-1,1.5), xlim=c(0,1),
     main="Model 3 vs. Model 1",
     type="l",
     lty=1,
     col="blue",
     xlab="g (proportion of PGC)", ylab="log2-ratio (Model 3 / Model 1)")
lines(m3$g[m3$p=="0.1"], m3$log2_c_ratio_GRC_m3[m3$p=="0.1"], lty=1, col="red")
lines(m3$g[m3$p=="0.2"], m3$log2_c_ratio_RG_m3[m3$p=="0.2"], lty=2, col="blue")
lines(m3$g[m3$p=="0.2"], m3$log2_c_ratio_GRC_m3[m3$p=="0.2"], lty=2, col="red")
lines(m3$g[m3$p=="0.4"], m3$log2_c_ratio_RG_m3[m3$p=="0.4"], lty=3, col="blue")
lines(m3$g[m3$p=="0.4"], m3$log2_c_ratio_GRC_m3[m3$p=="0.4"], lty=3, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.0,1.5,legend=c("RG-region, p=0.1",
                        "GRC-region, p=0.1",
                         "RG-region, p=0.2",
                        "GRC-region, p=0.2",
                         "RG-region, p=0.4",
                        "GRC-region, p=0.4"
                         ),
       col=c("blue",
             "red",
             "blue",
             "red",
             "blue",
             "red"),
       lty=c(1,
             1,
             2,
             2,
             3,
             3), ncol=1)


##### Model 3 vs. Model 1 #####

p <- seq(0.01,1,0.01)
g <- seq(0.01,1,0.01)
read_len <- 150
G_R <- 1.0*10^6
p_value <- 0.3
d_value <- 5
N_mA <- 4.0*10^5
N_mB <- 4.0*10^5

m3 <- expand.grid(p = p, g = g)

m4$log2_c_ratio_RG_m3 <- log2_c_ratio_RG_m3(N_mA = N_mA, N_mB = N_mB, m3$p, m3$g)
m4$log2_c_ratio_GRC_m3 <- log2_c_ratio_GRC_m3(N_mA = N_mA, N_mB = N_mB, m3$p, m3$g)


# log2-ratio

plot(m3$g[m3$p=="0.1"], m3$log2_c_ratio_RG_m3[m3$p=="0.1"], ylim=c(-1,1.5), xlim=c(0,1),
     main="Model 3 vs. Model 1",
     type="l",
     lty=1,
     col="blue",
     xlab="g (proportion of PGC)", ylab="log2-ratio (Model 3 / Model 1)")
lines(m3$g[m3$p=="0.1"], m3$log2_c_ratio_GRC_m3[m3$p=="0.1"], lty=1, col="red")
lines(m3$g[m3$p=="0.2"], m3$log2_c_ratio_RG_m3[m3$p=="0.2"], lty=2, col="blue")
lines(m3$g[m3$p=="0.2"], m3$log2_c_ratio_GRC_m3[m3$p=="0.2"], lty=2, col="red")
lines(m3$g[m3$p=="0.4"], m3$log2_c_ratio_RG_m3[m3$p=="0.4"], lty=3, col="blue")
lines(m3$g[m3$p=="0.4"], m3$log2_c_ratio_GRC_m3[m3$p=="0.4"], lty=3, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.0,1.5,legend=c("RG-region, p=0.1",
                        "GRC-region, p=0.1",
                        "RG-region, p=0.2",
                        "GRC-region, p=0.2",
                        "RG-region, p=0.4",
                        "GRC-region, p=0.4"
),
col=c("blue",
      "red",
      "blue",
      "red",
      "blue",
      "red"),
lty=c(1,
      1,
      2,
      2,
      3,
      3), ncol=1)




##### Model 4 vs. Model 1 #####

p <- seq(0.01,1,0.01)
d <- seq(1,10,1)
read_len <- 150
G_R <- 1.0*10^6
N_mA <- 4.0*10^5
N_mB <- 4.0*10^5

m4 <- expand.grid(p = p, d = d)

m4$log2_c_ratio_RG_m4 <- log2_c_ratio_RG_m4(N_mA = N_mA, N_mB = N_mB, m4$p)
m4$log2_c_ratio_GRC_m4 <- log2_c_ratio_GRC_m4(N_mA = N_mA, N_mB = N_mB, m4$p, m4$d)


# log2-ratio

plot(m4$d[m4$p=="0.1"], m4$log2_c_ratio_RG_m4[m4$p=="0.1"], ylim=c(-0.5,4.5), xlim=c(1,10),
     main="Model 4 vs. Model 1",
     type="l",
     lty=1,
     col="blue",
     xlab="d (number of times GRC maps)", ylab="log2-ratio (Model 4 / Model 1)")
lines(m4$d[m4$p=="0.1"], m4$log2_c_ratio_GRC_m4[m4$p=="0.1"], lty=1, col="red")
lines(m4$d[m4$p=="0.2"], m4$log2_c_ratio_RG_m4[m4$p=="0.2"], lty=2, col="blue")
lines(m4$d[m4$p=="0.2"], m4$log2_c_ratio_GRC_m4[m4$p=="0.2"], lty=2, col="red")
lines(m4$d[m4$p=="0.4"], m4$log2_c_ratio_RG_m4[m4$p=="0.4"], lty=3, col="blue")
lines(m4$d[m4$p=="0.4"], m4$log2_c_ratio_GRC_m4[m4$p=="0.4"], lty=3, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(1.0,4.55,legend=c("RG-region, p=0.1",
                        "GRC-region, p=0.1",
                        "RG-region, p=0.2",
                        "GRC-region, p=0.2",
                        "RG-region, p=0.4",
                        "GRC-region, p=0.4"
),
col=c("blue",
      "red",
      "blue",
      "red",
      "blue",
      "red"),
lty=c(1,
      1,
      2,
      2,
      3,
      3), ncol=1)



###### Preparing data for example calculations ######

samples <- read_excel(paste0(samples_path, "sample_summary_data.xlsx"),range = "A1:M22")

samples$sample <- as.factor(samples$sample)
samples$sample_type <- as.factor(samples$sample_type)

summary(samples)

# prepare variables
read_len <- 150
G_RG <- 1.5*10^9
p_value <- 0.2

N_m_s20 <- samples$num_reads_mapped[samples$sample=='s20']
N_m_s4 <- samples$num_reads_mapped[samples$sample=='s4']

# calculate coverage according to model 1
c_s20 <- coverage_model1(read_len, N_m_s20, G_RG)
c_s20

c_s4 <- coverage_model1(read_len, N_m_s4, G_RG)
c_s4

# calculate number of mapped reads according to model 2 
n_RG_mapped_to_RG_s20_m2 <- calc_n_RG_mapped_to_RG_m2(N_m_s20, p_value)
n_RG_mapped_to_RG_s20_m2

n_GRC_mapped_to_RG_s20_m2 <- calc_n_GRC_mapped_to_RG_m2(N_m_s20, p_value)
n_GRC_mapped_to_RG_s20_m2

# coverage for region where only RG sequences map to according to model 2
c_RG_mapped_to_RG_m2  <- coverage_RG_mapped_to_RG_m2(read_len, G_RG, N_m_s20, p_value)
c_RG_mapped_to_RG_m2


# coverage for region where RG and GRC sequences map to according to model 2
c_GRC_mapped_to_RG_m2  <- coverage_GRC_mapped_to_RG_m2(read_len, G_RG, N_m_s20, p_value)
c_GRC_mapped_to_RG_m2



# range of values for p (proportion of GRC of RG)
p <- seq(0.01,1,0.01)
# range of values for g (proportion of PGC in sample)
g <- seq(0.01,2,0.01)


#### Model 2

m3 <- data.frame(p)
m3$c_RG_m2 <- coverage_RG_mapped_to_RG_m2(L=read_len, G_RG=G_RG, N_m=N_m_s20, df$p)
m3$c_GRC_m2 <- coverage_GRC_mapped_to_RG_m2(L=read_len, G_RG=G_RG, N_m=N_m_s20, df$p)





#### plot expected coverage of s20

plot(df$p, df$c_RG, ylim=c(0,100),
     main="Sample s20 (Model 2)",
     type="l", col="blue",
     xlab="p (relative size of GRC)", ylab="Expected coverage")
lines(df$p, df$c_GRC, lty=2, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.7,98,legend=c("RG-region","GRC-region"), col=c("blue","red"),
       lty=c(1,2), ncol=1)


#### plot expected coverage log_2-ratio (s20/s4)

df$log2_c_ratio_RG_m2 <- log2_c_ratio_RG_m2(N_m_s20, N_m_s4, df$p)
df$log2_c_ratio_GRC_m2 <- log2_c_ratio_GRC_m2(N_m_s20, N_m_s4, df$p)

plot(df$p, df$log2_c_ratio_RG_m2, ylim=c(-2,2),
     main="s20 (Model 2) and s4 (Model 1)",
     type="l", col="blue",
     xlab="p (relative size of GRC)", ylab="log2-ratio (s20/s4)")
lines(df$p, df$log2_c_ratio_GRC_m2, lty=2, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.7,2,legend=c("RG-region","GRC-region"), col=c("blue","red"),
       lty=c(1,2), ncol=1)

#### Model 3

df <- expand.grid(p = p, g = g)
df$c_RG_m2 <- coverage_RG_mapped_to_RG_m2(L=read_len, G_RG=G_RG, N_m=N_m_s20, df$p)
df$c_GRC_m2 <- coverage_GRC_mapped_to_RG_m2(L=read_len, G_RG=G_RG, N_m=N_m_s20, df$p)

df$c_RG_m3 <- coverage_RG_mapped_to_RG_m3(L=read_len, G_RG=G_RG, N_m=N_m_s20, df$p, df$g)
df$c_GRC_m3 <- coverage_GRC_mapped_to_RG_m3(L=read_len, G_RG=G_RG, N_m=N_m_s20, df$p, df$g)

df$log2_c_ratio_RG_m2 <- log2_c_ratio_RG_m2(N_m_s20, N_m_s4, df$p)
df$log2_c_ratio_GRC_m2 <- log2_c_ratio_GRC_m2(N_m_s20, N_m_s4, df$p)

df$log2_c_ratio_RG_m3 <- log2_c_ratio_RG_m3(N_m_s20, N_m_s4, df$p, df$g)
df$log2_c_ratio_GRC_m3 <- log2_c_ratio_GRC_m3(N_m_s20, N_m_s4, df$p, df$g)


#### plot expected coverage of s20

plot(df$g[df$p=="0.1"], df$c_RG_m3[df$p=="0.1"], ylim=c(0,100), xlim=c(0,2),
     main="Sample s20 (Model 2)",
     type="l",
     lty=1,
     col="blue",
     xlab="g (proportion of PGC)", ylab="Expected coverage")
lines(df$g[df$p=="0.1"], df$c_GRC_m3[df$p=="0.1"], lty=1, col="red")
lines(df$g[df$p=="0.2"], df$c_RG_m3[df$p=="0.2"], lty=2, col="blue")
lines(df$g[df$p=="0.2"], df$c_GRC_m3[df$p=="0.2"], lty=2, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.6,103,legend=c("RG-region, p=0.1",
                       "RG-region, p=0.2",
                       "GRC-region, p=0.1",
                       "GRC-region, p=0.2"),
       col=c("blue",
             "blue",
             "red",
             "red"),
       lty=c(1,
             2,
             1,
             2), ncol=1)


# log2-ratio

plot(df$g[df$p=="0.1"], df$log2_c_ratio_RG_m3[df$p=="0.1"], ylim=c(-1,2), xlim=c(0,2),
     main="s20 (Model 2) and s4 (Model 1)",
     type="l",
     lty=1,
     col="blue",
     xlab="g (proportion of PGC)", ylab="log2-ratio (s20/s4)")
lines(df$g[df$p=="0.1"], df$log2_c_ratio_GRC_m3[df$p=="0.1"], lty=1, col="red")
lines(df$g[df$p=="0.2"], df$log2_c_ratio_RG_m3[df$p=="0.2"], lty=2, col="blue")
lines(df$g[df$p=="0.2"], df$log2_c_ratio_GRC_m3[df$p=="0.2"], lty=2, col="red")
lines(df$g[df$p=="0.4"], df$log2_c_ratio_RG_m3[df$p=="0.4"], lty=3, col="blue")
lines(df$g[df$p=="0.4"], df$log2_c_ratio_GRC_m3[df$p=="0.4"], lty=3, col="red")
lines(df$g[df$p=="0.6"], df$log2_c_ratio_RG_m3[df$p=="0.6"], lty=4, col="blue")
lines(df$g[df$p=="0.6"], df$log2_c_ratio_GRC_m3[df$p=="0.6"], lty=4, col="red")

#add a legend in top left corner of chart at (x, y) coordinates = (1, 19)
legend(0.6,2.05,legend=c("RG-region, p=0.1",
                        "RG-region, p=0.2",
                        "RG-region, p=0.4",
                        "RG-region, p=0.6",
                        "GRC-region, p=0.1",
                        "GRC-region, p=0.2",
                        "GRC-region, p=0.4",
                        "GRC-region, p=0.6"),
       col=c("blue",
             "blue",
             "blue",
             "blue",
             "red",
             "red",
             "red",
             "red"),
       lty=c(1,
             2,
             3,
             4,
             1,
             2,
             3,
             4), ncol=1)

library(manipulate)

# expected coverage
manipulate(
    plot(g, coverage_RG_mapped_to_RG_m3(L=read_len, G_RG=G_RG, N_m=N_m_s20, p, g),
         ylim=c(0, 60)),
    p = slider(0,1))

# expected log2-ratio
manipulate(
    plot(g, log2_c_ratio_GRC_m3(N_m_s20, N_m_s4, p, g),
         ylim=c(-1,1)),
    p = slider(0,1))



####### sdf ####







