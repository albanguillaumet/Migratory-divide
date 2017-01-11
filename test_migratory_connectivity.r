# This files contains the code necessary to replicate the simulations presented in Appendix 2 of the paper:
# Guillaumet, A., Dorr, B, Wang, G et al. 2011. Determinants of local and migratory movements of Great Lakes double-crested cormorants. Behavioral Ecology 22, 1096-1103.
# the function test_migratory_connectivity
# 3 simulated data sets and case studies

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(gtools) # needed to create one of the simulated data set (for odd / even)
dev.off() # to avoid possible interference with previous graph

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

test_migratory_connectivity = function(data) {

options(digits = 6)

#---

AICC = function(model) {
LL = logLik(model) ; df = as.numeric(attr(LL, "df")) ; n = attr(LL, "nobs")
aic = as.numeric(-2* LL) + 2*df  + (2*df *(df+1))/(n-df-1)
return(aic) } 

#---

set_flyways = function(data = data, var = data$long_w) {

d = data

n = length(d[,1])

ordination_col = sort(unique(d$long_col), decreasing = TRUE)
n_units = length(ordination_col)	
col_names = as.factor(1:n)
AICC_pattern_col = rbind(ordination_col, AICC = rep(NA, n_units))
colnames(AICC_pattern_col) = col_names 

fly_prov = d$flyway_col
cat = ordination_col 
nloop = nloop_col = length(cat) 

for(i in 1 : nloop) {
fly_prov[d$long_col == cat[1]] = 1 ; cat = cat[-1] 
AICC_pattern_col[2, i] = AICC(lm(var ~ fly_prov))
}

#-

# Fill in the data

cat = ordination_col 
nloop = which(AICC_pattern_col[2,]  == min(AICC_pattern_col[2,]))
last_divide_col = AICC_pattern_col[1, nloop]
divide_col = mean(c(last_divide_col, AICC_pattern_col[1, (nloop+1)]))	
for(i in 1 : nloop) {
d$flyway_col[d$long_col == cat[1]] = 1 ; cat = cat[-1] 
}

#---

# Figure

par(mar=c(4.3, 4, 2, 2) + 0.1, family = "serif", font.main = 1, cex.main = 1, font.sub = 1, font.lab = 1)
x_col = as.numeric(AICC_pattern_col[1,][-nloop_col])
y_col = as.numeric(AICC_pattern_col[2,][-nloop_col]) 
ry = range(y_col) ; rx = range(x_col)
col = rep(c("darkorange3", "darkgreen"), length(x_col))[1:length(x_col)]
for (i in 1: length(x_col)) {
if(x_col[i] < divide_col) col[i]  ="darkgreen"
if(x_col[i] > divide_col) col[i]  ="darkorange3"
}
ch_cex = 1
plot(rx, ry, type = "n", bty = "n", xlab = "Longitude Breeding", ylab = "AICc", 
main = "Position of the putative migratory divide",
cex.axis = ch_cex, cex.lab = ch_cex)
points(x_col, y_col, pch = 19, cex = ch_cex, col = col)
lines(x_col, y_col, col = "black") 
abline(v = divide_col, col = "black", lty = 2) 

#---

return(list(
d = d, 
last_divide_col = as.numeric(last_divide_col), divide_col = as.numeric(divide_col)
))

} # set_flyways 

#---

split.screen(c(2, 1))     
split.screen(c(1, 2), screen = 2) 

screen(1)
u = set_flyways(data = data, var = data$long_w)
data = u$d

col = rep(c("darkorange3", "darkgreen"), n)[1:n]
for (i in 1: n) ifelse(data$long_col[i] < u$divide_col, col[i] <-"darkgreen", col[i] <-"darkorange3")

par(mar=c(2, 2, 2, 2) + 0.1, family = "serif", font.main = 1, cex.main = 1, font.sub = 1, font.lab = 1)
screen(3)
plot(data$long_col, data$long_w, bty="n", las=1, col = col, pch = 19, cex = 1, cex.axis = 1, cex.lab = 1, xlab = "Longitude Breeding", ylab = "Longitude Winter", main = "Parallel migration")
screen(4)
plot(as.factor(data$flyway_col), data$long_w, bty="n", col = c("darkgreen", "darkorange3"), names = c("M", "A"), cex.axis = 1, cex.lab = 1, xlab = "Migratory Flyway", ylab = "Longitude Winter", main = "Flyway system")

m_par = lm(data$long_w ~ data$long_col)
m_fly = lm(data$long_w ~ data$flyway_col)
#summary(m_par)
#summary(m_fly)
aic = cbind(AICC(m_par), AICC(m_fly))

min_AICc = min(aic) ; delta = w = cbind(NA, NA); for(i in 1:2) delta[i] = aic[i] - min_AICc
s = sum( exp(-delta / 2) ) ; for(i in 1 : 2) { w[i] = exp(-delta[i] / 2) / s } 

res = as.data.frame(rbind(aic, w))
colnames(res) = c("Parallel", "Flyway")
rownames(res) = c("AICc", "Akaike Weight")
res

} # test_migratory_connectivity

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 1 - Simulate and test a  data set consistent with a parallel migration system

n = 100; long_col = 1: n
long_w = rep(NA, n) ; for(i in 1 : n) long_w[i] = long_col[i] + rnorm(1, sd = 5)
data_f = as.data.frame(cbind(long_col = long_col, long_w = long_w, flyway_col = rep(0, n)))

test_migratory_connectivity(data = data_f); close.screen(all = TRUE)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 2 - Simulate and test a  data set consistent with a flyway migration system:
	# central migratory divide
	# balanced sample size on each flyway

n = 100; long_col = 1: n
long_w = rep(NA, n) ; for(i in 1 : n) ifelse(i <= 50, long_w[i] <- 25 + rnorm(1, sd = 5), long_w[i] <- 75 + rnorm(1, sd = 5))
data_f = as.data.frame(cbind(long_col = long_col, long_w = long_w, flyway_col = rep(0, n)))

test_migratory_connectivity(data = data_f); close.screen(all = TRUE)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 3 - Simulate and test a more complex data set consistent with a flyway migration system:
	# off-center migratory divide (Longitude Breeding at 3/4 of range)
	#  sample size for flyways M and A = 1/4 and 3/4 of total sample size, respectively
	# Some overlap between flyways due to the proximity of winter locations (long_w)

n = 100; long_col = rep(NA, n); x = 1:75 ; b = rep(c(FALSE, FALSE, TRUE), 25); x = x[b]; for(i in 1 : 25) long_col[i] = x[i]
y1 = 76:100; y2 = y1-1/3; y3 = y1-2/3; z1 = y2[odd(1:50)]; z2 = y3[odd(1:50)]; a = sort(c(y1, y2, y3)); for(i in 26 : 100) long_col[i] = a[i-25]
long_w = rep(NA, n) ; for(i in 1 : n) ifelse(long_col[i] <= 75, long_w[i] <- 46 + rnorm(1, sd = 5), long_w[i] <- 54 + rnorm(1, sd = 5))
data_f = as.data.frame(cbind(long_col = long_col, long_w = long_w, flyway_col = rep(0, n)))

test_migratory_connectivity(data = data_f); close.screen(all = TRUE)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
