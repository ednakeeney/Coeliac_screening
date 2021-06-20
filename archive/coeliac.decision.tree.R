# Coeliac screening decision tree

# Number of PSA samples
n.samples<-100

# Number and names of strategies
n.treat<-3
t.names<-c("Test","Test + biopsy", "Double test")


#prevalence of coeliac disease
p.cd <- 0.5 

# Accuracy of test or test + biopsy
tp <- fn <- fp <- tn <- matrix(nrow=n.samples, ncol=n.treat)

#IgA EMA sensitivity in adults: 88.0 (75.2, 94.7)
Sens_test <- 0.88
SensSE <- (0.947 - 0.752)/3.92
SensAlpha <- (Sens^2*(1-Sens)/SensSE^2)-Sens
SensBeta <- (SensAlpha/Sens)-SensAlpha
#IgA EMA specificity in adults: 99.6 (92.3, 100.0)
Spec_test <- 0.996
SpecSE <- (1 - 0.923)/3.92 
SpecAlpha <- (Spec^2*(1-Spec)/SpecSE^2)-Spec
SpecBeta <- (SpecAlpha/Spec)-SpecAlpha

Sens_testbiopsy <- 1
Spec_testbiopsy <- 1
Sens_doubletest <- 1
Spec_doubletest <- 1

# Probabilities for test 
tp[,1]<- (n.samples*p.cd*rbeta(n=n.samples, shape1 = SensAlpha, shape2 = SensBeta))/n.samples
fn[,1]<- ((n.samples*p.cd) - tp[,1])/n.samples  
tn[,1] <- (n.samples*p.cd*rbeta(n=n.samples, shape1 = SpecAlpha, shape2 = SpecBeta))/n.samples
fp[,1] <- ((n.samples*p.cd) - tn[,1])/n.samples

# Probabilities for test + biopsy
tp[,2] <- (n.samples*p.cd*Sens_testbiopsy)/n.samples
fn[,2]<- ((n.samples*p.cd) - tp[,2])/n.samples
tn[,2] <- (n.samples*p.cd*Spec_testbiopsy)/n.samples
fp[,2] <- ((n.samples*p.cd) - tn[,2])/n.samples

# Probabilities for Double test
tp[,3] <- (n.samples*p.cd*Sens_doubletest)/n.samples
fn[,3]<- ((n.samples*p.cd) - tp[,3])/n.samples
tn[,3] <- (n.samples*p.cd*Spec_doubletest)/n.samples
fp[,3] <- ((n.samples*p.cd) - tn[,3])/n.samples

# probability start on gfd when diagnosed (adherence)
p.gfd <- 0.5

# Costs for FP - Cost of GFP for one year
c.fp<-rnorm(n=n.samples, mean=1000, sd=50)


# Cost of test, test + biopsy or double test
c.test<-t(matrix(rep(c(300,500, 400),n.samples),ncol=n.samples,nrow=3))


