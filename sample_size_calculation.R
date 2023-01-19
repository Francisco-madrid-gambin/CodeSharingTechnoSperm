# 
# example: Can height, age, and time spent at the gym, predict
# weight in adult males?
#   • H0=0, H1≠0
# 
# • You don’t have background info, so you guess that there is a
# medium effect size
# 
# • For f2-tests:
#   0.02=small, 0.15=medium, and 0.35 large effect sizes
# • Numerator degrees of freedom is the number of predictor
# variables (3)
# 
# • Output will be denominator degrees of freedom rather than
# sample size; will need to round up and add the total number
# of variables (4)
# 
# R Code: pwer -> pwr.f2.test

# pwr.f2.test(u =, v= , f2=, sig.level =, power = )
# • u=numerator degrees of freedom
# • v=denominator degrees of freedom
# • f2=effect size
# • sig.level=significant level
# • power=power of test

library(pwr)
pwr.f2.test(
  u = 1, #degrees of freedom for numerator (how many predictors)
  # v = , # degrees of freedom for denominator (en univarite linear se redondea y se suman ¿+2? y esto es el sample size)
  f2 = 0.35, # effect size  (raiz cuadrada de R2 in a model sqrt(x))
  sig.level = 0.05, 
  power = 0.6 #Power of test (1 minus Type II error probability)
)

