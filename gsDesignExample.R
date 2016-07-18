# Group sequential design with Pocock and O'Brien-Fleming boundaries


library(gsDesign)

# two-sided Pocock design
# k is the number of interim analysis
x <- gsDesign(k = 5, test.type=2, sfu="Pocock")
# print the sample size ratio
x$n.I

# two-sided O'Brien-Fleming design
gsDesign(k = 5, test.type=2, sfu="OF", alpha = 0.025, beta = 0.05)

# calculate the sample size for non-group sequential design
n.fix <- nBinomial(p1=0.3, p2=0.15, alpha = 0.025, beta = 0.1)
# print the fixed sample size
n.fix

# calculate group sequential sample sizes at each interim
n.GS <- n.fix*x$n.I
# print the sample sizes
n.GS

# direct calculation of GS design
x.new <- gsDesign(k = 5, test.type=2, sfu="OF", n.fix=n.fix)
# print the design
x.new
class(x.new)

# Make the boundary plot
print(plot(x.new, plottype=1))

print(plot(x.new, plottype=2))
print(plot(x.new, plottype=3))
print(plot(x.new, plottype=4))
print(plot(x.new, plottype=6))
