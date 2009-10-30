\name{someKFWER-package}
\Rdversion{1.1}
\alias{someKfwer-package}
\docType{package}
\title{
Controlling the k-FWER (Generalized Familywise Error Rate)
}
\description{
This library collects some procedures controlling the Generalized Familywise Error Rate: Lehmann and Romano (2005), Guo and Romano (2007) (single step and stepdown), Finos and Farcomeni (2009).}
\details{
\tabular{ll}{
Package: \tab kfwe\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-10-30\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
	
}
\author{
L. Finos and A. Farcomeni

Maintainer: <livio@stat.unipd.it>
}
\references{
L. Finos and A. Farcomeni (2009) k-FWER control without multiplicity correction, with application to detection
of genetic determinants of multiple sclerosis in Italian twins. University of Padua, Dept Statistical Science. Working Paper #7
}
\keyword{ package }
\examples{
set.seed(13)
y=matrix(rnorm(3000),3,1000)+2  #create toy data
p=apply(y,2,function(y) t.test(y)$p.value)  #compute p-values
M2=apply(y^2,2,mean)   #compute ordering criterion

kord=kfweOrd(p,k=5,ord=M2)  #ordinal procedure
klr=kfweLR(p,k=5)           #Lehaman and Romano
kgr=kfweGR(p,k=5)           #Guo and Romano
}

\keyword{ htest }
