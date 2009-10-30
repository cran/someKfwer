\name{Procedures controlling the k-FWER (Generalized Familywise Error Rate)}
\Rdversion{1.1}
\alias{kfweOrd}
\alias{kfweLR}
\alias{kfweGR}
\alias{getAlpha}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Controlling the Generalized Familywise Error Rate 
}
\description{
This library collects some procedures controlling the Generalized Familywise Error Rate: Lehmannn and Romano (2005), Guo and Romano (2007) (single step and stepdown), Finos and Farcomeni (2009).
}
\usage{
kfweLR(p, k = 1, alpha = 0.01, disp = TRUE)
kfweGR(p, k = 1, alpha = 0.01, disp = TRUE, SD=TRUE, const = 10, alpha.prime = getAlpha(k = k, s = length(p), alpha = alpha, const = const))
kfweOrd(p, k = 1, alpha = 0.01, ord = NULL, alpha.prime = alpha, J = qnbinom(alpha,k,alpha.prime), disp = TRUE, GD=FALSE)
                                                  
getAlpha (s, k = 1, alpha = 0.01, const = 10)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of p-values of length s
%%     ~~Describe \code{p} here~~
}
  \item{s}{ number of p-values (i.e. hypotheses)
%%     ~~Describe \code{p} here~~
}
  \item{k}{ number of allowed errors in kFWE controls
%%     ~~Describe \code{k} here~~
}
  \item{alpha}{ global significance level
%%     ~~Describe \code{alpha} here~~
}
  \item{ord}{ the vector of values based on which the p-values have to be ordered
%%     ~~Describe \code{ord} here~~
}
  \item{const}{ Bigger is better (more precise but slower)
%%     ~~Describe \code{const} here~~
}
  \item{J}{ number of allowed jumps befor stopping
%%     ~~Describe \code{J} here~~
}
  \item{disp}{ diplay output? TRUE/FALSE
%%     ~~Describe \code{disp} here~~
}
  \item{SD}{ Step-down version of the procedure? (TRUE/FALSE) the step-down version is uniformly more powerful than the single step one.
%%     ~~Describe \code{disp} here~~
}
  \item{alpha.prime}{ univariate alpha for single step Guo and Romano procedure
%%     ~~Describe \code{alpha} here~~
}
  \item{GD}{ Logic value. Should the correction for general dependence be applied? (See reference below for further details)
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
kfweOrd, kfweLR, kfweGR, kfweGR.SD return a vector of kFWE-adjusted p-values. It respect the order of input vector of p-values p.

getAlpha returns the alpha for Guo and Romano procedure.

}
\references{
For Lehmann and Romano procedure see: 

Lehmann and Romano (2005) Generalizations of the Familywise Error Rate, Annals of Statistics, 33, 1138-1154.

For Guo and Romano procedure see:

Guo and Romano (2007) A Generalized Sidak-Holm procedure and control of genralized error rates under independence, Statistical Applications in Genetics and Molecular Biology, 
6, 3.


For Ordinal procedure see:

Finos and Farcomeni (2010) k-FWER control without multiplicity correction, with application to detection
of genetic determinants of multiple sclerosis in Italian twins. Biometrics  (Articles online in advance of print: DOI 10.1111/j.1541-0420.2010.01443.x)
%% ~put references to the literature/web site here ~
}
\author{ L. Finos and A. Farcomeni
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(13)
y=matrix(rnorm(3000),3,1000)+2                      #create toy data
p=apply(y,2,function(y) t.test(y)$p.value)          #compute p-values
M2=apply(y^2,2,mean)                                #compute ordering criterion

kord=kfweOrd(p,k=5,ord=M2)                          #ordinal procedure
kgr=kfweGR(p,k=5)                                   #Guo and Romano

kord=kfweOrd(p,k=5,ord=M2,GD=TRUE)                  #ordinal procedure (any dependence)
klr=kfweLR(p,k=5)                                   #Lehaman and Romano (any dependence)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
