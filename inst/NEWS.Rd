\name{NEWS}
\title{palaeoSig News}
\encoding{UTF-8}

\section{Changes in version 1.1-3}{

  \subsection{GENERAL}{
    \itemize{
      \item Bug fixes and minor improvement, merged package with age-depth functions and R implementation of Minchin's COMPAS. 
      \item Changes to the namespace mean that the rioja package needs to be loaded explicitly before some of the function will work properly.
      
    }
  } % end GENERAL
  
  \subsection{NEW FUNCTIONS}{
    \itemize{      
      \item \code{agelme}, \code{predict.agelme} and \code{plot.fittedAgelme} age-depth modelling functions from Heegaard et al (2005)
      \item \code{coverage.plot} Diagnostic plot for reconstructions
      \item \code{centipede.plot} Plot WA optima and tolerances
      \item \code{make.env}, \code{species}, \code{make.env} Simulate species-environment relationships based on Minchin's (1983) compas.
    }
  } % end NEW FUNCTIONS

  \subsection{REMOVED FUNCTIONS}{
    \itemize{      
      \item \code{rotate} No longer needed as gstat now takes geodesic distance.
      \item \code{simulate.spatial} Was a wrapper for krige(), now easier to use krige() directly.
    }
  } % end NEW FUNCTIONS

  \subsection{MODIFIED FUNCTIONS}{
    \itemize{      
      \item \code{randomTF} and \code{randomTFmm} bug when partialling reconstructions out fixed
      
    }
  } % end NEW FUNCTIONS

  \subsection{BUG FIXES}{
    \itemize{      
      \item \code{obs.cor} bug when species not in same order in calibration set and fossil data fixed.
    }
  } % end NEW FUNCTIONS

}% end VERSTION 1.1-3
  

\section{Changes in version 1.1-2}{

  \subsection{GENERAL}{
    \itemize{
      
      \item Bug fixes and minor improvement.
      
    }
  } % end GENERAL
  
  \subsection{NEW FUNCTIONS}{
    \itemize{      
      \item \code{ModelMaker} and \code{randomTFmm} allows models for randomTF to be fitted once and used for several reconstructions. This can be much faster, but does not work with MAT.
      
    }
  } % end NEW FUNCTIONS

  \subsection{MODIFIED FUNCTIONS}{
    \itemize{      
      \item \code{plot.palaeoSig} p-value highlighted can be selected.
      \item \code{plot.obscor} p-value highlighted can be selected.
      \item \code{randomTF} Modified to accept predictions in a vector. 
      
      
      
    }
  } % end NEW FUNCTIONS

  \subsection{BUG FIXES}{
    \itemize{      
      \item \code{obs.cor} bug when species not in same order in calibration set and fossil data fixed.
    }
  } % end NEW FUNCTIONS

}% end VERSTION 1.1-2
  


\section{Changes in version 1.1-1}{

  \subsection{GENERAL}{
    \itemize{
      
      \item The package \pkg{autocorTF} has been merged with \pkg{palaeoSig} for ease of maintenance.
      
    }
  } % end GENERAL
  
  \subsection{NEW FUNCTIONS}{
    \itemize{      
      \item \code{identify.obscor} allows species names to be added to plots interactively.
      
      
    }
  } % end NEW FUNCTIONS

  \subsection{MODIFIED FUNCTIONS}{
    \itemize{      
      \item \code{obs.cor} Now includes correlations with several different species weights
      \item \code{plot.obscor} upgraded to allow the a choice of which abundance weighting is used, and the code for scaling points has been improved.
      
      
    }
  } % end NEW FUNCTIONS

  \subsection{BUG FIXES}{
    \itemize{      
      \item \code{randomTF} bug when partialling out reconstruction other than MAT now fixed.
    }
  } % end NEW FUNCTIONS
  
  

}% end VERSTION 1.1-1




\section{Changes in version 1.1}{

  \subsection{GENERAL}{
    \itemize{
      
      \item The package \pkg{autocorTF} has been merged with \pkg{palaeoSig} for ease of maintenance.
      
    }
  } % end GENERAL
  
  \subsection{NEW and REWRITTEN FUNCTIONS}{
    \itemize{
      
      \item \code{RNE} replaces function \code{mat.rne}, \code{mat.he} and \code{mat.rd}. \code{RNE} allows any of the transfer function methods in \pkg{rioja} to be used to find the dependence of transfer function performance on spatially close observations.
      \item \code{plot.RNE} replaces \code{plot.rne} to work with the output of \code{RNE}.
      \item \code{obs.cor} replaces the previous version of \code{obs.cor}, \code{sim.cor} and \code{obscor.sig} to make the function easier to use and more similar to \code{randomTF}
      \item \code{plot.obscor} replaces the previous version of \code{plot.obscor} and \code{plot.simcor} to use the new output of \code{obs.cor}. The two plot types can be selected with the argument \code{which}.
      \item \code{jointsig}: tests if two environmental variables have joint control on fossil assemblage composition. 
      \item \code{plot.js}: a \code{plot} function for \code{jointsig}
      \item \code{Hill.N2.core} now calculates the Minimum, first quartile and median effective number of species for all fossil observations. This makes much more sense than the previous version, which conflated diversity in individual levels and turnover between levels.
      \item \code{plot.palaeoSig} has been improved to give neater figures.
      
    }
  } % end NEW FUNCTIONS
  
  \subsection{REMOVED FUNCTIONS}{
    
    \itemize{
      
      \item \code{wajack} was removed as \code{WA} is now implemented in \pkg{rioja}. 
      
    }
  } % end REMOVED FUNCTIONS
  
  \subsection{FORTHCOMMING CHANGES}{
    \itemize{
      
      \item \code{mat.h} will be removed in the next release, when code for h-block resampling should be implimented in \pkg{rioja}.
      
    }
  } % end FORTHCOMMING CHANGES
}% end VERSTION 1.1
   	  