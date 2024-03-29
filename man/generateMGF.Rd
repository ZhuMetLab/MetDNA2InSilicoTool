% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunTools.R
\name{generateMGF}
\alias{generateMGF}
\title{generateMGF}
\usage{
generateMGF(
  file_name = "./test.mgf",
  title = NULL,
  precusormz = NULL,
  mslevel = c("2", "1"),
  charge = NULL,
  rt = NULL,
  spec = NULL
)
}
\arguments{
\item{file_name}{Default: 'test.mgf'}

\item{title}{feature name}

\item{precusormz}{numeric. Default: NULL}

\item{charge}{numeric. Default: NULL}

\item{rt}{numeric. Default: NULL}

\item{spec}{matrix}
}
\description{
generateMGF
}
\examples{
test <- ImmsTools::readMSP(file = 'F:/01 MetIMMS/00 data processing/190920_allccs_metabolite_id_demenstration/topscience_msms_pos_20v_190913.msp',
                           mode = 'all')
GenerateMGF(file_name = 'test.mgf',
title = test[[1]]$info$NAME,
precusormz = test[[1]]$info$PRECURSORMZ,
mslevel = "2",
charge = '1+',
spec = test[[1]]$spec)
}
\author{
Zhiwei Zhou
}
