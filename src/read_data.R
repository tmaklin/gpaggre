## BSD 3-Clause License
## 
## Copyright (c) 2023, Tommi Mäklin (tommi `at' maklin.fi)
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
## 
## 1. Redistributions of source code must retain the above copyright notice, this
##    list of conditions and the following disclaimer.
## 
## 2. Redistributions in binary form must reproduce the above copyright notice,
##    this list of conditions and the following disclaimer in the documentation
##    and/or other materials provided with the distribution.
## 
## 3. Neither the name of the copyright holder nor the names of its
##    contributors may be used to endorse or promote products derived from
##    this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
## FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Functions for processing the input data

ReadGallups <- function(path, start.date, time.unit="days") {
    gallups <- read.table(path, sep='\t', header=TRUE)
    gallups$Sample_size <- as.numeric(gsub(" ", "", gallups$Sample_size)) ## TODO fix format
    gallups$LIIK <- as.numeric(gsub("— ", NA, gallups$LIIK)) ## TODO set as empty

    ## Convere Finnish format dates into time from first entry
    gallups$Date <- as.POSIXlt(gallups$Date, format = "%d.%m.%Y")
    gallups$days_since_election <- difftime(time1=gallups$Date, time2=as.POSIXlt(start.date, format="%d.%m.%Y"), tz="GMT", units=time.unit)

    ## Remove rows that have missing values
    ## This could be modelled in the GP kernel, see
    ## the "Kernel specification in practice" subheader in:
    ## https://www.nature.com/articles/s41467-019-09785-8
    gallups <- gallups[!apply(gallups, 1, function(x) { any(is.na(x)) }), ]

    ## Kantar TNS changed name to Kantar Public at some point
    gallups$Pollster <- gsub("Public", "TNS", gallups$Pollster)

    ## Reorder so that the oldest poll comes first
    gallups <- gallups[order(gallups$days_since_election), ]

    gallups
}
