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

## Functions for plotting the model

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
    ## Create text with a `bg` (default: black) border
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo,
             labels, col=bg, ... )
        
    }
    text(xy$x, xy$y, labels, col=col, ... )
}

sigfig <- function(vec, digits){
    ## Format `vec` numbers nicely with significant `digits`
    return(gsub("\\.$", "", formatC(signif(vec,digits=digits), digits=digits, format="fg", flag="#")))
}

PlotFit <- function(samples, gallups) {
    ## Plot `samples` from the GP fit on `gallups`.

    ## Extract the GP predictions
    support.samples <- extract(samples, "f")

    ## Calculate 80% credible interval and mean
    support.average <- apply(support.samples[[1]], c(2, 3), mean)

    ## Extract the bias
    bias <- mean(unlist(extract(samples, "bias")))

    ## Transform the unbounded values to range (0, 100)
    support.mean <- t(apply(support.average, 1, function(x) exp(bias + x)/(1 + sum(exp(bias + x)))))
    ## Party colors
    cols <- c("#f54b4b", "lightblue", "#00517d", "#3aad2e", "#284735", "#f00A64", "#e1ad01", "#7851a9", "#b42079")

    ## legend in right panel, support values in left
    layout(matrix(1:2, 1, 2), widths=c(0.875, 0.125))

    par(mar=c(3, 4.5, 1, 4))
    plot(',', xlim=c(0, nrow(support.mean)), ylim=c(0, 30), ylab="Kannatus (%)", xlab="", bty='n', xaxt='n', cex.lab=1.2)
    other.parties.support <- 1 - rowSums(support.mean)
    lines(x=1:nrow(support.mean), y=other.parties.support*100, type='l', col="gray", lwd=2)
    lines(x=1:nrow(gallups), y=gallups[, 3 + ncol(gallups) - 3 - 1], type='p', col="gray")
    axis(side=1, at=c(seq(1, nrow(gallups), by=15), 90), labels=gsub("[-][0-9][0-9]$", "", as.character(c(gallups$Date, as.POSIXlt("2.4.2023", format="%d.%m.%Y"))))[c(seq(1, nrow(gallups), by=15), 90)])
    ## Plot the estimated party suppports
    for (i in 1:(ncol(gallups) - 3 - 1 - 1)) {
        lines(x=1:nrow(support.mean), y=support.mean[, i]*100, type='l', col=cols[i], lwd=2)
        lines(x=1:nrow(gallups), y=gallups[, 3 + i], type='p', col=cols[i])
    }
    ## Plot the constrained variable (other parties support in polling)

    abline(v=nrow(support.mean) - 1, col="gray50", lty="dashed")
    ## Add text which has the support on the requested prediction date
    y.vals <- c(18, 20, 22, 12, 10, 8, 6, 4, 0, 2)
    y.labels <- c(support.mean[nrow(support.mean) - 1, ], other.parties.support[nrow(support.mean) - 1])*100
    x.pos <- rep(nrow(support.mean) + 8, length(y.vals))
    x.pos[5:10] <- x.pos[5:10] + 1
    shadowtext(x=x.pos, y=y.vals, labels=paste(c(sigfig(y.labels[1:4], 3), sigfig(y.labels[5:10], 2)), "%"), xpd=TRUE, col=c(cols, "gray"), cex=1.5)

    ## Legend
    par(mar=c(4,0.25,0,0))
    plot(0, 0, xlab='', ylab='', xaxt='n', yaxt='n', pty=" ", bty='n', col="white")
    cols <- c("#00517d", "lightblue", "#f54b4b", "#3aad2e", "#284735", "#f00A64", "#e1ad01", "#7851a9", "#b42079")
    legend("left", legend=c("KOK", "PS", "SDP", "KESK", "VIHR", "VAS", "RKP", "KD", "LIIK", "Muut"), fill=c(cols, "gray"), bty='n')
##}
}
