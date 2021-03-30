context("Test make vertical profile")
prof <- verticalProfileLUF20(acousticsQA::echosounderSkrei2019, 31)

expect_true(!any(is.na(prof$log)))

plotStretch(prof)

trawl <- extractTrawlsBiotic(acousticsQA::bioticSkrei2019)
expect_true(!any(is.na(trawl$log)))
plotMap(prof, trawl, header="torsk")
