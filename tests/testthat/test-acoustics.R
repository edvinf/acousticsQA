context("Test make vertical profile")
prof <- verticalProfileLUF20(acousticsQA::echosounderSkrei2019, 31)
profAll <- verticalProfileLUF20(acousticsQA::echosounderSkrei2019, unique(acousticsQA::echosounderSkrei2019$acocat$acocat))

expect_true(!any(is.na(prof$log)))



context("Test extract trawls")
trawl <- extractTrawlsBiotic(acousticsQA::bioticSkrei2019)
expect_true(!any(is.na(trawl$log)))

context("Test plots")
plotStretch(prof)
plotMap(prof, trawl, header="torsk")
plotStretch(prof, trawls=trawl)
plotPurity(prof, profAll)




