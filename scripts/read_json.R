#load and clean nextstrain json file output for variant frequency

packages = c("jsonlite", "purrr", "stringr", "rlist")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

json_import = fromJSON(file = "data_input/CT-SARS-CoV-2_connecticut.json")



tmp <- tempfile()
url1 <- "https://raw.githubusercontent.com/rebecca-earnest/ncov/master/auspice/ncov_nextstrain-test_connecticut_tip-frequencies.json"
url2 <- "https://github.com/rebecca-earnest/ncov/blob/master/auspice/ncov_nextstrain-test_connecticut.json"
download.file(url2, destfile =tmp,quiet = FALSE, mode = "wb")
json_import <- jsonlite::fromJSON(tmp)

json = list.subset(json_import, 
                   str_detect(names(json_import), "USA/CT")
                   )

json$`USA/CT-BVAMC-749193/2021`

#ASK KANE WHY THIS ISNT WORKING
json1 = map(names(json_import), ~ keep(.x, .p = str_detect(names(.x), "USA/CT", negate = T)))
