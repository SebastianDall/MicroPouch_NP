env <- new.env()
source("../../R/functions/metadata_functions.R", local = env)
attach(env, name = "sourced_scripts")



### Loading data
test_that("Test correct metadata is loaded",{
  #Dim
  expect_equivalent(dim(load_metadata()), c(1277,20))
  #colnames
  expect_setequal(
    colnames(load_metadata()),
    c("id","sample_barcode","LibID","library_plate","ext_conc","lib_conc","pool_ng","stage","project","group","donor","batch_1","batch_2","batch_3","fecal_donation_number","fecal_batch_date","pdai_score","paired_reads","mean_quality","non_human_reads")
  )

})




test_that("Metaphlan3 loads correctly",{
  expect_equivalent(dim(load_metaphlan("species")), c(505,444))
  expect_equivalent(dim(load_metaphlan()), c(186,444))
  expect_equivalent(dim(load_metaphlan("species")), c(505,444))
})


### metadata functions

test_that("selectMetadata selects 13 columns",{
  expect_equivalent(dim(selectMetadata(metadata)), c(1277,13))
})