
source("../../R/functions/metadata_functions.R")



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

metadata <- load_metadata()


test_that("Metaphlan3 loads correctly",{
  expect_equivalent(dim(load_metaphlan("species")), c(505,444))
  expect_equivalent(dim(load_metaphlan()), c(186,444))
  expect_equivalent(dim(load_metaphlan("species")), c(505,444))
})

metaphlan <- load_metaphlan("species")




### metadata functions

test_that("selectMetadata selects 13 columns",{
  expect_equivalent(dim(selectMetadata(metadata)), c(1277,13))
})




test_that("test if isolateDonorBatchMetadata isolates donor metadata",{
  donorbatches <- isolateDonorBatchMetadata(metadata)

  expect_equal(all(donorbatches$project == "donor_batch"), TRUE)
})




## Metaphlan functions

test_that("calculateSpeciesRichness has the correct output",{
  expect_equivalent(dim(calculateSpeciesRichness(metaphlan)), c(442,2))
  expect_equivalent(colnames(calculateSpeciesRichness(metaphlan)), c("LibID", "richness"))
})


t_metaphlan <- transposeMetaphlan(metaphlan)

test_that("transpose Metaphlan transposes metaphlan",{
 expect_equivalent(sort(colnames(metaphlan)[-c(1,2)]), rownames(t_metaphlan))
})