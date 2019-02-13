context("Test spillover_ng")

control_path <- system.file("extdata", "compdata", "data",
                            package="flowCore")
matchfile_path <- system.file("extdata", "compdata", "comp_match",
                              package="flowCore")
controls <- lapply(dir(control_path, full.names=TRUE), 
                   read.FCS)
filenames <- sapply(controls, keyword, "$FIL")
names(controls) <- c("UNSTAINED", "FL1-H", "FL2-H", "FL4-H", "FL3-H")
controls <- as(controls, "flowSet")
# Scramble the columns (particularly fsc, ssc) to make sure
# the methods still handle it appropriately
controls <- controls[,c(3,1,4,5,2,6,7)]

test_that("spillover_ng: Using path to dir with files", {
  ref_file <- system.file("extdata", "compdata", "compref4", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  comp <- spillover_ng(path=control_path, fsc="FSC-H", ssc="SSC-H", 
                       patt="-H", matchfile = matchfile_path, pregate = FALSE)
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})

test_that("spillover_ng: Using preconstructed flowSet with filenames as sample names", {
  sampleNames(controls) <- filenames
  ref_file <- system.file("extdata", "compdata", "compref4", package="flowCore")
  comp_ref <- as.matrix(read.table(ref_file, check.names = FALSE))
  matched <- spillover_match(controls, fsc="FSC-H", ssc="SSC-H", matchfile = matchfile_path)
  comp <- spillover_ng(path=control_path, fsc="FSC-H", ssc="SSC-H", 
                       patt="-H", matchfile = matchfile_path, pregate = FALSE)
  expect_equal(colnames(comp), colnames(comp_ref))
  expect_equal(rownames(comp), rownames(comp_ref))
  expect_equivalent(comp, comp_ref, tolerance=3e-08)
})



