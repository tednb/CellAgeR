# Load the objects from the original Clocks.Rda
env <- new.env()
load("/mnt/local-disk/data/guoxiaolong/Renv/EhancerClock/Clocks/Clocks.Rda", envir = env)
clock_coefficients <- as.list(env)

# Save as internal data for the package
usethis::use_data(clock_coefficients, internal = TRUE, overwrite = TRUE)
