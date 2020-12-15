suppressPackageStartupMessages(library(renv))

# initialize a new project-local environment with a private R library
renv::init()

# save the state of your project to renv.lock
renv::snapshot()

# list previous revisions of your lockfile
renv::history()

# restore the state of your project from renv.lock
renv::restore()

