######################
### SWIRL TRAINING ###
######################

# The swirl R package makes it fun and easy 
# to learn R programming and data science.

# You find everything you need to know about SWIRL under:
# https://swirlstats.com/students.html 

### Install swirl ### ---------------------------
# Install swirl if not already installed
if(!requireNamespace("swirl", quietly = TRUE)) install.packages("swirl")

# update packages - only if needed
# update.packages(ask = F)

# load swirl
require("swirl")
#####################

### Getting swirl courses ### -------------------------------
# you have various options to install a swirl course.

## Option 1) courses from swirlstats
# The easiest way is to run: 
# > install_course("Course Name Here")
# Note that course names are case sensitive!!!
# i.e. if you like to install the "R-Grundlagen" course you simply run:
install_course("R-Grundlagen")

# You can find a list of available courses you can install like this under:
# https://swirlstats.com/scn/title.html
# Here is an example how you can install various courses from swirlstats in a loop:

# vector list of courses you wish to install 
var <- c("R-Grundlagen","A_(very)_short_introduction_to_R","Daten_einlesen_und_kennenlernen",
         "Einfuehrung_in_Datenaufbereitung_mit_tidyR",
         "Daten_visualisieren_mit_ggplot2","R Programming",
         "Advanced R Programming","Exploratory Data Analysis"
         )
# run loop
for(i in var) {
  swirl::install_course(i)
}

## Option 2) courses from github
# However this doesn't apply to all courses and there are more available on github
# https://github.com/swirldev/swirl_courses#swirl-courses
# If a course is on github but not swirlstat, install_course will not work.
# For example the course Mathematical_Biostatistics_Boot_Camp is only stored
# at github. So how do we install it then?

# Swirl offers you the options to install courses from various sources such as
# local files (ziped or unziped), dropbox, google drive or lgithub
# Ask R for help about the possible options.
?InstallCourses()

# If you wish to install a course from that github repo the easiest way is to
# download the course files locally and then install those you want to.
download.file(
  url = "https://github.com/swirldev/swirl_courses/archive/master.zip",
  destfile = "swirl_courses.zip"
  )
# Now the course files are stored in your current working directory.
# First we must unzip them: 
unzip(zipfile = "swirl_courses.zip")

# Have a look the available courses in the "swirl_course" directory
list.files("swirl_courses-master/")

# Now we can simply install i.e. the Mathematical_Biostatistics_Boot_Camp
# course running: 
install_course_directory("swirl_courses-master/Mathematical_Biostatistics_Boot_Camp/")
#############################

# That's it! Now you know two ways how to install swirl courses. 
# To run a course just start swirl and select the course you want to train with.
swirl()
# To leave swirl type: 0

# Have fun! :)