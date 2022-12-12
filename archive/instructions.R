
### Uncomment the following lines and run them one by one
### in the RStudio Terminal:

# which git
## -> should output a folder/directory on 
## your computer like "/usr/local/bin/git"
## IF NOT, you need to download GIT before proceeding

# git --version
## -> should output something like "git version 2.37.2"

# git config --global user.name "Your Name"
# git config --global user.email "your@github_email.com"

### Install "usethis" then run the create_from_github command
# install.packages('usethis')
# usethis::create_from_github(
#   "sheffardj/bayes-project", #LEAVE THIS!!!
#   destdir = "~/Documents/UvicTerm1/S564/", ### CHANGE the destination directory!!!
#   fork=FALSE
# )


### Change the name and email, this strictly sets the local user in
### the file .git/config (a hidden file!!!)
# git config --global user.name "Dayten Sheffar"
# git config --global user.email "sheffardj@gmail.com"
# git config credential.useHttpPath true

## If that all worked as expected (unlikely) then
## you should be able to pull and push changes via: 

## RStudio > "Tools" > "Version Control" > 
#### "Commit": opens a window to select your changes and commit 
#### with a descriptive message

#### "Push": after you commit your changes, you still need to push them
#### for other people to pull them!

#### "Pull": this pulls the latest state of the branch that reflects other
#### users having "Pushed" changes to github. This gets you up to date!


## more here:
## https://jennybc.github.io/2014-05-12-ubc/ubc-r/session03_git.html