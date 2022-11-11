print("Hello Hanan!")

### Uncomment the following lines and run them one by one
### in the RStudio Terminal:

# which git
## -> should output a folder/directory on 
## your computer like "/usr/local/bin/git"
## IF NOT, you need to download GIT before proceeding

# git --version
## -> should output something like "git version 2.37.2"

# git config  user.name "Your Name"
# git config  user.email "your@github_email.com"

usethis::create_from_github(
"https://github.com/sheffardj/bayes-project.git",
destdir = "~/Documents/UvicTerm1/S564/git-repo"
)


## more here:
## https://jennybc.github.io/2014-05-12-ubc/ubc-r/session03_git.html