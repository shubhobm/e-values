# install CRAN packages
packages = readLines("./requirements.txt")
print(packages)
for(i in 1:length(packages)){
    install.packages(packages[i])
}