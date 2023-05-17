# General_bioinformatics_tutorials
The general idea of this tutorial is to get you familar with general bioinformatics. Much of the examples are pertinent to the Texas A&M University (TAMU) High Performance Research Compututing (HPRC) resources.

## Basics
### Logging in 
There are a few easy ways to access the TAMU HPRC.    
**Important:You must be physically on campus or logged into TAMU VPN  to access the clusters**

1. Access using SSH. Open terminal on linux or macOS, cmd or powershell on Windows.     
```
# terminal command
ssh NetID@terra.tamu.edu
ssh NetID@grace.tamu.edu
ssh NetID@faster.hprc.edu
```

2. You can use MobaXterm on windows to access the terminal easier just reduces having to type in password and it has a file manager like sidebar that allows for FTP.

3. You can also use the HPRC Open OnDemand portal https://portal.hprc.tamu.edu/

Side note: I like to use a combination of all three depending on what you are doing and ease of procedure.


Let's login to GRACE because I do most of my work on GRACE.

Moving around in linux
```
pwd    # this tells you your current directory
cd     # change directory
ls     # list what is in current directory add -l to see list view
mv     # move
mkdir  # make directory
rm     # remove BE CAREFUL with this one

cd ..  # go one folder up
ls ..  # look at contents one folder structure up

```




In order to load a program we need for analysis or any step we should run the following lines.

```
module spider R 
ml spider R

module spider R/4.2.0
```