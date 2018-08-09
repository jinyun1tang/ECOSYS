# ECOSYS
The ECOSYS model for terrestrial ecosystem biogeochemistry

The code is developed by Professor Robert Grant at University of Alberta.

Jinyun Tang (jinyuntang@lbl.gov) helps maintaining the code on github.


To download the code
git clone https://github.com/jinyun1tang/ECOSYS.git

You can also make a fork using the fork button on the upper right corner.

To build the code, type

./build_ecosys make

Then follow the output information to locate the executable.

To clean the build directory, type

./build_ecosys clean

The following are several useful git commands

git log --graph

which will give you a history of the commits.

git status

will list files that are currently tracked or not tracked by git.

git checkout -b your_branch

will make you a branch with name "your_branch"

git push -u origin "your_branch"

will push your branch to remote for the first time. 

You can also create a new branch using the github web tools.

git add files_you_changed

will add the "files_you_changed" to the tracking list.

git commit

will ask you to write a message about what you've done and what you want to others know about (why) the changes are made.

git push origin "your_branch"

will push your changes to your branch.

git pull --ff-only

will do fast forward merge in case your local branch is not synchronized with your remote branch. 

git stash

will foget all your changes made in local directory. 



