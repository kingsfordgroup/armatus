For members of the Kingsford Group committing to armatus, this is how releases and binaries are made (instructions and procedures thanks to Rob!):

The auto build happens on every commit, and each time it gets pushed to the same location (https://github.com/kingsfordgroup/armatus/releases/tag/TravisCIBinary).  That means that each time you commit to the armatus repository, the commit gets uploaded to Travis CI (unless the commit message contains “[ci skip]” or “[skip ci]” somewhere).  If the build is successful (https://travis-ci.org/kingsfordgroup/armatus), then the resulting binary, built on TravisCI’s VM (currently ubuntu 12.04), gets uploaded to the TravisCIBinary git tag mentioned above.  This is done via the GitHub API.   You can use this binary for a release you want to make (i.e. download it, rename it, and upload it for a release).  Before you create a release make sure to update the version in the the Armatus help output.

In order to make a new release, you should do the following.

1. Create an annotated tag with git (on your local machine)

   git tag -a v2.1 -m "Armatus version X.X"

2. This will tag the *most recent commit*.  You can then push these tags upstream with

   git push origin —tags 

3. Go to GitHub and, under the "Releases" tab (https://github.com/kingsfordgroup/armatus/releases) click "Draft a new release"

4. From the box at the top left that says "Tag version", select the tag you just pushed, and fill out the rest of the release info as desired

5. GitHub will automatically make a source tarball for this release, and there is a box at the bottom of this page where you can drag and drop
     any binaries you want to associate with the release (if you forget to do that now, you can always edit the release later from the same page).

6. Click Publish release.

7. Profit.

The files that are associating with creating the binary using the Travis CI system are: .travis.yml and scripts/push_binary.sh
