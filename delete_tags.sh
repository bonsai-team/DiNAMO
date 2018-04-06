#Delete local tags.
git tag -d $(git tag -l)
#Fetch remote tags.
git fetch
#Delete remote tags.
git push origin --delete $(git tag -l) # Pushing once should be faster than multiple times
#Delete local tags.
git tag -d $(git tag -l)
