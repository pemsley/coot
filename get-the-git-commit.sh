
if git status > /dev/null ; then
    git rev-parse HEAD
else
    echo missing-git-commit
fi
