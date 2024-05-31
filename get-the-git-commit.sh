
if command -v git > /dev/null ; then
    if git status > /dev/null 2>&1 ; then
        git rev-parse HEAD
    else
        echo git-commit-unavailable
    fi
else
    echo git-commit-unavailable
fi
