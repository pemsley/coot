if git status > /dev/null ; then
    git rev-list --count HEAD 2>/dev/null
else
    echo 0
fi

