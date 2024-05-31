if command -v git > /dev/null ; then
    if git status > /dev/null 2>&1 ; then
        git rev-list --count HEAD 2>/dev/null
    else
        echo 0
    fi
else
    echo 0
fi

