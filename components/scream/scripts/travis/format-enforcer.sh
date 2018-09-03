#!/usr/bin/env bash

# Echo the commands and their output
set -x

commit_msg="Automatic code formatting commit"
scream_dir="${TRAVIS_BUILD_DIR}/components/scream"

test "${TRAVIS_COMMIT_MESSAGE}" != "${commit_msg}"
if [ "$?" -eq 0 ]; then
    git_hash=$(git log --pretty=format:'%H' -n 1)
    git checkout ${TRAVIS_BRANCH}
    for f in $(find ${scream_dir} -iname "*.[ch]pp" -not -path ${scream_dir}/extern); do
	clang-format -i ${f}
	git add ${f}
    done

    git commit -m "${commit_msg}"
    git remote add origin-oauth https://${TRAVIS_GITHUB_TOKEN}@github.com/E3SM-Project/scream.git
    echo WHILE THIS PR IS IN PROGRESS, WE ARE NOT PUSHING COMMITS TO MASTER. UNCOMMENT WHEN READY.
    #git push -u origin-oauth ${TRAVIS_BRANCH}

    git checkout ${git_hash}
fi
