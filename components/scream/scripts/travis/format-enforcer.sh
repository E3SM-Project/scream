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
    git push -u origin ${TRAVIS_BRANCH}

    git checkout ${git_hash}
fi
