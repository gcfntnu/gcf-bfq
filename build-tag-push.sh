#!/bin/bash

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: ./build-tag-push.sh [prod|test|base] [docker-tag] [optional args]

Test image options:
  -t BRANCH       gcf-tools branch (default: bfq-dev)
  -w BRANCH       gcf-workflows branch (default: bfq-dev)
  -b IMAGE        base image (default: gcfntnu/bfq:base-260925)
  -q BRANCH       legacy MultiQC argument (accepted; currently unused)

Set BFQ_DOCKER to a Docker executable or wrapper path to replace "sudo docker".
For an account already allowed to use Docker: BFQ_DOCKER=docker ./build-tag-push.sh ...
The script builds and immediately pushes the requested tag.
USAGE
}

if (( $# < 2 )); then
    usage >&2
    exit 1
fi

mode=$1
tag=$2
shift 2
case "$mode" in
    prod|test|base) ;;
    *) usage >&2; exit 1 ;;
esac
if [[ ! "$tag" =~ ^[a-zA-Z0-9_][a-zA-Z0-9_.-]{0,127}$ ]]; then
    echo "Invalid Docker tag: $tag" >&2
    exit 1
fi

gcf_tools=bfq-dev
gcf_workflows=bfq-dev
base_image=gcfntnu/bfq:base-260925
test_options=false
while getopts ':t:q:w:b:' flag; do
    test_options=true
    case "$flag" in
        t) gcf_tools=$OPTARG ;;
        q) echo 'The -q MultiQC branch option is retained for compatibility but is unused.' >&2 ;;
        w) gcf_workflows=$OPTARG ;;
        b) base_image=$OPTARG ;;
        *) usage >&2; exit 1 ;;
    esac
done
shift "$((OPTIND - 1))"
if (( $# )) || { [[ "$mode" != test ]] && "$test_options"; }; then
    echo 'Branch and base-image options apply only to test builds.' >&2
    usage >&2
    exit 1
fi

# The companion branches are public. Resolving them over Git avoids GitHub's
# unauthenticated REST rate limit and pins each build to one branch snapshot.
branch_revision() {
    local repository=$1 branch=$2 remote_line revision remote_ref
    git check-ref-format "refs/heads/$branch" >&2 || return 1
    remote_line=$(git ls-remote --exit-code "https://github.com/gcfntnu/$repository.git" "refs/heads/$branch") || return 1
    read -r revision remote_ref <<< "$remote_line"
    if [[ ! "$revision" =~ ^[0-9a-f]{40}$ || "$remote_ref" != "refs/heads/$branch" ]]; then
        echo "Cannot resolve $repository branch $branch to a commit." >&2
        return 1
    fi
    printf '%s\n' "$revision"
}

build_args=()
if [[ "$mode" == test ]]; then
    tools_revision=$(branch_revision gcf-tools "$gcf_tools")
    workflows_revision=$(branch_revision gcf-workflows "$gcf_workflows")
    bfq_revision=$(git rev-parse HEAD)
    build_args=(
        --build-arg "BASE_IMAGE=$base_image"
        --build-arg "BFQ_REVISION=$bfq_revision"
        --build-arg "GCF_TOOLS_BRANCH=$gcf_tools"
        --build-arg "GCF_TOOLS_REV=$tools_revision"
        --build-arg "GCF_WORKFLOWS_BRANCH=$gcf_workflows"
        --build-arg "GCF_WORKFLOWS_REV=$workflows_revision"
    )
fi

docker_command=(sudo docker)
if [[ -n "${BFQ_DOCKER:-}" ]]; then
    docker_command=("$BFQ_DOCKER")
fi
image="gcfntnu/bfq:$tag"
dockerfile="dockerfile-$mode"
printf 'Building %s with %s\n' "$image" "$dockerfile"
printf '%q ' "${docker_command[@]}" build -t "$image" . -f "$dockerfile" "${build_args[@]}"
printf '\n'
"${docker_command[@]}" build -t "$image" . -f "$dockerfile" "${build_args[@]}"
"${docker_command[@]}" push "$image"
