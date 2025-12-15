#!/usr/bin/env bash
# Create the GitHub repository and push this directory to it.
# Run from inside the hsds-proportions folder, once, after checking LICENSE
# and the USERNAME placeholders in pyproject.toml, CITATION.cff and README.md.
set -euo pipefail

REPO="${1:-hsds-proportions}"
VISIBILITY="${2:---public}"

command -v gh >/dev/null || { echo "gh is not installed: brew install gh"; exit 1; }
gh auth status >/dev/null 2>&1 || { echo "run 'gh auth login' first"; exit 1; }

USER=$(gh api user --jq .login)
echo "Creating ${USER}/${REPO} (${VISIBILITY#--})"

# Fill in the repository URL placeholders before the first commit.
for f in pyproject.toml CITATION.cff README.md; do
  [ -f "$f" ] && sed -i '' "s|USERNAME/hsds-proportions|${USER}/${REPO}|g" "$f" 2>/dev/null \
              || sed -i "s|USERNAME/hsds-proportions|${USER}/${REPO}|g" "$f"
done

git init -b main
git add -A
git commit -m "Rework the hsdS allele proportion script into a configurable tool

The original script hard-coded the spacer length, the motif families, the
exclusion contexts and the output layout, and picked up its inputs by globbing
the working directory. Everything is now a command line option and the code is
split into modules that can be tested without pysam or the external tools.

Defaults reproduce the published WW2842 numbers, including the 4.9961% error
rate at Q99 and a probability floor of 255, which is covered by a test."

gh repo create "$REPO" "$VISIBILITY" --source=. --remote=origin --push
echo "Done: https://github.com/${USER}/${REPO}"
