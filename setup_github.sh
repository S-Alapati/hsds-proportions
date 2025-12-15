#!/usr/bin/env bash
# One-time bootstrap: create the GitHub repository and push this directory to it.
# After the first push, use git add / commit / push as normal.
set -euo pipefail

REPO="${1:-hsds-proportions}"
VISIBILITY="${2:---public}"

command -v gh >/dev/null || { echo "gh is not installed: brew install gh"; exit 1; }
gh auth status >/dev/null 2>&1 || { echo "run 'gh auth login' first"; exit 1; }

ACCOUNT=$(gh api user --jq .login)
echo "Creating ${ACCOUNT}/${REPO} (${VISIBILITY#--})"

for f in pyproject.toml CITATION.cff README.md; do
  [ -f "$f" ] || continue
  sed -i '' "s|USERNAME/hsds-proportions|${ACCOUNT}/${REPO}|g" "$f" 2>/dev/null \
    || sed -i "s|USERNAME/hsds-proportions|${ACCOUNT}/${REPO}|g" "$f"
done

[ -d .git ] || git init -b main
git add -A
git diff --cached --quiet || git commit -m "hsds-proportions: assign long reads to hsdS alleles from their 6mA motifs

Motifs, spacer lengths, exclusion contexts and allele counts come from a JSON
definition supplied by the user, so any bacterium with a Type I
restriction-modification shufflon can be scored. IUPAC codes are supported,
reverse motifs are derived as reverse complements unless given explicitly, and
each allele may set its own spacer."

gh repo create "$REPO" "$VISIBILITY" --source=. --remote=origin --push
echo "Done: https://github.com/${ACCOUNT}/${REPO}"
