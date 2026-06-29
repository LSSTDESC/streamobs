#!/bin/bash
# Script to build documentation and deploy to gh-pages branch

set -e  # Exit on error

# Cleanup function to ensure we return to the original branch and restore stashed files
cleanup() {
    local exit_code=$?
    echo ""
    echo "=========================================="
    echo "Cleaning up..."
    echo "=========================================="
    
    # Remove temporary directory if it exists
    if [ -n "$TMP_DIR" ] && [ -d "$TMP_DIR" ]; then
        rm -rf "$TMP_DIR"
        echo "✓ Temporary directory cleaned up"
    fi
    
    # Switch back to original branch if we're not already there
    CURRENT=$(git branch --show-current)
    if [ "$CURRENT" != "$MUST_BRANCH" ]; then
        echo "Switching back to $MUST_BRANCH branch..."
        git checkout "$MUST_BRANCH"
    fi
    
    # Restore stashed files if they exist
    if git stash list | grep -q "Temporary stash for gh-pages deployment"; then
        echo "Restoring stashed untracked files..."
        git stash pop
        echo "✓ Untracked files restored"
    fi
    
    if [ $exit_code -ne 0 ]; then
        echo ""
        echo "⚠ Script exited with errors, but cleanup completed successfully"
    fi
    
    exit $exit_code
}

# Set trap to run cleanup on exit (whether successful or not)
trap cleanup EXIT

echo "=========================================="
echo "Building and Deploying Documentation"
echo "=========================================="

# Check we're on the MUST_BRANCH branch
MUST_BRANCH="main"
CURRENT_BRANCH=$(git branch --show-current)
if [ "$CURRENT_BRANCH" != "$MUST_BRANCH" ]; then
    echo "Error: This script must be run from the '$MUST_BRANCH' branch"
    echo "Current branch: $CURRENT_BRANCH"
    exit 1
fi

# Check for uncommitted changes
if ! git diff-index --quiet HEAD --; then
    echo "Error: You have uncommitted changes. Please commit or stash them first."
    exit 1
fi

#echo "✓ On $MUST_BRANCH branch with no uncommitted changes"

# Build the documentation
echo ""
echo "Building documentation..."
make clean
make html

echo "✓ Documentation built successfully"

# Save the HTML to a temporary location
TMP_DIR=$(mktemp -d)
cp -r build/html/* "$TMP_DIR/"
echo "✓ HTML files copied to temporary location: $TMP_DIR"

# Go back to project root
cd ..

# Stash any untracked and ignored files to prevent them from being carried over to gh-pages
echo ""
echo "Stashing untracked and ignored files..."
git stash push --all --keep-index -m "Temporary stash for gh-pages deployment" || true

# Switch to gh-pages branch
echo ""
echo "Switching to gh-pages branch..."
git checkout gh-pages

# Remove old files (except .git and .nojekyll)
echo "Removing old documentation files..."
find . -maxdepth 1 ! -name '.git' ! -name '.nojekyll' ! -name '.' -exec rm -rf {} \; 2>/dev/null || true

# Copy new files
echo "Copying new documentation..."
cp -r "$TMP_DIR"/* .
touch .nojekyll

# Add and commit
echo ""
echo "Committing changes..."
git add -A
if git diff --staged --quiet; then
    echo "No changes to commit"
else
    git commit -m "Update documentation ($(date '+%Y-%m-%d %H:%M:%S'))"
    echo "✓ Changes committed"
    
    # Push to remote (allow this to fail without stopping the script)
    echo ""
    echo "Pushing to GitHub..."
    #if git push origin gh-pages; then
    if git push origin gh-pages; then
        echo "✓ Documentation deployed successfully!"
    else
        echo "⚠ Warning: Push to GitHub failed. You may need to push manually."
        echo "  This can happen if the commit is too large."
        echo "  Try: git push origin gh-pages --no-verify"
    fi
fi

echo ""
echo "=========================================="
echo "✓ Documentation deployment complete!"
echo "=========================================="
echo ""
echo "Your documentation should be available at:"
echo "https://lsstdesc.github.io/streamobs/"
echo ""
echo "Note: It may take a few minutes for GitHub Pages to update."
