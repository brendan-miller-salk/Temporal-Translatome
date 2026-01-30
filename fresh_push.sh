#!/bin/bash
#===============================================================================
# Fresh Push: Reset Git History and Push as Initial Commit
#===============================================================================
# This script removes all git history and pushes the current state as a single
# fresh "initial commit". Use this when you want to publish a clean repo.
#
# Usage: ./fresh_push.sh "Your commit message"
#        ./fresh_push.sh  (uses default message)
#
# WARNING: This overwrites remote history! Only use on repos you control.
#===============================================================================

set -e  # Exit on error

# Default commit message
COMMIT_MSG="${1:-Initial commit: Multi-Tissue Atlas of smORFs and Muscle Microproteins}"

# Remote URL (update this for your repo)
REMOTE_URL="https://github.com/brendan-miller-salk/Temporal-Translatome.git"

echo "🗑️  Removing existing git history..."
rm -rf .git

echo "📦 Initializing fresh repository..."
git init

echo "➕ Staging all files..."
git add .

echo "💾 Creating commit: $COMMIT_MSG"
git commit -m "$COMMIT_MSG"

echo "🔀 Setting branch to main..."
git branch -M main

echo "🔗 Adding remote origin..."
git remote add origin "$REMOTE_URL"

echo "🚀 Force pushing to remote..."
git push -f origin main

echo ""
echo "✅ Done! Repository pushed as fresh initial commit."
echo "   View at: ${REMOTE_URL%.git}"
