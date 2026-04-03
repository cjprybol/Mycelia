#!/usr/bin/env bash
set -euo pipefail

echo "Installing fleet dispatch tools..."

# Claude Code
npm install -g @anthropic-ai/claude-code@latest

# Codex CLI
npm install -g @openai/codex@latest

# Gemini CLI
npm install -g @google/gemini-cli@latest

# Configure Codex for container use (bubblewrap fails in Docker)
CODEX_DIR="/home/vscode/.codex"
mkdir -p "${CODEX_DIR}"
cat > "${CODEX_DIR}/config.toml" << 'TOML'
sandbox_mode = "danger-full-access"
approval_policy = "never"
TOML

# Claude Code onboarding bypass (prevents interactive prompt on first run)
CLAUDE_DIR="/home/vscode/.claude"
mkdir -p "${CLAUDE_DIR}"
cat > "/home/vscode/.claude.json" << 'JSON'
{
  "hasCompletedOnboarding": true
}
JSON

# Ensure vscode user owns all config dirs
chown -R vscode:vscode "${CODEX_DIR}" "${CLAUDE_DIR}" "/home/vscode/.claude.json" 2>/dev/null || true

echo "Fleet dispatch tools installed."
