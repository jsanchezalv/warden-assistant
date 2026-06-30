# WARDEN Assistant Skill

An AI coding assistant skill for building, validating, and optimizing health economic discrete event simulation (DES) cost-effectiveness analysis (CEA) models using the [WARDEN R package](https://github.com/jsanchezalv/WARDEN) (v2.0+).

## Features

- **Build** — Create complete WARDEN models from clinical pathway descriptions and parameter inputs (Excel, CSV, or plain text)
- **Validate** — Two-pass static + runtime validation covering API misuse, DES logic errors, health economics issues, PSA/DSA problems, and HTA submission risks
- **Optimize** — Profile-first efficiency review with 9 optimization patterns and before/after benchmarking

## Installation

### Automatic (recommended)

Works with Claude Code, OpenAI Codex, Cursor, Windsurf, Cline, GitHub Copilot, and more:

```bash
npx skills@latest add jsanchezalv/WARDEN-assistant_skill
```

The installer auto-detects which agents you have and configures each one.

### Manual installation by agent

If you prefer manual setup, clone this repository and follow the instructions for your agent:

| Agent | Where to place the skill |
|-------|--------------------------|
| Claude Code | Copy or symlink `skills/warden-assistant/` into your project's `.claude/skills/` or global `~/.claude/skills/` |
| OpenAI Codex | Copy the content of `skills/warden-assistant/SKILL.md` into your project's `AGENTS.md` |
| Cursor | Copy the content into `.cursor/rules/warden-assistant.md` |
| Windsurf | Append the content to `.windsurfrules` |
| GitHub Copilot | Copy the content into `.github/copilot-instructions.md` |
| Other agents | Paste the content of `SKILL.md` into the agent's system prompt or instruction file |

For agents that use a single instruction file (Codex, Windsurf, Copilot), copy both `SKILL.md` and the relevant `references/` files your workflow needs.

## Multi-project reuse

To use this skill across multiple projects and receive updates automatically:

```bash
# 1. Clone once to a permanent location
git clone https://github.com/jsanchezalv/WARDEN-assistant_skill.git ~/skills/warden-assistant
```

Then symlink into each project:

```bash
# Linux / macOS
ln -s ~/skills/warden-assistant/skills/warden-assistant .claude/skills/warden-assistant

# Windows (PowerShell as admin)
New-Item -ItemType SymbolicLink -Path .claude\skills\warden-assistant -Target $HOME\skills\warden-assistant\skills\warden-assistant
```

### Updating

Pull once, all symlinked projects update immediately:

```bash
cd ~/skills/warden-assistant && git pull
```

## What's included

```
skills/warden-assistant/
├── SKILL.md                    # Main skill (model creation, validation, optimization workflows)
└── references/
    ├── api_reference.md        # Complete WARDEN v2.0+ API documentation
    ├── cea_patterns.md         # 11 standard CEA model archetypes with worked examples
    ├── efficiency_patterns.md  # 9 optimization patterns with expected speedups
    ├── input_block_guide.md    # Parameter management for PSA, DSA, and scenarios
    └── validation_rules.md     # 7-category validation framework (40+ checks)
```

## Requirements

- R 4.2+
- [WARDEN package](https://github.com/jsanchezalv/WARDEN) v2.0+

## License

GPL-3.0
