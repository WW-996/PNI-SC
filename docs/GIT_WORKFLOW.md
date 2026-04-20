# Git Workflow (PNI SC)

## 1) Repository setup baseline

- Default branch: `main`
- Keep `main` always runnable.
- Use short-lived feature branches for focused work.

## 2) Commit rhythm (record completely)

- Commit every meaningful change set (not by day, by intent).
- One commit should answer one question: "what changed and why?"
- Prefer small commits over large mixed commits.

Recommended commit types:

- `feat:` new analysis capability or output
- `fix:` bug fix or logic correction
- `refactor:` structure cleanup with no behavior change
- `docs:` log or documentation update
- `chore:` tooling / environment / maintenance
- `test:` add or update tests

Examples:

- `feat(integration): add tumor-only PNI UMAP export`
- `fix(dotplot): align marker keys to canonical broad labels`
- `docs(log): add 2026-04-21 progress checkpoint`

## 3) Branch habit

- Branch naming: `feat/<topic>`, `fix/<topic>`, `docs/<topic>`
- Merge to `main` only after:
  - tests pass
  - notebook outputs are reproducible for touched steps
  - `docs/logs/` updated if analysis behavior changed

## 4) Notebook + log discipline

- For each notebook behavior change:
  - update relevant test(s) under `tests/`
  - add/append one dev log file in `docs/logs/`
- Log template (minimal):
  - Date
  - Scope
  - Files changed
  - Verification commands + result
  - Next step

## 5) Suggested daily command checklist

```powershell
git status
git add -p
git commit -m "type(scope): message"
git push
```

## 6) Tags for milestones

- Add tags for stable analysis milestones:
  - `v0.1-qc`
  - `v0.2-integration`
  - `v0.3-subclustering`

