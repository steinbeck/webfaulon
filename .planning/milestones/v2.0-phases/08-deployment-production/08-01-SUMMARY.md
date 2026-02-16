---
phase: 08-deployment-production
plan: 01
subsystem: backend, frontend, deployment
tags: [deployment, configuration, environment, render]
dependency-graph:
  requires: [phase-07]
  provides: [production-config, requirements-txt, env-variables]
  affects: [backend-deps, api-client, sse-connection]
tech-stack:
  added: [rdkit-dependency]
  patterns: [vite-env-variables, build-time-injection]
key-files:
  created:
    - backend/requirements.txt
    - .env.production
  modified:
    - backend/pyproject.toml
    - src/api/client.ts
    - src/api/sse.ts
    - .gitignore
decisions:
  - Manual requirements.txt creation (not poetry export) for clean >= constraints
  - VITE_API_URL environment variable for build-time backend URL injection
  - .env.production committed to repo (Vite convention), .env excluded via gitignore
  - Conditional URL construction with fallback to relative paths for dev workflow
metrics:
  duration: 102s
  completed: 2026-02-16T17:27:26Z
  tasks: 2
  files: 6
---

# Phase 08 Plan 01: Production Configuration Summary

Production deployment configuration files created with environment-driven backend URL injection via Vite build system.

## Tasks Completed

### Task 1: Create backend requirements.txt and add rdkit dependency

**Status:** Complete
**Commit:** f9892f2
**Duration:** ~50s

Created production-ready requirements.txt with all 6 dependencies and added rdkit to pyproject.toml.

**Changes:**
- Added `rdkit = "^2025.9.4"` to `backend/pyproject.toml` [tool.poetry.dependencies]
- Created `backend/requirements.txt` with 6 production packages (fastapi, uvicorn, pydantic, pydantic-settings, sse-starlette, rdkit)
- Updated `.gitignore` to exclude `.env`, `.env.local`, and `backend/.env` files from version control

**Verification:**
```bash
$ cat backend/requirements.txt
fastapi>=0.129.0
uvicorn[standard]>=0.40.0
pydantic>=2.12.5
pydantic-settings>=2.13.0
sse-starlette>=2.0
rdkit>=2025.9.4

$ grep rdkit backend/pyproject.toml
rdkit = "^2025.9.4"

$ grep -c "env" .gitignore
3
```

**Key decision:** Used manual requirements.txt creation with `>=` constraints (not `poetry export`) to avoid dev dependencies and overly pinned versions. This ensures clean, production-only dependencies with minimum version requirements suitable for Render deployment.

### Task 2: Add VITE_API_URL support to frontend and GitHub Actions workflow

**Status:** Complete
**Commit:** da75bdc
**Duration:** ~52s

Updated frontend API client and SSE connection to use environment-driven backend URL with fallback to development proxy.

**Changes:**
- Created `.env.production` with placeholder Render URL: `VITE_API_URL=https://webfaulon-backend.onrender.com`
- Updated `src/api/client.ts` baseURL to use conditional:
  ```typescript
  private baseURL = import.meta.env.VITE_API_URL
    ? `${import.meta.env.VITE_API_URL}/api/sa`
    : '/api/sa';
  ```
- Updated `src/api/sse.ts` EventSource URL construction:
  ```typescript
  const apiBase = import.meta.env.VITE_API_URL || '';
  this.eventSource = new EventSource(`${apiBase}/api/sa/${sessionId}/stream`);
  ```

**Verification:**
```bash
$ npm run build
✓ built in 437ms

$ grep -r "onrender.com" dist/
# Found production URL embedded in JavaScript bundle

$ grep "import.meta.env.VITE_API_URL" src/api/client.ts src/api/sse.ts
# Both files use environment variable with fallback
```

**Key decision:** Used conditional URL construction with fallback to relative paths. This ensures:
- Production: Full Render URL injected at build time from `.env.production`
- Development: Relative `/api` path uses Vite proxy to localhost:8000
- Zero changes to dev workflow — proxy continues working when VITE_API_URL is unset

**GitHub Actions workflow:** No changes required. Vite automatically reads `.env.production` during `npm run build` when NODE_ENV=production (default for `vite build`). The workflow already runs `npm run build`, which will inject the production URL.

## Overall Verification

All success criteria met:

1. **backend/requirements.txt** — Valid pip requirements file with all 6 production dependencies ✓
2. **.env.production** — Contains VITE_API_URL pointing to Render backend ✓
3. **client.ts and sse.ts** — Both use import.meta.env.VITE_API_URL with fallback to relative URL ✓
4. **Production build** — `npm run build` injects backend URL into JavaScript bundle ✓
5. **Development workflow** — Vite proxy still works unchanged (verified by conditional logic) ✓

```bash
# Verification summary
$ cat backend/requirements.txt  # 6 dependencies
$ grep "import.meta.env.VITE_API_URL" src/api/*.ts  # 3 occurrences
$ npm run build && grep -r "onrender.com" dist/  # URL embedded in bundle
```

## Deviations from Plan

None — plan executed exactly as written.

## Architecture Notes

**Vite Environment Variables:**
- Vite injects `import.meta.env.*` variables at build time (NOT runtime)
- `.env.production` is read automatically by `vite build`
- Variables prefixed with `VITE_` are exposed to client code
- Conditional expressions allow graceful fallback for development

**Deployment Flow:**
1. GitHub Actions runs `npm run build`
2. Vite reads `.env.production` from repo
3. `VITE_API_URL` is injected into JavaScript bundle during build
4. Built `dist/` folder contains hardcoded production URL
5. Backend deployed separately to Render with requirements.txt

**Dev/Prod Separation:**
- Dev: `npm run dev` → no VITE_API_URL → relative `/api` → Vite proxy → localhost:8000
- Prod: `npm run build` → reads `.env.production` → full URL → direct Render connection

## Next Steps

Plan 08-02 will:
1. Create Render account and deploy backend service using requirements.txt
2. Update `.env.production` with actual Render URL
3. Deploy frontend to GitHub Pages
4. Verify end-to-end production functionality

## Self-Check: PASSED

**Files created:**
```bash
$ [ -f "backend/requirements.txt" ] && echo "FOUND: backend/requirements.txt"
FOUND: backend/requirements.txt

$ [ -f ".env.production" ] && echo "FOUND: .env.production"
FOUND: .env.production
```

**Files modified:**
```bash
$ git log --oneline -2
da75bdc feat(08-01): add VITE_API_URL support for production backend configuration
f9892f2 chore(08-01): add rdkit dependency and create requirements.txt for production deployment
```

**Commits exist:**
```bash
$ git log --oneline --all | grep -q "f9892f2" && echo "FOUND: f9892f2"
FOUND: f9892f2

$ git log --oneline --all | grep -q "da75bdc" && echo "FOUND: da75bdc"
FOUND: da75bdc
```

All files created, all commits recorded. Self-check passed.
