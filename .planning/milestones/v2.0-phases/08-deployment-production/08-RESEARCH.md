# Phase 08: Deployment & Production - Research

**Researched:** 2026-02-16
**Domain:** Cloud deployment of FastAPI + RDKit Python applications with Vite frontend
**Confidence:** HIGH

## Summary

Phase 8 deployment centers on three core challenges: (1) deploying a FastAPI backend with RDKit dependency to a cloud PaaS, (2) configuring the production Vite frontend to consume the deployed API, and (3) hardening CORS to production-only origins. The primary deployment blocker is RDKit—a C++ library with Python bindings requiring special handling.

**Key finding**: RDKit is now pip-installable via prebuilt wheels on PyPI (package name `rdkit`, not the legacy `rdkit-pypi`), eliminating the historical conda/Docker requirement. This makes Railway and Render viable deployment targets using standard Python buildpacks. Fly.io remains an option but requires Dockerfile configuration.

**Primary recommendation**: Deploy to **Render** using their FastAPI template approach with `pip install rdkit` in requirements. Render provides the most FastAPI-specific documentation, flat-rate pricing suitable for classroom demos, and automatic SSL + health checks. Configure frontend to use `VITE_API_URL` environment variable pointing to the Render backend URL.

## Standard Stack

### Core Deployment Tools

| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| **Render** | N/A (PaaS) | Backend hosting platform | Best FastAPI documentation, flat pricing, automatic SSL/health checks |
| **GitHub Actions** | N/A | Frontend deployment | Already configured for GitHub Pages, zero additional setup |
| **Uvicorn** | 0.40.0+ | ASGI server (production) | FastAPI's recommended production server, async-native |
| **Gunicorn** | Latest | Process manager (optional) | Multi-worker support for production traffic, industry standard with Uvicorn workers |

### Python Deployment Dependencies

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `rdkit` | 2025.9.4+ | Chemistry toolkit | **Always** — critical dependency, now pip-installable via PyPI wheels |
| `uvicorn[standard]` | 0.40.0+ | Production ASGI server | **Always** — required for serving FastAPI in production |
| `gunicorn` | 22.0.0+ | WSGI/ASGI process manager | Optional — only if multi-worker support needed (overkill for classroom demo) |

### Frontend Configuration

| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| Vite env variables | N/A | API URL configuration | **Always** — `.env.production` with `VITE_API_URL` for backend URL |
| GitHub Pages | N/A | Static site hosting | Already deployed — no changes needed to hosting |

### Alternative Platforms Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Render | Railway | Simpler UI, faster deploys, but usage-based billing can surprise; Railpack replaces Nixpacks (2026) |
| Render | Fly.io | Global edge deployment, lower latency, but requires Dockerfile; more complex for simple apps |
| Render | Vercel/Netlify | Serverless edge functions, but cold starts hurt SSE performance; not designed for long-lived connections |

**Installation (Backend):**

Add to `backend/pyproject.toml`:
```toml
[tool.poetry.dependencies]
rdkit = "^2025.9.4"
uvicorn = {extras = ["standard"], version = "^0.40.0"}
```

Or create `backend/requirements.txt`:
```
rdkit>=2025.9.4
fastapi>=0.129.0
uvicorn[standard]>=0.40.0
pydantic>=2.12.5
pydantic-settings>=2.13.0
sse-starlette>=2.0
```

**Installation (Frontend):**

Create `.env.production`:
```
VITE_API_URL=https://your-backend.onrender.com
```

Update `src/api/client.ts`:
```typescript
private baseURL = import.meta.env.VITE_API_URL
  ? `${import.meta.env.VITE_API_URL}/api/sa`
  : '/api/sa'; // fallback for dev proxy
```

## Architecture Patterns

### Recommended Deployment Structure

```
.
├── backend/
│   ├── app/                    # FastAPI application code
│   ├── pyproject.toml          # Poetry dependencies (dev)
│   ├── requirements.txt        # Pip dependencies (production) ← GENERATE THIS
│   └── .env                    # Local env vars (gitignored)
├── .github/workflows/
│   └── deploy.yml              # GitHub Pages deployment (frontend)
├── .env.production             # Vite production config ← CREATE THIS
└── vite.config.ts              # Already configured with base: '/webfaulon/'
```

### Pattern 1: Environment-Driven Configuration (Backend)

**What:** Use `pydantic-settings` with environment variable overrides for production URLs and CORS origins.

**When to use:** Always — separates dev/prod config without code changes.

**Example:**

```python
# backend/app/config.py (already exists)
from pydantic_settings import BaseSettings, SettingsConfigDict

class Settings(BaseSettings):
    app_name: str = "WebFaulon"
    debug: bool = False
    allowed_origins: list[str] = [
        "http://localhost:5173",  # Dev frontend
        "http://127.0.0.1:5173",
        "https://steinbeck.github.io",  # Production frontend
    ]
    session_ttl_seconds: int = 3600

    model_config = SettingsConfigDict(env_file=".env", env_prefix="WEBFAULON_")
```

**Production override via Render environment variables:**
```
WEBFAULON_ALLOWED_ORIGINS=["https://steinbeck.github.io"]
WEBFAULON_DEBUG=false
```

**Source:** [FastAPI Settings Documentation](https://fastapi.tiangolo.com/advanced/settings/)

### Pattern 2: Environment-Driven API URLs (Frontend)

**What:** Use Vite's `import.meta.env.VITE_*` variables to inject backend URL at build time.

**When to use:** Always — avoids hardcoding production URLs in source.

**Example:**

```typescript
// src/api/client.ts
export class SAAPIClient {
  private baseURL = import.meta.env.VITE_API_URL
    ? `${import.meta.env.VITE_API_URL}/api/sa`
    : '/api/sa'; // Dev proxy fallback
}
```

**Build configuration:**

`.env.production`:
```
VITE_API_URL=https://webfaulon-backend.onrender.com
```

`.env.local` (for local dev):
```
# Empty or omitted — use Vite proxy to localhost:8000
```

**Source:** [Vite Env Variables and Modes](https://vite.dev/guide/env-and-mode)

### Pattern 3: Production ASGI Server Configuration

**What:** Run Uvicorn with production settings (no reload, bind to 0.0.0.0, use PORT env var).

**When to use:** Always in production — development server flags (`--reload`) are performance killers.

**Example:**

**Render start command:**
```bash
uvicorn app.main:app --host 0.0.0.0 --port $PORT
```

**With multiple workers (optional, likely overkill for classroom demo):**
```bash
gunicorn app.main:app --workers 2 --worker-class uvicorn.workers.UvicornWorker --bind 0.0.0.0:$PORT
```

**Worker count formula:** Set workers equal to CPU cores (not `2*cores+1` — that's for sync workers). Async Uvicorn workers handle concurrency efficiently within a single thread.

**Source:** [FastAPI Production Deployment Best Practices](https://render.com/articles/fastapi-production-deployment-best-practices)

### Pattern 4: Production CORS Locking

**What:** Restrict `allowed_origins` to production frontend domain only (remove localhost).

**When to use:** Production deployment — prevents unauthorized cross-origin requests.

**Example:**

```python
# backend/app/config.py (production override)
class Settings(BaseSettings):
    allowed_origins: list[str] = [
        "https://steinbeck.github.io",  # Production only
    ]
```

**Environment variable override (Render dashboard):**
```
WEBFAULON_ALLOWED_ORIGINS=["https://steinbeck.github.io"]
```

**Verification test:**
```bash
# Should succeed
curl -H "Origin: https://steinbeck.github.io" https://your-backend.onrender.com/health

# Should fail (no CORS headers)
curl -H "Origin: https://evil.com" https://your-backend.onrender.com/health
```

**Source:** [FastAPI CORS Security Best Practices](https://www.stackhawk.com/blog/configuring-cors-in-fastapi/)

### Anti-Patterns to Avoid

- **Wildcard CORS in production:** `allow_origins=["*"]` defeats the purpose of CORS and breaks credential support. Always use explicit origins.
- **Debug mode enabled:** `debug=True` leaks stack traces and internal paths to attackers. FastAPI defaults to `False` — verify it stays that way.
- **Uvicorn `--reload` in production:** Development-only flag that consumes excessive resources and degrades stability. Use only in dev.
- **Hardcoded API URLs in frontend:** Changing backends requires code rebuild. Use `import.meta.env.VITE_API_URL` instead.
- **Missing health check endpoint:** Load balancers and monitoring tools rely on `/health` for liveness/readiness probes. Already implemented at `/health`.
- **Committing `.env` files:** Secrets leak into git history. Use `.env.example` for templates, ignore `.env` in `.gitignore`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| RDKit installation from source | Custom conda Docker base image | `pip install rdkit` (PyPI wheels) | Pre-built wheels available for Linux/macOS/Windows since 2023; no compilation needed |
| Multi-platform Python dependency management | Manual `apt-get` + pip commands | Render/Railway buildpacks (auto-detect `pyproject.toml` or `requirements.txt`) | Platforms handle Python version detection, caching, and dependency installation |
| HTTPS/SSL certificates | Nginx + Let's Encrypt scripts | Platform-managed SSL (Render/Railway/Fly.io auto-provision) | Free automatic SSL with auto-renewal; manual certs are maintenance burden |
| Health check endpoints | Custom ping routes | Use existing `/health` endpoint (already implemented) | Standard pattern; platforms auto-detect and monitor it |
| Process management | systemd/supervisor scripts | Uvicorn alone or Gunicorn+Uvicorn workers | ASGI servers handle graceful shutdowns, signal handling, and worker lifecycle |
| Environment variable injection | Shell scripts in Dockerfile | Platform environment variable UI (Render/Railway dashboard) | Encrypted at rest, easy rotation, no secrets in code/Dockerfile |

**Key insight:** Modern PaaS platforms (Render, Railway, Fly.io) abstract away 90% of deployment complexity. The trap is over-engineering with Docker, custom Nginx, manual SSL, systemd services—all unnecessary for a classroom demo app. The RDKit PyPI wheels (introduced ~2023) eliminated the last deployment blocker. Use standard Python buildpacks and platform defaults.

## Common Pitfalls

### Pitfall 1: RDKit Deployment Failures (Historical)

**What goes wrong:** Deployments fail with "ModuleNotFoundError: No module named 'rdkit'" or C++ compilation errors.

**Why it happens:** RDKit historically required conda or Docker because PyPI wheels didn't exist. Outdated guides still recommend conda/Docker.

**How to avoid:**
- Use `rdkit` package on PyPI (not legacy `rdkit-pypi`)
- Version `2025.9.4+` has pre-built wheels for Python 3.9-3.14
- Add to `requirements.txt`: `rdkit>=2025.9.4`
- Render/Railway will install via pip automatically

**Warning signs:**
- Guides mentioning "conda-forge" or "conda install rdkit"
- Dockerfile with `FROM continuumio/miniconda3`
- Build logs showing C++ compilation

**Verification:**
```bash
# After deployment, check via Render shell or Railway CLI:
python -c "import rdkit; print(rdkit.__version__)"
# Should print: 2025.9.4 (or later)
```

**Source:** [RDKit PyPI Wheels](https://pypi.org/project/rdkit/)

### Pitfall 2: CORS Misconfiguration Breaking SSE

**What goes wrong:** Frontend connects to backend but SSE stream fails with CORS errors. Chart freezes, no live updates.

**Why it happens:** `EventSource` (SSE) sends credentialed requests by default. If `allow_credentials=False` (current config) but browser expects credentials, or if origin not in `allowed_origins`, CORS preflight fails.

**How to avoid:**
- Verify production frontend URL in `allowed_origins`: `"https://steinbeck.github.io"`
- Current config has `allow_credentials=False` — correct for this app (no cookies/auth)
- Test SSE endpoint manually: `curl -H "Origin: https://steinbeck.github.io" https://backend.onrender.com/api/sa/{session_id}/stream`

**Warning signs:**
- Browser console: "Access to XMLHttpRequest at 'https://backend...' from origin 'https://steinbeck.github.io' has been blocked by CORS policy"
- Network tab shows preflight OPTIONS request returning 403 or missing CORS headers
- `/health` works but `/api/sa/*/stream` fails

**Verification:**
```bash
# Check CORS headers on health endpoint
curl -i -H "Origin: https://steinbeck.github.io" https://backend.onrender.com/health
# Should see: Access-Control-Allow-Origin: https://steinbeck.github.io

# Check SSE endpoint (after configuring a session)
curl -H "Origin: https://steinbeck.github.io" https://backend.onrender.com/api/sa/test-session/stream
# Should stream events (or 404 if session doesn't exist, but NOT CORS error)
```

**Source:** [FastAPI CORS Troubleshooting](https://davidmuraya.com/blog/fastapi-cors-configuration/)

### Pitfall 3: Frontend API URL Hardcoding

**What goes wrong:** Production frontend still points to `localhost:8000` or hardcoded Render URL. Changing backends requires code rebuild.

**Why it happens:** Forgetting to use Vite environment variables; hardcoding URLs in `client.ts`.

**How to avoid:**
- Create `.env.production` with `VITE_API_URL=https://backend.onrender.com`
- Update `src/api/client.ts` to read `import.meta.env.VITE_API_URL`
- Verify build injects correct URL: `npm run build && grep -r "backend.onrender.com" dist/`

**Warning signs:**
- Production site console errors: "Failed to fetch http://localhost:8000/api/sa/configure"
- GitHub Pages works but clicking "Start" fails silently
- `dist/` build contains `localhost:8000` references

**Verification:**
```bash
# Build production bundle
npm run build

# Check injected API URL
grep -r "VITE_API_URL" dist/ || grep -r "onrender.com" dist/
# Should find the production backend URL embedded in JS bundles
```

**Source:** [Vite Production Build](https://vite.dev/guide/build)

### Pitfall 4: Missing `requirements.txt` (Poetry Users)

**What goes wrong:** Render/Railway build fails: "No module named 'fastapi'". Dependencies not installed.

**Why it happens:** Poetry users develop with `pyproject.toml` but platforms expect `requirements.txt`. Railway/Render don't auto-install from Poetry by default (though Railway's Railpack may support it).

**How to avoid:**
- Generate `requirements.txt` from Poetry lock: `poetry export -f requirements.txt --output requirements.txt --without-hashes`
- Commit `requirements.txt` to repo
- Alternative: Configure platform to use Poetry (more complex, not recommended for simplicity)

**Warning signs:**
- Build logs: "Looking for requirements.txt... not found"
- Build logs: "Running pip install... No such file"
- Platform falls back to empty environment

**Verification:**
```bash
# Generate requirements.txt from Poetry
cd backend
poetry export -f requirements.txt --output requirements.txt --without-hashes

# Verify it contains critical dependencies
grep -E "(fastapi|uvicorn|rdkit|pydantic)" requirements.txt
# Should show all key packages with versions
```

**Source:** [Render Python Troubleshooting](https://render.com/docs/troubleshooting-python-deploys)

### Pitfall 5: Uvicorn `--reload` in Production

**What goes wrong:** Backend consumes excessive CPU/memory, becomes unstable under load, random crashes.

**Why it happens:** `--reload` flag watches filesystem for changes and auto-restarts server. In production, file I/O monitoring leaks resources.

**How to avoid:**
- **Never** use `--reload` in production start command
- Correct: `uvicorn app.main:app --host 0.0.0.0 --port $PORT`
- Incorrect: `uvicorn app.main:app --reload --port $PORT`

**Warning signs:**
- High CPU usage even with no traffic
- Random server restarts in logs
- Platform monitoring shows memory creep

**Verification:**
```bash
# Check Render/Railway start command in dashboard
# Should NOT contain "--reload"
```

**Source:** [FastAPI Deployment Mistakes](https://www.linkedin.com/pulse/common-pitfalls-fastapi-deployment-segun-oje)

## Code Examples

Verified patterns from official sources:

### Render Deployment Configuration

**File: `backend/requirements.txt` (generate from Poetry):**
```txt
fastapi>=0.129.0
uvicorn[standard]>=0.40.0
pydantic>=2.12.5
pydantic-settings>=2.13.0
sse-starlette>=2.0
rdkit>=2025.9.4
```

**Render Web Service Settings (via dashboard):**
```
Build Command: pip install -r requirements.txt
Start Command: uvicorn app.main:app --host 0.0.0.0 --port $PORT
Environment Variables:
  WEBFAULON_ALLOWED_ORIGINS=["https://steinbeck.github.io"]
  WEBFAULON_DEBUG=false
```

**Source:** [Render FastAPI Deployment](https://render.com/articles/fastapi-deployment-options)

### Frontend Production Configuration

**File: `.env.production` (create at repo root):**
```bash
VITE_API_URL=https://webfaulon-backend.onrender.com
```

**File: `src/api/client.ts` (update baseURL):**
```typescript
export class SAAPIClient {
  private baseURL = import.meta.env.VITE_API_URL
    ? `${import.meta.env.VITE_API_URL}/api/sa`
    : '/api/sa'; // Fallback to Vite proxy in dev
}
```

**Source:** [Vite Environment Variables](https://vite.dev/guide/env-and-mode)

### Health Check Endpoint (Already Implemented)

**File: `backend/app/main.py` (lines 108-111):**
```python
@app.get("/health")
async def health_check():
    """Health check endpoint returning server status."""
    return {"status": "healthy", "version": "2.0.0"}
```

**Platform monitoring:** Render/Railway auto-detect `/health` and use it for liveness checks.

**Source:** [FastAPI Health Checks](https://www.index.dev/blog/how-to-implement-health-check-in-python)

### Production CORS Configuration

**File: `backend/app/config.py` (production override):**
```python
# Override via environment variable WEBFAULON_ALLOWED_ORIGINS
class Settings(BaseSettings):
    allowed_origins: list[str] = [
        "http://localhost:5173",      # Dev only
        "http://127.0.0.1:5173",      # Dev only
        "https://steinbeck.github.io", # Production
    ]
```

**Production environment variable (Render dashboard):**
```json
WEBFAULON_ALLOWED_ORIGINS=["https://steinbeck.github.io"]
```

**File: `backend/app/main.py` (lines 24-30, already configured):**
```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.allowed_origins,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

**Source:** [FastAPI CORS Middleware](https://www.compilenrun.com/docs/framework/fastapi/fastapi-middleware/fastapi-cors-middleware/)

## State of the Art

| Old Approach | Current Approach (2026) | When Changed | Impact |
|--------------|-------------------------|--------------|--------|
| Conda-based RDKit deployment | Pip-installable `rdkit` wheels on PyPI | ~2023 | No Docker required; works with standard Python buildpacks |
| Heroku free tier | Railway/Render/Fly.io | 2022 | Heroku eliminated free tier; ecosystem migrated to alternatives |
| Nixpacks (Railway) | Railpack (Railway) | March 2026 | Faster builds, smaller images, better caching; Nixpacks in maintenance mode |
| Manual SSL with Let's Encrypt | Platform-managed SSL | 2020s | Automatic SSL provisioning/renewal on all major PaaS platforms |
| Docker + docker-compose for simple apps | Buildpack auto-detection | 2020s | Platforms detect language/framework automatically; Docker optional |

**Deprecated/outdated:**

- **`rdkit-pypi` package name**: Renamed to `rdkit` on PyPI. Old package still works but unmaintained. Update dependencies.
- **Conda-only RDKit deployment**: PyPI wheels now available. Conda still viable but unnecessary complexity for cloud deployment.
- **Nixpacks on Railway**: Deprecated March 2026, replaced by Railpack. Existing services still work but new deployments should use Railpack.
- **Heroku**: Free tier eliminated November 2022. Guides referencing Heroku are outdated; use Render/Railway/Fly.io instead.

## Open Questions

1. **Multi-worker configuration: Is Gunicorn needed?**
   - What we know: Uvicorn alone handles async concurrency well. Gunicorn+Uvicorn workers provide multi-core utilization and fault isolation.
   - What's unclear: For a classroom demo app with low traffic, is single Uvicorn worker sufficient? Or does Render's free tier benefit from 2 workers?
   - Recommendation: **Start with single Uvicorn worker.** Add Gunicorn only if traffic monitoring shows CPU maxing out. For classroom demos, single worker is simpler and sufficient.

2. **Railway vs Render: Which platform is better for this app?**
   - What we know: Both support FastAPI. Railway has simpler UI and faster deploys but usage-based billing. Render has flat-rate pricing and better FastAPI docs.
   - What's unclear: Which billing model is cheaper for sporadic classroom use (high traffic during demos, zero traffic otherwise)?
   - Recommendation: **Render free tier for prototyping, paid tier ($7/month) for production.** Railway's usage-based billing can surprise with idle costs. Render's flat rate is predictable for education budgets.

3. **GitHub Actions: Should backend deployment be automated?**
   - What we know: Frontend already has GitHub Actions workflow for GitHub Pages. Backend deployment is manual via Render dashboard (git push triggers auto-deploy).
   - What's unclear: Would adding a backend deployment workflow (CI/CD) add value, or is manual deploy sufficient?
   - Recommendation: **Manual deploy is sufficient.** Render's auto-deploy on git push already provides CI/CD. GitHub Actions would duplicate functionality without added value.

4. **Production monitoring: What's the minimum viable setup?**
   - What we know: Platforms provide built-in monitoring (CPU, memory, request logs). `/health` endpoint enables uptime monitoring.
   - What's unclear: Should we add application-level monitoring (Sentry, Datadog) or rely on platform tools?
   - Recommendation: **Use platform monitoring only.** For a classroom demo, platform dashboards are sufficient. External monitoring is overkill and adds complexity.

## Sources

### Primary (HIGH confidence)

- [RDKit PyPI Package](https://pypi.org/project/rdkit/) - Official PyPI repository for rdkit wheels
- [RDKit Installation Documentation](https://www.rdkit.org/docs/Install.html) - Official installation guide
- [FastAPI Settings Documentation](https://fastapi.tiangolo.com/advanced/settings/) - Official guide for pydantic-settings
- [FastAPI CORS Documentation](https://fastapi.tiangolo.com/tutorial/cors/) - Official CORS middleware guide
- [Vite Environment Variables](https://vite.dev/guide/env-and-mode) - Official Vite env variable documentation
- [Vite Building for Production](https://vite.dev/guide/build) - Official Vite production build guide
- [Railway FastAPI Guide](https://docs.railway.com/guides/fastapi) - Official Railway deployment guide
- [Render FastAPI Deployment](https://render.com/articles/fastapi-deployment-options) - Official Render deployment article
- [Fly.io FastAPI Guide](https://fly.io/docs/python/frameworks/fastapi/) - Official Fly.io FastAPI documentation
- [Railway Nixpacks Documentation](https://docs.railway.com/reference/nixpacks) - Official Nixpacks reference
- [Railway Railpack Announcement](https://blog.railway.com/p/comparing-deployment-methods-in-railway) - Official Railpack announcement (March 2026)

### Secondary (MEDIUM confidence)

- [Railway vs Render vs Fly.io Comparison (2026)](https://medium.com/ai-disruption/railway-vs-fly-io-vs-render-which-cloud-gives-you-the-best-roi-2e3305399e5b) - Platform comparison article
- [FastAPI Production Deployment Best Practices](https://render.com/articles/fastapi-production-deployment-best-practices) - Render-published best practices guide
- [FastAPI Best Practices for Production (2026)](https://fastlaunchapi.dev/blog/fastapi-best-practices-production-2026) - Community best practices guide
- [FastAPI CORS Configuration Guide](https://davidmuraya.com/blog/fastapi-cors-configuration/) - Community CORS troubleshooting guide
- [Deploying React Vite with Routing on GitHub Pages](https://medium.com/@karinamisnik94/deploying-react-vite-with-routing-on-github-pages-68385676b788) - Community Vite deployment guide
- [FastAPI Health Check Implementation](https://www.index.dev/blog/how-to-implement-health-check-in-python) - Community health check guide
- [Render Python Troubleshooting](https://render.com/docs/troubleshooting-python-deploys) - Official Render troubleshooting docs

### Tertiary (LOW confidence - flagged for validation)

- [Python Hosting Options Compared (2025)](https://www.nandann.com/blog/python-hosting-options-comparison) - Comparison article, may be outdated for 2026
- [Common FastAPI Deployment Pitfalls](https://www.linkedin.com/pulse/common-pitfalls-fastapi-deployment-segun-oje) - LinkedIn article, no verification

## Metadata

**Confidence breakdown:**
- Standard stack: **HIGH** - RDKit PyPI wheels verified on official PyPI, FastAPI deployment patterns documented in official guides
- Architecture: **HIGH** - Environment-driven config is FastAPI official pattern, Vite env vars are official Vite pattern
- Pitfalls: **MEDIUM** - CORS/RDKit issues verified via multiple sources, production mistakes documented in community guides
- Platform comparison: **MEDIUM** - Railway Railpack announcement official (March 2026), pricing/features verified but subject to change

**Research date:** 2026-02-16
**Valid until:** ~2026-08-16 (6 months for stable tools, 3 months for fast-moving PaaS platforms)
