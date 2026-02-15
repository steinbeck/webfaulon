# Technology Stack

**Project:** WebFaulon v2.0 (Python Backend + TypeScript Frontend)
**Researched:** 2026-02-15 (Backend additions)
**Confidence:** HIGH (Backend), MEDIUM (Frontend — from 2026-02-14)

## Overview

**v1.0 (Browser-only):** Vite + TypeScript + RDKit.js + Web Workers
**v2.0 (Backend migration):** FastAPI + RDKit (Python) backend with SSE streaming + Vite frontend visualization layer

This document focuses on **backend stack additions** for v2.0. Frontend stack (Vite, Alpine.js, Chart.js) is **retained** from v1.0 but transitions from computation to visualization-only role.

---

## Backend Stack (NEW for v2.0)

### Core Backend Framework

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **FastAPI** | 0.129.0+ | ASGI web framework for Python backend | Industry standard for async Python APIs in 2026. Native async/await support, automatic OpenAPI docs, built on Starlette/Pydantic. Requires Python 3.10+ (dropped 3.9 in latest version). 5-50x faster than traditional WSGI frameworks due to async architecture. |
| **RDKit** | 2025.09.5+ | Molecular informatics library | The definitive open-source toolkit for cheminformatics. Provides native Mol objects, graph operations (Wiener index via `Chem.GetDistanceMatrix()`), SMILES parsing, and SVG rendering via `rdMolDraw2D.MolDraw2DSVG`. Essential for Faulon displacement operations on molecular graphs. |
| **Pydantic** | v2.x | Data validation and serialization | Required by FastAPI 0.128.0+ (v1 support removed). 5-50x faster than v1 due to Rust-based validation core. Provides runtime type checking using Python type hints. Use `BaseModel` for request/response schemas, `Field` for constraints. |
| **Uvicorn** | 0.34.0+ | ASGI server | Lightning-fast ASGI server for running FastAPI. Single-process works for development and low-traffic production. For production, run multiple Uvicorn workers under Gunicorn for process management and load balancing. |
| **sse-starlette** | 3.2.0+ | Server-Sent Events for FastAPI | Production-ready SSE implementation following W3C spec. Released Jan 17, 2026. Provides `EventSourceResponse` for streaming SA progress updates to frontend. Native FastAPI integration, automatic disconnect detection, multi-threaded event loop safety. Requires Python 3.9+. |

### Supporting Backend Libraries

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| **httpx** | 0.27.0+ | Async HTTP client | For FastAPI async testing with `AsyncClient`. Modern replacement for requests. Use in tests to make async calls to FastAPI endpoints without network overhead. |
| **python-multipart** | 0.0.9+ | Form data parsing | Required if accepting file uploads or form-encoded data. FastAPI dependency for multipart form handling. |
| **python-dotenv** | 1.0.0+ | Environment configuration | Load environment variables from `.env` files. Use for local development configuration (API keys, database URLs). |
| **pydantic-settings** | 2.x | Settings management | Extends Pydantic for environment-based configuration. Use for structured config from `.env` with validation. |

### Development Tools (Backend)

| Tool | Purpose | Notes |
|------|---------|-------|
| **Poetry** | Dependency and project management | Recommended over pip+requirements.txt for 2026. Uses `pyproject.toml` for config, generates `poetry.lock` for reproducible builds. Run `poetry config virtualenvs.in-project true` to keep `.venv/` in repo. Install deps with `poetry add fastapi uvicorn rdkit`, dev deps with `poetry add --group dev pytest pytest-asyncio pytest-cov`. |
| **pytest** | Testing framework | Standard for Python testing in 2026. Use with `pytest-asyncio` for async endpoint tests. TestClient (sync) for simple tests, AsyncClient for async/lifespan testing. |
| **pytest-asyncio** | Async test support | Enables `async def` test functions. Set `asyncio_default_fixture_loop_scope = "function"` in `pyproject.toml` for test isolation. Required for testing FastAPI async endpoints. |
| **pytest-cov** | Code coverage reporting | Built on coverage.py 7.13.4 (supports Python 3.10-3.15). Run `pytest --cov=app --cov-report=html` for detailed coverage reports. 2026 updates improved async code coverage tracking. |
| **Ruff** | Linting and formatting | Fast Rust-based linter/formatter replacing Black + isort + Flake8. Configure in `pyproject.toml`. Run `ruff check .` for linting, `ruff format .` for formatting. |
| **mypy** | Static type checking | Verify type hints match Pydantic models. Configure strict mode in `pyproject.toml`. Catches type errors before runtime. |

---

## Frontend Stack (RETAINED from v1.0)

### Core Frontend Framework

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **Vite** | ^6.x | Build tool and dev server | Industry standard for modern browser apps in 2026. Provides instant HMR, esbuild-powered TypeScript transpilation (20-30x faster than traditional tools). Zero-config dev experience with production-ready optimizations. |
| **TypeScript** | ^5.x | Type-safe development | Mandatory for maintainable scientific applications. Vite's esbuild integration provides instantaneous type stripping. Use strict mode for early error detection in API client code. |

### UI and Visualization (Frontend)

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **Chart.js** | ^4.x | Real-time SA progress charting | Retained from v1.0. Lightweight canvas-based charting for live SA progress updates via SSE. Receives data from backend, no computation. Alternative: Apache ECharts for >10K data points, but Chart.js sufficient for typical SA runs. |
| **Alpine.js** | ^3.x | Lightweight UI interactivity | 7.1KB framework perfect for parameter controls and UI state without SPA overhead. Declarative syntax in HTML reduces build complexity. Ideal for educational demos where students read source. Alternative: vanilla JS for minimal footprint. |

### Supporting Frontend Libraries

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| **EventSource API** | Native | SSE client for streaming SA progress | Browser-native API for consuming Server-Sent Events from FastAPI `/api/sa/stream` endpoint. No library needed. Falls back to polyfill for older browsers. |

---

## Installation

### Backend Setup (Poetry)

```bash
# Navigate to backend directory
cd backend/

# Initialize Poetry project
poetry init --name webfaulon-backend --python "^3.10"

# Core backend framework
poetry add fastapi uvicorn[standard] pydantic pydantic-settings

# RDKit (molecular operations)
# OPTION 1: Via Poetry (may have binary dependency issues)
poetry add rdkit

# OPTION 2: Via conda in Poetry venv (RECOMMENDED)
# Create conda env first, then use Poetry with system Python
# See "RDKit Installation Methods" section below

# SSE streaming
poetry add sse-starlette

# Supporting libraries
poetry add httpx python-multipart python-dotenv

# Dev dependencies
poetry add --group dev pytest pytest-asyncio pytest-cov httpx ruff mypy

# Install all dependencies
poetry install
```

### Frontend Setup (npm)

```bash
# Navigate to frontend directory (or root if monorepo)
cd frontend/  # or stay in root

# Core (already installed from v1.0, but shown for reference)
npm install alpinejs chart.js

# Build tooling (already installed from v1.0)
npm install -D vite typescript

# Dev dependencies (already installed from v1.0)
npm install -D vitest eslint prettier @typescript-eslint/parser @typescript-eslint/eslint-plugin
```

### RDKit Installation Methods

**Conda (Recommended):**
```bash
conda install -c conda-forge rdkit python=3.12
```
- Simpler binary dependency management
- Better for development environments
- Works seamlessly with C++ dependencies

**Pip (Alternative):**
```bash
pip install rdkit
```
- Pre-built wheels for Linux, Windows, macOS
- Faster for CI/CD if conda isn't available
- May have issues with binary dependencies on some platforms

**Recommendation:** Use conda for development, conda-forge Docker images for production. If CI/CD is pip-based, use `pip install rdkit` with pre-built wheels.

---

## Alternatives Considered

### Backend Alternatives

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| **FastAPI** | Flask + Flask-RESTX | If team is deeply familiar with Flask and doesn't need async. FastAPI is faster and has better async support, but Flask has larger ecosystem. |
| **sse-starlette** | Custom SSE implementation | Never. sse-starlette is production-ready and follows W3C spec. Custom implementations often have event loop binding issues. |
| **Poetry** | pip + requirements.txt | If team workflow is pip-based and Poetry adds friction. Poetry is superior for dependency management (lockfiles, separation of prod/dev deps) but has learning curve. |
| **Uvicorn + Gunicorn** | Hypercorn, Daphne | Uvicorn is standard. Hypercorn for HTTP/2 support, Daphne for Django Channels compatibility. For basic FastAPI, stick with Uvicorn. |
| **RDKit pip install** | RDKit conda install | **Conda preferred for RDKit**. Pip wheels exist but conda has better binary dependency management. For CI/CD, use conda-forge Docker images or mamba. |
| **pytest** | unittest | pytest is 2026 standard. unittest is Python stdlib but more verbose. Use pytest unless stdlib-only constraint exists. |

### Frontend Alternatives (from v1.0)

| Recommended | Alternative | Why Not |
|-------------|-------------|---------|
| **Alpine.js** | React | React adds 40KB+ overhead unnecessary for parameter form + visualization. Alpine's in-HTML syntax more educational (students see interactivity in source). No build complexity. |
| **Alpine.js** | Vanilla JS | Viable alternative. Choose Alpine if you want reactive data binding without manual DOM updates. Choose vanilla if targeting zero dependencies. |
| **Chart.js** | Apache ECharts | ECharts (canvas rendering) is 10x faster for >10K points. For typical SA runs (<5K steps), Chart.js is simpler and sufficient. Upgrade to ECharts if performance issues. |

---

## What NOT to Use

### Backend Anti-Patterns

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| **Pydantic v1** | Deprecated as of FastAPI 0.128.0. No longer supported. 5-50x slower than v2. | **Pydantic v2** with Rust-based validation. Migration guide at docs.pydantic.dev. |
| **Python 3.9 or below** | FastAPI 0.129.0 dropped Python 3.9 support. Pydantic v2 and modern type hints require 3.10+. | **Python 3.10, 3.11, or 3.12**. 3.12 is fastest for async workloads. |
| **Synchronous route handlers** | Blocks event loop, kills FastAPI performance. SSE streaming requires async. | **async def** route handlers. Use `await` for I/O. For CPU-bound work (SA iterations), offload to threadpool or process pool. |
| **requests library** | Synchronous, blocks event loop. Not compatible with async FastAPI. | **httpx** for async HTTP calls. Drop-in replacement with `httpx.AsyncClient()`. |
| **WebSockets for SA streaming** | Bidirectional overhead for one-way server-to-client stream. Requires more complex frontend code. | **SSE (Server-Sent Events)** via sse-starlette. Simpler, auto-reconnect, works over HTTP, firewall-friendly. |
| **Flask** | WSGI-based (synchronous). Async support is bolted-on, not native. Slower for async workloads. | **FastAPI** with native async/await and ASGI. |
| **rdkit-pypi package** | Renamed. Old PyPI package name. | **rdkit** (new PyPI name as of 2019.03+). Update dependencies from `rdkit-pypi` to `rdkit`. |

### Frontend Anti-Patterns (from v1.0, retained for v2.0)

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| **RDKit.js (v2.0 context)** | Replaced by Python RDKit backend. Remove from frontend dependencies. | Python RDKit in FastAPI backend for all molecular operations. |
| **Web Workers (v2.0 context)** | No longer needed. Computation moved to backend. Remove Comlink and worker infrastructure. | Backend async tasks. Frontend is thin visualization layer. |
| **CDN imports for dependencies** | Vite bundler ensures correct resolution and optimization. CDN adds CORS complexity. | `npm install` + Vite bundler. |

---

## Stack Patterns by Variant

### Development Mode

**Backend:**
- Use single Uvicorn process: `uvicorn app.main:app --reload --port 8000`
- Enable CORS for Vite dev server (http://localhost:5173): `app.add_middleware(CORSMiddleware, allow_origins=["http://localhost:5173"], ...)`
- Use `.env` file with `python-dotenv` for config
- RDKit via conda for easier dependency management

**Frontend:**
- Use Vite dev server: `npm run dev` (port 5173)
- API calls to http://localhost:8000 (FastAPI backend)
- SSE connection to `/api/sa/stream` endpoint

### Production Mode

**Backend:**
- Use Gunicorn with Uvicorn workers: `gunicorn app.main:app -k uvicorn.workers.UvicornWorker -w 4 --bind 0.0.0.0:8000`
- Workers = 2-4 × CPU cores (balance parallelism vs memory)
- Docker with conda-forge base image: `FROM continuumio/miniconda3`
- Reverse proxy (Nginx) for SSL termination and static file serving
- CORS locked to production frontend domain only

**Frontend:**
- Vite production build: `npm run build`
- Deploy static assets to GitHub Pages, Netlify, or Vercel
- API calls to production FastAPI domain (e.g., https://api.webfaulon.com)
- SSE connection to production `/api/sa/stream` endpoint

### Testing

**Backend:**
- Use `TestClient` for synchronous tests (simple routes, no lifespan events)
- Use `AsyncClient` with `ASGITransport` for async tests (SSE endpoints, lifespan events)
- Mock RDKit operations for fast unit tests (test SA logic without molecular graph overhead)
- Integration tests with real RDKit Mol objects for end-to-end validation

**Frontend:**
- Vitest for unit tests (Alpine.js components, Chart.js integration)
- Mock fetch/EventSource for API tests (no backend dependency)
- E2E tests with Playwright for full SSE streaming validation

---

## Version Compatibility

### Backend

| Package A | Compatible With | Notes |
|-----------|-----------------|-------|
| FastAPI 0.129.0 | Python 3.10+ | **Breaking:** Dropped Python 3.9 support in 0.129.0. Pin Python ≥3.10. |
| FastAPI 0.128.0+ | Pydantic v2 only | **Breaking:** Removed Pydantic v1 support in 0.128.0. Use `pydantic>=2.0`. |
| sse-starlette 3.2.0 | Python 3.9+ | Requires Python ≥3.9. Compatible with FastAPI 0.129.0 (both use Starlette). |
| RDKit 2025.09.5 | Python 3.10-3.12 | Works with modern Python. **Conda install recommended** for binary dependencies. |
| pytest-cov 7.x | coverage.py 7.13.4 | Supports Python 3.10-3.15. 2026 updates improved async coverage tracking. |
| Uvicorn 0.34.0 | Python 3.10+ | Use `uvicorn[standard]` for HTTP/2 and better performance. |
| httpx 0.27.0 | Python 3.10+ | Modern async HTTP client. Required for AsyncClient testing in FastAPI. |

### Frontend

| Package A | Compatible With | Notes |
|-----------|-----------------|-------|
| Vite ^6.x | Node.js 18+ | Breaking changes in Vite 6 around worker imports (irrelevant for v2.0 since workers removed). |
| Chart.js ^4.x | TypeScript ^5.x | Types included. Modern TS has native support. |
| Alpine.js ^3.x | Any modern browser | No build step needed. Works with ES6 modules via CDN or npm. |

---

## CORS Configuration for Vite Frontend

FastAPI backend must allow Vite dev server (http://localhost:5173) during development:

```python
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

# Development CORS (allow Vite dev server)
origins = [
    "http://localhost:5173",  # Vite default dev server
    "http://localhost:5174",  # Vite preview server
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,  # Lockdown in production to deployment domain
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

**Production:** Replace `origins` with production domain. Never use `"*"` with `allow_credentials=True`.

---

## Backend Project Structure

Recommended structure for Poetry-based FastAPI + RDKit backend:

```
backend/
├── pyproject.toml          # Poetry config, dependencies, dev tools
├── poetry.lock             # Locked dependency versions
├── .env                    # Local environment variables (gitignored)
├── app/
│   ├── __init__.py
│   ├── main.py             # FastAPI app, CORS, startup/shutdown
│   ├── config.py           # Pydantic settings from .env
│   ├── models.py           # Pydantic request/response models
│   ├── api/
│   │   ├── __init__.py
│   │   ├── routes.py       # REST endpoints
│   │   └── sse.py          # SSE streaming endpoints
│   ├── services/
│   │   ├── __init__.py
│   │   ├── sa_engine.py    # SA algorithm (Faulon displacement)
│   │   └── rdkit_ops.py    # RDKit Mol operations (Wiener index, SVG rendering)
│   └── target_functions/
│       ├── __init__.py
│       ├── base.py         # Abstract base class for target functions
│       ├── wiener.py       # Wiener index target function
│       └── registry.py     # Pluggable target function registry
├── tests/
│   ├── __init__.py
│   ├── conftest.py         # Pytest fixtures
│   ├── test_api.py         # API endpoint tests
│   ├── test_sse.py         # SSE streaming tests
│   ├── test_sa_engine.py   # SA algorithm unit tests
│   └── test_rdkit_ops.py   # RDKit operations tests
└── README.md
```

---

## Confidence Assessment

### Backend Stack

| Technology | Confidence | Source | Notes |
|------------|------------|--------|-------|
| FastAPI | **HIGH** | Official release notes, PyPI | Version 0.129.0 verified. Python 3.10+ requirement confirmed. |
| RDKit | **HIGH** | Official docs, PyPI | Version 2025.09.5 verified. Conda vs pip installation methods confirmed. |
| sse-starlette | **HIGH** | PyPI, official docs | Version 3.2.0 released Jan 17, 2026. W3C spec compliance verified. |
| Pydantic v2 | **HIGH** | Official docs, migration guides | Rust-based validation, 5-50x performance improvement verified. |
| Uvicorn + Gunicorn | **HIGH** | Multiple 2026 production guides | Industry standard for FastAPI production deployment. |
| Poetry | **HIGH** | 2026 best practices articles | De facto standard for Python dependency management in 2026. |
| pytest + pytest-asyncio | **HIGH** | Official FastAPI testing docs | Standard for async FastAPI testing. coverage.py 7.13.4 supports Python 3.10-3.15. |

### Frontend Stack (from v1.0)

| Technology | Confidence | Source | Notes |
|------------|------------|--------|-------|
| Vite | **HIGH** | Official documentation, 2026 production setup guides | De facto standard for TypeScript browser apps in 2026. |
| Alpine.js | **MEDIUM** | Framework comparison articles 2026 | Lightweight category recommendation consistent across sources. Not as mainstream as React but established for progressive enhancement use cases. |
| Chart.js | **MEDIUM** | npm-compare, official docs | Sufficient for typical SA runs. Upgrade to ECharts if >10K points. |

---

## Sources

### Backend Stack (HIGH Confidence)

**Official Documentation:**
- [FastAPI Release Notes](https://fastapi.tiangolo.com/release-notes/) — Latest version 0.129.0, Python 3.10+ requirement
- [RDKit Installation](https://www.rdkit.org/docs/Install.html) — Conda vs pip, Python version support
- [sse-starlette PyPI](https://pypi.org/project/sse-starlette/) — Version 3.2.0, SSE implementation details
- [RDKit Draw Module](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html) — SVG rendering with MolDraw2DSVG
- [FastAPI CORS](https://fastapi.tiangolo.com/tutorial/cors/) — CORSMiddleware configuration
- [FastAPI Async Tests](https://fastapi.tiangolo.com/advanced/async-tests/) — AsyncClient testing patterns
- [Pydantic Welcome](https://docs.pydantic.dev/latest/) — Pydantic v2 features

**Community Resources (MEDIUM-HIGH Confidence):**
- [FastAPI Best Practices Production 2026](https://fastlaunchapi.dev/blog/fastapi-best-practices-production-2026) — Production deployment patterns
- [Uvicorn vs Gunicorn FastAPI 2026](https://www.geeksforgeeks.org/python/fast-api-gunicorn-vs-uvicorn/) — Worker configuration for production
- [Poetry for Python Projects 2026](https://thelinuxcode.com/poetry-for-python-projects-overview-benefits-and-practical-workflow-2026/) — Modern dependency management
- [Server-Sent Events with FastAPI](https://medium.com/@nandagopal05/server-sent-events-with-python-fastapi-f1960e0c8e4b) — SSE implementation patterns
- [Pydantic Complete Guide 2026](https://devtoolbox.dedyn.io/blog/pydantic-complete-guide) — Pydantic v2 features and performance
- [Testing FastAPI Applications](https://pythoneo.com/testing-fastapi-applications/) — Pytest, coverage, best practices
- [The Complete FastAPI × pytest Guide](https://blog.greeden.me/en/2026/01/06/the-complete-fastapi-x-pytest-guide-building-fearless-to-change-apis-with-unit-tests-api-tests-integration-tests-and-mocking-strategies/) — Comprehensive testing strategies (Feb 2026)

**RDKit Cheminformatics (MEDIUM Confidence):**
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html) — Wiener index calculation examples
- [Revisiting Wiener Index](https://bertiewooster.github.io/2023/03/10/Revisiting-a-Classic-Cheminformatics-Paper-The-Wiener-Index.html) — Implementation patterns with GetDistanceMatrix

### Frontend Stack (from 2026-02-14 research)

**RDKit.js:**
- [GitHub - rdkit/rdkit-js](https://github.com/rdkit/rdkit-js) — Official repository
- [@rdkit/rdkit on npm](https://libraries.io/npm/@rdkit%2Frdkit) — Version verification (2025.3.4-1.0.0)

**Charting Libraries:**
- [JavaScript Chart Libraries in 2026 - Luzmo](https://www.luzmo.com/blog/javascript-chart-libraries) — Chart.js vs ECharts performance
- [npm-compare: chart.js vs echarts](https://npm-compare.com/chart.js,echarts) — Package statistics

**Frontend Frameworks:**
- [Alpine.js: The Minimalist JavaScript Framework - Medium](https://medium.com/@zulfikarditya/alpine-js-the-minimalist-javascript-framework-for-modern-web-development-839382997988) — Size and approach

**Vite and TypeScript:**
- [Vite Official Documentation](https://vite.dev/guide/) — Official source

---
*Stack research for: Python RDKit backend with FastAPI and SSE streaming + Vite frontend*
*Backend researched: 2026-02-15*
*Frontend researched: 2026-02-14*
