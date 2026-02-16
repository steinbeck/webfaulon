"""Main FastAPI application with CORS, health check, and error handling."""

import asyncio

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as StarletteHTTPException

from app.config import Settings
from app.models.errors import ErrorResponse
from app.dependencies import session_manager
from app.api.sa_configure import router as configure_router
from app.api.sa_control import router as control_router
from app.api.sa_status import router as status_router
from app.api.sa_stream import router as stream_router

settings = Settings()

app = FastAPI(title="WebFaulon API", version="2.0.0")

# Add CORS middleware FIRST (before any other middleware)
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.allowed_origins,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Exception handlers
@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle validation errors with structured ErrorResponse."""
    return JSONResponse(
        status_code=422,
        content=ErrorResponse(
            error="Validation Error",
            detail=str(exc),
            status_code=422,
        ).model_dump(),
    )


@app.exception_handler(StarletteHTTPException)
async def http_exception_handler(request: Request, exc: StarletteHTTPException):
    """Handle HTTP exceptions with structured ErrorResponse (catches default 404s)."""
    return JSONResponse(
        status_code=exc.status_code,
        content=ErrorResponse(
            error=exc.detail or "HTTP Error",
            detail=None,
            status_code=exc.status_code,
        ).model_dump(),
    )


@app.exception_handler(HTTPException)
async def fastapi_http_exception_handler(request: Request, exc: HTTPException):
    """Handle FastAPI HTTP exceptions with structured ErrorResponse."""
    return JSONResponse(
        status_code=exc.status_code,
        content=ErrorResponse(
            error=exc.detail or "HTTP Error",
            detail=None,
            status_code=exc.status_code,
        ).model_dump(),
    )


@app.exception_handler(Exception)
async def generic_exception_handler(request: Request, exc: Exception):
    """Handle unexpected errors with structured ErrorResponse."""
    detail = str(exc) if settings.debug else None
    return JSONResponse(
        status_code=500,
        content=ErrorResponse(
            error="Internal Server Error",
            detail=detail,
            status_code=500,
        ).model_dump(),
    )


# Register API routers
app.include_router(configure_router)
app.include_router(control_router)
app.include_router(status_router)
app.include_router(stream_router)


# Startup event for periodic session cleanup
@app.on_event("startup")
async def startup_event():
    """Start periodic session cleanup."""

    async def periodic_cleanup():
        while True:
            await asyncio.sleep(300)  # Every 5 minutes
            session_manager.cleanup_expired()

    asyncio.create_task(periodic_cleanup())


# Health check endpoint
@app.get("/health")
async def health_check():
    """Health check endpoint returning server status."""
    return {"status": "healthy", "version": "2.0.0"}
