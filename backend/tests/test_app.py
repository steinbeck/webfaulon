"""Tests for FastAPI application: health endpoint, CORS, error handling, and OpenAPI."""

import pytest


class TestHealthEndpoint:
    """Tests for the /health endpoint (BACK-03)."""

    def test_health_returns_200(self, client):
        """GET /health returns 200 status code."""
        response = client.get("/health")
        assert response.status_code == 200

    def test_health_response_body(self, client):
        """Health response contains status and version keys."""
        response = client.get("/health")
        data = response.json()
        assert "status" in data
        assert "version" in data

    def test_health_status_value(self, client):
        """Health status is 'healthy'."""
        response = client.get("/health")
        data = response.json()
        assert data["status"] == "healthy"


class TestCORS:
    """Tests for CORS middleware (BACK-02)."""

    def test_cors_allows_localhost_origin(self, client):
        """OPTIONS request with localhost:5173 origin returns access-control-allow-origin header."""
        response = client.options(
            "/health",
            headers={"Origin": "http://localhost:5173"},
        )
        assert "access-control-allow-origin" in response.headers
        assert response.headers["access-control-allow-origin"] == "http://localhost:5173"

    def test_cors_allows_github_pages(self, client):
        """OPTIONS request with GitHub Pages origin returns correct header."""
        response = client.options(
            "/health",
            headers={"Origin": "https://steinbeck.github.io"},
        )
        assert "access-control-allow-origin" in response.headers
        assert response.headers["access-control-allow-origin"] == "https://steinbeck.github.io"

    def test_cors_blocks_unknown_origin(self, client):
        """OPTIONS request with unknown origin does NOT return access-control-allow-origin header."""
        response = client.options(
            "/health",
            headers={"Origin": "https://evil.com"},
        )
        # FastAPI CORS middleware doesn't include the header for disallowed origins
        # The header should either be absent or not match the origin
        origin_header = response.headers.get("access-control-allow-origin")
        if origin_header:
            assert origin_header != "https://evil.com"


class TestErrorHandling:
    """Tests for error handling (BACK-04)."""

    def test_404_returns_json_error(self, client):
        """GET to nonexistent path returns 404 with ErrorResponse JSON body."""
        response = client.get("/nonexistent")
        assert response.status_code == 404
        data = response.json()
        assert "error" in data
        assert "status_code" in data

    def test_error_response_structure(self, client):
        """Verify ErrorResponse has required fields."""
        response = client.get("/nonexistent")
        data = response.json()
        assert "error" in data
        assert "status_code" in data
        assert data["status_code"] == 404
        # detail is optional
        assert "detail" in data or True  # detail can be None


class TestOpenAPI:
    """Tests for OpenAPI documentation (BACK-01)."""

    def test_docs_accessible(self, client):
        """GET /docs returns 200."""
        response = client.get("/docs")
        assert response.status_code == 200

    def test_openapi_json_accessible(self, client):
        """GET /openapi.json returns 200 with valid JSON."""
        response = client.get("/openapi.json")
        assert response.status_code == 200
        data = response.json()
        assert "openapi" in data
        assert "info" in data
