"""Application configuration using Pydantic settings."""

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    app_name: str = "WebFaulon"
    debug: bool = False
    allowed_origins: list[str] = [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
        "https://steinbeck.github.io",
    ]

    model_config = SettingsConfigDict(env_file=".env", env_prefix="WEBFAULON_")
