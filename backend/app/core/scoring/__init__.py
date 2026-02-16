"""Scoring component framework - Protocol interface and component registry."""
from app.core.scoring.protocol import ScoringComponent
from app.core.scoring.registry import ComponentRegistry, get_registry

__all__ = ['ScoringComponent', 'ComponentRegistry', 'get_registry']
