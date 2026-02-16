"""Component registry for scoring components."""
from app.core.scoring.protocol import ScoringComponent
from app.core.scoring.wiener import WienerIndexComponent
from app.core.scoring.logp import LogPComponent


class ComponentRegistry:
    """
    Registry for scoring components.

    Provides registration, retrieval, and listing of scoring components
    with validation to prevent duplicate names and provide helpful error messages.
    """

    def __init__(self):
        """Initialize empty registry."""
        self._components: dict[str, ScoringComponent] = {}

    def register(self, component: ScoringComponent) -> None:
        """
        Register a scoring component.

        Args:
            component: ScoringComponent to register

        Raises:
            ValueError: If component with same name already registered
        """
        if component.name in self._components:
            raise ValueError(f"Component '{component.name}' already registered")
        self._components[component.name] = component

    def get(self, name: str) -> ScoringComponent:
        """
        Get scoring component by name.

        Args:
            name: Component name to retrieve

        Returns:
            ScoringComponent with given name

        Raises:
            KeyError: If component not found (with helpful message listing available components)
        """
        if name not in self._components:
            available = ', '.join(sorted(self._components.keys()))
            raise KeyError(f"Unknown component '{name}'. Available: {available}")
        return self._components[name]

    def list_components(self) -> list[str]:
        """
        List all registered component names.

        Returns:
            List of component names
        """
        return list(self._components.keys())


# Create default registry with built-in components
_registry = ComponentRegistry()
_registry.register(WienerIndexComponent())
_registry.register(LogPComponent())


def get_registry() -> ComponentRegistry:
    """
    Get the default component registry.

    Returns:
        Default ComponentRegistry singleton
    """
    return _registry
