"""Tests for scoring components, protocol compliance, and registry."""
import pytest
from app.core.scoring.wiener import WienerIndexComponent
from app.core.scoring.logp import LogPComponent
from app.core.scoring.registry import ComponentRegistry, get_registry
from app.core.wiener import compute_wiener_index
from app.core.molecule import MoleculeGraph


class TestWienerIndexComponent:
    """Test WienerIndexComponent satisfies protocol and computes correctly"""

    def test_name_is_wiener_index(self):
        """Component name should be 'wiener_index'"""
        component = WienerIndexComponent()
        assert component.name == "wiener_index"

    def test_linear_hexane_wiener_35(self):
        """Linear hexane (C6) should have Wiener Index 35"""
        graph = MoleculeGraph.create_linear_alkane(6)
        component = WienerIndexComponent()
        assert component.compute(graph) == 35

    def test_isobutane_wiener_9(self):
        """Isobutane should have Wiener Index 9"""
        graph = MoleculeGraph.create_branched('isobutane')
        component = WienerIndexComponent()
        assert component.compute(graph) == 9

    def test_neopentane_wiener_16(self):
        """Neopentane should have Wiener Index 16"""
        graph = MoleculeGraph.create_branched('neopentane')
        component = WienerIndexComponent()
        assert component.compute(graph) == 16

    def test_cyclohexane_wiener_27(self):
        """Cyclohexane should have Wiener Index 27"""
        graph = MoleculeGraph.create_cyclohexane()
        component = WienerIndexComponent()
        assert component.compute(graph) == 27

    def test_matches_standalone_function(self):
        """Component should return identical values to standalone compute_wiener_index()"""
        graph = MoleculeGraph.create_linear_alkane(6)
        component = WienerIndexComponent()
        standalone_result = compute_wiener_index(graph)
        component_result = component.compute(graph)
        assert component_result == standalone_result


class TestLogPComponent:
    """Test LogPComponent satisfies protocol and computes correctly"""

    def test_name_is_logp(self):
        """Component name should be 'logp'"""
        component = LogPComponent()
        assert component.name == "logp"

    def test_linear_vs_branched_logp_differs(self):
        """LogP should differ between linear and branched isomers"""
        linear = MoleculeGraph.create_linear_alkane(6)
        branched = MoleculeGraph.create_branched('isobutane')
        component = LogPComponent()
        linear_logp = component.compute(linear)
        branched_logp = component.compute(branched)
        assert linear_logp != branched_logp, "LogP should vary across constitutional isomers"

    def test_hexane_logp_positive(self):
        """Hexane (hydrocarbon) should have positive LogP"""
        graph = MoleculeGraph.create_linear_alkane(6)
        component = LogPComponent()
        logp = component.compute(graph)
        assert logp > 0, f"Expected positive LogP for hydrocarbon, got {logp}"


class TestProtocolCompliance:
    """Test that components satisfy the ScoringComponent protocol"""

    def test_wiener_has_name_and_compute(self):
        """WienerIndexComponent should have name attribute and compute method"""
        component = WienerIndexComponent()
        assert hasattr(component, 'name')
        assert callable(getattr(component, 'compute'))

    def test_logp_has_name_and_compute(self):
        """LogPComponent should have name attribute and compute method"""
        component = LogPComponent()
        assert hasattr(component, 'name')
        assert callable(getattr(component, 'compute'))


class TestComponentRegistry:
    """Test ComponentRegistry registration, retrieval, and validation"""

    def test_default_registry_has_wiener(self):
        """Default registry should have wiener_index component"""
        registry = get_registry()
        component = registry.get("wiener_index")
        assert component is not None

    def test_default_registry_has_logp(self):
        """Default registry should have logp component"""
        registry = get_registry()
        component = registry.get("logp")
        assert component is not None

    def test_list_components(self):
        """list_components() should return both component names"""
        registry = get_registry()
        components = registry.list_components()
        assert "wiener_index" in components
        assert "logp" in components

    def test_unknown_component_raises_keyerror(self):
        """Getting unknown component should raise KeyError"""
        registry = get_registry()
        with pytest.raises(KeyError):
            registry.get("nonexistent")

    def test_unknown_component_error_lists_available(self):
        """KeyError message should list available components"""
        registry = get_registry()
        try:
            registry.get("nonexistent")
            assert False, "Should have raised KeyError"
        except KeyError as e:
            error_msg = str(e)
            assert "wiener_index" in error_msg
            assert "logp" in error_msg

    def test_duplicate_registration_raises(self):
        """Registering duplicate component name should raise ValueError"""
        registry = ComponentRegistry()
        component1 = WienerIndexComponent()
        component2 = WienerIndexComponent()

        registry.register(component1)
        with pytest.raises(ValueError):
            registry.register(component2)
