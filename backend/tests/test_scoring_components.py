"""Tests for scoring components, protocol compliance, and registry."""
import pytest
from app.core.scoring.wiener import WienerIndexComponent
from app.core.scoring.molecular_weight import MolecularWeightComponent
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


class TestMolecularWeightComponent:
    """Test MolecularWeightComponent satisfies protocol and computes correctly"""

    def test_name_is_molecular_weight(self):
        """Component name should be 'molecular_weight'"""
        component = MolecularWeightComponent()
        assert component.name == "molecular_weight"

    def test_hexane_mw(self):
        """Hexane (C6H14) should have molecular weight ~86.18 Da"""
        graph = MoleculeGraph.create_linear_alkane(6)
        component = MolecularWeightComponent()
        mw = component.compute(graph)
        assert 85.0 < mw < 87.0, f"Expected ~86.18, got {mw}"

    def test_methane_mw(self):
        """Methane (CH4) should have molecular weight ~16.04 Da"""
        graph = MoleculeGraph.create_linear_alkane(1)
        component = MolecularWeightComponent()
        mw = component.compute(graph)
        assert 15.0 < mw < 17.0, f"Expected ~16.04, got {mw}"


class TestProtocolCompliance:
    """Test that components satisfy the ScoringComponent protocol"""

    def test_wiener_has_name_and_compute(self):
        """WienerIndexComponent should have name attribute and compute method"""
        component = WienerIndexComponent()
        assert hasattr(component, 'name')
        assert callable(getattr(component, 'compute'))

    def test_molecular_weight_has_name_and_compute(self):
        """MolecularWeightComponent should have name attribute and compute method"""
        component = MolecularWeightComponent()
        assert hasattr(component, 'name')
        assert callable(getattr(component, 'compute'))


class TestComponentRegistry:
    """Test ComponentRegistry registration, retrieval, and validation"""

    def test_default_registry_has_wiener(self):
        """Default registry should have wiener_index component"""
        registry = get_registry()
        component = registry.get("wiener_index")
        assert component is not None

    def test_default_registry_has_molecular_weight(self):
        """Default registry should have molecular_weight component"""
        registry = get_registry()
        component = registry.get("molecular_weight")
        assert component is not None

    def test_list_components(self):
        """list_components() should return both component names"""
        registry = get_registry()
        components = registry.list_components()
        assert "wiener_index" in components
        assert "molecular_weight" in components

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
            assert "molecular_weight" in error_msg

    def test_duplicate_registration_raises(self):
        """Registering duplicate component name should raise ValueError"""
        registry = ComponentRegistry()
        component1 = WienerIndexComponent()
        component2 = WienerIndexComponent()

        registry.register(component1)
        with pytest.raises(ValueError):
            registry.register(component2)
