import pytest

from core.universe_simulator import tau_from_redshift
from island.islands_module.entropy_model import (
    IslandToyParameters,
    evaporation_redshift_from_time,
    no_island_entropy,
    run_island_scan,
    tau_macro_from_time,
)


def test_island_scan_has_page_transition():
    sim = run_island_scan(IslandToyParameters(n_time=120, n_depth=120))

    assert sim.page_time is not None
    assert sim.page_time == pytest.approx(0.0968, abs=0.01)

    page_row = next(row for row in sim.rows if row.island_active)
    assert 0.15 < page_row.island_depth < 0.30
    assert page_row.island_margin > 0.0

    late_row = sim.rows[-1]
    assert late_row.island_active
    assert late_row.island_depth > 0.90


def test_hawking_branch_uses_native_area_law():
    params = IslandToyParameters(t_final=1.0)

    assert no_island_entropy(0.0, params) == pytest.approx(0.0)
    assert no_island_entropy(1.0, params) == pytest.approx(params.hawking_rate)
    assert no_island_entropy(0.5, params) == pytest.approx(params.hawking_rate * 0.75)


def test_macro_tau_uses_log_horizon_map():
    params = IslandToyParameters(t_final=1.0)
    t = 0.25

    z = evaporation_redshift_from_time(t, params)
    assert tau_macro_from_time(t, params) == pytest.approx(tau_from_redshift(z))
