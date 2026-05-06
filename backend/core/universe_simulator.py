from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import json
import math
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Axioms:
    pi: float = math.pi
    e: float = math.e
    i: complex = 1j
    dim_boundary: int = 2
    dim_volume: int = 3
    tau_macro_fixed: float = 1.0
    tau_micro_fixed: float = -1.0


@dataclass(frozen=True)
class PhysicalInputs:
    c: float = 299_792_458.0
    h: float = 6.626_070_15e-34
    hbar: float = 1.054_571_817e-34
    k_B: float = 1.380_649e-23
    alpha: float = 1 / 137.035_999_084
    m_e_kg: float = 9.109_383_7139e-31
    m_p_MeV: float = 938.272_088_16
    m_n_MeV: float = 939.565_420_52
    T_rec_K: float = 3000.0
    z_rec: float = 1100.0
    T_CMB_observed_K: float = 2.7255
    m_W_observed_GeV: float = 80.377
    m_Z_observed_GeV: float = 91.1876
    m_s_PDG_central_MeV: float = 95.0
    m_s_PDG_min_MeV: float = 83.0
    m_s_PDG_max_MeV: float = 117.0
    Omega_Lambda_PDG_2024: float = 0.689
    Omega_DM_Planck_derived: float = 0.2626
    Casimir_parallel_ratio_reference: float | None = None
    Boyer_spherical_Casimir_reference: float = 0.09235
    R_star_reference_fm: float = 2.197
    H0_SI: float | None = None
    rho_Lambda_observed_SI: float | None = None


@dataclass(frozen=True)
class DerivedConstants:
    C0_axis: float
    C1: float
    Omega_Lambda_fp: float
    Omega_DM_base: float
    Omega_DM_fp: float
    f_sph: float
    C_sphere: float
    F0_D2_to_D3: float
    F_D3_to_D2: float
    lambda_capacity: float
    kappa_geom: float
    Schwinger_ratio: float


@dataclass(frozen=True)
class CosmicEpoch:
    tau: float
    z: float
    scale_factor: float
    temperature_K: float
    R_dimensionless: float
    l_UV_dimensionless: float
    photon_energy_J: float
    photon_energy_eV: float


@dataclass(frozen=True)
class BranchProjection:
    branch_name: str
    epoch_tau: float
    tau_R: float
    tau_I: float
    tau_micro: float | None
    phase_lead: float | None
    C_value: complex
    C_abs: float
    C_ratio_to_C1: float
    Z_coherence: float
    quadrant: str
    zone: str
    physical_regime: str
    event_label: str | None = None


@dataclass(frozen=True)
class ProjectedState:
    epoch: CosmicEpoch
    projection: BranchProjection


@dataclass(frozen=True)
class SimulationEvent:
    name: str
    tau: float
    z: float | None
    temperature_K: float | None
    description: str
    first_principles_status: str


@dataclass(frozen=True)
class ParticleEvent:
    name: str
    branch_name: str
    epoch_tau: float
    trigger: str
    rest_mass_eV: float | None
    rest_mass_source: str
    thermal_kBT_eV: float
    kinetic_energy_proxy_eV: float
    kinetic_energy_proxy_source: str
    temperature_K: float
    z: float
    C_ratio_to_C1: float
    Z_coherence: float
    quadrant: str
    zone: str
    first_principles_status: str


@dataclass(frozen=True)
class QuarkMassEstimate:
    flavor: str
    model_mass_MeV: float
    model_mass_eV: float
    model_definition: str
    reference_mass_MeV: float | None
    reference_definition: str
    difference_MeV: float | None
    percent_error: float | None
    status: str


@dataclass(frozen=True)
class ForceValidation:
    name: str
    sector: str
    model_value: float
    reference_value: float | None
    unit: str
    difference: float | None
    percent_error: float | None
    formula: str
    reference: str
    status: str


@dataclass(frozen=True)
class DMRelaxationPoint:
    label: str
    z: float
    tau_desitter: float
    tau_process: float
    C_ratio_to_C1: float
    delta_C_ratio_to_C1: float
    eta_suppression_kernel: float
    f_sigma8_suppression_fraction: float
    omega_dm_base: float
    omega_dm_bias: float
    omega_dm_observer_inferred: float
    process_note: str


@dataclass(frozen=True)
class QCaptureEstimate:
    particle: str
    tau_event: float
    rest_mass_eV: float
    lP_eff_m: float
    E_UV_eV: float
    q_capture: float
    minus_log_Q: float
    nearest_alpha_power: int
    alpha_power_value: float
    ratio_to_alpha_power: float
    nearest_C1_power: int
    C1_power_value: float
    ratio_to_C1_power: float
    status: str


@dataclass(frozen=True)
class AlphaMassRelation:
    name: str
    electron_mass_MeV: float
    light_quark_effective_mass_MeV: float
    alpha_from_mass_ratio: float
    alpha_reference: float
    m_s_reference_MeV: float
    m_s_percent_error: float
    within_reference_range: bool
    formula: str
    status: str


@dataclass(frozen=True)
class EffectiveParticleMass:
    sector: str
    event_time_definition: str
    tau_event: float
    C_ratio_to_C1: float
    E_UV_eV: float | None
    effective_mass_MeV: float | None
    reference_mass_MeV: float | None
    percent_error: float | None
    q_capture: float | None
    minus_log_Q: float | None
    formula: str
    parameter_layer: str
    status: str


def derive_C0_axis(pi: float = math.pi) -> float:
    return pi**2 / 16


def derive_C1(pi: float = math.pi, e: float = math.e) -> float:
    return pi**2 / (16 * e**2)


def derive_omega_lambda_fp(pi: float = math.pi, e: float = math.e) -> float:
    return pi**3 / (6 * e**2)


def derive_f_sph(pi: float = math.pi, e: float = math.e) -> float:
    return math.sqrt(4 * pi) / e


def derive_projection_factors(pi: float = math.pi, e: float = math.e) -> tuple[float, float]:
    return pi / e, e / pi


def derive_schwinger_ratio(e: float = math.e, alpha: float = PhysicalInputs().alpha) -> float:
    return e ** (2 * (1 - alpha))


def derive_omega_dm(pi: float = math.pi, e: float = math.e) -> float:
    return pi**4 * (math.sqrt(4 * pi) - e) / (6 * e**4)


def derive_omega_dm_fp(pi: float = math.pi, e: float = math.e) -> float:
    return derive_kappa_geom(pi, e) * derive_omega_dm(pi, e)


def derive_lambda_capacity(pi: float = math.pi, e: float = math.e) -> float:
    return 1 / (derive_f_sph(pi, e) * e**2)


def derive_kappa_geom(pi: float = math.pi, e: float = math.e) -> float:
    lam = derive_lambda_capacity(pi, e)
    return 1 / math.sqrt(1 - lam)


def derive_constants(axioms: Axioms, inputs: PhysicalInputs) -> DerivedConstants:
    c0_axis = derive_C0_axis(axioms.pi)
    c1 = derive_C1(axioms.pi, axioms.e)
    f_sph = derive_f_sph(axioms.pi, axioms.e)
    f_d2_d3, f_d3_d2 = derive_projection_factors(axioms.pi, axioms.e)
    return DerivedConstants(
        C0_axis=c0_axis,
        C1=c1,
        Omega_Lambda_fp=derive_omega_lambda_fp(axioms.pi, axioms.e),
        Omega_DM_base=derive_omega_dm(axioms.pi, axioms.e),
        Omega_DM_fp=derive_omega_dm_fp(axioms.pi, axioms.e),
        f_sph=f_sph,
        C_sphere=c1 * f_sph,
        F0_D2_to_D3=f_d2_d3,
        F_D3_to_D2=f_d3_d2,
        lambda_capacity=derive_lambda_capacity(axioms.pi, axioms.e),
        kappa_geom=derive_kappa_geom(axioms.pi, axioms.e),
        Schwinger_ratio=derive_schwinger_ratio(axioms.e, inputs.alpha),
    )


def tau_from_redshift(z: float) -> float:
    if z <= -1:
        raise ValueError("redshift must be greater than -1")
    return 1 - math.log1p(z)


def redshift_from_tau(tau: float) -> float:
    return math.exp(1 - tau) - 1


def scale_factor_from_tau(tau: float) -> float:
    return math.exp(tau - 1)


def e_weighted_micro_tau(tau_macro: float, axioms: Axioms = Axioms()) -> float:
    delta_macro = axioms.tau_macro_fixed - tau_macro
    return axioms.tau_micro_fixed + delta_macro / axioms.e


def phase_lead_from_tau(tau_macro: float, axioms: Axioms = Axioms()) -> float:
    delta_macro = axioms.tau_macro_fixed - tau_macro
    delta_micro = delta_macro / axioms.e
    return 2 * (delta_macro - delta_micro)


def quadrant_of(tau_R: float, tau_I: float) -> str:
    if tau_R == 0 or tau_I == 0:
        return "axis / singular cut"
    if tau_R > 0 and tau_I > 0:
        return "QI stable matter"
    if tau_R < 0 and tau_I > 0:
        return "QII Z0 saddle"
    if tau_R < 0 and tau_I < 0:
        return "QIII W- saddle"
    return "QIV radiation background"


def C_base(R: float, l_UV: float, pi: float = math.pi) -> float:
    if R <= 0 or l_UV <= 0:
        raise ValueError("R and l_UV must be positive")
    return (pi**2 / 16) * (l_UV / R) ** 2


def C_complex(R: float, tau_R: float, tau_I: float, l_UV: float, pi: float = math.pi) -> complex:
    phase = complex(math.cos(-2 * tau_I), math.sin(-2 * tau_I))
    return C_base(R, l_UV, pi) * math.exp(-2 * tau_R) * phase


def C_conjugate_partner(R: float, tau_R: float, tau_I: float, l_UV: float) -> complex:
    return C_complex(R, -tau_R, tau_I, l_UV)


def conservation_product(R: float, tau_R: float, tau_I: float, l_UV: float) -> complex:
    return C_complex(R, tau_R, tau_I, l_UV) * C_complex(R, -tau_R, -tau_I, l_UV)


def coherence_Z(R: float, tau_R: float, tau_I: float, l_UV: float) -> float:
    return abs(C_complex(R, tau_R, tau_I, l_UV) + C_conjugate_partner(R, tau_R, tau_I, l_UV)) ** 2


def temperature_from_tau(T_initial: float, tau_initial: float, tau: float) -> float:
    return T_initial * math.exp(-(tau - tau_initial))


def T_CMB_from_recombination(T_rec: float, z_rec: float) -> float:
    tau_rec = tau_from_redshift(z_rec)
    return temperature_from_tau(T_rec, tau_rec, 1.0)


def rho_lambda_fp_from_C1(C1: float, hbar: float, c: float, l_P: float, R_H: float) -> float:
    return C1 * hbar * c / (l_P**2 * R_H**2)


def G_from_holographic_identity_energy_density(C1: float, H0: float, rho_lambda_energy: float, c: float) -> float:
    return C1 * c**2 * H0**2 / rho_lambda_energy


def effective_lP_from_tau(lP_ref: float, tau: float, tau_ref: float = 1.0) -> float:
    return lP_ref * math.exp(-(tau - tau_ref))


def causal_scale_from_tau(R_ref: float, tau: float, tau_ref: float = 1.0) -> float:
    return R_ref * math.exp(tau - tau_ref)


def photon_energy_from_R(R: float, inputs: PhysicalInputs) -> float:
    return inputs.h * inputs.c / (4 * R)


def default_horizon_radius_m(inputs: PhysicalInputs) -> float:
    # Temporary SI bridge used only for photon-energy tracking in the MVP.
    h0_s_inv = 70_000 / 3.085_677_581_491_367e22
    return inputs.c / h0_s_inv


def omega_r_over_omega_m_from_QIV(qiv_peak_ratio: float, z_rec: float) -> float:
    return qiv_peak_ratio / (1 + z_rec)


def wz_ratio_from_stability_inversion(tau_W: float = 0.9468, tau_Z: float = 0.9125) -> float:
    return math.exp(-4 * (tau_W - tau_Z))


def joule_to_eV(energy_j: float) -> float:
    return energy_j / 1.602_176_634e-19


def kg_to_eV_mass(mass_kg: float, inputs: PhysicalInputs) -> float:
    return joule_to_eV(mass_kg * inputs.c**2)


def thermal_kBT_eV(temperature_K: float, inputs: PhysicalInputs) -> float:
    return joule_to_eV(inputs.k_B * temperature_K)


def thermal_kinetic_proxy_eV(temperature_K: float, inputs: PhysicalInputs, particle_kind: str) -> float:
    kBT = thermal_kBT_eV(temperature_K, inputs)
    if particle_kind == "photon":
        return 2.701 * kBT
    return 1.5 * kBT


def classify_zone(c_ratio: float, schwinger_ratio: float, tol: float = 1e-9) -> str:
    if c_ratio < 1 - tol:
        return "under-saturated virtual"
    if abs(c_ratio - 1) <= tol:
        return "boundary-critical"
    if c_ratio < 4 - tol:
        return "stable storage"
    if abs(c_ratio - 4) <= tol:
        return "UV cutoff boundary"
    if c_ratio < schwinger_ratio - tol:
        return "forced shared-coding"
    if abs(c_ratio - schwinger_ratio) <= tol:
        return "Schwinger / Zone-I boundary"
    return "forced entanglement / confinement"


def build_epoch(
    tau: float,
    inputs: PhysicalInputs,
    tau_rec: float,
    R_ref_dimensionless: float = 1.0,
    l_UV_ref_dimensionless: float = 1.0,
) -> CosmicEpoch:
    photon_R_m = causal_scale_from_tau(default_horizon_radius_m(inputs), tau)
    photon_energy_j = photon_energy_from_R(photon_R_m, inputs)
    return CosmicEpoch(
        tau=tau,
        z=redshift_from_tau(tau),
        scale_factor=scale_factor_from_tau(tau),
        temperature_K=temperature_from_tau(inputs.T_rec_K, tau_rec, tau),
        R_dimensionless=causal_scale_from_tau(R_ref_dimensionless, tau),
        l_UV_dimensionless=effective_lP_from_tau(l_UV_ref_dimensionless, tau),
        photon_energy_J=photon_energy_j,
        photon_energy_eV=photon_energy_j / 1.602_176_634e-19,
    )


def branch_coordinates(epoch_tau: float, branch_name: str, axioms: Axioms) -> tuple[float, float, str]:
    """Map a cosmological epoch onto a complex-tau physics branch."""
    if branch_name == "micro_mirror":
        tau_r = epoch_tau
        tau_i = e_weighted_micro_tau(tau_r, axioms)
        return tau_r, tau_i, "canonical e-weighted macro/micro mirror"
    if branch_name == "radiation_qiv":
        tau_r = max(abs(epoch_tau), 1e-6)
        tau_i = -abs(e_weighted_micro_tau(tau_r, axioms))
        return tau_r, tau_i, "QIV cosmological radiation background projection"
    if branch_name == "matter_qi":
        tau_r = max(abs(epoch_tau), 1e-6)
        tau_i = abs(e_weighted_micro_tau(tau_r, axioms))
        return tau_r, tau_i, "QI stable matter / atomic binding projection"
    if branch_name == "weak_z_qii":
        return -0.9125, 0.9125, "QII Z0 diagonal saddle reference"
    if branch_name == "weak_w_qiii":
        return -0.9468, -1e-12, "QIII W- charged saddle reference"
    raise ValueError(f"unknown branch: {branch_name}")


def project_epoch(
    epoch: CosmicEpoch,
    branch_name: str,
    axioms: Axioms,
    constants: DerivedConstants,
    R_ref_dimensionless: float = 1.0,
    l_UV_ref_dimensionless: float = 1.0,
) -> BranchProjection:
    tau_r, tau_i, regime = branch_coordinates(epoch.tau, branch_name, axioms)
    C = C_complex(R_ref_dimensionless, tau_r, tau_i, l_UV_ref_dimensionless, axioms.pi)
    c_ratio = abs(C) / constants.C1
    tau_micro = e_weighted_micro_tau(tau_r, axioms) if branch_name == "micro_mirror" else None
    phase_lead = phase_lead_from_tau(tau_r, axioms) if branch_name == "micro_mirror" else None
    return BranchProjection(
        branch_name=branch_name,
        epoch_tau=epoch.tau,
        tau_R=tau_r,
        tau_I=tau_i,
        tau_micro=tau_micro,
        phase_lead=phase_lead,
        C_value=C,
        C_abs=abs(C),
        C_ratio_to_C1=c_ratio,
        Z_coherence=coherence_Z(R_ref_dimensionless, tau_r, tau_i, l_UV_ref_dimensionless),
        quadrant=quadrant_of(tau_r, tau_i),
        zone=classify_zone(c_ratio, constants.Schwinger_ratio),
        physical_regime=regime,
    )


def event_from_tau(
    name: str,
    tau: float,
    inputs: PhysicalInputs,
    description: str,
    first_principles_status: str,
) -> SimulationEvent:
    return SimulationEvent(
        name=name,
        tau=tau,
        z=redshift_from_tau(tau),
        temperature_K=temperature_from_tau(inputs.T_rec_K, tau_from_redshift(inputs.z_rec), tau),
        description=description,
        first_principles_status=first_principles_status,
    )


def build_events(inputs: PhysicalInputs, constants: DerivedConstants) -> list[SimulationEvent]:
    tau_rec = tau_from_redshift(inputs.z_rec)
    tau_first_qiv = 1e-6
    return [
        event_from_tau(
            "origin / early boundary state",
            -10.0,
            inputs,
            "Start of the finite MVP tau grid, not a resolved singular origin.",
            "provisional grid boundary; first-principles singular-origin event is not yet closed",
        ),
        event_from_tau(
            "first photon event",
            tau_first_qiv,
            inputs,
            "Earliest positive-tau QIV radiation-background marker with photon energy tracked by E=h c/(4R).",
            "provisional event rule; photon-formation threshold is not yet derived from first principles",
        ),
        event_from_tau(
            "first quark-level encoding",
            -10.0,
            inputs,
            "Earliest sampled state already above the Zone-I threshold C/C1 = exp(2(1-alpha)).",
            "provisional grid-boundary event; full quark event requires an internal confinement/node criterion",
        ),
        event_from_tau(
            "recombination",
            tau_rec,
            inputs,
            "Hydrogen recombination bridge event using z_rec and T_rec.",
            "temporary external atomic/cosmological bridge input",
        ),
        event_from_tau(
            "first hydrogen atom",
            tau_rec,
            inputs,
            "First stable hydrogen event currently coincident with recombination in the MVP.",
            "temporary external atomic bridge input; Saha/binding trigger not yet internally closed",
        ),
        event_from_tau(
            "today / de Sitter reference",
            1.0,
            inputs,
            "Reference state tau=1.",
            "model reference fixed point",
        ),
    ]


def epoch_timeline(tau_start: float, tau_end: float, steps: int, inputs: PhysicalInputs) -> list[CosmicEpoch]:
    tau_rec = tau_from_redshift(inputs.z_rec)
    return [
        build_epoch(tau_start + (tau_end - tau_start) * i / (steps - 1), inputs, tau_rec)
        for i in range(steps)
    ]


def project_timeline(
    epochs: list[CosmicEpoch],
    branches: list[str],
    axioms: Axioms,
    constants: DerivedConstants,
) -> list[ProjectedState]:
    return [
        ProjectedState(epoch=epoch, projection=project_epoch(epoch, branch, axioms, constants))
        for epoch in epochs
        for branch in branches
    ]


def json_ready(value: Any) -> Any:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if hasattr(value, "__dataclass_fields__"):
        return {key: json_ready(item) for key, item in asdict(value).items()}
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    return value


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(json_ready(payload), indent=2), encoding="utf-8")


def write_timeline_csv(path: Path, rows: list[ProjectedState]) -> None:
    fields = [
        "epoch_tau",
        "branch_name",
        "tau_R",
        "tau_I",
        "tau_micro",
        "phase_lead",
        "z",
        "scale_factor",
        "temperature_K",
        "R_dimensionless",
        "l_UV_dimensionless",
        "photon_energy_J",
        "photon_energy_eV",
        "C_abs",
        "C_ratio_to_C1",
        "Z_coherence",
        "quadrant",
        "zone",
        "physical_regime",
        "event_label",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            epoch = row.epoch
            projection = row.projection
            writer.writerow(
                {
                    "epoch_tau": epoch.tau,
                    "branch_name": projection.branch_name,
                    "tau_R": projection.tau_R,
                    "tau_I": projection.tau_I,
                    "tau_micro": projection.tau_micro,
                    "phase_lead": projection.phase_lead,
                    "z": epoch.z,
                    "scale_factor": epoch.scale_factor,
                    "temperature_K": epoch.temperature_K,
                    "R_dimensionless": epoch.R_dimensionless,
                    "l_UV_dimensionless": epoch.l_UV_dimensionless,
                    "photon_energy_J": epoch.photon_energy_J,
                    "photon_energy_eV": epoch.photon_energy_eV,
                    "C_abs": projection.C_abs,
                    "C_ratio_to_C1": projection.C_ratio_to_C1,
                    "Z_coherence": projection.Z_coherence,
                    "quadrant": projection.quadrant,
                    "zone": projection.zone,
                    "physical_regime": projection.physical_regime,
                    "event_label": projection.event_label,
                }
            )


def write_photon_spectrum_csv(path: Path, rows: list[CosmicEpoch]) -> None:
    fields = ["tau", "z", "scale_factor", "R_dimensionless", "photon_energy_J", "photon_energy_eV", "temperature_K"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fields})


def write_event_table(path: Path, events: list[SimulationEvent]) -> None:
    lines = [
        "# Universe Simulator Event Table",
        "",
        "| Event | tau | z | temperature_K | first-principles status | Description |",
        "|---|---:|---:|---:|---|---|",
    ]
    for event in events:
        z = "" if event.z is None else f"{event.z:.8g}"
        temp = "" if event.temperature_K is None else f"{event.temperature_K:.8g}"
        lines.append(
            f"| {event.name} | {event.tau:.8g} | {z} | {temp} | "
            f"{event.first_principles_status} | {event.description} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def event_payload(
    event: SimulationEvent,
    branch_name: str,
    inputs: PhysicalInputs,
    axioms: Axioms,
    constants: DerivedConstants,
) -> dict[str, Any]:
    epoch = build_epoch(event.tau, inputs, tau_from_redshift(inputs.z_rec))
    projection = project_epoch(epoch, branch_name, axioms, constants)
    return {
        "event": event,
        "selected_branch": branch_name,
        "epoch": epoch,
        "projection": projection,
    }


def build_particle_event(
    name: str,
    branch_name: str,
    epoch_tau: float,
    trigger: str,
    rest_mass_eV: float | None,
    rest_mass_source: str,
    particle_kind: str,
    first_principles_status: str,
    inputs: PhysicalInputs,
    axioms: Axioms,
    constants: DerivedConstants,
) -> ParticleEvent:
    epoch = build_epoch(epoch_tau, inputs, tau_from_redshift(inputs.z_rec))
    projection = project_epoch(epoch, branch_name, axioms, constants)
    return ParticleEvent(
        name=name,
        branch_name=branch_name,
        epoch_tau=epoch_tau,
        trigger=trigger,
        rest_mass_eV=rest_mass_eV,
        rest_mass_source=rest_mass_source,
        thermal_kBT_eV=thermal_kBT_eV(epoch.temperature_K, inputs),
        kinetic_energy_proxy_eV=thermal_kinetic_proxy_eV(epoch.temperature_K, inputs, particle_kind),
        kinetic_energy_proxy_source=(
            "2.701 k_B T photon thermal proxy"
            if particle_kind == "photon"
            else "1.5 k_B T nonrelativistic thermal proxy"
        ),
        temperature_K=epoch.temperature_K,
        z=epoch.z,
        C_ratio_to_C1=projection.C_ratio_to_C1,
        Z_coherence=projection.Z_coherence,
        quadrant=projection.quadrant,
        zone=projection.zone,
        first_principles_status=first_principles_status,
    )


def build_particle_events(inputs: PhysicalInputs, axioms: Axioms, constants: DerivedConstants) -> list[ParticleEvent]:
    tau_schwinger = 1 - 0.5 * math.log(constants.Schwinger_ratio)
    tau_sphere = 1 - 0.5 * math.log(constants.f_sph)
    tau_stable_storage_upper = 1 - 0.5 * math.log(4)
    tau_z = 0.9125
    tau_w = 0.9468
    m_z_eV = inputs.m_Z_observed_GeV * 1e9
    m_w_predicted_eV = wz_ratio_from_stability_inversion(tau_w, tau_z) * m_z_eV
    m_electron_eV = kg_to_eV_mass(inputs.m_e_kg, inputs)
    return [
        build_particle_event(
            name="QIV radiation / CMB precursor",
            branch_name="radiation_qiv",
            epoch_tau=-10.0,
            trigger="finite-grid earliest QIV radiation-background sample",
            rest_mass_eV=0.0,
            rest_mass_source="photon rest mass set to zero",
            particle_kind="photon",
            first_principles_status="provisional grid boundary; true first-photon threshold not closed",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
        build_particle_event(
            name="QIV Schwinger saddle marker",
            branch_name="radiation_qiv",
            epoch_tau=tau_schwinger,
            trigger="C/C1 = exp(2(1-alpha)); QIV reaches Schwinger / Zone-I boundary",
            rest_mass_eV=0.0,
            rest_mass_source="photon rest mass set to zero",
            particle_kind="photon",
            first_principles_status="model-endogenous threshold; interpretation as CMB maximum/saddle remains provisional",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
        build_particle_event(
            name="electron charged-coherence seed",
            branch_name="matter_qi",
            epoch_tau=tau_sphere,
            trigger="C(tau)=C_sphere, equivalently C/C1=f_sph; coherent charged D2-D3 projection begins",
            rest_mass_eV=m_electron_eV,
            rest_mass_source="external SI electron mass bridge converted to eV; only absolute micro-mass scale input",
            particle_kind="massive",
            first_principles_status="coherent charged precursor is model-endogenous; electron absolute mass remains Layer-3 bridge until QI island integral closes",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
        build_particle_event(
            name="first quark confinement-antinode candidate",
            branch_name="matter_qi",
            epoch_tau=tau_schwinger,
            trigger="tau=alpha, equivalently C/C1=exp(2(1-alpha)); stable quark antinode encoding threshold",
            rest_mass_eV=None,
            rest_mass_source="effective d/s masses are generated in effective_particle_masses; u remains unclosed",
            particle_kind="massive",
            first_principles_status="P4 node/antinode rule; d and s ratios closed relative to m_e bridge, u sector not yet closed",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
        build_particle_event(
            name="QI coherent plane-wave critical",
            branch_name="matter_qi",
            epoch_tau=tau_stable_storage_upper,
            trigger="C/C1 = 4 stable-storage upper boundary",
            rest_mass_eV=None,
            rest_mass_source="not a particle rest-mass event",
            particle_kind="massive",
            first_principles_status="model-endogenous zone crossing; standing-wave trigger not closed",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
        build_particle_event(
            name="Z0 saddle candidate",
            branch_name="weak_z_qii",
            epoch_tau=tau_schwinger,
            trigger="post-QIV Schwinger transition into fixed QII saddle coordinates tau_R=-0.9125, tau_I=+0.9125",
            rest_mass_eV=m_z_eV,
            rest_mass_source="external observed Z0 mass; used as absolute anchor for current table",
            particle_kind="massive",
            first_principles_status="event epoch tied to tau=alpha transition; saddle location model-specified; absolute mass not yet derived internally",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
        build_particle_event(
            name="W- saddle candidate",
            branch_name="weak_w_qiii",
            epoch_tau=tau_schwinger,
            trigger="post-QIV Schwinger transition into fixed QIII W- saddle coordinate tau_R=-0.9468",
            rest_mass_eV=m_w_predicted_eV,
            rest_mass_source="predicted from P7 W/Z ratio using observed Z0 mass as anchor",
            particle_kind="massive",
            first_principles_status="event epoch tied to tau=alpha transition; ratio is model-endogenous; absolute mass uses external Z0 anchor",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
        build_particle_event(
            name="W+ conjugate bridge candidate",
            branch_name="matter_qi",
            epoch_tau=tau_schwinger,
            trigger="charged weak bridge paired with W- after QIV Schwinger transition; QI-side conjugate projection in MVP",
            rest_mass_eV=m_w_predicted_eV,
            rest_mass_source="same predicted W mass as W- using observed Z0 anchor",
            particle_kind="massive",
            first_principles_status="event epoch tied to tau=alpha transition; charge-conjugate branch assignment is provisional",
            inputs=inputs,
            axioms=axioms,
            constants=constants,
        ),
    ]


def write_particle_events_csv(path: Path, events: list[ParticleEvent]) -> None:
    fields = list(asdict(events[0]).keys()) if events else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for event in events:
            writer.writerow(asdict(event))


def write_particle_events_md(path: Path, events: list[ParticleEvent]) -> None:
    lines = [
        "# Particle Event Candidate Table",
        "",
        "| Event | branch | tau | z | T_K | rest_mass_eV | kinetic_proxy_eV | C/C1 | quadrant | status |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for event in events:
        rest = "" if event.rest_mass_eV is None else f"{event.rest_mass_eV:.8g}"
        lines.append(
            f"| {event.name} | {event.branch_name} | {event.epoch_tau:.8g} | {event.z:.8g} | "
            f"{event.temperature_K:.8g} | {rest} | {event.kinetic_energy_proxy_eV:.8g} | "
            f"{event.C_ratio_to_C1:.8g} | {event.quadrant} | {event.first_principles_status} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_quark_mass_estimates(inputs: PhysicalInputs) -> list[QuarkMassEstimate]:
    """Estimate first-generation constituent quark masses from nucleon closure.

    This is a confined/coherence-node mass estimate, not an MSbar current-quark
    mass derivation. The simple closure equations are:
        proton ~= 2 u + d
        neutron ~= u + 2 d
    """
    u_mev = (2 * inputs.m_p_MeV - inputs.m_n_MeV) / 3
    d_mev = (2 * inputs.m_n_MeV - inputs.m_p_MeV) / 3
    q_avg_mev = ((inputs.m_p_MeV + inputs.m_n_MeV) / 2) / 3
    references = {
        "u": (2.16, "PDG 2024 MSbar current mass at 2 GeV"),
        "d": (4.70, "PDG 2024 MSbar current mass at 2 GeV"),
        "q_avg": (3.49, "PDG 2024 average current mass (m_u + m_d)/2 at 2 GeV"),
    }
    estimates = [
        ("u", u_mev, "constituent/coherence-node mass from 2u+d=proton, u+2d=neutron"),
        ("d", d_mev, "constituent/coherence-node mass from 2u+d=proton, u+2d=neutron"),
        ("q_avg", q_avg_mev, "average first-generation constituent/coherence-node mass m_N_avg/3"),
    ]
    rows = []
    for flavor, model_mev, model_definition in estimates:
        ref_mev, ref_definition = references[flavor]
        diff = model_mev - ref_mev
        rows.append(
            QuarkMassEstimate(
                flavor=flavor,
                model_mass_MeV=model_mev,
                model_mass_eV=model_mev * 1e6,
                model_definition=model_definition,
                reference_mass_MeV=ref_mev,
                reference_definition=ref_definition,
                difference_MeV=diff,
                percent_error=diff / ref_mev * 100,
                status="not directly comparable: constituent/confined mass vs MSbar current-quark mass",
            )
        )
    return rows


def write_quark_mass_estimates_csv(path: Path, estimates: list[QuarkMassEstimate]) -> None:
    fields = list(asdict(estimates[0]).keys()) if estimates else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for estimate in estimates:
            writer.writerow(asdict(estimate))


def write_quark_mass_estimates_md(path: Path, estimates: list[QuarkMassEstimate]) -> None:
    lines = [
        "# Quark Mass Estimate Table",
        "",
        "| flavor | model_mass_MeV | reference_mass_MeV | percent_error | status |",
        "|---|---:|---:|---:|---|",
    ]
    for estimate in estimates:
        ref = "" if estimate.reference_mass_MeV is None else f"{estimate.reference_mass_MeV:.8g}"
        pct = "" if estimate.percent_error is None else f"{estimate.percent_error:.8g}"
        lines.append(
            f"| {estimate.flavor} | {estimate.model_mass_MeV:.8g} | {ref} | {pct} | {estimate.status} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validation_row(
    name: str,
    sector: str,
    model_value: float,
    reference_value: float | None,
    unit: str,
    formula: str,
    reference: str,
    status: str,
) -> ForceValidation:
    difference = None if reference_value is None else model_value - reference_value
    percent_error = None if reference_value in (None, 0) else difference / reference_value * 100
    return ForceValidation(
        name=name,
        sector=sector,
        model_value=model_value,
        reference_value=reference_value,
        unit=unit,
        difference=difference,
        percent_error=percent_error,
        formula=formula,
        reference=reference,
        status=status,
    )


def build_force_validations(inputs: PhysicalInputs, axioms: Axioms, constants: DerivedConstants) -> list[ForceValidation]:
    re_m = inputs.alpha * inputs.hbar / (inputs.m_e_kg * inputs.c)
    re_fm = re_m / 1e-15
    r_star_fm = re_fm * (axioms.pi / 4) * math.exp(-inputs.alpha)
    c_parallel_ratio = math.e**2 / 15
    c_parallel_direct_ratio = (axioms.pi**2 / 240) / constants.C1
    qcd_tau_s = 1 - math.log(2)
    qcd_fixed_point = math.exp(-2 * (qcd_tau_s - 1))
    return [
        validation_row(
            "Omega_Lambda_fp",
            "gravity / dark energy",
            constants.Omega_Lambda_fp,
            inputs.Omega_Lambda_PDG_2024,
            "dimensionless",
            "pi^3/(6 e^2)",
            "PDG 2024 cosmological parameters, Planck+BAO observed rho_Lambda density fraction ~= 0.689 +/- 0.006",
            "fixed-point dark-energy density fraction; rho_Lambda observed is validation-only",
        ),
        validation_row(
            "Omega_DM_fp",
            "dark matter / D2-to-D3 overflow",
            constants.Omega_DM_fp,
            inputs.Omega_DM_Planck_derived,
            "dimensionless",
            "Omega_DM_fp = kappa*pi^4*(sqrt(4*pi)-e)/(6 e^4)",
            "derived from PDG 2024 Planck+BAO Omega_m and Omega_b h^2/h^2",
            "kappa-corrected fixed-point dark-matter density; observed value is validation-only",
        ),
        validation_row(
            "Casimir parallel-plate ratio",
            "static EM / Casimir",
            c_parallel_ratio,
            c_parallel_direct_ratio,
            "dimensionless",
            "C_parallel/C1 = e^2/15 = (pi^2/240)/(pi^2/(16e^2))",
            "exact coefficient identity from Paper III / overview",
            "internal exact coefficient check",
        ),
        validation_row(
            "spherical Casimir coefficient",
            "static EM / spherical Casimir",
            constants.C_sphere,
            inputs.Boyer_spherical_Casimir_reference,
            "dimensionless",
            "C_sphere = C1 * f_sph",
            "Boyer ideal conducting spherical shell coefficient ~= 0.09235",
            "roughly comparable but geometry/sign conventions require care",
        ),
        validation_row(
            "entanglement critical scale R_star",
            "EM / entanglement",
            r_star_fm,
            inputs.R_star_reference_fm,
            "fm",
            "R_star = r_e*(pi/4)*exp(-alpha)",
            "Paper VI / overview reference value ~= 2.197 fm; future EIC target, not yet experimental",
            "formula closure check; not yet experimental validation",
        ),
        validation_row(
            "electron classical radius closure",
            "EM",
            re_fm,
            2.817_940_3205,
            "fm",
            "r_e = alpha*hbar/(m_e*c)",
            "CODATA 2022 classical electron radius",
            "derived from external constants; closure check, not independent prediction",
        ),
        validation_row(
            "QCD fixed-point coupling",
            "strong",
            qcd_fixed_point,
            4.0,
            "dimensionless",
            "C_QCD/C1 = exp[-2*(tau_s - 1)], tau_s = 1 - ln(2)",
            "model closure: QCD fixed point identified with stable-storage upper boundary C/C1=4",
            "closed model identity; experimental comparison requires mapping to a running alpha_s convention",
        ),
        validation_row(
            "kappa geometric compression",
            "dark matter / D2-to-D3 overflow",
            constants.kappa_geom,
            inputs.Omega_DM_Planck_derived / constants.Omega_DM_base,
            "dimensionless",
            "kappa = 1/sqrt(1 - lambda), lambda = 1/(f_sph e^2)",
            "P4 derived kappa vs Planck-inferred diagnostic kappa_obs=Omega_DM_obs/Omega_DM_base",
            "model-endogenous compression factor; observed kappa is validation-only",
        ),
        validation_row(
            "QCD fixed-point tau_s",
            "strong",
            qcd_tau_s,
            0.3,
            "dimensionless tau",
            "tau_s = 1 - ln(2)",
            "overview structural alpha_s ~= 0.3 reference",
            "closed model value; comparison to alpha_s is interpretive, not a direct scheme-independent observable",
        ),
    ]


def write_force_validations_csv(path: Path, rows: list[ForceValidation]) -> None:
    fields = list(asdict(rows[0]).keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_force_validations_md(path: Path, rows: list[ForceValidation]) -> None:
    lines = [
        "# Force Equivalence Validation Table",
        "",
        "| name | sector | model | reference | unit | percent_error | status |",
        "|---|---|---:|---:|---|---:|---|",
    ]
    for row in rows:
        ref = "" if row.reference_value is None else f"{row.reference_value:.8g}"
        pct = "" if row.percent_error is None else f"{row.percent_error:.8g}"
        lines.append(
            f"| {row.name} | {row.sector} | {row.model_value:.8g} | {ref} | {row.unit} | {pct} | {row.status} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def dm_relaxation_c_ratio(z: float) -> float:
    """Matter-corrected supersaturation shape used for the MVP process trace.

    P5.1 gives the landmark window and peak: C/C1=1 at z=0.001 and z=2.846,
    C/C1=1.467 at z=0.898. This smooth log-Gaussian interpolant is used only
    to record the process around the stated landmarks until the full ODE from
    Paper IV is implemented.
    """
    z_left = 0.001
    z_peak = 0.898
    c_peak = 1.467
    sigma = (math.log1p(z_peak) - math.log1p(z_left)) / math.sqrt(2 * math.log(c_peak))
    x = math.log1p(z) - math.log1p(z_peak)
    return 1 + (c_peak - 1) * math.exp(-(x * x) / (2 * sigma * sigma))


def dm_relaxation_eta(z: float, c_ratio: float) -> float:
    if c_ratio <= 1:
        return 1.0
    z_peak = 0.898
    tau_peak = 0.808
    tau_sat = 1.0
    c_peak_ratio = 1.467
    tau_z = tau_peak + (tau_sat - tau_peak) * min(1.0, abs(z - z_peak) / max(1.721 - z_peak, z_peak - 0.345))
    delta_c_ratio = max(0.0, c_ratio - 1.0)
    loss_fraction = delta_c_ratio * (1 - math.exp(-2 * (tau_sat - tau_z))) / (2 * c_peak_ratio)
    return max(0.0, 1 - loss_fraction)


def build_dm_relaxation_process(axioms: Axioms) -> list[DMRelaxationPoint]:
    omega_dm_base = derive_omega_dm(axioms.pi, axioms.e)
    peak_loss_fraction = 0.0508
    omega_dm_bias = peak_loss_fraction * omega_dm_base
    omega_dm_inferred = omega_dm_base + omega_dm_bias
    landmarks = [
        ("weak supersaturation begins", 0.001, None),
        ("strong supersaturation begins", 0.345, None),
        ("relaxation peak / triple coincidence", 0.898, 0.808),
        ("H(z) DESI imprint", 0.93, None),
        ("strong supersaturation ends", 1.721, None),
        ("weak supersaturation ends", 2.846, None),
    ]
    rows = []
    for label, z, tau_process_override in landmarks:
        c_ratio = 1.467 if "peak" in label else dm_relaxation_c_ratio(z)
        if "begins" in label or "ends" in label:
            if "weak" in label:
                c_ratio = 1.0
            elif "strong" in label:
                c_ratio = derive_f_sph(axioms.pi, axioms.e)
        eta = dm_relaxation_eta(z, c_ratio)
        suppression = 1 - eta
        if "peak" in label:
            suppression = peak_loss_fraction
            eta = 1 - suppression
        tau_desitter = tau_from_redshift(z)
        tau_process = tau_process_override if tau_process_override is not None else tau_desitter
        rows.append(
            DMRelaxationPoint(
                label=label,
                z=z,
                tau_desitter=tau_desitter,
                tau_process=tau_process,
                C_ratio_to_C1=c_ratio,
                delta_C_ratio_to_C1=max(0.0, c_ratio - 1.0),
                eta_suppression_kernel=eta,
                f_sigma8_suppression_fraction=suppression,
                omega_dm_base=omega_dm_base,
                omega_dm_bias=omega_dm_bias,
                omega_dm_observer_inferred=omega_dm_inferred,
                process_note="P5.1 relaxation: dDeltaC/dtau=-2DeltaC-2C1; dOmegaDM/dtau=(pi/e)DeltaC eta(tau)",
            )
        )
    return rows


def write_dm_relaxation_csv(path: Path, rows: list[DMRelaxationPoint]) -> None:
    fields = list(asdict(rows[0]).keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_dm_relaxation_md(path: Path, rows: list[DMRelaxationPoint]) -> None:
    lines = [
        "# Dark Matter Relaxation Process",
        "",
        "P5.1 process equations:",
        "",
        "- `dDeltaC/dtau = -2 DeltaC - 2 C1`",
        "- `dOmegaDM/dtau = (pi/e) DeltaC eta(tau)`",
        "- strong supersaturation zone: `z in [0.345, 1.721]`",
        "- peak/triple coincidence: `z ~= 0.898`, `tau ~= 0.808`, `C/C1 ~= 1.467`",
        "- relaxation loss / observational bias: `5.08% * Omega_DM_base`",
        "",
        "| label | z | tau_desitter | tau_process | C/C1 | suppression | Omega_DM_base | bias | observer_inferred |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row.label} | {row.z:.8g} | {row.tau_desitter:.8g} | {row.tau_process:.8g} | {row.C_ratio_to_C1:.8g} | "
            f"{row.f_sigma8_suppression_fraction:.8g} | {row.omega_dm_base:.8g} | "
            f"{row.omega_dm_bias:.8g} | {row.omega_dm_observer_inferred:.8g} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def planck_length_from_inputs(inputs: PhysicalInputs, G_SI: float = 6.674_30e-11) -> float:
    return math.sqrt(inputs.hbar * G_SI / inputs.c**3)


def uv_energy_eV_from_tau(tau: float, inputs: PhysicalInputs) -> tuple[float, float]:
    lP_today = planck_length_from_inputs(inputs)
    lP_eff = effective_lP_from_tau(lP_today, tau)
    energy_j = inputs.hbar * inputs.c / lP_eff
    return lP_eff, joule_to_eV(energy_j)


def nearest_power(value: float, base: float, max_power: int = 80) -> tuple[int, float, float]:
    best_n = 0
    best_value = 1.0
    best_score = float("inf")
    for n in range(max_power + 1):
        candidate = base**n
        score = abs(math.log(value / candidate))
        if score < best_score:
            best_n = n
            best_value = candidate
            best_score = score
    return best_n, best_value, value / best_value


def build_q_capture_estimates(inputs: PhysicalInputs, constants: DerivedConstants) -> list[QCaptureEstimate]:
    m_z_eV = inputs.m_Z_observed_GeV * 1e9
    m_w_eV = wz_ratio_from_stability_inversion() * m_z_eV
    m_e_eV = kg_to_eV_mass(inputs.m_e_kg, inputs)
    tau_sphere = 1 - 0.5 * math.log(constants.f_sph)
    m_s_effective_eV = light_quark_mass_from_alpha_relation(m_e_eV / 1e6, inputs.alpha) * 1e6
    m_q_constituent_eV = ((inputs.m_p_MeV + inputs.m_n_MeV) / 2 / 3) * 1e6
    particles = [
        ("electron", tau_sphere, m_e_eV, "external mass at C_sphere node; Q is posterior"),
        ("light_quark_effective_s", tau_sphere, m_s_effective_eV, "P4 alpha relation target at C_sphere antinode; Q is posterior"),
        ("constituent_q_avg", 1 - math.log(2), m_q_constituent_eV, "nucleon-closure constituent mass; Q is posterior"),
        ("Z0", 0.9125, m_z_eV, "external Z mass anchor; Q is posterior"),
        ("W", 0.9468, m_w_eV, "W from model W/Z ratio using Z anchor; Q is posterior"),
    ]
    rows = []
    for particle, tau, mass_eV, status in particles:
        lP_eff, e_uv = uv_energy_eV_from_tau(tau, inputs)
        q_capture = mass_eV / e_uv
        alpha_n, alpha_value, alpha_ratio = nearest_power(q_capture, inputs.alpha)
        c1_n, c1_value, c1_ratio = nearest_power(q_capture, constants.C1)
        rows.append(
            QCaptureEstimate(
                particle=particle,
                tau_event=tau,
                rest_mass_eV=mass_eV,
                lP_eff_m=lP_eff,
                E_UV_eV=e_uv,
                q_capture=q_capture,
                minus_log_Q=-math.log(q_capture),
                nearest_alpha_power=alpha_n,
                alpha_power_value=alpha_value,
                ratio_to_alpha_power=alpha_ratio,
                nearest_C1_power=c1_n,
                C1_power_value=c1_value,
                ratio_to_C1_power=c1_ratio,
                status=status,
            )
        )
    return rows


def write_q_capture_csv(path: Path, rows: list[QCaptureEstimate]) -> None:
    fields = list(asdict(rows[0]).keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_q_capture_md(path: Path, rows: list[QCaptureEstimate]) -> None:
    lines = [
        "# Q Capture Estimate Table",
        "",
        "`q_capture = rest_mass_energy / E_UV(tau)` with `E_UV = hbar c / lP_eff(tau)`.",
        "",
        "| particle | tau | rest_mass_eV | E_UV_eV | q_capture | -ln(Q) | nearest alpha^n | ratio | nearest C1^n | ratio |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row.particle} | {row.tau_event:.8g} | {row.rest_mass_eV:.8g} | {row.E_UV_eV:.8g} | "
            f"{row.q_capture:.8e} | {row.minus_log_Q:.8g} | {row.nearest_alpha_power} | "
            f"{row.ratio_to_alpha_power:.8g} | {row.nearest_C1_power} | {row.ratio_to_C1_power:.8g} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_alpha_mass_relation_md(path: Path, relation: AlphaMassRelation) -> None:
    lines = [
        "# Alpha Mass Relation",
        "",
        "`alpha = (4/3) * (m_e / m_s)` is recorded as a P4 model relation.",
        "",
        "| name | m_e_MeV | m_s_effective_MeV | m_s_reference_MeV | m_s_percent_error | in_reference_range | alpha_from_ratio | status |",
        "|---|---:|---:|---:|---:|---|---:|---|",
        (
            f"| {relation.name} | {relation.electron_mass_MeV:.8g} | "
            f"{relation.light_quark_effective_mass_MeV:.8g} | {relation.m_s_reference_MeV:.8g} | "
            f"{relation.m_s_percent_error:.8g} | {relation.within_reference_range} | "
            f"{relation.alpha_from_mass_ratio:.8g} | {relation.status} |"
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def alpha_from_electron_light_quark_ratio(m_e_MeV: float, m_s_MeV: float) -> float:
    return (4 / 3) * (m_e_MeV / m_s_MeV)


def light_quark_mass_from_alpha_relation(m_e_MeV: float, alpha: float) -> float:
    return (4 / 3) * m_e_MeV / alpha


def down_quark_mass_from_lambda_relation(m_e_MeV: float, lambda_capacity: float) -> float:
    return m_e_MeV / lambda_capacity


def build_alpha_mass_relation(inputs: PhysicalInputs) -> AlphaMassRelation:
    m_e_mev = kg_to_eV_mass(inputs.m_e_kg, inputs) / 1e6
    m_s_effective_mev = light_quark_mass_from_alpha_relation(m_e_mev, inputs.alpha)
    alpha_from_ratio = alpha_from_electron_light_quark_ratio(m_e_mev, m_s_effective_mev)
    m_s_percent_error = (m_s_effective_mev - inputs.m_s_PDG_central_MeV) / inputs.m_s_PDG_central_MeV * 100
    return AlphaMassRelation(
        name="alpha electron-light-quark mass relation",
        electron_mass_MeV=m_e_mev,
        light_quark_effective_mass_MeV=m_s_effective_mev,
        alpha_from_mass_ratio=alpha_from_ratio,
        alpha_reference=inputs.alpha,
        m_s_reference_MeV=inputs.m_s_PDG_central_MeV,
        m_s_percent_error=m_s_percent_error,
        within_reference_range=inputs.m_s_PDG_min_MeV <= m_s_effective_mev <= inputs.m_s_PDG_max_MeV,
        formula="alpha = (4/3)*(m_e/m_s); m_s = (4/3)*m_e/alpha",
        status="P4 model relation; m_s is the QI antinode effective light-quark mass target, not a fitted PDG current mass",
    )


def effective_mass_row(
    sector: str,
    event_time_definition: str,
    tau_event: float,
    C_ratio_to_C1: float,
    E_UV_eV: float | None,
    effective_mass_MeV: float | None,
    reference_mass_MeV: float | None,
    formula: str,
    parameter_layer: str,
    status: str,
) -> EffectiveParticleMass:
    percent_error = None
    q_capture = None
    minus_log_q = None
    if effective_mass_MeV is not None and reference_mass_MeV not in (None, 0):
        percent_error = (effective_mass_MeV - reference_mass_MeV) / reference_mass_MeV * 100
    if effective_mass_MeV is not None and E_UV_eV not in (None, 0):
        q_capture = effective_mass_MeV * 1e6 / E_UV_eV
        minus_log_q = -math.log(q_capture)
    return EffectiveParticleMass(
        sector=sector,
        event_time_definition=event_time_definition,
        tau_event=tau_event,
        C_ratio_to_C1=C_ratio_to_C1,
        E_UV_eV=E_UV_eV,
        effective_mass_MeV=effective_mass_MeV,
        reference_mass_MeV=reference_mass_MeV,
        percent_error=percent_error,
        q_capture=q_capture,
        minus_log_Q=minus_log_q,
        formula=formula,
        parameter_layer=parameter_layer,
        status=status,
    )


def build_effective_particle_masses(inputs: PhysicalInputs, constants: DerivedConstants) -> list[EffectiveParticleMass]:
    m_e_mev = kg_to_eV_mass(inputs.m_e_kg, inputs) / 1e6
    tau_sphere = 1 - 0.5 * math.log(constants.f_sph)
    tau_alpha = inputs.alpha
    c_alpha_ratio = constants.Schwinger_ratio
    _, e_uv_sphere = uv_energy_eV_from_tau(tau_sphere, inputs)
    _, e_uv_alpha = uv_energy_eV_from_tau(tau_alpha, inputs)
    m_d_mev = down_quark_mass_from_lambda_relation(m_e_mev, constants.lambda_capacity)
    m_s_mev = light_quark_mass_from_alpha_relation(m_e_mev, inputs.alpha)
    return [
        effective_mass_row(
            sector="electron charged-coherence seed",
            event_time_definition="C=C_sphere",
            tau_event=tau_sphere,
            C_ratio_to_C1=constants.f_sph,
            E_UV_eV=e_uv_sphere,
            effective_mass_MeV=m_e_mev,
            reference_mass_MeV=0.510_998_95,
            formula="m_e supplied as the only Layer-3 absolute micro-mass bridge",
            parameter_layer="Layer 3 mass bridge",
            status="not an internal mass prediction; coherent charge onset is model-endogenous",
        ),
        effective_mass_row(
            sector="d effective",
            event_time_definition="tau=alpha",
            tau_event=tau_alpha,
            C_ratio_to_C1=c_alpha_ratio,
            E_UV_eV=e_uv_alpha,
            effective_mass_MeV=m_d_mev,
            reference_mass_MeV=4.70,
            formula="m_d = m_e/lambda, lambda = 1/(f_sph e^2)",
            parameter_layer="Layer 2 ratio using Layer 3 m_e bridge",
            status="P4 lambda interpretation; comparable to PDG current-mass scale with caution",
        ),
        effective_mass_row(
            sector="s effective",
            event_time_definition="tau=alpha",
            tau_event=tau_alpha,
            C_ratio_to_C1=c_alpha_ratio,
            E_UV_eV=e_uv_alpha,
            effective_mass_MeV=m_s_mev,
            reference_mass_MeV=inputs.m_s_PDG_central_MeV,
            formula="m_s = (4/3)*m_e/alpha",
            parameter_layer="Layer 2 P4 alpha relation using Layer 3 m_e bridge",
            status="P4 alpha geometry relation; inside PDG strange-mass range",
        ),
        effective_mass_row(
            sector="u effective",
            event_time_definition="tau=alpha",
            tau_event=tau_alpha,
            C_ratio_to_C1=c_alpha_ratio,
            E_UV_eV=e_uv_alpha,
            effective_mass_MeV=None,
            reference_mass_MeV=2.16,
            formula="missing same-level internal relation",
            parameter_layer="not closed",
            status="not computed; no free parameter added",
        ),
    ]


def write_effective_particle_masses_csv(path: Path, rows: list[EffectiveParticleMass]) -> None:
    fields = list(asdict(rows[0]).keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def write_effective_particle_masses_md(path: Path, rows: list[EffectiveParticleMass]) -> None:
    lines = [
        "# Effective Particle Masses",
        "",
        "Rules:",
        "",
        "- electron charged-coherence seed is evaluated at `C=C_sphere`;",
        "- quark effective masses are evaluated at `tau=alpha`; ",
        "- `m_e` is the only Layer-3 absolute micro-mass bridge;",
        "- alpha is treated as the P4 geometric relation, not a fitted knob;",
        "- `u` remains unclosed until a same-level internal relation is found.",
        "",
        "| sector | event | tau | C/C1 | E_UV_eV | mass_MeV | reference_MeV | percent_error | q_capture | status |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        e_uv = "" if row.E_UV_eV is None else f"{row.E_UV_eV:.8g}"
        mass = "" if row.effective_mass_MeV is None else f"{row.effective_mass_MeV:.8g}"
        ref = "" if row.reference_mass_MeV is None else f"{row.reference_mass_MeV:.8g}"
        pct = "" if row.percent_error is None else f"{row.percent_error:.8g}"
        q_capture = "" if row.q_capture is None else f"{row.q_capture:.8e}"
        lines.append(
            f"| {row.sector} | {row.event_time_definition} | {row.tau_event:.8g} | "
            f"{row.C_ratio_to_C1:.8g} | {e_uv} | {mass} | {ref} | {pct} | {q_capture} | {row.status} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_me_bridge_validation_table(
    inputs: PhysicalInputs,
    axioms: Axioms,
    constants: DerivedConstants,
) -> list[ForceValidation]:
    """Validation rows under the current micro-mass convention.

    In this table, m_e is the only empirical absolute micro-mass bridge. Other
    experimental quantities are references only, not fitting inputs.
    """
    m_e_mev = kg_to_eV_mass(inputs.m_e_kg, inputs) / 1e6
    m_d_mev = down_quark_mass_from_lambda_relation(m_e_mev, constants.lambda_capacity)
    m_s_mev = light_quark_mass_from_alpha_relation(m_e_mev, inputs.alpha)
    wz_ratio = wz_ratio_from_stability_inversion()
    tau_s = 1 - math.log(2)
    return [
        validation_row(
            "C1 fixed-point coefficient",
            "geometry",
            constants.C1,
            None,
            "dimensionless",
            "C1 = pi^2/(16e^2)",
            "no experimental input",
            "model-endogenous",
        ),
        validation_row(
            "Omega_Lambda_fp",
            "gravity / dark energy",
            constants.Omega_Lambda_fp,
            inputs.Omega_Lambda_PDG_2024,
            "dimensionless",
            "Omega_Lambda_fp = pi^3/(6e^2)",
            "PDG/Planck observed rho_Lambda fraction; validation only",
            "model-endogenous prediction compared to observation",
        ),
        validation_row(
            "Omega_DM_fp",
            "dark matter",
            constants.Omega_DM_fp,
            inputs.Omega_DM_Planck_derived,
            "dimensionless",
            "Omega_DM_fp = kappa*pi^4*(sqrt(4pi)-e)/(6e^4)",
            "Planck-inferred dark matter abundance; validation only",
            "kappa-corrected model-endogenous fixed-point value",
        ),
        validation_row(
            "kappa_geom",
            "dark matter",
            constants.kappa_geom,
            inputs.Omega_DM_Planck_derived / constants.Omega_DM_base,
            "dimensionless",
            "kappa = 1/sqrt(1 - 1/(f_sph e^2))",
            "kappa_obs=Omega_DM_obs/Omega_DM_base; validation only",
            "model-endogenous compression factor",
        ),
        validation_row(
            "W/Z ratio",
            "weak",
            wz_ratio,
            inputs.m_W_observed_GeV / inputs.m_Z_observed_GeV,
            "dimensionless",
            "exp[-4(tau_W - tau_Z)]",
            "PDG mass ratio; validation only",
            "ratio prediction; no Z absolute-mass anchor used in this row",
        ),
        validation_row(
            "QCD fixed-point coupling",
            "strong",
            math.exp(-2 * (tau_s - 1)),
            4.0,
            "dimensionless",
            "C_QCD/C1 = exp[-2(tau_s - 1)], tau_s=1-ln2",
            "model closure target",
            "closed internal identity",
        ),
        validation_row(
            "d effective mass",
            "micro mass",
            m_d_mev,
            4.70,
            "MeV",
            "m_d = m_e/lambda, lambda=1/(f_sph e^2)",
            "PDG central current-quark mass; validation only",
            "uses only m_e as absolute mass bridge",
        ),
        validation_row(
            "s effective mass",
            "micro mass",
            m_s_mev,
            inputs.m_s_PDG_central_MeV,
            "MeV",
            "m_s = (4/3)*m_e/alpha",
            "PDG central strange-quark mass; validation only",
            "uses only m_e as absolute mass bridge and P4 alpha relation",
        ),
        validation_row(
            "electron seed mass bridge",
            "micro mass",
            m_e_mev,
            0.510_998_95,
            "MeV",
            "m_e supplied as Layer-3 bridge",
            "CODATA electron mass; this is the bridge itself",
            "not an independent prediction",
        ),
    ]


def write_me_bridge_validation_md(path: Path, rows: list[ForceValidation]) -> None:
    lines = [
        "# m_e-Bridge Validation Table",
        "",
        "Convention:",
        "",
        "- `m_e` is the only empirical absolute micro-mass bridge.",
        "- `alpha` is treated as the P4 geometric relation, not a fitted knob.",
        "- Reference values are used only for validation; they are not used to tune the model rows.",
        "- `u` remains unclosed and is intentionally absent from the prediction rows.",
        "",
        "| quantity | sector | model | reference | unit | percent_error | status |",
        "|---|---|---:|---:|---|---:|---|",
    ]
    for row in rows:
        ref = "" if row.reference_value is None else f"{row.reference_value:.8g}"
        pct = "" if row.percent_error is None else f"{row.percent_error:.8g}"
        lines.append(
            f"| {row.name} | {row.sector} | {row.model_value:.8g} | {ref} | "
            f"{row.unit} | {pct} | {row.status} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path = Path("outputs")) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    axioms = Axioms()
    inputs = PhysicalInputs()
    constants = derive_constants(axioms, inputs)
    events = build_events(inputs, constants)
    tau_rec = tau_from_redshift(inputs.z_rec)
    sample_taus = {row.tau for row in epoch_timeline(-10.0, 1.0, 1201, inputs)}
    sample_taus.update(event.tau for event in events)
    epochs = [build_epoch(tau, inputs, tau_rec) for tau in sorted(sample_taus)]

    event_lookup: dict[float, list[str]] = {}
    for event in events:
        event_lookup.setdefault(round(event.tau, 9), []).append(event.name)
    branches = ["micro_mirror", "radiation_qiv", "matter_qi"]
    labelled_rows: list[ProjectedState] = []
    for state in project_timeline(epochs, branches, axioms, constants):
        event_label = "; ".join(event_lookup.get(round(state.epoch.tau, 9), [])) or None
        labelled_rows.append(
            ProjectedState(
                epoch=state.epoch,
                projection=BranchProjection(**{**asdict(state.projection), "event_label": event_label}),
            )
        )

    cmb_predicted = T_CMB_from_recombination(inputs.T_rec_K, inputs.z_rec)
    wz_ratio = wz_ratio_from_stability_inversion()
    particle_events = build_particle_events(inputs, axioms, constants)
    quark_mass_estimates = build_quark_mass_estimates(inputs)
    force_validations = build_force_validations(inputs, axioms, constants)
    dm_relaxation_process = build_dm_relaxation_process(axioms)
    q_capture_estimates = build_q_capture_estimates(inputs, constants)
    alpha_mass_relation = build_alpha_mass_relation(inputs)
    effective_particle_masses = build_effective_particle_masses(inputs, constants)
    me_bridge_validations = build_me_bridge_validation_table(inputs, axioms, constants)
    constants_payload = {
        "axioms": axioms,
        "physical_inputs": inputs,
        "derived_constants": constants,
        "wz_ratio_prediction": {
            "m_W_over_m_Z": wz_ratio,
            "observed_ratio": inputs.m_W_observed_GeV / inputs.m_Z_observed_GeV,
        },
        "cmb_prediction": {
            "T_CMB_K": cmb_predicted,
            "observed_T_CMB_K": inputs.T_CMB_observed_K,
            "deviation_fraction": (cmb_predicted - inputs.T_CMB_observed_K) / inputs.T_CMB_observed_K,
        },
        "alpha_mass_relation": alpha_mass_relation,
        "micro_mass_scale_status": {
            "only_layer3_absolute_micro_mass_bridge": "m_e",
            "alpha_status": "P4 geometric relation alpha=(4/3)*(m_e/m_s); not a fitted knob",
            "u_sector_status": "not closed; no free parameter added",
        },
    }

    write_json(out_dir / "simulation_constants.json", constants_payload)
    write_timeline_csv(out_dir / "universe_timeline.csv", labelled_rows)
    write_photon_spectrum_csv(out_dir / "photon_spectrum.csv", epochs)
    write_event_table(out_dir / "event_table.md", events)
    write_json(out_dir / "first_photon_event.json", event_payload(events[1], "radiation_qiv", inputs, axioms, constants))
    write_json(out_dir / "first_quark_event.json", event_payload(events[2], "micro_mirror", inputs, axioms, constants))
    write_json(out_dir / "cmb_prediction.json", constants_payload["cmb_prediction"])
    write_json(out_dir / "first_hydrogen_event.json", event_payload(events[4], "matter_qi", inputs, axioms, constants))
    write_json(out_dir / "particle_events.json", particle_events)
    write_particle_events_csv(out_dir / "particle_events.csv", particle_events)
    write_particle_events_md(out_dir / "particle_event_table.md", particle_events)
    write_json(out_dir / "quark_mass_estimates.json", quark_mass_estimates)
    write_quark_mass_estimates_csv(out_dir / "quark_mass_estimates.csv", quark_mass_estimates)
    write_quark_mass_estimates_md(out_dir / "quark_mass_estimates.md", quark_mass_estimates)
    write_json(out_dir / "force_validations.json", force_validations)
    write_force_validations_csv(out_dir / "force_validations.csv", force_validations)
    write_force_validations_md(out_dir / "force_validation_table.md", force_validations)
    write_json(out_dir / "dm_relaxation_process.json", dm_relaxation_process)
    write_dm_relaxation_csv(out_dir / "dm_relaxation_process.csv", dm_relaxation_process)
    write_dm_relaxation_md(out_dir / "dm_relaxation_process.md", dm_relaxation_process)
    write_json(out_dir / "q_capture_estimates.json", q_capture_estimates)
    write_q_capture_csv(out_dir / "q_capture_estimates.csv", q_capture_estimates)
    write_q_capture_md(out_dir / "q_capture_estimates.md", q_capture_estimates)
    write_json(out_dir / "alpha_mass_relation.json", alpha_mass_relation)
    write_alpha_mass_relation_md(out_dir / "alpha_mass_relation.md", alpha_mass_relation)
    write_json(out_dir / "effective_particle_masses.json", effective_particle_masses)
    write_effective_particle_masses_csv(out_dir / "effective_particle_masses.csv", effective_particle_masses)
    write_effective_particle_masses_md(out_dir / "effective_particle_masses.md", effective_particle_masses)
    write_json(out_dir / "me_bridge_validations.json", me_bridge_validations)
    write_force_validations_csv(out_dir / "me_bridge_validations.csv", me_bridge_validations)
    write_me_bridge_validation_md(out_dir / "me_bridge_validation_table.md", me_bridge_validations)


if __name__ == "__main__":
    write_outputs()
