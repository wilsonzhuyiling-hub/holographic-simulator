from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import math
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
RENDERED = OUTPUTS / "rendered"
AGE_TODAY_LCDM_GYR = 13.8


@dataclass(frozen=True)
class ParticleEvent:
    name: str
    branch_name: str
    tau: float
    z: float
    temperature_K: float
    c_ratio: float
    quadrant: str
    status: str


@dataclass(frozen=True)
class DMPoint:
    label: str
    z: float
    tau_desitter: float
    tau_process: float
    c_ratio: float
    suppression: float


@dataclass(frozen=True)
class Validation:
    name: str
    model_value: float
    reference_value: float | None
    percent_error: float | None
    status: str


PALETTE = {
    "paper": (248, 247, 242),
    "ink": (29, 38, 48),
    "muted": (93, 106, 116),
    "line": (204, 199, 189),
    "radiation": (45, 112, 161),
    "matter": (50, 137, 97),
    "weak_z": (139, 88, 164),
    "weak_w": (177, 83, 92),
    "dm": (43, 49, 64),
    "gold": (216, 152, 48),
    "soft_blue": (223, 236, 244),
    "soft_green": (224, 239, 230),
    "soft_red": (242, 225, 226),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def load_particle_events(path: Path = OUTPUTS / "particle_events.csv") -> list[ParticleEvent]:
    rows = []
    for row in read_csv(path):
        rows.append(
            ParticleEvent(
                name=row["name"],
                branch_name=row["branch_name"],
                tau=float(row["epoch_tau"]),
                z=float(row["z"]),
                temperature_K=float(row["temperature_K"]),
                c_ratio=float(row["C_ratio_to_C1"]),
                quadrant=row["quadrant"],
                status=row["first_principles_status"],
            )
        )
    return rows


def load_dm_points(path: Path = OUTPUTS / "dm_relaxation_process.csv") -> list[DMPoint]:
    rows = []
    for row in read_csv(path):
        rows.append(
            DMPoint(
                label=row["label"],
                z=float(row["z"]),
                tau_desitter=float(row["tau_desitter"]),
                tau_process=float(row["tau_process"]),
                c_ratio=float(row["C_ratio_to_C1"]),
                suppression=float(row["f_sigma8_suppression_fraction"]),
            )
        )
    return rows


def load_validations(path: Path = OUTPUTS / "force_validations.csv") -> list[Validation]:
    rows = []
    for row in read_csv(path):
        rows.append(
            Validation(
                name=row["name"],
                model_value=float(row["model_value"]),
                reference_value=as_float(row["reference_value"]),
                percent_error=as_float(row["percent_error"]),
                status=row["status"],
            )
        )
    return rows


def load_constants(path: Path = OUTPUTS / "simulation_constants.json") -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_fonts() -> tuple[ImageFont.ImageFont, ImageFont.ImageFont, ImageFont.ImageFont, ImageFont.ImageFont]:
    candidates = [
        "C:/Windows/Fonts/segoeui.ttf",
        "C:/Windows/Fonts/arial.ttf",
    ]
    bold_candidates = [
        "C:/Windows/Fonts/segoeuib.ttf",
        "C:/Windows/Fonts/arialbd.ttf",
    ]

    def truetype(paths: list[str], size: int) -> ImageFont.ImageFont:
        for path in paths:
            try:
                return ImageFont.truetype(path, size)
            except OSError:
                continue
        return ImageFont.load_default()

    return (
        truetype(bold_candidates, 32),
        truetype(bold_candidates, 20),
        truetype(candidates, 17),
        truetype(candidates, 14),
    )


def tau_to_x(tau: float, left: int, right: int) -> float:
    """Piecewise compression: keep early-universe and late-time events both readable."""
    width = right - left
    if tau <= -6:
        t = max(0.0, min(1.0, (tau + 10) / 4))
        return left + width * (0.04 + 0.20 * t)
    if tau <= 0:
        t = (tau + 6) / 6
        return left + width * (0.24 + 0.33 * t)
    t = min(1.0, tau)
    return left + width * (0.57 + 0.39 * t)


def age_model_gyr_from_tau(tau: float) -> float:
    """Model-anchored age: today's LCDM age is the only external time anchor."""
    return AGE_TODAY_LCDM_GYR * math.exp(tau - 1)


def age_to_x(age_gyr: float, left: int, right: int) -> float:
    """Chronological display axis with extra room around the 5-7.5 Gyr model era."""
    min_age_gyr = 1e-4
    max_age_gyr = AGE_TODAY_LCDM_GYR
    age = max(min_age_gyr, min(max_age_gyr, age_gyr))
    if age <= 1:
        low = math.log10(min_age_gyr)
        high = 0.0
        t = 0.28 * (math.log10(age) - low) / (high - low)
    elif age <= 5:
        t = 0.28 + 0.20 * (age - 1) / 4
    elif age <= 7.5:
        t = 0.48 + 0.30 * (age - 5) / 2.5
    else:
        t = 0.78 + 0.22 * (age - 7.5) / (max_age_gyr - 7.5)
    return left + (right - left) * t


def age_label(age_gyr: float) -> str:
    if age_gyr < 1e-6:
        return f"{age_gyr * 1e9:.0f} yr"
    if age_gyr < 1e-3:
        return f"{age_gyr * 1e6:.0f} kyr"
    if age_gyr < 1:
        return f"{age_gyr * 1e3:.2g} Myr"
    return f"{age_gyr:.2f} Gyr"


def draw_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    text: str,
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    anchor: str | None = None,
) -> None:
    draw.text((int(xy[0]), int(xy[1])), text, font=font, fill=fill, anchor=anchor)


def text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=font)
    return box[2] - box[0], box[3] - box[1]


def draw_chip(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    text: str,
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    outline: tuple[int, int, int],
    text_fill: tuple[int, int, int] = PALETTE["ink"],
) -> None:
    w, h = text_size(draw, text, font)
    x, y = xy
    box = (int(x), int(y), int(x + w + 18), int(y + h + 10))
    draw.rounded_rectangle(box, radius=6, fill=fill, outline=outline, width=1)
    draw_text(draw, (x + 9, y + 5), text, font, text_fill)


def draw_axis(draw: ImageDraw.ImageDraw, left: int, right: int, y: int, font: ImageFont.ImageFont) -> None:
    draw.line((left, y, right, y), fill=PALETTE["line"], width=2)
    ticks = [
        (-10.0, "early grid"),
        (-6.0039741, "recomb."),
        (-0.3470336, "DM weak end"),
        (0.00729735, "alpha"),
        (0.3068528, "QCD"),
        (0.808, "horizon overflow"),
        (1.0, "today"),
    ]
    for idx, (tau, label) in enumerate(ticks):
        age = age_model_gyr_from_tau(tau)
        x = age_to_x(age, left, right)
        draw.line((x, y - 8, x, y + 8), fill=PALETTE["muted"], width=2)
        y_offset = -46 if idx % 2 == 0 else -70
        draw_text(draw, (x, y + y_offset), label, font, PALETTE["muted"], anchor="ma")
        prefix = "Age today" if label == "today" else "Age_model"
        draw_text(draw, (x, y + y_offset + 20), f"{prefix} {age_label(age)}", font, PALETTE["muted"], anchor="ma")


def draw_lane(
    draw: ImageDraw.ImageDraw,
    left: int,
    right: int,
    y: int,
    color: tuple[int, int, int],
    label: str,
    font: ImageFont.ImageFont,
) -> None:
    draw.line((left, y, right, y), fill=(221, 218, 210), width=1)
    draw.line((left, y, right, y), fill=color, width=5)
    draw_text(draw, (left - 16, y), label, font, PALETTE["ink"], anchor="rm")


def event_color(event: ParticleEvent) -> tuple[int, int, int]:
    if event.branch_name == "radiation_qiv":
        return PALETTE["radiation"]
    if event.branch_name == "matter_qi":
        return PALETTE["matter"]
    if event.branch_name == "weak_z_qii":
        return PALETTE["weak_z"]
    if event.branch_name == "weak_w_qiii":
        return PALETTE["weak_w"]
    return PALETTE["gold"]


def event_lane_y(event: ParticleEvent, lanes: dict[str, int]) -> int:
    if event.branch_name == "radiation_qiv":
        return lanes["radiation"]
    if event.branch_name == "matter_qi":
        return lanes["matter"]
    if event.branch_name == "weak_z_qii":
        return lanes["weak_z"]
    if event.branch_name == "weak_w_qiii":
        return lanes["weak_w"]
    return lanes["matter"]


def compact_event_name(name: str) -> str:
    replacements = {
        "QIV radiation / CMB precursor": "CMB precursor",
        "QIV Schwinger saddle marker": "QIV Schwinger",
        "electron charged-coherence seed": "charge seed",
        "first quark confinement-antinode candidate": "quark antinode",
        "QI coherent plane-wave critical": "QI plane-wave",
        "Z0 saddle candidate": "Z0 saddle",
        "W- saddle candidate": "W- saddle",
        "W+ conjugate bridge candidate": "W+ bridge",
    }
    return replacements.get(name, name)


def draw_events(
    draw: ImageDraw.ImageDraw,
    events: Iterable[ParticleEvent],
    left: int,
    right: int,
    lanes: dict[str, int],
    font: ImageFont.ImageFont,
    small: ImageFont.ImageFont,
) -> None:
    placed: list[tuple[float, float, float, float]] = []
    manual_offsets = {
        "QIV radiation / CMB precursor": (0, 22),
        "QIV Schwinger saddle marker": (-64, -44),
        "electron charged-coherence seed": (-60, -118),
        "first quark confinement-antinode candidate": (-58, -62),
        "QI coherent plane-wave critical": (-66, -158),
        "Z0 saddle candidate": (-56, -92),
        "W- saddle candidate": (-50, -76),
        "W+ conjugate bridge candidate": (-54, -96),
    }
    for i, event in enumerate(sorted(events, key=lambda row: (row.tau, row.branch_name))):
        age = age_model_gyr_from_tau(event.tau)
        x = age_to_x(age, left, right)
        y = event_lane_y(event, lanes)
        color = event_color(event)
        draw.ellipse((x - 7, y - 7, x + 7, y + 7), fill=color, outline=(255, 255, 255), width=2)
        draw.line((x, y - 7, x, y - 48), fill=color, width=1)

        label = compact_event_name(event.name)
        w, h = text_size(draw, label, font)
        dx, dy = manual_offsets.get(event.name, (-w / 2 - 10, -82 - (i % 3) * 24))
        label_x = max(left + 4, min(right - w - 20, x + dx))
        label_y = y + dy
        box = (label_x, label_y, label_x + w + 20, label_y + h + 10)
        while any(not (box[2] < p[0] or box[0] > p[2] or box[3] < p[1] or box[1] > p[3]) for p in placed):
            box = (box[0], box[1] - 28, box[2], box[3] - 28)
        placed.append(box)
        draw.rounded_rectangle(tuple(int(v) for v in box), radius=6, fill=(255, 255, 255), outline=color, width=1)
        draw_text(draw, (box[0] + 10, box[1] + 5), label, font, PALETTE["ink"])
        draw_text(draw, (box[0] + 10, box[1] + h + 5), f"Age_model {age_label(age)}", small, PALETTE["muted"])


def draw_dm_overlay(
    draw: ImageDraw.ImageDraw,
    points: list[DMPoint],
    left: int,
    right: int,
    y: int,
    font: ImageFont.ImageFont,
    small: ImageFont.ImageFont,
) -> None:
    weak = [p for p in points if "weak" in p.label]
    strong = [p for p in points if "strong" in p.label]
    if len(weak) >= 2:
        x0 = age_to_x(age_model_gyr_from_tau(min(p.tau_desitter for p in weak)), left, right)
        x1 = age_to_x(age_model_gyr_from_tau(max(p.tau_desitter for p in weak)), left, right)
        draw.rounded_rectangle((x0, y - 34, x1, y + 34), radius=6, fill=(236, 234, 226), outline=(177, 174, 164), width=1)
        draw_text(draw, (x0 + 10, y - 56), "weak relaxation window", small, PALETTE["muted"])
    if len(strong) >= 2:
        x0 = age_to_x(age_model_gyr_from_tau(min(p.tau_desitter for p in strong)), left, right)
        x1 = age_to_x(age_model_gyr_from_tau(max(p.tau_desitter for p in strong)), left, right)
        draw.rounded_rectangle((x0, y - 24, x1, y + 24), radius=6, fill=(215, 222, 231), outline=PALETTE["dm"], width=2)
        draw_text(draw, (x0 + 10, y - 82), "strong overflow window", small, PALETTE["dm"])

    draw.line((left, y, right, y), fill=PALETTE["dm"], width=2)
    for point in points:
        tau_for_display = point.tau_process if "peak" in point.label else point.tau_desitter
        age = age_model_gyr_from_tau(tau_for_display)
        x = age_to_x(age, left, right)
        radius = 5 + int(point.suppression * 110)
        draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=PALETTE["dm"], outline=(255, 255, 255), width=2)
        if "peak" in point.label:
            draw.line((x, y - 112, x, y + 12), fill=PALETTE["gold"], width=2)
            draw.ellipse((x - 9, y - 9, x + 9, y + 9), fill=PALETTE["gold"], outline=(255, 255, 255), width=2)
            draw_text(draw, (x - 122, y - 134), "horizon overflow", small, PALETTE["ink"])
            draw_text(draw, (x - 122, y - 112), f"tau_process={point.tau_process:.3f}", small, PALETTE["muted"])
            draw_text(draw, (x - 122, y - 90), f"Age_model {age_label(age)}, loss={point.suppression:.2%}", small, PALETTE["muted"])
        elif "DESI" in point.label:
            draw_text(draw, (x + 18, y + 22), "H(z) DESI imprint", small, PALETTE["ink"])
            draw_text(draw, (x + 18, y + 44), f"Age_model {age_label(age)}", small, PALETTE["muted"])

    draw_text(draw, (left - 16, y), "DM relaxation", font, PALETTE["ink"], anchor="rm")


def draw_validation_panel(
    draw: ImageDraw.ImageDraw,
    validations: list[Validation],
    constants: dict,
    x: int,
    y: int,
    width: int,
    title_font: ImageFont.ImageFont,
    font: ImageFont.ImageFont,
    small: ImageFont.ImageFont,
) -> None:
    draw.rounded_rectangle((x, y, x + width, y + 330), radius=8, fill=(255, 255, 255), outline=PALETTE["line"], width=1)
    draw_text(draw, (x + 18, y + 16), "Validation Snapshot", title_font, PALETTE["ink"])

    cmb = constants["cmb_prediction"]
    wz = constants["wz_ratio_prediction"]
    rows = [
        ("T_CMB", f"{cmb['T_CMB_K']:.6f} K", f"{cmb['deviation_fraction'] * 100:+.3f}%"),
        ("W/Z", f"{wz['m_W_over_m_Z']:.6f}", f"{(wz['m_W_over_m_Z'] - wz['observed_ratio']) / wz['observed_ratio'] * 100:+.3f}%"),
    ]
    names = {"Omega_Lambda_fp", "Omega_DM_fp", "QCD fixed-point coupling", "entanglement critical scale R_star"}
    for item in validations:
        if item.name in names:
            error = "" if item.percent_error is None else f"{item.percent_error:+.3f}%"
            label = item.name.replace("QCD fixed-point coupling", "QCD C/C1")
            label = label.replace("Omega_Lambda_fp", "Omega_Lambda_fp")
            label = label.replace("entanglement critical scale R_star", "R_star")
            rows.append((label, f"{item.model_value:.6g}", error))

    yy = y + 58
    for name, model, error in rows[:7]:
        draw_text(draw, (x + 18, yy), name, font, PALETTE["ink"])
        draw_text(draw, (x + width - 118, yy), model, small, PALETTE["muted"], anchor="ra")
        draw_text(draw, (x + width - 18, yy), error, small, PALETTE["gold"], anchor="ra")
        yy += 34

    draw.line((x + 18, yy + 3, x + width - 18, yy + 3), fill=PALETTE["line"], width=1)
    yy += 20
    draw_text(draw, (x + 18, yy), "Status flags", font, PALETTE["ink"])
    draw_text(draw, (x + 18, yy + 28), "m_e only absolute micro-mass bridge", small, PALETTE["muted"])
    draw_text(draw, (x + 18, yy + 50), "Z mass anchor external", small, PALETTE["muted"])
    draw_text(draw, (x + 18, yy + 72), "charge seed precedes quark alpha threshold", small, PALETTE["muted"])


def render_compressed_timeline(out_path: Path = RENDERED / "compressed_universe_timeline.png") -> Path:
    events = load_particle_events()
    dm_points = load_dm_points()
    validations = load_validations()
    constants = load_constants()

    RENDERED.mkdir(parents=True, exist_ok=True)
    width, height = 1680, 1120
    left, right = 210, 1115
    top = 112
    title_font, heading_font, font, small = load_fonts()

    img = Image.new("RGB", (width, height), PALETTE["paper"])
    draw = ImageDraw.Draw(img)

    draw_text(draw, (48, 34), "Universe Simulator Rendering Layer", title_font, PALETTE["ink"])
    draw_text(
        draw,
        (50, 76),
        "Chronological Age_model timeline anchored only by today's LCDM age; tau remains the model coordinate behind each marker.",
        font,
        PALETTE["muted"],
    )

    lanes = {
        "radiation": top + 170,
        "matter": top + 370,
        "weak_z": top + 560,
        "weak_w": top + 670,
        "dm": top + 850,
    }
    draw_lane(draw, left, right, lanes["radiation"], PALETTE["radiation"], "QIV radiation", heading_font)
    draw_lane(draw, left, right, lanes["matter"], PALETTE["matter"], "QI matter", heading_font)
    draw_lane(draw, left, right, lanes["weak_z"], PALETTE["weak_z"], "QII Z0 saddle", heading_font)
    draw_lane(draw, left, right, lanes["weak_w"], PALETTE["weak_w"], "QIII W saddle", heading_font)

    draw_axis(draw, left, right, top + 55, small)
    draw_events(draw, events, left, right, lanes, font, small)
    draw_dm_overlay(draw, dm_points, left, right, lanes["dm"], heading_font, small)

    draw_chip(draw, (left, top + 90), "log/linear Age_model axis", small, (255, 255, 255), PALETTE["line"], PALETTE["muted"])
    draw_chip(draw, (left + 194, top + 90), "read-only renderer", small, (255, 255, 255), PALETTE["line"], PALETTE["muted"])
    draw_chip(draw, (left + 352, top + 90), "MVP thresholds, not final event ordering", small, PALETTE["soft_red"], (210, 170, 174), PALETTE["ink"])

    draw_validation_panel(draw, validations, constants, 1210, 152, 410, heading_font, font, small)

    note_y = 1015
    draw.line((50, note_y - 18, width - 50, note_y - 18), fill=PALETTE["line"], width=1)
    draw_text(draw, (50, note_y), "Rendering contract", heading_font, PALETTE["ink"])
    draw_text(
        draw,
        (50, note_y + 32),
        "This image consumes outputs/*.csv and outputs/*.json. Ages are display-only: Age_model = 13.8 Gyr * exp(tau - 1).",
        font,
        PALETTE["muted"],
    )
    draw_text(
        draw,
        (50, note_y + 60),
        "Current caveat: W/Z saddle coordinates are projection parameters; their event epoch is tied to the post-QIV tau=alpha transition.",
        font,
        PALETTE["muted"],
    )

    img.save(out_path)
    return out_path


def main() -> None:
    path = render_compressed_timeline()
    print(path)


if __name__ == "__main__":
    main()
