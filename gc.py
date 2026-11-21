#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Combina nueve series de contenido GC:
DreamTaq, Klenow y RPA × evgA, ygjV y phnH.

Ejecución:
    python figuras_gc_paper.py

Entrada:
    series.tsv
    Nueve TSV con gen, GC_pred_percent y DP_sum.

Salida:
    resultados/
        Figure_1.pdf / .png / .tif
        Figure_2.pdf / .png / .tif
        estadisticas_gc.tsv
        datos_figuras.tsv

Las filas de cada serie se ordenan por su identificador original
y se asignan a generaciones experimentales consecutivas 1...n.

Los intervalos y valores p corresponden al modelo convencional WLS.
"""

from pathlib import Path
import io
import string

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, FuncFormatter

from PIL import Image
from scipy.stats import spearmanr
import statsmodels.api as sm


# ============================================================
# 1. CONFIGURACIÓN
# ============================================================

BASE_DIR = Path(__file__).resolve().parent

INPUT_MANIFEST = BASE_DIR / "series.tsv"
OUTPUT_DIR = BASE_DIR / "resultados"

# None: buscar Arial o Helvetica instaladas.
# También puedes escribir una ruta:
# FONT_FILE = "/ruta/a/Arial.ttf"
FONT_FILE = None

USE_WLS_WEIGHTS = True

FIRST_FIGURE_NUMBER = 1

FIGURE_WIDTH_MM = 174
FIGURE_HEIGHT_MM = 200

DPI = 600

GENES = ["evgA", "ygjV", "phnH"]
SYSTEMS = ["DreamTaq", "Klenow", "RPA"]

SYSTEM_LABELS = {
    "DreamTaq": "DreamTaq",
    "Klenow": "Klenow fragment",
    "RPA": "RPA",
}

COLORS = {
    "DreamTaq": "#194F73",
    "Klenow": "#B76528",
    "RPA": "#287B68",
}

BAND_COLORS = {
    "DreamTaq": "#BCD1E2",
    "Klenow": "#E8C9AC",
    "RPA": "#B9D9CE",
}


# ============================================================
# 2. FUENTE Y ESTILO
# ============================================================

def configure_font(font_file=None):
    if font_file:
        path = Path(font_file).expanduser()

        if not path.is_absolute():
            path = BASE_DIR / path

        if not path.is_file():
            raise FileNotFoundError(
                f"No se encontró el archivo de fuente: {path}"
            )

        font_manager.fontManager.addfont(str(path))
        family = font_manager.FontProperties(
            fname=str(path)
        ).get_name()

        if family not in ("Arial", "Helvetica"):
            raise ValueError(
                f"La fuente seleccionada es '{family}'. "
                "Se requiere Arial o Helvetica."
            )

    else:
        family = None

        for candidate in ("Arial", "Helvetica"):
            try:
                font_manager.findfont(
                    candidate,
                    fallback_to_default=False,
                )
                family = candidate
                break
            except ValueError:
                continue

        if family is None:
            raise ValueError(
                "No se encontró Arial ni Helvetica.\n"
                "Instala una de estas fuentes o configura FONT_FILE "
                "con la ruta de un archivo Arial.ttf o Helvetica.ttf."
            )

    plt.rcParams.update({
        "font.family": family,
        "font.size": 8,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.7,
        "axes.edgecolor": "#333333",
        "text.color": "#111111",
        "axes.labelcolor": "#111111",
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "lines.linewidth": 1.1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "figure.autolayout": False,
        "figure.constrained_layout.use": False,
        "savefig.bbox": None,
        "axes.unicode_minus": True,
    })

    print(f"Fuente utilizada: {family}")


# ============================================================
# 3. FORMATOS NUMÉRICOS
# ============================================================

def format_four_decimals(value, position=None):
    """Etiquetas del eje Y con cuatro decimales."""
    if abs(value) < 0.00005:
        value = 0.0

    return f"{value:.4f}"


def format_beta(value):
    """Tres cifras significativas, sin notación científica."""
    value = float(value)

    if not np.isfinite(value):
        return "NA"

    if value == 0:
        return "0.00"

    exponent = int(np.floor(np.log10(abs(value))))
    decimals = max(0, 2 - exponent)

    rounded = round(value, decimals)

    if rounded != 0:
        exponent = int(np.floor(np.log10(abs(rounded))))
        decimals = max(0, 2 - exponent)

    return f"{rounded:.{decimals}f}"


def format_rho(value):
    if not np.isfinite(value):
        return "NA"

    return f"{value:.3f}"


# ============================================================
# 4. CINCO MARCAS EN Y Y MENOS ESPACIO VACÍO
# ============================================================

def five_y_ticks(values, reserve_statistics=False):
    """
    Genera cinco marcas equidistantes.

    Los intervalos son múltiplos de 0.0001 para evitar etiquetas
    repetidas cuando se muestran cuatro decimales.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        raise ValueError(
            "No hay valores finitos para calcular el eje Y."
        )

    low = float(values.min())
    high = float(values.max())

    span = max(high - low, 0.0004)
    center = (low + high) / 2.0

    data_low = center - span / 2.0
    data_high = center + span / 2.0

    target_low = data_low - 0.06 * span

    if reserve_statistics:
        target_high = data_high + 0.30 * span
    else:
        target_high = data_high + 0.08 * span

    quantum = 0.0001

    start_units = int(np.floor(target_low / quantum))
    end_units = int(np.ceil(target_high / quantum))

    step_units = max(
        1,
        int(np.ceil((end_units - start_units) / 4.0)),
    )

    ticks = (
        start_units
        + np.arange(5, dtype=np.int64) * step_units
    ) * quantum

    ticks = np.round(ticks, 4)

    margin = 0.015 * (ticks[-1] - ticks[0])

    return (
        ticks,
        float(ticks[0] - margin),
        float(ticks[-1] + margin),
    )


# ============================================================
# 5. LECTURA DE LAS NUEVE SERIES
# ============================================================

def read_series(manifest_path):
    manifest_path = Path(manifest_path)

    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"No se encontró:\n{manifest_path}\n"
            "Coloca series.tsv en la misma carpeta que el script."
        )

    config = pd.read_csv(
        manifest_path,
        sep="\t",
        dtype=str,
    )

    config.columns = config.columns.str.strip()

    required = {"system", "gene", "tsv"}
    missing = required - set(config.columns)

    if missing:
        raise ValueError(
            f"Faltan columnas en series.tsv: {sorted(missing)}"
        )

    for column in required:
        if config[column].isna().any():
            raise ValueError(
                f"Hay valores vacíos en la columna '{column}'."
            )

        config[column] = config[column].str.strip()

        if config[column].eq("").any():
            raise ValueError(
                f"Hay valores vacíos en la columna '{column}'."
            )

    if config.duplicated(["system", "gene"]).any():
        raise ValueError(
            "Hay combinaciones system/gene duplicadas en series.tsv."
        )

    expected = {
        (system, gene)
        for system in SYSTEMS
        for gene in GENES
    }

    observed = set(zip(config["system"], config["gene"]))

    if observed != expected:
        raise ValueError(
            "series.tsv debe contener exactamente las nueve "
            "combinaciones de DreamTaq, Klenow y RPA con "
            "evgA, ygjV y phnH.\n"
            f"Faltan: {sorted(expected - observed)}\n"
            f"No reconocidas: {sorted(observed - expected)}"
        )

    if "baseline_gen" in config.columns:
        print(
            "Nota: baseline_gen se ignora. "
            "La primera observación de cada serie corresponde "
            "a la generación experimental 1."
        )

    data = {}

    for row in config.itertuples(index=False):
        path = Path(row.tsv).expanduser()

        if not path.is_absolute():
            path = manifest_path.parent / path

        if not path.is_file():
            raise FileNotFoundError(
                f"No se encontró el archivo de "
                f"{row.system}/{row.gene}:\n{path}"
            )

        df = pd.read_csv(path, sep="\t")
        df.columns = df.columns.str.strip()

        required_columns = {"gen", "GC_pred_percent"}

        if USE_WLS_WEIGHTS:
            required_columns.add("DP_sum")

        missing_columns = required_columns - set(df.columns)

        if missing_columns:
            raise ValueError(
                f"{path.name}: faltan columnas "
                f"{sorted(missing_columns)}"
            )

        if len(df) < 3:
            raise ValueError(
                f"{row.system}/{row.gene}: se necesitan "
                "al menos tres generaciones."
            )

        for column in required_columns:
            df[column] = pd.to_numeric(
                df[column],
                errors="raise",
            )

            if not np.isfinite(
                df[column].to_numpy(dtype=float)
            ).all():
                raise ValueError(
                    f"{path.name}: '{column}' contiene "
                    "valores ausentes o no finitos."
                )

        if (df["gen"] % 1 != 0).any():
            raise ValueError(
                f"{path.name}: gen debe contener enteros."
            )

        if df["gen"].duplicated().any():
            raise ValueError(
                f"{path.name}: hay generaciones originales repetidas."
            )

        if not df["GC_pred_percent"].between(0, 100).all():
            raise ValueError(
                f"{path.name}: GC_pred_percent debe estar entre 0 y 100."
            )

        if USE_WLS_WEIGHTS and (df["DP_sum"] <= 0).any():
            raise ValueError(
                f"{path.name}: DP_sum debe ser mayor que cero."
            )

        df["gen"] = df["gen"].astype(int)

        df = (
            df.sort_values("gen")
            .reset_index(drop=True)
            .copy()
        )

        # Conserva el identificador original de las muestras.
        df["source_gen"] = df["gen"]

        # El usuario confirma que el orden corresponde a generaciones
        # experimentales consecutivas, aunque cambien los barcodes.
        df["generation"] = np.arange(
            1,
            len(df) + 1,
            dtype=int,
        )

        df["gen_ord"] = df["generation"]

        df["system"] = row.system
        df["gene"] = row.gene

        if (
            "n_muts" in df.columns
            and df["n_muts"].nunique(dropna=False) > 1
        ):
            print(
                f"REVISAR {row.system}/{row.gene}: "
                "n_muts cambia entre generaciones."
            )

        data[(row.system, row.gene)] = df

        print(
            f"Leído {row.system}/{row.gene}: "
            f"{len(df)} observaciones; "
            f"identificadores originales "
            f"{df['source_gen'].iloc[0]}–"
            f"{df['source_gen'].iloc[-1]}; "
            f"generaciones experimentales 1–{len(df)}."
        )

    return data


# ============================================================
# 6. SPEARMAN, WLS, IC95 % Y R²
# ============================================================

def calculate_statistics(df):
    df = df.copy()

    x = df["generation"].to_numpy(dtype=float)
    y = df["GC_pred_percent"].to_numpy(dtype=float)

    if USE_WLS_WEIGHTS:
        raw_weights = df["DP_sum"].to_numpy(dtype=float)
        model_name = "WLS_DP_sum"
    else:
        raw_weights = np.ones_like(y)
        model_name = "OLS"

    # Reescalar todos los pesos por una constante no cambia el ajuste.
    weights = raw_weights / raw_weights.mean()

    df["relative_weight"] = weights

    baseline_gc = float(y[0])
    df["baseline_gc_percent"] = baseline_gc
    df["change_pp"] = y - baseline_gc

    constant_series = bool(np.all(y == y[0]))

    if constant_series:
        fitted = np.full_like(y, y[0])

        ci_low = np.full_like(y, np.nan)
        ci_high = np.full_like(y, np.nan)

        slope = 0.0
        slope_ci_low = np.nan
        slope_ci_high = np.nan
        p_slope = np.nan

        rho = np.nan
        p_rho = np.nan
        r2_weighted = np.nan

        status = "constant_series"

    else:
        rho, p_rho = spearmanr(x, y)

        X = sm.add_constant(
            x,
            has_constant="add",
        )

        fit = sm.WLS(
            y,
            X,
            weights=weights,
        ).fit()

        prediction = fit.get_prediction(X).summary_frame(
            alpha=0.05
        )

        fitted = prediction["mean"].to_numpy(dtype=float)

        ci_low = prediction[
            "mean_ci_lower"
        ].to_numpy(dtype=float)

        ci_high = prediction[
            "mean_ci_upper"
        ].to_numpy(dtype=float)

        slope = float(fit.params[1])
        p_slope = float(fit.pvalues[1])

        confidence_intervals = np.asarray(
            fit.conf_int(alpha=0.05),
            dtype=float,
        )

        slope_ci_low = float(confidence_intervals[1, 0])
        slope_ci_high = float(confidence_intervals[1, 1])

        weighted_mean = np.average(y, weights=weights)

        sse = np.sum(weights * (y - fitted) ** 2)
        sst = np.sum(weights * (y - weighted_mean) ** 2)

        r2_weighted = (
            float(1.0 - sse / sst)
            if sst > 0
            else np.nan
        )

        status = "ok"

    df["GC_fit"] = fitted
    df["GC_ci_low"] = ci_low
    df["GC_ci_high"] = ci_high
    df["residual_GC_pp"] = y - fitted

    # Figura 2: se resta el mismo GC observado en G1.
    # Esta banda es la banda del ajuste absoluto trasladada.
    # No incorpora una estimación adicional de la incertidumbre de G1.
    df["change_fit"] = fitted - baseline_gc
    df["change_ci_low"] = ci_low - baseline_gc
    df["change_ci_high"] = ci_high - baseline_gc

    statistics = {
        "system": str(df["system"].iloc[0]),
        "gene": str(df["gene"].iloc[0]),
        "model": model_name,
        "status": status,
        "n_generations": len(df),
        "first_source_gen": int(df["source_gen"].iloc[0]),
        "last_source_gen": int(df["source_gen"].iloc[-1]),
        "GC_first_percent": float(y[0]),
        "GC_last_percent": float(y[-1]),
        "change_first_to_last_pp": float(y[-1] - y[0]),
        "slope_pp_per_generation": slope,
        "slope_ci95_low_model": slope_ci_low,
        "slope_ci95_high_model": slope_ci_high,
        "p_slope_model": p_slope,
        "R2_weighted": r2_weighted,
        "spearman_rho": float(rho),
        "p_spearman_standard": float(p_rho),
    }

    return df, statistics


# ============================================================
# 7. EXPORTACIÓN: PDF, PNG Y TIFF
# ============================================================

def export_figure(fig, stem):
    stem = Path(stem)

    # Sin recorte automático: conserva el tamaño físico.
    fig.savefig(
        str(stem) + ".pdf",
        facecolor="white",
        bbox_inches=None,
    )

    buffer = io.BytesIO()

    fig.savefig(
        buffer,
        format="png",
        dpi=DPI,
        facecolor="white",
        transparent=False,
        bbox_inches=None,
    )

    buffer.seek(0)

    with Image.open(buffer) as source:
        # RGB, ocho bits por canal, sin alfa.
        rgb = source.convert("RGB")

        rgb.save(
            str(stem) + ".png",
            dpi=(DPI, DPI),
        )

        rgb.save(
            str(stem) + ".tif",
            dpi=(DPI, DPI),
            compression="tiff_lzw",
        )

    buffer.close()
    plt.close(fig)

    print(f"Guardado: {stem.name} [PDF, PNG, TIFF]")


# ============================================================
# 8. FIGURAS DE NUEVE PANELES
# ============================================================

def make_figure(data, statistics_by_series, metric, number):
    absolute = metric == "GC_pred_percent"

    if absolute:
        observed_column = "GC_pred_percent"
        fitted_column = "GC_fit"
        lower_column = "GC_ci_low"
        upper_column = "GC_ci_high"

        y_label = "Estimated GC content (%)"

    else:
        observed_column = "change_pp"
        fitted_column = "change_fit"
        lower_column = "change_ci_low"
        upper_column = "change_ci_high"

        y_label = "Change in estimated GC content (pp)"

    fig, axes = plt.subplots(
        nrows=3,
        ncols=3,
        figsize=(
            FIGURE_WIDTH_MM / 25.4,
            FIGURE_HEIGHT_MM / 25.4,
        ),
        sharex=False,
        sharey=False,
    )

    # Espacio izquierdo reservado para las etiquetas.
    # Márgenes superior e inferior más equilibrados.
    fig.subplots_adjust(
        left=0.18,
        right=0.985,
        bottom=0.085,
        top=0.93,
        wspace=0.75,
        hspace=0.48,
    )

    for i, system in enumerate(SYSTEMS):
        for j, gene in enumerate(GENES):
            ax = axes[i, j]

            df = data[(system, gene)]
            stats = statistics_by_series[(system, gene)]

            x = df["generation"].to_numpy(dtype=float)

            y = df[observed_column].to_numpy(dtype=float)
            fitted = df[fitted_column].to_numpy(dtype=float)

            ci_low = df[lower_column].to_numpy(dtype=float)
            ci_high = df[upper_column].to_numpy(dtype=float)

            valid_band = (
                np.isfinite(ci_low)
                & np.isfinite(ci_high)
            )

            if valid_band.any():
                ax.fill_between(
                    x,
                    ci_low,
                    ci_high,
                    where=valid_band,
                    color=BAND_COLORS[system],
                    alpha=0.65,
                    linewidth=0,
                    zorder=1,
                )

            if not absolute:
                ax.axhline(
                    0.0,
                    color="#888888",
                    linewidth=0.65,
                    linestyle=":",
                    zorder=0,
                )

            ax.plot(
                x,
                fitted,
                color="#222222",
                linewidth=1.1,
                linestyle=(0, (4, 2.5)),
                zorder=2,
            )

            ax.scatter(
                x,
                y,
                s=21,
                color=COLORS[system],
                edgecolor="#253038",
                linewidth=0.55,
                zorder=3,
            )

            # Escala independiente por panel.
            axis_values = np.concatenate([
                y,
                fitted,
                ci_low,
                ci_high,
            ])

            if not absolute:
                axis_values = np.append(axis_values, 0.0)

            ticks, ymin, ymax = five_y_ticks(
                axis_values,
                reserve_statistics=absolute,
            )

            ax.yaxis.set_major_locator(FixedLocator(ticks))

            ax.yaxis.set_major_formatter(
                FuncFormatter(format_four_decimals)
            )

            ax.set_ylim(ymin, ymax)

            # Marcas X: generaciones experimentales enteras.
            n_generations = len(df)

            if n_generations <= 8:
                x_ticks = np.arange(
                    1,
                    n_generations + 1,
                    dtype=int,
                )
            else:
                x_ticks = np.unique(
                    np.rint(
                        np.linspace(1, n_generations, 6)
                    ).astype(int)
                )

            ax.xaxis.set_major_locator(
                FixedLocator(x_ticks)
            )

            ax.xaxis.set_major_formatter(
                FuncFormatter(
                    lambda value, position: f"{int(value)}"
                )
            )

            ax.set_xlim(
                0.7,
                n_generations + 0.3,
            )

            ax.tick_params(
                axis="both",
                which="major",
                direction="out",
                width=0.65,
                length=3,
                pad=2,
                labelsize=8,
            )

            ax.minorticks_off()
            ax.grid(False)

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            # Letras de panel a–i.
            panel_letter = string.ascii_lowercase[i * 3 + j]

            ax.text(
                -0.04,
                1.035,
                panel_letter,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                color="#111111",
                clip_on=False,
            )

            # Solo dos estadísticos, en la primera figura.
            if absolute:
                slope = stats["slope_pp_per_generation"]
                rho = stats["spearman_rho"]

                ax.text(
                    0.04,
                    0.97,
                    f"β = {format_beta(slope)}\n"
                    f"ρ = {format_rho(rho)}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=8,
                    linespacing=1.15,
                    color="#222222",
                    zorder=5,
                )

    # Encabezados de columnas centrados en cada panel.
    for j, gene in enumerate(GENES):
        position = axes[0, j].get_position()

        fig.text(
            (position.x0 + position.x1) / 2.0,
            0.972,
            gene,
            ha="center",
            va="center",
            fontsize=10,
            fontstyle="italic",
            color="#111111",
            clip_on=False,
        )

    # Etiquetas de sistemas centradas en cada fila.
    for i, system in enumerate(SYSTEMS):
        position = axes[i, 0].get_position()
        row_center = (position.y0 + position.y1) / 2.0

        fig.text(
            0.067,
            row_center,
            SYSTEM_LABELS[system],
            ha="center",
            va="center",
            rotation=90,
            rotation_mode="anchor",
            fontsize=9,
            fontweight="bold",
            color=COLORS[system],
            clip_on=False,
        )

    # Eje Y alineado verticalmente con la fila central.
    middle_position = axes[1, 0].get_position()
    center_y = (
        middle_position.y0 + middle_position.y1
    ) / 2.0

    fig.text(
        0.030,
        center_y,
        y_label,
        ha="center",
        va="center",
        rotation=90,
        rotation_mode="anchor",
        fontsize=9,
        color="#111111",
        clip_on=False,
    )

    # Eje X centrado en la imagen, con menor margen inferior.
    fig.text(
        0.50,
        0.025,
        "Generation",
        ha="center",
        va="center",
        fontsize=9,
        color="#111111",
        clip_on=False,
    )

    # Sin títulos, leyendas ni notas inferiores.
    # No usar tight_layout ni bbox_inches="tight".
    export_figure(
        fig,
        OUTPUT_DIR / f"Figure_{number}",
    )


# ============================================================
# 9. PROGRAMA PRINCIPAL
# ============================================================

def main():
    configure_font(FONT_FILE)

    data = read_series(INPUT_MANIFEST)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    statistics_by_series = {}
    statistics_rows = []

    for system in SYSTEMS:
        for gene in GENES:
            key = (system, gene)

            processed, statistics = calculate_statistics(
                data[key]
            )

            data[key] = processed
            statistics_by_series[key] = statistics
            statistics_rows.append(statistics)

    # Exporta los valores numéricos completos.
    # El redondeo de las figuras no se aplica a estos archivos.
    combined = pd.concat(
        [
            data[(system, gene)]
            for system in SYSTEMS
            for gene in GENES
        ],
        ignore_index=True,
    )

    combined.to_csv(
        OUTPUT_DIR / "datos_figuras.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    statistics_table = pd.DataFrame(statistics_rows)

    statistics_table.to_csv(
        OUTPUT_DIR / "estadisticas_gc.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    make_figure(
        data=data,
        statistics_by_series=statistics_by_series,
        metric="GC_pred_percent",
        number=FIRST_FIGURE_NUMBER,
    )

    make_figure(
        data=data,
        statistics_by_series=statistics_by_series,
        metric="change_pp",
        number=FIRST_FIGURE_NUMBER + 1,
    )

    print()
    print("Proceso completado.")
    print(f"Carpeta de salida: {OUTPUT_DIR.resolve()}")
    print(
        f"Tamaño de cada figura: "
        f"{FIGURE_WIDTH_MM} × {FIGURE_HEIGHT_MM} mm."
    )
    print(f"PNG y TIFF: RGB, 8 bits por canal, {DPI} ppp.")
    print("PDF: formato vectorial.")
    print("Estadísticas: estadisticas_gc.tsv")
    print("Datos y ajustes: datos_figuras.tsv")
    print(
        "Figura 2: cambio respecto al GC observado en G1; "
        "la banda corresponde al IC del ajuste absoluto "
        "desplazado por ese mismo valor."
    )


if __name__ == "__main__":
    main()
