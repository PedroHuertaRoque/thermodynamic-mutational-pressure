#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

# =========================
# CONFIGURACIÓN
# =========================
INFILE = "alt_dominante_strand_filteredQ20.tsv"
OUTFILE = "mutaciones_clasificadasQ20.tsv"

GEN_MIN = 1
GEN_MAX = 11

MIN_AF_START  = 0.0005   # AF mínima para considerar "presente"
ALLOW_GAP     = 0        # ausencias (AF < MIN_AF_START) permitidas desde aparición hasta GEN_MAX

# Persistencia + crecimiento + estabilidad
MIN_DELTA_AF  = 0.0003   # AF_final - AF_inicio mínimo
ALLOW_DROP    = 0.4     # fracción máxima de "bajadas" dentro del tail (AF[t] < AF[t-1])

MIN_PRESENCE_LEN = 2     # evita "aparece solo al final"

# =========================
# CARGA + PARSEO SEGURO
# =========================
df = pd.read_csv(INFILE, sep="\t", dtype=str)

df["gen"] = df["sample"].astype(str).str.extract(r"barcode(\d+)")[0]
df["gen"] = pd.to_numeric(df["gen"], errors="coerce")
df = df.dropna(subset=["gen"]).copy()
df["gen"] = df["gen"].astype(int)

df = df[(df["gen"] >= GEN_MIN) & (df["gen"] <= GEN_MAX)].copy()

df["AF_alt"] = df["AF_alt"].astype(str).str.replace(",", ".", regex=False)
df["AF_alt"] = pd.to_numeric(df["AF_alt"], errors="coerce")
df = df.dropna(subset=["AF_alt", "pos", "ref_base", "alt_base"]).copy()

df["ref_base"] = df["ref_base"].astype(str).str.upper().str.strip()
df["alt_base"] = df["alt_base"].astype(str).str.upper().str.strip()
df["mut_id"] = df["pos"].astype(str) + "_" + df["ref_base"] + ">" + df["alt_base"]

# =========================
# EVALUACIÓN + CLASIFICACIÓN
# =========================
gens_window = np.arange(GEN_MIN, GEN_MAX + 1, dtype=int)
results = []

for mut, d in df.groupby("mut_id", sort=False):

    # si hay duplicados por gen, toma el AF máximo
    g2af = d.groupby("gen")["AF_alt"].max()

    # serie completa (gens faltantes = 0)
    af_series = np.array([float(g2af.get(g, 0.0)) for g in gens_window], dtype=float)
    present = af_series >= MIN_AF_START

    # defaults por si cae temprano
    status = "PROBABLE_ERROR"
    motivo = ""

    if not present.any():
        # nunca aparece en la ventana
        motivo = "NO_APARECE_EN_VENTANA"
        # opcional: no guardarlo
        continue

    # primer gen donde aparece
    first_idx = int(np.argmax(present))
    gen_inicio = int(gens_window[first_idx])

    # tail (desde aparición hasta GEN_MAX)
    af_tail = af_series[first_idx:]
    present_tail = present[first_idx:]
    tail_len = int(GEN_MAX - gen_inicio + 1)

    # métricas básicas
    AF_inicio = float(af_tail[0])
    AF_final  = float(af_tail[-1])
    delta_AF  = float(AF_final - AF_inicio)

    # gaps (ausencias bajo umbral en el tail)
    gaps_observados = int((~present_tail).sum())

    # drops (bajadas paso a paso)
    diffs = np.diff(af_tail)
    n_drops = int(np.sum(diffs < 0))
    drop_frac = n_drops / max(len(diffs), 1)

    # =========================
    # REGLAS DE CLASIFICACIÓN
    # =========================
    if tail_len < MIN_PRESENCE_LEN:
        motivo = "APARECE_MUY_TARDE"
    elif gaps_observados > ALLOW_GAP:
        motivo = f"SE_PIERDE_EN_EL_TAIL(gaps={gaps_observados}>{ALLOW_GAP})"
    elif delta_AF < MIN_DELTA_AF:
        motivo = f"NO_CRECE(delta_AF={delta_AF:.6g}<{MIN_DELTA_AF})"
    elif drop_frac > ALLOW_DROP:
        motivo = f"DEMASIADAS_BAJADAS(drop_frac={drop_frac:.3f}>{ALLOW_DROP})"
    else:
        status = "SEÑAL_BIOLOGICA"
        motivo = "CUMPLE_PERSISTENCIA_Y_CRECE"

    results.append({
        "mut_id": mut,
        "pos": d["pos"].iloc[0],
        "ref": d["ref_base"].iloc[0],
        "alt": d["alt_base"].iloc[0],
        "gen_inicio": gen_inicio,
        "gen_final": GEN_MAX,
        "AF_inicio": AF_inicio,
        "AF_final": AF_final,
        "delta_AF": delta_AF,
        "n_gens_desde_aparicion": tail_len,
        "gaps_permitidos": int(ALLOW_GAP),
        "gaps_observados": gaps_observados,
        "drop_frac": float(drop_frac),
        "n_drops": n_drops,
        "AF_mean_tail": float(af_tail.mean()),
        "AF_min_tail": float(af_tail.min()),
        "AF_max_tail": float(af_tail.max()),
        "status": status,
        "motivo": motivo
    })

out = pd.DataFrame(results)
if out.empty:
    raise SystemExit("No quedó ninguna mutación evaluable en la ventana (revisa GEN_MIN/GEN_MAX).")

# Orden: primero señales, luego errores; dentro, las más fuertes arriba
out.sort_values(
    ["status", "gen_inicio", "delta_AF", "AF_final"],
    ascending=[True, True, False, False],
    inplace=True
)

out.to_csv(OUTFILE, sep="\t", index=False)
print(f"✔ Clasificación completada: {OUTFILE} | n={len(out)} | señales={(out['status']=='SEÑAL_BIOLOGICA').sum()}")

