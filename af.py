#!/usr/bin/env python3
"""Panel de sumas AF. Ejecutar: python figuras_balance_af.py"""
from pathlib import Path
import io
import string
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter
from PIL import Image

BASE_DIR = Path(__file__).resolve().parent
INPUT_MANIFEST = BASE_DIR / 'series_af.tsv'
OUTPUT_DIR = BASE_DIR / 'resultados_af'
FONT_FILE = None
# False exige Arial/Helvetica. True permite Liberation Sans, avisando.
ALLOW_FONT_FALLBACK = True
# Las series ausentes se muestran como "Data not provided".
# False exige los nueve archivos antes de exportar.
ALLOW_MISSING_SERIES = False
SHOW_LEGEND = True
SHOW_ANNOTATIONS = True  # Solo Δnet; sin correlación ni valores p
# 'panel': muestra cada trayectoria; 'global': comparación en escala común.
Y_SCALE = 'panel'
DPI = 600
WIDTH_MM, HEIGHT_MM = 174, 200
GENES = ['evgA', 'ygjV', 'phnH']
SYSTEMS = ['DreamTaq', 'Klenow', 'RPA']
LABELS = {'DreamTaq':'DreamTaq', 'Klenow':'Klenow fragment', 'RPA':'RPA'}
COLUMNS = ['sum_AF_Fav_GC', 'sum_AF_Fav_AT', 'net_AF_signed']
COLORS = ['#194F73', '#92278F', '#222222']
MARKERS = ['o', 's', '^']
LINESTYLES = ['-', '-', '--']
CURVE_LABELS = ['AT→GC', 'GC→AT', 'Net (GC − AT)']


def configure_font():
    family = None
    if FONT_FILE:
        path = Path(FONT_FILE).expanduser()
        if not path.is_absolute(): path = BASE_DIR / path
        font_manager.fontManager.addfont(str(path))
        family = font_manager.FontProperties(fname=str(path)).get_name()
        if family not in ('Arial', 'Helvetica'):
            raise ValueError('FONT_FILE debe corresponder a Arial o Helvetica.')
    else:
        for candidate in ['Arial', 'Helvetica']:
            try:
                font_manager.findfont(candidate, fallback_to_default=False)
                family = candidate
                break
            except ValueError: pass
    if family is None:
        if not ALLOW_FONT_FALLBACK:
            raise ValueError('Instala Arial/Helvetica o configura FONT_FILE.')
        for candidate in ['Liberation Sans', 'DejaVu Sans']:
            try:
                font_manager.findfont(candidate, fallback_to_default=False)
                family = candidate
                break
            except ValueError: pass
        if family is None: raise ValueError('No hay una fuente sans serif disponible.')
        warnings.warn(f'Se utiliza {family}: Arial/Helvetica no disponibles.')
    plt.rcParams.update({'font.family':family, 'font.size':8,
        'axes.labelsize':9, 'xtick.labelsize':8, 'ytick.labelsize':8,
        'axes.linewidth':.7, 'pdf.fonttype':42, 'ps.fonttype':42,
        'figure.autolayout':False, 'figure.constrained_layout.use':False,
        'savefig.bbox':None, 'savefig.facecolor':'white'})
    print('Fuente:', family)
    return family


def scientific_label(value, decimals=2):
    """Notación científica explícita; no cambia los datos exportados."""
    if value == 0:
        return "0"
    if not np.isfinite(value):
        return "NA"
    mantissa, exponent = f"{float(value):.{decimals}e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


def read_data():
    cfg = pd.read_csv(INPUT_MANIFEST, sep='\t', dtype=str)
    cfg.columns = cfg.columns.str.strip()
    if not {'system','gene','tsv'}.issubset(cfg):
        raise ValueError('series_af.tsv requiere system, gene, tsv.')
    for col in ['system','gene','tsv']:
        if cfg[col].isna().any(): raise ValueError(f'Valores vacíos en {col}.')
        cfg[col] = cfg[col].str.strip()
    expected = {(s,g) for s in SYSTEMS for g in GENES}
    if cfg.duplicated(['system','gene']).any() or set(zip(cfg.system,cfg.gene)) != expected:
        raise ValueError('El manifiesto requiere exactamente las nueve combinaciones.')
    data, missing = {}, []
    for row in cfg.itertuples(index=False):
        path = Path(row.tsv)
        if not path.is_absolute(): path = INPUT_MANIFEST.parent / path
        key = (row.system,row.gene)
        if not path.is_file():
            missing.append('/'.join(key))
            continue
        d = pd.read_csv(path,sep='\t')
        d.columns = d.columns.str.strip()
        req = ['gen', *COLUMNS, 'n_mut']
        if not set(req).issubset(d): raise ValueError(f'{path}: faltan {set(req)-set(d)}')
        if d.empty: raise ValueError(f'{path}: archivo vacío.')
        for col in req + (['gen_norm'] if 'gen_norm' in d else []):
            d[col] = pd.to_numeric(d[col],errors='raise')
            if not np.isfinite(d[col]).all(): raise ValueError(f'{path}: {col} no finito.')
        for col in ['gen','n_mut']:
            if (d[col]%1 != 0).any(): raise ValueError(f'{path}: {col} no entero.')
        if d.gen.duplicated().any(): raise ValueError(f'{path}: gen duplicado.')
        if (d[['sum_AF_Fav_GC','sum_AF_Fav_AT','n_mut']] < 0).any().any():
            raise ValueError(f'{path}: conteos o sumas AF negativos.')
        if ((d.sum_AF_Fav_GC+d.sum_AF_Fav_AT)>d.n_mut+1e-6).any():
            raise ValueError(f'{path}: sumas AF incompatibles con n_mut.')
        if not np.allclose(d.net_AF_signed,d.sum_AF_Fav_GC-d.sum_AF_Fav_AT,rtol=1e-6,atol=1e-8):
            raise ValueError(f'{path}: net_AF_signed no coincide con GC menos AT.')
        if 'gen_norm' in d:
            if (d.gen_norm%1 != 0).any() or d.gen_norm.duplicated().any() or (d.gen_norm<1).any():
                raise ValueError(f'{path}: gen_norm debe contener enteros positivos únicos.')
            d = d.sort_values('gen_norm').copy()
            d['generation'] = d.gen_norm.astype(int)
        else:
            d = d.sort_values('gen').copy()
            d['generation'] = np.arange(1,len(d)+1)
        d['system'],d['gene'] = key
        data[key] = d
        if d.n_mut.nunique()>1: print(f'REVISAR {key}: n_mut varía entre generaciones.')
    if missing:
        print('Series sin archivo:', ', '.join(missing))
        if not ALLOW_MISSING_SERIES: raise FileNotFoundError('Faltan series del panel.')
    if not data: raise ValueError('No hay ninguna serie disponible.')
    return data, missing


def five_ticks(values):
    """Cinco marcas con pasos sencillos, incluyendo cero y valores negativos."""
    low, high = min(0,float(np.min(values))), max(0,float(np.max(values)))
    span = max(high-low, .0004)
    target_high = high + (.16 if SHOW_ANNOTATIONS else .08)*span
    target_low = low - (.04*span if low<0 else 0)
    required = (target_high-target_low)/4
    exponent = int(np.floor(np.log10(max(required,1e-4))))
    for e in range(exponent-1,exponent+5):
        for m in [1,1.5,2,2.5,3,4,5,6,8,10]:
            step=m*10.**e
            if step<max(required,1e-4): continue
            if not np.isclose(step/1e-4,round(step/1e-4)): continue
            start=np.floor(target_low/step+1e-10)*step
            if start+4*step>=target_high-1e-12:
                ticks=start+np.arange(5)*step
                return ticks, (start-.025*4*step,start+4.025*step)
    raise ValueError('No se pudieron determinar las marcas Y.')


def export(fig, stem):
    fig.savefig(str(stem)+'.pdf',bbox_inches=None)
    buf=io.BytesIO()
    fig.savefig(buf,format='png',dpi=DPI,bbox_inches=None,facecolor='white')
    buf.seek(0)
    with Image.open(buf) as src:
        rgb=src.convert('RGB')
        rgb.save(str(stem)+'.png',dpi=(DPI,DPI))
        rgb.save(str(stem)+'.tif',dpi=(DPI,DPI),compression='tiff_lzw')
    buf.close(); plt.close(fig)


def make_panel(data, missing):
    fig,axes=plt.subplots(3,3,figsize=(WIDTH_MM/25.4,HEIGHT_MM/25.4))
    fig.subplots_adjust(left=.20,right=.985,bottom=.13 if SHOW_LEGEND else .085,
                        top=.93,wspace=.75,hspace=.48)
    global_values=np.concatenate([d[COLUMNS].to_numpy().ravel() for d in data.values()])
    summaries=[]
    for i,s in enumerate(SYSTEMS):
        for j,g in enumerate(GENES):
            ax=axes[i,j]
            ax.text(-.04,1.035,string.ascii_lowercase[i*3+j],transform=ax.transAxes,
                    fontsize=10,fontweight='bold',ha='left',va='bottom')
            if (s,g) not in data:
                ax.set_xticks([]);ax.set_yticks([])
                for spine in ax.spines.values(): spine.set_visible(False)
                ax.text(.5,.5,'Data not provided',transform=ax.transAxes,
                        ha='center',va='center',fontsize=8,color='.5')
                continue
            d=data[(s,g)];x=d.generation.to_numpy()
            # Interrumpe las líneas si gen_norm tiene generaciones ausentes.
            plot=d.set_index('generation').reindex(range(int(x.min()),int(x.max())+1))
            for col,color,marker,ls in zip(COLUMNS,COLORS,MARKERS,LINESTYLES):
                ax.plot(plot.index,plot[col],color=color,marker=marker,linestyle=ls,
                        linewidth=1.05,markersize=3.1,markeredgewidth=.45,
                        markerfacecolor='white' if col=='net_AF_signed' else color,
                        markeredgecolor=color,zorder=3)
            values=global_values if Y_SCALE=='global' else d[COLUMNS].to_numpy().ravel()
            ticks,limits=five_ticks(values)
            ax.set_yticks(ticks);ax.set_ylim(*limits)
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p: scientific_label(v, 2)))
            xt=x if len(x)<=8 else x[np.unique(np.rint(np.linspace(0,len(x)-1,6)).astype(int))]
            ax.xaxis.set_major_locator(FixedLocator(xt))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{int(v)}'))
            ax.set_xlim(x.min()-.3,x.max()+.3)
            ax.axhline(0,color='.7',linewidth=.55,zorder=0)
            ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
            ax.tick_params(width=.65,length=3,pad=2)
            net=d.net_AF_signed.to_numpy()
            change=float(net[-1]-net[0])
            if SHOW_ANNOTATIONS:
                ax.text(.04,.97,'Δnet = ' + scientific_label(change, 2),transform=ax.transAxes,
                        ha='left',va='top',fontsize=8)
            summaries.append(dict(system=s,gene=g,n_generations=len(d),
                first_generation=int(x[0]),last_generation=int(x[-1]),
                net_initial=net[0],net_final=net[-1],change_net_AF=change,
                n_mut_min=int(d.n_mut.min()),n_mut_max=int(d.n_mut.max())))
    for j,g in enumerate(GENES):
        p=axes[0,j].get_position()
        fig.text((p.x0+p.x1)/2,.972,g,ha='center',va='center',fontstyle='italic',fontsize=10)
    for i,s in enumerate(SYSTEMS):
        p=axes[i,0].get_position()
        fig.text(.060,(p.y0+p.y1)/2,LABELS[s],rotation=90,ha='center',va='center',fontsize=9,fontweight='bold')
    p=axes[1,0].get_position()
    fig.text(.025,(p.y0+p.y1)/2,'Summed allele frequency (ΣAF)',rotation=90,ha='center',va='center',fontsize=9)
    fig.text(.5,.071 if SHOW_LEGEND else .025,'Generation',ha='center',va='center',fontsize=9)
    if SHOW_LEGEND:
        handles=[Line2D([],[],color=c,marker=m,linestyle=ls,markersize=3.5,linewidth=1.1,
                       label=label) for c,m,ls,label in zip(COLORS,MARKERS,LINESTYLES,CURVE_LABELS)]
        fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.55,.012),
                   ncol=3,frameon=False,fontsize=8,handlelength=2,columnspacing=1.4)
    stem='Panel_AF_PREVIEW_incompleto' if missing else 'Figure_balance_AF'
    export(fig,OUTPUT_DIR/stem)
    pd.DataFrame(summaries).to_csv(OUTPUT_DIR/'resumen_balance_AF.tsv',sep='\t',index=False,na_rep='NA')


def main():
    if Y_SCALE not in ('panel','global'): raise ValueError('Y_SCALE: panel o global.')
    font=configure_font();data,missing=read_data()
    OUTPUT_DIR.mkdir(parents=True,exist_ok=True)
    pd.concat(data.values(),ignore_index=True).to_csv(OUTPUT_DIR/'datos_panel_AF.tsv',sep='\t',index=False)
    make_panel(data,missing)
    (OUTPUT_DIR/'informacion_exportacion.txt').write_text(
        f'Fuente: {font}\nTamaño: {WIDTH_MM} × {HEIGHT_MM} mm\nRGB 8 bits/canal; {DPI} ppp\n'
        f'Escalas Y: {Y_SCALE}\nSeries disponibles: {len(data)}/9\n'
        'Δnet = net_AF_signed final − inicial. Unidad: suma AF, no pp de GC.\n'
        'No se calcula Spearman, valores p, regresiones ni intervalos.\n'
        'Las líneas conectan observaciones; no son regresiones. No se añade G0.\n',encoding='utf-8')
    print('Salida:',OUTPUT_DIR)

if __name__=='__main__': main()
