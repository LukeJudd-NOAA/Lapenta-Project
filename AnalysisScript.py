########################################################################
#                       IMPORTS + VARIABLES                            #
########################################################################
from __future__ import annotations
from pathlib import Path
from scipy import stats
from scipy.signal import get_window
from scipy.spatial import cKDTree
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.weightstats import ttost_paired
import json, pathlib, urllib.request, csv, shutil
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
import os, re, numpy as np, pandas as pd
import pyproj
import xarray as xr

#plotly imports
pio.templates["nos_white"] = pio.templates["plotly_white"].update({
     "layout": {
         "font": {"family": "Times New Roman", "size": 12},
         "xaxis": {"showline": True, "linewidth": 1, "linecolor": "black",
                   "mirror": True, "ticks": "outside", "gridcolor": "#cccccc"},
         "yaxis": {"showline": True, "linewidth": 1, "linecolor": "black",
                   "mirror": True, "ticks": "outside", "gridcolor": "#cccccc"},
         "legend": {"orientation": "h", "yanchor": "bottom",
                    "y": 1.02, "xanchor": "center", "x": 0.5},
         "margin": {"l": 80, "r": 40, "t": 90, "b": 60}
     }
 })
pio.templates.default = "nos_white"
BLUE    = "#1f77b4"       
ORANGE  = "#ff7f0e"        

STYLE_BLUE_SOLID   = dict(line=dict(color=BLUE,   width=2, dash="solid"),
                          mode="lines")
STYLE_ORANGE_DOT   = dict(line=dict(color=ORANGE, width=2, dash="dot"),
                          mode="lines")
STYLE_ORANGE_DASH  = dict(line=dict(color=ORANGE, width=2, dash="dash"),
                          mode="lines")  
#

DATUM_VAR = {
    "MHHW": "mhhwtomsl",
    "MHW" : "mhwtomsl",
    "MTL" : "mtltomsl",
    "MSL" : None,         
    "DTL" : "dtltomsl",
    "MLW" : "mlwtomsl",
    "MLLW": "mllwtomsl",
    "NAVD": "navd88tomsl",
    "STND": "stndtomsl",
}

#overwritten by user input
OFS         = "cbofs"           
OFS_ABBR    = OFS[:2].lower()   
YEAR        = 2024            
DATE_START  = pd.Timestamp("2024-01-01 00:00:00")
DATE_END    = pd.Timestamp("2024-07-01 00:00:00")

VAR_FOLDERS = {
    'CU':   'currents',
    'Salt': 'salinity',
    'Temp': 'water_temperature',
    'WL':   'water_level',
}

OBS_MOD = {
    'currents':          ('obs_speed',         'mod_speed'),
    'water_level':       ('obs_water_depth',   'mod_water_depth'),
    'water_temperature': ('obs_water_temperature', 'mod_water_temperature'),
    'salinity':          ('obs_salinity',      'mod_salinity'),
}

UNIT_FACTOR = {
    'currents':          0.01,
    'water_level':       0.01,
    'water_temperature': 1.0,
    'salinity':          1.0,
}

CSV_DIR, OUTPUT_FILE, ALPHA = ".", "compare.xlsx", 0.05

ABS_MARGINS = {
    ("water_temperature", "bias"): 1.00,    
    ("water_temperature", "rmse"): 1.25,    
    ("water_level", "bias"):       0.12,     
    ("water_level", "rmse"):       0.18,  
}

REL_MARGINS = {
    "rmse":                 0.15,   
    "bias":                 0.15,
    "bias_perc":            0.15,  
    "bias_standard_dev":    0.15,
    "central_freq":         0.05,   
    "pos_outlier_freq":     0.50,   
    "neg_outlier_freq":     0.50,
}
ABS_MARGIN_R = .05

VARS = ["currents","salinity","water_level","water_temperature"]

NUM_METRICS = ["rmse","r","bias","bias_perc","central_freq",
               "pos_outlier_freq","neg_outlier_freq","bias_standard_dev"]

BIN_METRICS = ["central_freq_pass_fail",
               "pos_outlier_freq_pass_fail","neg_outlier_freq_pass_fail"]

ID_RE = re.compile(r"(\d{5,8})")

ROOT      = Path(__file__).resolve().parent
PY_ROOT   = ROOT / "python run"
FR_ROOT   = ROOT / "fortran run"

FR_WL_RAW     = FR_ROOT / "cbofs 2024 WL"  
FR_WL_SHIFTED = FR_ROOT / "wl MLLW shifted"

PLOTS_ROOT        = ROOT / "plots"
BODE_ROOT         = ROOT / "bode"      

OFFSET_CSV    = ROOT / "station_offsets.csv"

VAR_INFO = {
    "wl"  : ( r"wl",    1, "Water-level", "m",  "C0" ),
    "temp": ( r"temp",  1, "Temperature", "°C", "C1" ),
    "salt": ( r"salt",  1, "Salinity",    "PSU","C2" ),
    "cu"  : ( r"cu",    4, "Currents",    "m s⁻¹ / °", None ),
}

RUNTYPES = ("nowcast", "forecast")






########################################################################
#                           USER INPUT HELPERS                         #
########################################################################
def configure(
    ofs: str,
    start_str: str,
    end_str: str,
    datum: str
) -> None:
    global OFS, OFS_ABBR, YEAR, DATE_START, DATE_END
    global DATUM, FR_WL_RAW, FR_WL_SHIFTED, OUTPUT_FILE, VDAT_FILE

    OFS   = ofs.lower()
    DATUM = datum.upper()
    if DATUM not in DATUM_VAR:
        raise ValueError(f"Datum {DATUM} not in {list(DATUM_VAR)}")

    FR_WL_RAW     = FR_ROOT / f"{OFS} {YEAR} WL"
    FR_WL_SHIFTED = FR_ROOT / f"wl_{DATUM.lower()}_shifted"
    FR_WL_SHIFTED.mkdir(parents=True, exist_ok=True)

    OFS_ABBR   = OFS[:2].lower()
    DATE_START = pd.to_datetime(start_str)
    DATE_END   = pd.to_datetime(end_str)
    YEAR       = DATE_START.year

    OUTPUT_FILE   = f"compare_{OFS}_{YEAR}.xlsx"

    VDAT_FILE = ROOT / f"{OFS}_vdatums.nc"
    _fetch_vdatum(OFS, VDAT_FILE)     

def ofs_folder(suffix: str) -> str:
    return f"{OFS} {YEAR} {suffix}"

def ofs_skill_tag(var: str, tag: str) -> str:
    return f"skill_{OFS}_{var}_{tag}.csv"







########################################################################
#                            WL SHIFT MODULE                           #
########################################################################
def _fetch_vdatum(ofs: str, dest: Path) -> None:
    if dest.exists():
        print(f"✓ {dest.name} already present")
        return

    url = f"https://noaa-nos-ofs-pds.s3.amazonaws.com/OFS_Grid_Datum/{ofs}_vdatums.nc"
    print(f"Downloading {url} → {dest} …")
    try:
        urllib.request.urlretrieve(url, dest)
        print("✓ download complete")
    except Exception as e:
        raise RuntimeError(f"Could not fetch vdatum file: {e}") from e
    
def build_station_coords():
    root   = pathlib.Path(__file__).resolve().parent
    fr_csv = root / "fortran_csv"
    py_dir = root / "python run"
    ids    = set()

    for p in fr_csv.glob(f"for_{OFS}_*_nowcast.csv"):
        try:
            ids |= set(pd.read_csv(p, usecols=["ID"])["ID"].astype(str))
        except Exception:
            pass   

    for p in py_dir.glob("**/data/skill/stats/**/*.csv"):
        try:
            ids |= set(pd.read_csv(p, usecols=["ID"])["ID"].astype(str))
        except Exception:
            pass

    id_re = re.compile(r"(\d{7,8})")
    for dat in (root / "fortran run").rglob("*_obsandnow.dat"):
        m = id_re.search(dat.name)
        if m:
            ids.add(m.group(1))

    if not ids:
        print("Still didn’t find any station IDs – double-check your CSV paths.")
        return

    print(f"Found {len(ids)} unique station IDs – querying CO-OPS …")
    rows, missing = [], []
    def pick(d, *keys):
        for k in keys:
            if k in d:
                return float(d[k])
        raise KeyError

    for st in sorted(ids):
        url = f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{st}.json"
        try:
            meta = json.load(urllib.request.urlopen(url))["stations"][0]
            lon  = pick(meta, "lng", "longitude")
            lat  = pick(meta, "lat", "latitude")
            rows.append((st, lon, lat))
            print(f"✓ {st}  lon={lon:.4f}  lat={lat:.4f}")
        except Exception as e:
            print(f"✗ {st} ({type(e).__name__}) – add manually if needed")
            missing.append(st)

    out_csv = root / "station_coords.csv"
    with open(out_csv, "w", newline="") as f:
        csv.writer(f).writerows([("id", "lon", "lat"), *rows])

    print(f"\nWrote {len(rows)} coordinates to {out_csv.name}.")
    if missing:
        print("Missing stations:", ", ".join(missing))

def make_station_offsets() -> None:
    if DATUM == "MSL":
        Path("station_offsets.csv").write_text("id,offset_m\n")
        print("Datum = MSL → no offsets necessary.")
        return

    var_name = DATUM_VAR[DATUM]                 
    if var_name is None:
        raise RuntimeError("Internal error: DATUM_VAR mapping missing!")

    VDAT_FILE = ROOT / f"{OFS}_vdatums.nc"
    STN_FILE  = ROOT / "station_coords.csv"
    OUT_FILE  = ROOT / "station_offsets.csv"
    FILLVALUE = -9999.9

    print(f"[offsets] opening {VDAT_FILE.name} …")
    ds        = xr.open_dataset(VDAT_FILE)
    offs      = ds[var_name]                      
    lon_grid  = ds["longitude"]
    lat_grid  = ds["latitude"]

    print("[offsets] building KD-tree of wet cells …")
    wet_mask        = offs != FILLVALUE
    j_idx, i_idx    = np.where(wet_mask)
    proj            = pyproj.Proj("epsg:3857")
    gx, gy          = proj(lon_grid.values[j_idx, i_idx],
                           lat_grid.values[j_idx, i_idx])
    kdt             = cKDTree(np.column_stack([gx, gy]))

    stns            = pd.read_csv(STN_FILE)
    rows            = []

    for sid, slon, slat in stns.itertuples(index=False):
        sx, sy        = proj(slon, slat)
        dist, idx     = kdt.query([sx, sy])
        j, i          = int(j_idx[idx]), int(i_idx[idx])
        off           = float(offs[j, i])

        if off == FILLVALUE or np.isnan(off):
            print(f"! {sid}: no wet grid cell nearby – skipped")
            continue

        rows.append({"id": str(sid), "offset_m": off})
        print(f"✓ {sid}: offset {off:+.3f} m  (distance {dist/1000:.1f} km)")

    pd.DataFrame(rows).to_csv(OUT_FILE, index=False)
    print(f"✓ {OUT_FILE.name} written  ({len(rows)} stations)")

def load_station_offsets() -> dict[str,float]:
    df = pd.read_csv(OFFSET_CSV)
    return dict(zip(df.iloc[:,0].astype(str), df["offset_m"].astype(float)))

#Main shift function, applies to fortran wl time series
def shift_fortran_wl(offsets: dict[str, float]) -> None:
    if DATUM == "MSL":
        print("WL already referenced to MSL – skipping datum shift.")
        return

    if FR_WL_SHIFTED.exists():
        shutil.rmtree(FR_WL_SHIFTED)
    FR_WL_SHIFTED.mkdir(parents=True, exist_ok=True)

    pat_ok = re.compile(rf"_{OFS_ABBR.upper()}_(nowcast|forecast)\.dat$", re.I)
    pat_skip = re.compile(r"(timeseries|obsand|alldata|table|extreme)", re.I)

    for fp in FR_WL_RAW.glob("*.dat"):
        name = fp.name
        if pat_skip.search(name) or not pat_ok.search(name):
            continue

        station, _, _ = extract_station_var_run(name, "fortran")
        offset = offsets.get(station)
        if offset is None or np.isnan(offset):
            print(f"[shift] offset missing for station {station} – skipped")
            continue

        suffix   = f"_{DATUM.lower()}.dat"
        out_path = FR_WL_SHIFTED / name.replace(".dat", suffix)
        with fp.open() as fin, out_path.open("w") as fout:
            for line in fin:
                toks = line.strip().split()
                if len(toks) < 7:
                    continue
                try:
                    elev = float(toks[-1]) - offset   
                    toks[-1] = f"{elev:.5f}"
                    fout.write(" ".join(toks) + "\n")
                except ValueError:
                    fout.write(line)
        print(f"[shift] {name} → {out_path.name}  (offset {offset:+.3f} m)")



########################################################################
#                        CONVERT FORTRAN TABLES                        #
########################################################################
#grab the metrics row (U-u, H-h …)
def parse_fortran_table(lines, agg_tag, fc_tag, uf, var):
    now_m = {}
    fc_m  = {}
    for L in lines:
        tok = L.split()
        if not tok:
            continue
        tag = tok[0]
        if tag not in (agg_tag, fc_tag):
            continue
        try:
            thr = float(tok[1]) * uf
        except Exception:
            thr = ''
        vals = tok[5:]
        nums = []
        for v in vals:
            try:
                nums.append(float(v))
            except Exception:
                pass

        if var == 'water_level':
            # WL has an extra WOF column at nums[8]
            if len(nums) >= 10:
                metrics = {
                    'threshold': thr,
                    'SM':   nums[0],
                    'RMSE': nums[1],
                    'SD':   nums[2],
                    'NOF':  nums[3],
                    'CF':   nums[4],
                    'POF':  nums[5],
                    'MDNO': nums[6],
                    'MDPO': nums[7],
                    'CORR': nums[9],
                }
            else:
                metrics = {}
        else:
            if len(nums) >= 9:
                metrics = {
                    'threshold': thr,
                    'SM':   nums[0],
                    'RMSE': nums[1],
                    'SD':   nums[2],
                    'NOF':  nums[3],
                    'CF':   nums[4],
                    'POF':  nums[5],
                    'MDNO': nums[6],
                    'MDPO': nums[7],
                    'CORR': nums[8],
                }
            else:
                metrics = {}

        if tag == agg_tag:
            now_m = metrics
        else:
            fc_m = metrics
            break

    return now_m, fc_m

def parse_phase_bias(lines):
    bias_now = bias_fc = ''
    for L in lines:
        parts = L.split()
        if len(parts) < 6:
            continue
        tag = parts[0]
        if tag == 'D-d':
            bias_now = float(parts[5])
        if tag == 'D006-d006':          
            bias_fc = float(parts[5])
    return bias_now, bias_fc

#grab the observation mean (row 'u', 'h', 't', 's')
def parse_obs_mean(lines, obs_tag):
    for L in lines:
        tok = L.split()
        if tok and tok[0] == obs_tag and len(tok) >= 3:
            try:
                return float(tok[2])
            except Exception:
                pass
    return ''  

#Main conversion function, generates csvs identical to python formatting
def convert_var(var, out_dir):
    suffix  = next(s for s, v in VAR_FOLDERS.items() if v == var)
    folder  = ofs_folder(suffix)

    agg_tag = {'currents': 'U-u',  'water_level': 'H-h',
               'water_temperature': 'T-t', 'salinity': 'S-s'}[var]

    fc_tag  = {'currents': 'U006-u006',  'water_level': 'H006-h006',
               'water_temperature': 'T006-t006', 'salinity': 'S006-s006'}[var]

    obs_tag = {'currents': 'u', 'water_level': 'h',
               'water_temperature': 't', 'salinity': 's'}[var]

    uf      = UNIT_FACTOR[var]
    obs_c, mod_c = OBS_MOD[var]

    now_rows = []
    fc_rows  = []

    files  = os.listdir(folder)
    tables = sorted(f for f in files if f.endswith("_table.out") and "phase" not in f)
    phases = {f for f in files if f.endswith("phase_table.out")}

    for tbl in tables:
        stn   = tbl.replace("_table.out", "")
        path  = os.path.join(folder, tbl)
        lines = open(path).read().splitlines()

        now_m, fc_m = parse_fortran_table(lines, agg_tag, fc_tag, uf, var)

        obs_mean = parse_obs_mean(lines, obs_tag)

        bias_now = bias_fc = ''
        if var == 'currents':
            ph = tbl.replace("_table.out", "phase_table.out")
            if ph in phases:
                plines = open(os.path.join(folder, ph)).read().splitlines()
                bias_now, bias_fc = parse_phase_bias(plines)

        template = {
            'ID': stn, 'NODE': '',
            obs_c: '', mod_c: '',
            'datum': '', 'Y': '', 'X': '',
            'start_date': '', 'end_date': ''
        }

        nc = dict(template)
        nc_bias = now_m.get('SM', '')
        nc.update({
            'rmse':  now_m.get('RMSE', ''),
            'r':     now_m.get('CORR', ''),
            'bias':  nc_bias,
            'bias_perc': '' if obs_mean in ('', 0) or nc_bias == '' else (nc_bias / obs_mean) * 100,
            'bias_dir': bias_now,
            'central_freq': now_m.get('CF', ''),
            'central_freq_pass_fail': 'pass' if now_m.get('CF', 0) >= 90 else 'fail',
            'pos_outlier_freq': now_m.get('POF', ''),
            'pos_outlier_freq_pass_fail': 'pass' if now_m.get('POF', 0) <= 1 else 'fail',
            'neg_outlier_freq': now_m.get('NOF', ''),
            'neg_outlier_freq_pass_fail': 'pass' if now_m.get('NOF', 0) <= 1 else 'fail',
            'bias_standard_dev': now_m.get('SD', ''),
            'target_error_range': now_m.get('threshold', ''),
        })
        now_rows.append(nc)

        fc = dict(template)
        fc_bias = fc_m.get('SM', '')
        fc.update({
            'rmse':  fc_m.get('RMSE', ''),
            'r':     fc_m.get('CORR', ''),
            'bias':  fc_bias,
            'bias_perc': '' if obs_mean in ('', 0) or fc_bias == '' else (fc_bias / obs_mean) * 100,
            'bias_dir': bias_fc,
            'central_freq': fc_m.get('CF', ''),
            'central_freq_pass_fail': 'pass' if fc_m.get('CF', 0) >= 90 else 'fail',
            'pos_outlier_freq': fc_m.get('POF', ''),
            'pos_outlier_freq_pass_fail': 'pass' if fc_m.get('POF', 0) <= 1 else 'fail',
            'neg_outlier_freq': fc_m.get('NOF', ''),
            'neg_outlier_freq_pass_fail': 'pass' if fc_m.get('NOF', 0) <= 1 else 'fail',
            'bias_standard_dev': fc_m.get('SD', ''),
            'target_error_range': fc_m.get('threshold', ''),
        })
        fc_rows.append(fc)

    cols = [
        'ID', 'NODE', obs_c, mod_c,
        'rmse', 'r', 'bias', 'bias_perc', 'bias_dir',
        'central_freq', 'central_freq_pass_fail',
        'pos_outlier_freq', 'pos_outlier_freq_pass_fail',
        'neg_outlier_freq', 'neg_outlier_freq_pass_fail',
        'bias_standard_dev', 'target_error_range',
        'datum', 'Y', 'X', 'start_date', 'end_date'
    ]

    now_name = f"for_{OFS}_{var}_nowcast.csv"
    fc_name  = f"for_{OFS}_{var}_forecast_b.csv"

    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame(now_rows, columns=cols)\
        .to_csv(os.path.join(out_dir, now_name), index=False)
    pd.DataFrame(fc_rows, columns=cols)\
        .to_csv(os.path.join(out_dir, fc_name), index=False)

    print(f"[{var}] → csv/{now_name} & csv/{fc_name}")


########################################################################
#                        COMPARISON MODULE                             #
########################################################################
def csv_path(name:str)->Path:
    root = Path(__file__).resolve().parent
    direct = root / name
    if direct.exists():
        return direct
    
    if name.startswith(f"for_{OFS}"):
        p = root / "fortran_csv" / name
        if p.exists():
            return p
        
    m = re.match(rf'^(?:pyt_|skill_){OFS}_(.+?)_(.+?)\.csv$', name, re.I)
    if m:
        var, tag = m.group(1), m.group(2)
        for scen in ("nowcast","forecast"):
            p = root / "python run" / scen / "data" / "skill" / "stats" / ofs_skill_tag(var, tag)
            if p.exists():
                return p
    
    for p in (root / "python run").rglob(name):
        return p
    raise FileNotFoundError(name)

def norm_id(raw:str)->str:
    m=ID_RE.match(str(raw).strip())
    return m.group(1) if m else str(raw)

def extract_tost_pvals(res):
    if hasattr(res,"pvalue_low"):
        return float(res.pvalue_low),float(res.pvalue_high)
    tup=tuple(res)
    if len(tup)==4: _,plo,_,phi=tup
    elif len(tup)==2: plo,phi=tup
    elif len(tup)==3 and isinstance(tup[1],(tuple,list,np.ndarray)):
        plo,phi=tup[1][1],tup[2][1]
    else:
        raise ValueError("Unknown ttost return format")
    return float(plo),float(phi)

def analyse(var:str,scenario:str):
    tag="nowcast" if scenario=="nowcast" else "forecast_b"
    df_ft = pd.read_csv(csv_path(f"for_{OFS}_{var}_{tag}.csv"), dtype=str)
    df_py = pd.read_csv(csv_path(f"pyt_{OFS}_{var}_{tag}.csv"), dtype=str)
    if df_py.columns[0].startswith("Unnamed"): df_py=df_py.iloc[:,1:]
    df_ft["sid"]=df_ft["ID"].map(norm_id); df_py["sid"]=df_py["ID"].map(norm_id)
    df=pd.merge(df_ft,df_py,on="sid",suffixes=("_ft","_py"))

    tableA,tableB=[],[]
    mean_pass=mean_tot=stats_pass=stats_tot=bin_pass=bin_tot=0

    for m in NUM_METRICS:
        a, b = f"{m}_ft", f"{m}_py"
        if a not in df or b not in df:
            continue

        x = pd.to_numeric(df[a], errors="coerce")
        y = pd.to_numeric(df[b], errors="coerce")
        ok = x.notna() & y.notna()
        if ok.sum() == 0:         
            continue
        elif ok.sum() == 1:        #stats tests need 2+ stations
            sid = df.loc[ok, "sid"].iloc[0]

            margin = ABS_MARGINS.get((var, m),
                                      ABS_MARGIN_R if m == "r"
                                      else REL_MARGINS[m] * x.iloc[0])

            pass_ratio = float(abs(y.iloc[0] - x.iloc[0]) <= margin)
            eq_mean    = "YES" if pass_ratio == 1 else "NO"

            st_cols = {f"{sid}_ft": float(x.iloc[0]),
                       f"{sid}_py": float(y.iloc[0])}

            tableA.append(dict(variable=var, metric=m, n=1, **st_cols,
                               pass_ratio=pass_ratio,
                               overall_equiv=eq_mean,
                               mean_ft=float(x.iloc[0]),
                               mean_py=float(y.iloc[0])))

            tableB.append(dict(variable=var, metric=m, n=1,
                               p_ttest=np.nan, p_wilcoxon=np.nan, p_ks=np.nan,
                               cohen_dz=np.nan, p_tost=np.nan,
                               equivalent="Not Enough Data",
                               agreement=np.nan, p_mcnemar=np.nan))
            
            mean_tot  += 1
            mean_pass += (eq_mean == "YES")
            continue

        x, y, sids = x[ok], y[ok], df.loc[ok, "sid"]

        mean_ft, mean_py = x.mean(), y.mean()
        p_tt  = stats.ttest_rel(y, x).pvalue
        try:
            p_wl = stats.wilcoxon(y - x).pvalue
        except ValueError:
            p_wl = 1.0
        p_ks  = stats.ks_2samp(x, y).pvalue
        d_z   = (y - x).mean() / (y - x).std(ddof=1)

        margin = ABS_MARGINS.get((var, m),
                                 ABS_MARGIN_R if m == "r"
                                 else REL_MARGINS[m] * mean_ft)
        p_lo, p_hi = extract_tost_pvals(
            ttost_paired(x, y, -margin, margin)
        )
        p_tost = max(p_lo, p_hi)         

        b_tost   = p_tost <  ALPHA         
        b_ttest  = p_tt   >  ALPHA       
        b_wilcox = p_wl   >  ALPHA
        b_ks     = p_ks   >  ALPHA
        b_dz     = abs(d_z) < 0.3     

        score = (3 * b_tost) + b_ttest + b_wilcox + b_ks + b_dz
        eq_stats = "YES" if score >= 4 else "NO"   

        stats_tot += 1
        stats_pass += (eq_stats == "YES")

        pass_ratio = (np.abs(y - x) <= margin).mean()
        eq_mean = "YES" if pass_ratio >= 0.5 else "NO"
        mean_tot += 1
        mean_pass += (eq_mean == "YES")

        st_cols = {f"{sid}_ft": float(xf) for sid, xf in zip(sids, x)}
        st_cols.update({f"{sid}_py": float(yf) for sid, yf in zip(sids, y)})

        tableA.append(dict(
            variable=var, metric=m, n=len(x), **st_cols,
            pass_ratio=pass_ratio, overall_equiv=eq_mean,
            mean_ft=mean_ft, mean_py=mean_py
        ))
        tableB.append(dict(
            variable=var, metric=m, n=len(x),
            p_ttest=p_tt, p_wilcoxon=p_wl, p_ks=p_ks,
            cohen_dz=d_z, p_tost=p_tost, equivalent=eq_stats,
            agreement=np.nan, p_mcnemar=np.nan
        ))

    for m in BIN_METRICS:
        a,b=f"{m}_ft",f"{m}_py"
        if a not in df or b not in df: continue
        ft_pass=df[a].str.lower()=="pass"; py_pass=df[b].str.lower()=="pass"
        agree=(ft_pass==py_pass).mean()
        tbl=[[((~ft_pass)&(~py_pass)).sum(),((~ft_pass)&py_pass).sum()],
             [(ft_pass&(~py_pass)).sum(),(ft_pass&py_pass).sum()]]
        try: pm=mcnemar(tbl,exact=True).pvalue
        except Exception: pm=np.nan
        eq_bin="YES" if (not np.isnan(pm) and pm>ALPHA) else "NO"
        bin_tot+=1; bin_pass+=(eq_bin=="YES")
        tableB.append(dict(variable=var,metric=m,n=len(df),
                           p_ttest=np.nan,p_wilcoxon=np.nan,p_ks=np.nan,
                           cohen_dz=np.nan,p_tost=np.nan,equivalent=eq_bin,
                           agreement=agree,p_mcnemar=pm))

    dfA=pd.DataFrame(tableA)
    sid_cols=sorted({c[:-3] for c in dfA.columns if c.endswith("_ft")})
    ordered=[f"{sid}_{sfx}" for sid in sid_cols for sfx in ("ft","py")]
    dfA=dfA.reindex(columns=["variable","metric","n"]+ordered+
                    ["pass_ratio","overall_equiv","mean_ft","mean_py"])
    dfB=pd.DataFrame(tableB)
    return dfA,dfB,mean_pass,mean_tot,stats_pass,stats_tot,bin_pass,bin_tot

def add_accuracy(df, pcol, tcol, acol):
    df[acol]=df[pcol]/df[tcol].replace(0,np.nan)
    return df

def build_summary(base, pcol, tcol, acol, label):
    base = base.copy()
    base = add_accuracy(base, pcol, tcol, acol)

    scen_tot = base.groupby("scenario")[[pcol,tcol]].sum().reset_index()
    scen_tot["variable"] = "ALL"
    scen_tot = add_accuracy(scen_tot, pcol, tcol, acol)

    grand_pass = base[pcol].sum()
    grand_tot  = base[tcol].sum()
    grand_acc  = grand_pass / grand_tot if grand_tot else np.nan

    total_row = pd.DataFrame({
        "scenario":[""],
        "variable":[label],
        pcol:[grand_pass],
        tcol:[grand_tot],
        acol:[round(grand_acc,3)]
    })

    df = pd.concat([base, scen_tot, total_row], ignore_index=True)
    return df

#Main comparison function, builds csv with summary statistics
def comparison_main():
    writer=pd.ExcelWriter(OUTPUT_FILE,engine="xlsxwriter")
    mean_rows,stats_rows,bin_rows=[],[],[]

    for scen in ("nowcast","forecast"):
        ptr=0
        for var in VARS:
            dfA,dfB,mp,mt,sp,st,bp,bt = analyse(var,scen)
            dfA.to_excel(writer,scen,index=False,startrow=ptr); ptr+=len(dfA)+1
            dfB.to_excel(writer,scen,index=False,startrow=ptr); ptr+=len(dfB)+2
            mean_rows.append(dict(scenario=scen,variable=var,
                                  mean_pass=mp,mean_total=mt))
            stats_rows.append(dict(scenario=scen,variable=var,
                                   stats_pass=sp,stats_total=st))
            bin_rows.append(dict(scenario=scen,variable=var,
                                 bin_pass=bp,bin_total=bt))

    mean_base  = pd.DataFrame(mean_rows)
    stats_base = pd.DataFrame(stats_rows)
    bin_base   = pd.DataFrame(bin_rows)

    df_mean  = build_summary(mean_base ,"mean_pass" ,"mean_total" ,"mean_accuracy" ,"Total Mean Accuracy:")
    df_stats = build_summary(stats_base,"stats_pass","stats_total","stats_accuracy","Total Stats Accuracy:")
    df_bin   = build_summary(bin_base ,"bin_pass" ,"bin_total" ,"bin_accuracy" ,"Total Bin Accuracy:")

    total_accuracy = (mean_base.mean_pass.sum()+stats_base.stats_pass.sum()+bin_base.bin_pass.sum()) / (
                     mean_base.mean_total.sum()+stats_base.stats_total.sum()+bin_base.bin_total.sum())
    overall_row = pd.DataFrame({"scenario":[""],"variable":["Total Accuracy:"],"mean_accuracy":[round(total_accuracy,3)]})

    start=0
    for df in (df_mean,df_stats,df_bin):
        df.to_excel(writer,"summary",index=False,startrow=start); start+=len(df)+2
    overall_row.to_excel(writer,"summary",index=False,startrow=start)

    writer.close(); print(f"✓ {OUTPUT_FILE} written")





########################################################################
#                          PLOTTERS MODULE                             #
########################################################################
def julian_to_datetime(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, unit="D", origin="julian")

def doyfrac_to_datetime(year: pd.Series, doy: pd.Series) -> pd.Series:
    start_year = pd.to_datetime(year.astype(int).astype(str) + "-01-01")
    days       = doy - 1.0
    return start_year + pd.to_timedelta(days, unit="D")

def parse_fixed_width(
        path: Path, n_val: int, use_julian: bool
) -> pd.DataFrame:
    expected = 6 + n_val                        

    try:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            comment="#",
            engine="python",
            na_values=["nan", "NaN", "NA", ""],
            dtype=float,
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    df = df[df.columns[:expected]]
    df = df.dropna(how="all")

    if df.shape[1] < expected:
        raise ValueError(f"{path.name}: expected ≥{expected} numeric columns "
                         f"but found {df.shape[1]}")

    df.columns = (["col0", "year", "mo", "dy", "hr", "mn"] +
                  [f"v{i}" for i in range(n_val)])
    
    if use_julian:
        dt = julian_to_datetime(df["col0"])
    else:
        dt = doyfrac_to_datetime(df["year"], df["col0"])
    df.index = dt.dt.tz_localize("UTC")

    return df[[c for c in df.columns if c.startswith("v")]]

def extract_station_var_run(filename: str, source: str) -> tuple[str, str, str]:
    bn  = filename.lower()
    run = "forecast" if "forecast" in bn else "nowcast"

    m = re.search(r"(\d{5,7})", bn)                          
    if m:
        station = m.group(1)
    else:
        m = re.search(rf"{OFS_ABBR}(\d{{4,5}})", bn)      
        station = m.group(1) if m else "unknown"

    m = re.match(rf"\d{{5,7}}([wst])_{OFS_ABBR}_", bn, re.I)
    if m:
        code = m.group(1).lower()
        return station, {"w": "wl", "s": "salt", "t": "temp"}[code], run

    if re.search(rf"^{OFS_ABBR}\d{{4,5}}_(speed|dir)_", bn):
        return station, "cu", run
    
    if re.search(rf"_{OFS}_cu_", bn, re.I) or re.search(r"_cu_\d+_", bn):
        return station, "cu", run
    
    if re.match(rf"{OFS_ABBR}\d{{4,5}}_", bn):
        return station, "cu", run

    for vk, (frag, *_rest) in VAR_INFO.items():
        if vk == "cu":
            continue
        if re.search(rf"_{frag}(?:_|[0-9])", bn):
            return station, vk, run

    raise ValueError(f"Cannot identify variable from {filename!s}")

#Main function to locate all time series
def build_file_index() -> pd.DataFrame:
    recs: list[dict] = []

    print("\n=== BUILDING FILE INDEX ===")

    for run in RUNTYPES:
        model_dir = PY_ROOT / run / "data" / "model" / "1d_node"
        obs_dir   = PY_ROOT / run / "data" / "observations" / "1d_station"

        print(f"\nScanning Python {run} directories:")
        print(f"  Model dir: {model_dir}")
        print(f"  Obs dir:   {obs_dir}")

        if model_dir.exists():
            prd_files = list(model_dir.glob("*.prd"))
            print(f"  Found {len(prd_files)} .prd files")
            for fp in prd_files:
                try:
                    station, var, _ = extract_station_var_run(fp.name, "python_model")
                    recs.append(dict(
                        src="py_model",
                        run=run,
                        var=var,
                        station=station,
                        path=fp
                    ))
                except ValueError as e:
                    print(f"    ERROR parsing {fp.name}: {e}")
        else:
            print("  Model dir does not exist!")

        if obs_dir.exists():
            obs_files = list(obs_dir.glob("*.obs"))
            print(f"  Found {len(obs_files)} .obs files")
            for fp in obs_files:
                try:
                    station, var, _ = extract_station_var_run(fp.name, "obs")
                    recs.append(dict(
                        src="obs",
                        run=run,
                        var=var,
                        station=station,
                        path=fp
                    ))
                except ValueError as e:
                    print(f"    ERROR parsing {fp.name}: {e}")
        else:
            print("  Obs dir does not exist!")

    print("\nScanning Fortran directories:")

    pat_ok = re.compile(rf"_{OFS_ABBR.upper()}_(nowcast|forecast)(?:_{DATUM.lower()})?\.dat$", re.I)
    pat_skip = re.compile(r"(timeseries|obsand|alldata|table|extreme)", re.I)

    folders = {
    "cu"  : FR_ROOT / ofs_folder("CU"),
    "salt": FR_ROOT / ofs_folder("Salt"),
    "temp": FR_ROOT / ofs_folder("Temp"),
    "wl"  : FR_WL_SHIFTED,
}

    for var_key, folder in folders.items():
        print(f"\n  {var_key.upper()} folder: {folder}")
        if not folder.exists():
            print("    Folder does not exist!")
            continue

        dat_files = list(folder.glob("*.dat"))
        print(f"    Found {len(dat_files)} .dat files")

        for fp in dat_files:
            name = fp.name

            if pat_skip.search(name):
                continue

            if var_key == "cu":
                if not re.search(r"(nowcast|forecast)\.dat$", name, re.I):
                    print(f"    SKIPPED CU: {name} (no nowcast/forecast)")
                    continue
            else:
                if not pat_ok.search(name):
                    continue

            try:
                station, _, run_tag = extract_station_var_run(name, "fortran")
            except ValueError as e:
                print(f"    ERROR parsing Fortran file {name}: {e}")
                continue

            recs.append(dict(
                src="fr_model",
                run=run_tag,
                var=var_key,
                station=station,
                path=fp
            ))

            if var_key == "cu":
                print(f"    CU FR: {name} → station={station}, run={run_tag}")

    df = pd.DataFrame.from_records(recs)

    print("\n=== FILE INDEX SUMMARY ===")
    print(f"Total files found: {len(df)}")
    for var in VAR_INFO.keys():
        count = len(df[df["var"] == var])
        print(f"  {var.upper()}: {count} files")
        if var == "cu":
            src_counts = df[df["var"] == "cu"]["src"].value_counts().to_dict()
            print(f"    Sources: {src_counts}")

    return df

#Plotly helpers
def _save_plotly(fig: go.Figure, out_path: Path) -> None:
    html_path = out_path.with_suffix(".html")
    fig.write_html(html_path, include_plotlyjs="cdn", full_html=True)
    try:
        fig.write_image(out_path.with_suffix(".png"), scale=2)
    except Exception:
        pass  
    print(f"  → {html_path.name}")

def _plot_triplet_plotly(
        station: str, run: str, tag: str, units: str,
        py_ser: pd.Series | None,
        fr_ser: pd.Series | None,
        obs_ser: pd.Series | None) -> None:
    
    if any(s is None for s in (py_ser, fr_ser, obs_ser)):
        return

    fig = make_subplots(
        rows=3, cols=1, shared_xaxes=True,
        vertical_spacing=0.10,                 
        subplot_titles=("Fortran vs Obs",
                        "Python vs Obs",
                        "Python vs Fortran")
    )

    def _add(row: int, y, style: dict, name: str):
        fig.add_trace(
            go.Scatter(x=y.index, y=y, name=name,
                       legend=f"legend{row}" if row > 1 else None,
                       showlegend=True, **style),
            row=row, col=1
        )

    _add(1, fr_ser,  STYLE_BLUE_SOLID,  "Fortran model")
    _add(1, obs_ser, STYLE_ORANGE_DOT,  "Obs")

    _add(2, py_ser,  STYLE_BLUE_SOLID,  "Python model")
    _add(2, obs_ser, STYLE_ORANGE_DOT,  "Obs")

    _add(3, py_ser,  STYLE_BLUE_SOLID,  "Python model")
    _add(3, fr_ser,  STYLE_ORANGE_DOT,  "Fortran model")

    for r in (1, 2, 3):
        fig.update_yaxes(title_text=units, row=r, col=1)
    fig.update_xaxes(title_text="UTC time", title_standoff=45, row=3, col=1)

    fig.update_layout(
        height=1150,
        title=dict(text=f"{station} — {tag} ({run})", x=0.5, y=0.99),
        margin=dict(l=80, r=40, t=110, b=180),
        legend=dict(  
            orientation="h", x=0.5, xanchor="center",
            y=0.72, yanchor="top", bordercolor="black", borderwidth=1),
        legend2=dict( 
            orientation="h", x=0.5, xanchor="center",
            y=0.35, yanchor="top", bordercolor="black", borderwidth=1),
        legend3=dict( 
            orientation="h", x=0.5, xanchor="center",
            y=-.03, yanchor="top", bordercolor="black", borderwidth=1)
    )

    out_dir = PLOTS_ROOT / run / tag
    out_dir.mkdir(parents=True, exist_ok=True)
    _save_plotly(fig, out_dir / f"{station}_{tag}_{run}")


#forecast horizon sorting helper
def _clean_series(ser: pd.Series, src: str) -> pd.Series:
    ser = ser.sort_index(kind="stable")          

    if src == "fr_model":
        starts = ser.index.to_series().diff().le(pd.Timedelta(0))
        ser = ser[~starts]
        ser = ser[~ser.index.duplicated(keep="last")]
        ser = ser.groupby(ser.index.floor("6min"), sort=False).last()

    elif src == "obs":
        ser = ser[~ser.index.duplicated(keep="first")]
    return ser


def make_plots(df_index: pd.DataFrame) -> None:
    print(f"\n=== MAKING PLOTS ===")
        
    def _load(paths: dict, src: str, var_key: str) -> pd.Series | None:
        if src not in paths:
            return None

        n_val      = VAR_INFO[var_key][1]
        use_julian = (src != "fr_model")                
        df = parse_fixed_width(paths[src], n_val, use_julian)
        if df.empty or "v0" not in df.columns:
            return None

        ser = df["v0"]

        if ser.index.tz is not None:
            ser.index = ser.index.tz_convert("UTC").tz_localize(None)

        if not (src == "fr_model" and run == "nowcast"):
            ser = _clean_series(ser, src)


        mask = (ser.index >= DATE_START) & (ser.index < DATE_END)
        if not mask.any():
            return None
        return ser.loc[mask]


    def _load_pair(paths: dict, src: str, idx: int) -> pd.Series | None:
        if src not in paths:
            return None

        use_julian = (src != "fr_model")
        df = parse_fixed_width(paths[src], VAR_INFO["cu"][1], use_julian)
        if df.empty or f"v{idx}" not in df.columns:
            return None

        ser = df[f"v{idx}"]

        if ser.index.tz is not None:
            ser.index = ser.index.tz_convert("UTC").tz_localize(None)

        if not (src == "fr_model" and run == "nowcast"):
            ser = _clean_series(ser, src)

        mask = (ser.index >= DATE_START) & (ser.index < DATE_END)
        if not mask.any():
            return None
        return ser.loc[mask]


    for run in RUNTYPES:
        for var_key in VAR_INFO:
            sub = df_index[(df_index["run"] == run) & (df_index["var"] == var_key)]
            
            if len(sub) == 0:
                print(f"[plot] {run:<8} {var_key.upper():<4} - NO FILES FOUND")
                continue

            for station in sorted(sub["station"].unique()):
                print(f"[plot] {run:<8} {var_key.upper():<4} {station}")

                paths = {
                    r["src"]: r["path"]
                    for _, r in sub[sub["station"] == station].iterrows()
                }

                if var_key == "cu":        
                    todos = [
                        ("CU_speed", "m s⁻¹",
                         _load_pair(paths, "py_model", 0),
                         _load_pair(paths, "fr_model", 0),
                         _load_pair(paths, "obs",      0)),
                        ("CU_dir", "°",
                         _load_pair(paths, "py_model", 1),
                         _load_pair(paths, "fr_model", 1),
                         _load_pair(paths, "obs",      1)),
                    ]
                else:                      
                    todos = [
                        (var_key.upper(), VAR_INFO[var_key][3],
                         _load(paths, "py_model", var_key),
                         _load(paths, "fr_model", var_key),
                         _load(paths, "obs",      var_key)),
                    ]

                for tag, units, py_ser, fr_ser, obs_ser in todos:
                    _plot_triplet_plotly(station, run, tag, units,
                                  py_ser, fr_ser, obs_ser)


#bode helpers
def _log_bin(freq, y, bins_per_dec: int = 20):
    mask = freq > 0
    f, y = freq[mask], y[mask]

    logf     = np.log10(f)
    n_bins   = int((logf.max() - logf.min()) * bins_per_dec)
    edges    = np.linspace(logf.min(), logf.max(), n_bins)
    idxs     = np.digitize(logf, edges)

    fout, yout = [], []
    for i in range(1, len(edges)):
        sel = idxs == i
        if sel.any():
            fout.append(f[sel].mean())
            yout.append(y[sel].mean())
    return np.asarray(fout), np.asarray(yout)

def _station_bode_traces(obs_fft, model_fft, freqs, name, style):
    gain  = np.abs(model_fft) / np.abs(obs_fft)
    phase = np.angle(model_fft) - np.angle(obs_fft)

    f_sm, g_sm = _log_bin(freqs, 20 * np.log10(gain))
    _,   p_sm  = _log_bin(freqs, np.rad2deg(phase))

    period = 1.0 / f_sm / 3600.0          
    ord_   = np.argsort(period)           
    period, g_sm, p_sm = period[ord_], g_sm[ord_], p_sm[ord_]

    if name == "Fortran":
        g_style = STYLE_BLUE_SOLID            
        p_style = STYLE_BLUE_SOLID
    else:                                   
        g_style = STYLE_ORANGE_DASH          
        p_style = STYLE_ORANGE_DASH

    tr_gain = go.Scatter(x=period, y=g_sm,
                         name=f"{name} – gain",
                         showlegend=True, **g_style)
    tr_phase = go.Scatter(x=period, y=p_sm,
                          name=f"{name} – phase",
                          showlegend=True, **p_style)
    return tr_gain, tr_phase

#Main bode function for plotting time series
def make_bode_all(df_index: pd.DataFrame) -> None:
    root_out = BODE_ROOT
    if root_out.exists():
        shutil.rmtree(root_out)
    root_out.mkdir(parents=True, exist_ok=True)

    print(f"\n=== MAKING BODE PLOTS ")

    def _loader(path: Path, n_val: int, src: str) -> pd.Series | None:                           
        use_julian = (src != "fr_model")
        df = parse_fixed_width(path, n_val, use_julian)

        if df.empty or "v0" not in df.columns:
            return None
        ser = df["v0"]

        ser = ser[~ser.index.duplicated(keep="first")]
        if ser.index.tz is not None:
            ser = ser.tz_convert("UTC").tz_localize(None)
        return ser.sort_index()

    for var_key, info in VAR_INFO.items():
        n_val, nice = info[1], info[2]

        for run in RUNTYPES:
            sub = df_index[(df_index["run"] == run) &
                            (df_index["var"] == var_key)]
            if sub.empty:
                continue

            SUM_obs = SUM_fr = SUM_py = None  
            SUM_fs  = None                     
            stations_used = []                

            for station in sorted(sub["station"].unique()):

                print(f"[bode] "
                    f"{run}/{var_key} — {station}")
            
                paths = {r["src"]: r["path"]
                for _, r in sub[sub["station"] == station].iterrows()}

                if not all(k in paths for k in ("obs", "py_model", "fr_model")):
                    continue

                obs = _loader(paths["obs"],      n_val, "obs")
                py  = _loader(paths["py_model"], n_val, "py_model")
                fr  = _loader(paths["fr_model"], n_val, "fr_model")

                if (obs is None or py is None or fr is None or
                    obs.empty   or py.empty   or fr.empty):
                    continue

                common = obs.index.intersection(py.index).intersection(fr.index)
                if len(common) < 4:
                    continue

                y_obs = obs.loc[common].values
                y_py  =  py.loc[common].values
                y_fr  =  fr.loc[common].values
                N     = len(common)
                dt    = (common[1] - common[0]).total_seconds()
                fs    = 1.0 / dt
                win = get_window("hann", N)
                G_obs = np.fft.rfft(y_obs * win)
                G_fr  = np.fft.rfft(y_fr  * win)
                G_py  = np.fft.rfft(y_py  * win)
                freqs = np.fft.rfftfreq(N, d=1/fs)

                fig = make_subplots(rows=2, cols=1,
                                    shared_xaxes=True,
                                    subplot_titles=("Gain (dB)", "Phase (°)"))

                tr_g_fr, tr_p_fr = _station_bode_traces(G_obs, G_fr, freqs,
                                                        "Fortran",
                                                        {"line": {"dash": "solid"}})
                tr_g_py, tr_p_py = _station_bode_traces(G_obs, G_py, freqs,
                                                        "Python",
                                                        {"line": {"dash": "dash"}})
                
                fig.update_layout(template="nos_white")

                fig.add_trace(tr_g_fr, 1, 1)
                fig.add_trace(tr_g_py, 1, 1)
                fig.add_trace(tr_p_fr, 2, 1)
                fig.add_trace(tr_p_py, 2, 1)

                fig.update_xaxes(title_text="Period (hours)", type="log", row=2, col=1)
                fig.update_yaxes(title_text="dB", row=1, col=1)
                fig.update_yaxes(title_text="Phase °",  row=2, col=1)
                fig.update_layout(
                    height=600,                          
                    title=dict(text=f"{station} • {nice} ({run})", x=0.5, y=0.97),
                    legend=dict(orientation="h",
                                x=0.5, xanchor="center",
                                y=-0.16, yanchor="top",
                    bordercolor="black", borderwidth=1),
                    margin=dict(l=80, r=40, t=90, b=120)
                )
                fig.update_xaxes(type="log", title_text="Period (hours)",
                                showgrid=True, ticks="outside", row=2, col=1)
                fig.update_yaxes(showgrid=True, ticks="outside")

                out_dir = root_out / run / var_key
                out_dir.mkdir(parents=True, exist_ok=True)
                _save_plotly(fig, out_dir / f"{station}_{var_key}_{run}_bode")

                if SUM_obs is None:
                    SUM_obs, SUM_fr, SUM_py = (G_obs.copy(),
                                               G_fr.copy(),
                                               G_py.copy())
                    SUM_fs = fs
                else:
                    min_len = min(len(SUM_obs), len(G_obs))
                    SUM_obs = SUM_obs[:min_len];  SUM_fr  = SUM_fr[:min_len]
                    SUM_py  = SUM_py [:min_len]
                    G_obs   = G_obs [:min_len];   G_fr    = G_fr [:min_len]
                    G_py    = G_py  [:min_len]

                    SUM_obs += G_obs
                    SUM_fr  += G_fr
                    SUM_py  += G_py

                stations_used.append(station)

            if not stations_used:
                print(f"[bode] {nice} ({run}) – no usable stations")
                continue

            print(f"[bode-sum] {run}/{var_key}: "
                  f"{len(stations_used)} stations → "
                  f"{', '.join(stations_used)}")

            N_time = (SUM_obs.size - 1) * 2
            freqs  = np.fft.rfftfreq(N_time, d=1/SUM_fs)

            fig = make_subplots(rows=2, cols=1,
                                shared_xaxes=True,
                                subplot_titles=("Gain (dB)", "Phase (°)"))

            tr_g_fr, tr_p_fr = _station_bode_traces(SUM_obs, SUM_fr, freqs,
                                                    "Fortran",
                                                    {"line": {"dash": "solid"}})
            tr_g_py, tr_p_py = _station_bode_traces(SUM_obs, SUM_py, freqs,
                                                    "Python",
                                                    {"line": {"dash": "dash"}})

            fig.add_trace(tr_g_fr, 1, 1)
            fig.add_trace(tr_g_py, 1, 1)
            fig.add_trace(tr_p_fr, 2, 1)
            fig.add_trace(tr_p_py, 2, 1)

            fig.update_xaxes(title_text="Period (hours)", type="log", row=2, col=1)
            fig.update_yaxes(title_text="dB", row=1, col=1)
            fig.update_yaxes(title_text="Phase °",  row=2, col=1)
            fig.update_layout(
                height=600,                        
                title=dict(text=f"{station} • {nice} ({run})", x=0.5, y=0.97),
                legend=dict(orientation="h",
                            x=0.5, xanchor="center",
                            y=-0.16, yanchor="top",
                bordercolor="black", borderwidth=1),
                margin=dict(l=80, r=40, t=90, b=120)
            )

            fig.update_xaxes(type="log", title_text="Period (hours)",
                            showgrid=True, ticks="outside", row=2, col=1)
            fig.update_yaxes(showgrid=True, ticks="outside")

            ofs_dir = root_out / "OFS" / run
            ofs_dir.mkdir(parents=True, exist_ok=True)
            _save_plotly(fig, ofs_dir / f"{var_key}_{run}_bode")


########################################################################
#                            SCRIPT EXECUTION                          #
########################################################################
def run_all() -> None:
    ofs_in   = input("OFS code (cbofs / dbofs / sscofs …): ").strip().lower()
    start_in = input("Start UTC  (YYYY-MM-DD HH:MM): ").strip()
    end_in   = input("End   UTC  (YYYY-MM-DD HH:MM): ").strip()
    datum_in = input("WL datum used in Python [MHHW MHW MTL MSL DTL MLW MLLW NAVD STND]: ").strip().upper()
    configure(ofs_in, start_in, end_in, datum_in)

    root   = Path(__file__).resolve().parent
    fr_run = root / "fortran run"
    out_csv = root / "fortran_csv"
    out_csv.mkdir(exist_ok=True)

    print("\nConverting Fortran tables → CSV …")
    here = Path.cwd()
    os.chdir(fr_run)
    try:
        for var in VAR_FOLDERS.values():
            convert_var(var, out_csv)
    finally:
        os.chdir(here)

    print("\nBuilding station coordinates & datum offsets …")
    build_station_coords()

    if DATUM == "MSL":
        print("Datum = MSL → offset build & WL shift skipped.")
    else:
        make_station_offsets()
        offsets = load_station_offsets()

        run_shift = input("Apply WL datum shift now? [y/N] ").strip().lower() in {"y", "yes"}
        if run_shift:
            shift_fortran_wl(offsets)
        else:
            print("  ↳ Skipping WL shift (assumed done already).")

    print("\nBuilding Python vs Fortran comparison workbook …")
    global CSV_DIR
    CSV_DIR = str(root / "fortran_csv")
    comparison_main()

    run_ts   = input("Generate time-series plots?        [y/N] ").strip().lower() in {"y", "yes"}
    run_bode = input("Generate Bode plots?         [y/N] ").strip().lower() in {"y", "yes"}

    if any((run_ts, run_bode)):
        wipe = input("Wipe existing figure folders first? [y/N] ").strip().lower() in {"y", "yes"}
        if wipe:
            for p in (PLOTS_ROOT, BODE_ROOT):
                shutil.rmtree(p, ignore_errors=True)

        df_index = build_file_index()        

        if run_ts:
            make_plots(df_index)

        if run_bode:
            make_bode_all(df_index)

    print("\n✓ DONE – outputs written.\n")

if __name__ == "__main__":
    run_all()