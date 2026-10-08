#!/usr/bin/env python3
"""Builds js/data/catalog.js (the reference catalogue shown on each suite's Equations tab) from the project
specification documents.

    python3 tools/catalog.py <equations.docx> <ibcs.docx> <cvv.docx> <modules.docx>

Needs pandoc on the PATH. The catalogue lists, per suite: equation families grouped by module, hybrid
formulations, initial and boundary conditions, input and output data, and calibration / verification /
validation practice."""
import json, re, subprocess, sys, html, pathlib

def pandoc(path, fmt):
    return subprocess.run(['pandoc', '-f', 'docx', '-t', fmt, '--wrap=none', path], check=True, capture_output=True, text=True).stdout

def clean(s):
    s = re.sub(r'\\\\?\(?\\?', '', s) if False else s
    s = s.replace('\\\\k-\\epsilon\\\\', 'k–ε').replace('\\\\k-\\omega\\\\', 'k–ω').replace('\\(R^2\\)', 'R²').replace('\\(C_v/K_v\\)', 'Cv/Kv').replace('\\(C_v\\)', 'Cv').replace('\\(U\\)', 'U')
    s = re.sub(r'\\\((.*?)\\\)', r'\1', s)
    return re.sub(r'\s+', ' ', s.replace('**', '').replace('\\', '')).strip(' .;:')

def cap(s):
    return s[:1].upper() + s[1:] if s else s

def split_items(text):
    out, depth, cur = [], 0, ''
    for ch in text:
        if ch in '([': depth += 1
        if ch in ')]': depth -= 1
        if ch == ';' and depth <= 0:
            out.append(cur); cur = ''
        else: cur += ch
    out.append(cur)
    return [cap(clean(x)) for x in out if clean(x)]

def equations(md):
    suites, cur, pending = {}, None, None
    for line in md.splitlines():
        line = line.strip()
        m = re.match(r'^\*\*(\d)\. (.+?)\*\*$', line)
        if m: cur = suites.setdefault(m.group(1), {'groups': []}); pending = None; continue
        if cur is None or line.startswith('**How the seven modules'):
            if line.startswith('**How the seven modules'): cur = None
            continue
        m = re.match(r'^\*\*(.+?):\*\*\s*(.+)$', line)
        if m: cur['groups'].append({'title': clean(m.group(1)), 'items': split_items(m.group(2))}); pending = None; continue
        m = re.match(r'^\*\*(.+?)\*\*$', line)
        if m: pending = clean(m.group(1)); continue
        if pending and line and not line.startswith('--'):
            if '→' in line: # hybrid techno-economic couplings written as "A + B → purpose"
                g = next((x for x in cur['groups'] if x['title'] == pending), None)
                if not g: g = {'title': pending, 'items': []}; cur['groups'].append(g)
                g['items'].append(cap(clean(line.split('→')[0])))
            elif ';' in line: cur['groups'].append({'title': pending, 'items': split_items(line)}); pending = None
    for s in suites.values():
        hyb = [g for g in s['groups'] if re.search(r'hybrid', g['title'], re.I)]
        s['hybrid'] = [i for g in hyb for i in g['items']]
        s['classical'] = [i for g in s['groups'] if g not in hyb for i in g['items']]
        s['modules'] = [g['title'] for g in s['groups'] if g not in hyb]
    return suites

LEAD = r'(?:should (?:specify|define|begin with|prescribe|provide|include|also define|additionally define|distinguish|accommodate|conserve|impose)|should also (?:include|define|receive))'
def phrases(text):
    out = []
    for sent in re.split(r'(?<=[a-z\)])\.\s+', text):
        m = re.search(LEAD + r'\s+(.*)', sent)
        body = m.group(1) if m else None
        if not body: continue
        for part in re.split(r';|,\s+(?:and\s+)?(?![^()]*\))', body):
            p = clean(re.sub(r'^(?:and|the|an|a|any|all|every)\s+', '', part.strip(), flags=re.I))
            if 3 < len(p) < 110: out.append(cap(p))
    seen, uniq = set(), []
    for p in out:
        if p.lower() not in seen: seen.add(p.lower()); uniq.append(p)
    return uniq

def icbc(md):
    out, cur = {}, None
    for line in md.splitlines():
        line = line.strip()
        m = re.match(r'^\*\*(\d)\. (.+?)\*\*$', line)
        if m: cur = out.setdefault(m.group(1), {'initialText': '', 'boundaryText': ''}); continue
        if cur is None: continue
        m = re.match(r'^\*\*(Initial conditions|Boundary conditions):\*\*\s*(.+)$', line)
        if m: cur['initialText' if m.group(1).startswith('Initial') else 'boundaryText'] = clean(m.group(2)) + '.'
    for s in out.values():
        s['initial'] = phrases(s['initialText']); s['boundary'] = phrases(s['boundaryText'])
    return out

def cvv(txt):
    out, cur, key = {}, None, None
    for line in txt.splitlines():
        t = line.strip()
        m = re.match(r'^(\d)\. (Fluid|Geometry|Multiphase|Hydrate|Operations|Integrity|Economics)', t)
        if m: cur = out.setdefault(m.group(1), {'calibration': [], 'verification': [], 'validation': [], 'notes': []}); key = None; continue
        if cur is None: continue
        if t.startswith('Cross-Module Validation'): cur = None; continue
        if t in ('Calibration',): key = 'calibration'; continue
        if t in ('Verification', 'Numerical Verification'): key = 'verification'; continue
        if t == 'Validation': key = 'validation'; continue
        if t in ('Typical Inputs', 'Typical Outputs'): key = None; continue
        if t.startswith('- ') and key: cur[key].append(cap(clean(t[2:])))
        elif key and len(t) > 60 and not t.startswith('-') and not re.search(r'\bI would\b|\byour\b|\bYes\b|\bI am\b', t): cur['notes'].append(clean(t) + '.')
    return out

def modules(htmltext):
    out = {}
    for row in re.findall(r'<tr[^>]*>(.*?)</tr>', htmltext, re.S):
        cells = [clean(html.unescape(re.sub(r'<[^>]+>', '', c))) for c in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', row, re.S)]
        if len(cells) >= 4 and cells[0].isdigit(): out[cells[0]] = {'title': cells[1], 'inputs': split_items(cells[2]), 'outputs': split_items(cells[3])}
    return out

ECON_IO = {
    'inputs': 'CAPEX; installation cost; pipeline/riser cost; equipment cost; insulation cost; chemical-injection infrastructure; chemical cost; energy cost; maintenance cost; inspection cost; pigging cost; repair cost; replacement cost; production rates; commodity prices; production decline; downtime; hydrate/wax/scale intervention frequency; failure probabilities from Module 6; production deferment; shutdown duration; project life; inflation; discount rate; tax; royalties; depreciation; financing assumptions; abandonment/decommissioning costs',
    'outputs': 'CAPEX; OPEX; TOTEX/life-cycle cost; revenue; cash flow; discounted cash flow; NPV; IRR; MIRR; payback period; discounted payback; profitability index; breakeven commodity price; breakeven production; unit technical cost; cost per barrel/BOE; chemical cost; energy cost; maintenance cost; expected downtime cost; production-deferment cost; expected failure cost; expected annual loss; lifecycle intervention cost; risk-adjusted NPV; economic operating envelope; alternative-ranking metrics',
}
LEVELS = [
    ['Unit verification', 'individual equations, correlations and algorithms'], ['Solver verification', 'manufactured and analytical solutions, grid and time-step convergence, conservation'],
    ['Model calibration', 'uncertain physical parameters estimated from designated calibration data'], ['Component validation', 'each model validated independently against separate data'],
    ['Coupled-module verification', 'data transfer and conservation between suites'], ['Integrated validation', 'complete laboratory or field operating scenarios reproduced end to end'],
    ['Blind validation', 'data not used in development or calibration are predicted'], ['Uncertainty quantification', 'input uncertainty propagated through the complete chain'],
    ['Sensitivity and identifiability', 'screening, variance-based indices, parameter correlation'], ['Probabilistic validation', 'prediction intervals, coverage, reliability diagrams'],
    ['Field validation', 'start-up, normal production, slugging, shutdown, cooldown, restart and intervention data'], ['Independent validation', 'data and tests reproduced by people who did not calibrate the models'],
]

def main(eq, ib, cv, mo):
    E, I, V, M = equations(pandoc(eq, 'markdown_strict-raw_html')), icbc(pandoc(ib, 'markdown_strict-raw_html')), cvv(pandoc(cv, 'plain')), modules(pandoc(mo, 'html'))
    M['7'] = {'title': 'Economics, Techno-Economics & Decision Analysis', 'inputs': split_items(ECON_IO['inputs']), 'outputs': split_items(ECON_IO['outputs'])}
    cat = {}
    for k in '1234567':
        e, i, v, m = E[k], I[k], V[k], M[k]
        cat[k] = {'title': m['title'], 'groups': e['groups'], 'classical': e['classical'], 'hybrid': e['hybrid'], 'modules': e['modules'],
                  'initial': i['initial'], 'boundary': i['boundary'], 'initialText': i['initialText'], 'boundaryText': i['boundaryText'],
                  'inputs': m['inputs'], 'outputs': m['outputs'], 'calibration': v['calibration'], 'verification': v['verification'], 'validation': v['validation'], 'practice': v['notes']}
    body = json.dumps(cat, ensure_ascii=False, indent=1)
    head = ('// Reference catalogue for the 7 suites: equation families grouped by module, hybrid formulations, initial and\n'
            '// boundary conditions, input and output data, and calibration / verification / validation practice.\n'
            '// Generated by tools/catalog.py from the project specification.\n')
    pathlib.Path(__file__).resolve().parent.parent.joinpath('js/data/catalog.js').write_text(head + 'export const CATALOG = ' + body + ';\n\nexport const QUALIFICATION_LEVELS = ' + json.dumps(LEVELS, ensure_ascii=False) + ';\n', encoding='utf-8')
    for k, c in cat.items(): print(k, c['title'], '| groups', len(c['groups']), 'classical', len(c['classical']), 'hybrid', len(c['hybrid']), 'IC', len(c['initial']), 'BC', len(c['boundary']), 'in', len(c['inputs']), 'out', len(c['outputs']), 'cal', len(c['calibration']), 'ver', len(c['verification']), 'val', len(c['validation']))

if __name__ == '__main__':
    if len(sys.argv) != 5: sys.exit(__doc__)
    main(*sys.argv[1:])
