"""Build the static website (site/index.html) from data/*.json."""
import json, os
from . import model


def _window(data):
    """{site_key: {MM-DD: {...}}} from data/history_window/<site>/MM-DD.json"""
    import glob
    out = {}
    for f in sorted(glob.glob(os.path.join(data, 'history_window', '*', '*.json'))):
        d = json.load(open(f))
        out.setdefault(os.path.basename(os.path.dirname(f)), {})[d['date']] = d
    return out


def build(root, cfg):
    data = os.path.join(root, 'data')
    load = lambda n, d: json.load(open(os.path.join(data, n))) if os.path.exists(os.path.join(data, n)) else d
    payload = dict(latest=load('latest.json', None), history=load('history.json', {}),
                   climatology=_window(data),
                   config=dict(burst=cfg['burst_mode'], float=cfg['float_mode'], site=cfg['site'], sites=model.sites(cfg), flight=cfg['flight'],
                               avoid=cfg.get('avoid_zones') or [], map=cfg.get('map') or {}))
    tpl = open(os.path.join(root, 'forecast', 'template.html'), encoding='utf-8').read()
    html = tpl.replace('/*__DATA__*/null', json.dumps(payload, separators=(',', ':')).replace('</', '<\\/'))
    os.makedirs(os.path.join(root, 'site'), exist_ok=True)
    with open(os.path.join(root, 'site', 'index.html'), 'w', encoding='utf-8') as f:
        f.write(html)
    print('site/index.html written')
