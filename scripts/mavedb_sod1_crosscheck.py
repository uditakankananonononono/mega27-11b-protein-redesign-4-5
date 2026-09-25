"""Independent SOD1 variant-map comparison against the frozen 1SPD GNN scan.

Uses public MaveDB score-set CSVs; a distinct experimental functional assay,
not a direct ddG target. No fitting, selection, or score calibration.
"""
import csv, hashlib, json, pathlib, re, urllib.request
from scipy.stats import spearmanr

ROOT = pathlib.Path(__file__).resolve().parents[1]
API = 'https://api.mavedb.org/api/v1/score-sets/'
SETS = {'abundance':'urn:mavedb:00001217-a-3', 'activity':'urn:mavedb:00001217-a-4'}
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
              'Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val'.split()))


def analyze(scan, score_rows):
    scores = {r['hgvs_pro']: r for r in score_rows if r['score'] not in ('', 'NA')}
    pairs = []
    for row in scan['all_rows']:
        if row['ood_wt_pro']:
            continue
        # 1SPD mature SOD1 Ala1 = canonical precursor Ala2. Check this
        # offset against the independent SIFTS/ClinVar records in results/.
        key = f"p.{AA[row['wt']]}{row['position']+1}{AA[row['mut']]}"
        found = scores.get(key)
        if found:
            pairs.append((row['ddg_corr'], float(found['score'])))
    rho, p = spearmanr([x[0] for x in pairs], [x[1] for x in pairs])
    return {'n_matched':len(pairs), 'spearman':float(rho), 'p_unadjusted':float(p),
            'expected_relation':'positive if predicted stability tracks abundance/activity',
            'interpretation':'independent functional-map check only, not a ddG benchmark'}


def main():
    scan = json.loads((ROOT/'build/scan_1spd_whole_v2.json').read_text())
    out = {'scan': 'build/scan_1spd_whole_v2.json', 'mapping':'PDB mature position +1 to canonical HGVS', 'maps':{}}
    for label, urn in SETS.items():
        url = API + urn + '/scores'
        with urllib.request.urlopen(url, timeout=40) as response:
            raw = response.read()
        path = ROOT/'results'/('mavedb_sod1_'+label+'.csv')
        path.write_bytes(raw)
        rows = list(csv.DictReader(raw.decode('utf-8').splitlines()))
        out['maps'][label] = {'source_url':url,'sha256':hashlib.sha256(raw).hexdigest(),
                              'n_source_rows':len(rows), **analyze(scan,rows)}
    (ROOT/'results/mavedb_sod1_crosscheck.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__ == '__main__': main()
