"""
merge_avalanche_records.py — append coordinate-assigned CAIC records to the
named CAIC highway export, in the export's own column layout.

Input:  the named export (e.g. CAIC_HWY_avalanches_2009-01-01_2021-05-04.csv)
        and the output of assign_avalanche_paths.py (records with lat/lon,
        assigned to a path).
Rules:
  - only records assigned to a path are added;
  - dates are converted from UTC to the local (America/Denver) date;
  - a record is skipped as a duplicate when the named export already has a
    record for the same path on the same local date (2020-03-22 Sister #6,
    2020-12-16 Sister #6 and 2021-02-06 Sister #3 overlap this way);
  - added rows carry Source = 'coords-assigned' plus the latitude/longitude,
    the distance to the path corridor and the ambiguity flag.

Usage:
  python src/merge_avalanche_records.py NAMED.csv ASSIGNED.csv OUT.csv
"""
from __future__ import annotations

import sys

import pandas as pd

PATH_NAMES = {'Sister #1': 'Seven Sister #1', 'Sister #2': 'Seven Sister #2',
              'Sister #3': 'Seven Sister #3', 'Sister #4': 'Seven Sister #4',
              'Sister #6': 'Seven Sister #6'}


def main(named_csv: str, assigned_csv: str, out_csv: str) -> None:
    old = pd.read_csv(named_csv, dtype=str, keep_default_na=False)
    new = pd.read_csv(assigned_csv)
    new = new[new['assigned'].notna() & (new['assigned'] != '')].copy()
    new['local'] = pd.to_datetime(new['Date'], utc=True).dt.tz_convert('America/Denver')
    new['HW Path'] = new['assigned'].map(PATH_NAMES)

    old_day = pd.to_datetime(old['Date'], format='mixed', errors='coerce').dt.date
    have = set(zip(old_day, old['HW Path'].str.strip()))
    dup = [(d, p) in have for d, p in zip(new['local'].dt.date, new['HW Path'])]
    print(f"{len(new)} assigned records; {sum(dup)} already in the named export (skipped)")
    for (_, r), dp in zip(new.iterrows(), dup):
        if dp:
            print(f"  duplicate: {r['local'].date()} {r['HW Path']}")
    new = new[[not x for x in dup]]

    rows = pd.DataFrame({c: [''] * len(new) for c in old.columns})
    rows['Date'] = new['local'].dt.strftime('%m/%d/%Y').values
    rows['Time'] = new['local'].dt.strftime('%I:%M %p').values
    rows['Date Known'] = new['Date Known'].values
    rows['HW Zone'] = 'Loveland Pass'
    rows['BC Zone'] = 'Front Range'
    rows['HW Path'] = new['HW Path'].values
    rows['#'] = '1'
    rows['Elev'] = new['Elevation'].fillna('').replace('-', '').values
    rows['Asp'] = new['Aspect'].fillna('').values
    rows['Type'] = new['Type'].values
    rows['Trigger'] = new['Trigger'].values
    rows['Rsize'] = new['Relative Size'].fillna('').replace('-', '').values
    rows['Dsize'] = new['Destructive Size'].fillna('').replace('-', '').values
    rows['Source'] = 'coords-assigned'
    rows['Latitude'] = new['latitude'].values
    rows['Longitude'] = new['Longitude'].values
    rows['Corridor dist (m)'] = new['dist_m'].values
    rows['Assignment ambiguous'] = new['ambiguous'].values
    old['Source'] = 'CAIC named export'
    out = pd.concat([old, rows], ignore_index=True)
    out.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}: {len(old)} named + {len(rows)} added records")
    print(rows.groupby('HW Path').size().to_string())


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
