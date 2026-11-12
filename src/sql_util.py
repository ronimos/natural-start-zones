"""
sql_util.py — read CAIC station observations from the CAIC MySQL database.

Adapted from the avachain sql_util.py.  Differences: stations come from the
site's 'weather_stations' entry in pipeline.SITES (or --station on the CLI)
instead of a slope TOML, and a column requested from two stations is an
error rather than a silent duplicate.

The database runs on the same host as the pipeline, so this connects to it
directly (no SSH tunnel).

Stations are given as {sql_id: [columns]}.  ["*"] means every column; those
columns are then prefixed '<sql_id>.' so stations cannot collide.

The frame returned by fetch_weather() has a `time` index (UTC, naive, as
stored) and the raw CAIC columns in raw CAIC units (tenths of °F, tenths of
mph, ...), each taken from the station that lists it.

Credentials (db_user, db_password, db_name, local_host, db_port) are read,
first file found wins, from:
  1. $AVACHAIN_DB_ENV
  2. ~/.config/avachain/db.env
  3. the repository's .env
Variables already set in the environment take precedence over all of them.

CLI:
    python src/sql_util.py --site sisters [--start 2009-10-01] [--end 2021-06-01]
    python src/sql_util.py --station CALVP --station LBAC2 --columns   # list columns
"""

from __future__ import annotations

import argparse
import datetime
import logging
import os
from pathlib import Path
from typing import Optional

import pandas as pd
import sqlalchemy
from sqlalchemy import create_engine, text

logger = logging.getLogger(__name__)

ROOT_PATH = Path(__file__).resolve().parents[1]
TABLES = ["obsBattery", "obsHydro", "obsSnow", "obsSolar", "obsWX"]
ENV_VARS = ("db_user", "db_password", "db_name", "local_host")


# --------- Stations ---------

def stations_from_cfg(cfg: dict) -> dict:
    """{sql_id: [columns]} from a SITES entry's 'weather_stations'.

    Raises ValueError when the site lists no stations or a station no columns.
    """
    stations = cfg.get('weather_stations') or {}
    if not stations:
        raise ValueError(f"site {cfg.get('name')!r} has no 'weather_stations'; "
                         f"station data cannot be fetched")
    for sid, cols in stations.items():
        if not cols:
            raise ValueError(f"station {sid!r} lists no columns")
    return {sid: list(cols) for sid, cols in stations.items()}


# --------- Connection ---------

def _env_file() -> Optional[Path]:
    cands = [os.getenv("AVACHAIN_DB_ENV"),
             Path.home() / ".config" / "avachain" / "db.env",
             ROOT_PATH / ".env"]
    for c in cands:
        if c and Path(c).is_file():
            return Path(c)
    return None


def get_engine() -> sqlalchemy.engine.Engine:
    """Engine for the CAIC database; raises with the remedy when unconfigured."""
    from dotenv import load_dotenv
    env = _env_file()
    if env is not None:
        load_dotenv(env, override=False)
    missing = [k for k in ENV_VARS if not os.getenv(k)]
    if missing:
        raise RuntimeError(
            f"CAIC database credentials missing ({', '.join(missing)}). Put "
            f"db_user, db_password, db_name, local_host and db_port in "
            f"$AVACHAIN_DB_ENV or ~/.config/avachain/db.env (mode 600)"
            + (f"; read {env}" if env else "; no credentials file found"))
    url = sqlalchemy.engine.URL.create(
        "mysql+pymysql", username=os.getenv("db_user"),
        password=os.getenv("db_password"), host=os.getenv("local_host"),
        port=int(os.getenv("db_port", 3306)), database=os.getenv("db_name"))
    return create_engine(url, connect_args={"connect_timeout": 20})


# --------- Queries ---------

def _fetch_table_data(engine, tables, start_time, st_id, end_time=None):
    """One frame per table holding rows for st_id with start < time <= end.

    A table with no rows is skipped.  A connection or SQL error raises, so a
    database outage cannot look like a station with no data.
    """
    frames = []
    q = "SELECT * FROM {t} WHERE time > :start AND staname = :sid"
    params = {"start": start_time, "sid": st_id}
    if end_time is not None:
        q += " AND time <= :end"
        params["end"] = end_time
    for table in tables:
        df = pd.read_sql(text(q.format(t=table)), con=engine, params=params)
        if df.empty:
            logger.info(f"no rows in {table} for {st_id}")
            continue
        df["time"] = pd.to_datetime(df["time"])
        df = df.drop(columns=["staname"], errors="ignore").set_index("time")
        # A station can report the same timestamp twice; keep the last.
        frames.append(df[~df.index.duplicated(keep="last")])
    return frames


def get_weather_from_sql_db(st_id: str,
                            start_time: Optional[datetime.datetime] = None,
                            end_time: Optional[datetime.datetime] = None,
                            engine=None) -> pd.DataFrame:
    """All observation columns for one station, start < time <= end.

    start defaults to 7 days ago.  Times are the database's (UTC, naive).
    """
    if start_time is None:
        start_time = (datetime.datetime.now(datetime.timezone.utc)
                      .replace(tzinfo=None) - datetime.timedelta(hours=168))
    engine = engine or get_engine()
    frames = _fetch_table_data(engine, TABLES, _ts(start_time), st_id,
                               _ts(end_time) if end_time is not None else None)
    if not frames:
        logger.warning(f"no data retrieved for station {st_id}")
        return pd.DataFrame()
    df = pd.concat(frames, axis=1).sort_index()
    # Same column name in two tables (e.g. an id field): keep the first.
    return df.loc[:, ~df.columns.duplicated()]


def fetch_weather(stations: dict, start_time, end_time=None, engine=None) -> pd.DataFrame:
    """Station observations, one column per requested column.

    @param stations: {sql_id: [columns]}; ["*"] = all columns, prefixed '<sql_id>.'
    @returns: DataFrame indexed by `time`; empty when no rows
    @raises ValueError: a column is requested from more than one station
    """
    named = [c for cols in stations.values() if cols != ["*"] for c in cols]
    dupes = sorted({c for c in named if named.count(c) > 1})
    if dupes:
        raise ValueError(f"columns requested from more than one station: {dupes}")

    engine = engine or get_engine()
    parts = []
    for sid, cols in stations.items():
        df = get_weather_from_sql_db(sid, start_time, end_time, engine=engine)
        if df.empty:
            logger.warning(f"station {sid}: no rows between {start_time} and "
                           f"{end_time or 'now'}")
            continue
        if cols == ["*"]:
            parts.append(df.add_prefix(f"{sid}."))
            continue
        absent = [c for c in cols if c not in df.columns]
        if absent:
            logger.warning(f"station {sid}: columns not in the database: {absent}")
        parts.append(df[[c for c in cols if c in df.columns]])
    if not parts:
        return pd.DataFrame()
    out = pd.concat(parts, axis=1).sort_index()
    out.index.name = "time"
    return out


def station_columns(st_id: str, engine=None, start_time=None, end_time=None) -> pd.DataFrame:
    """Per column of one station: table, non-null count, first and last time."""
    engine = engine or get_engine()
    rows = []
    for table in TABLES:
        q = "SELECT * FROM {t} WHERE staname = :sid".format(t=table)
        p = {"sid": st_id}
        if start_time is not None:
            q += " AND time > :start"; p["start"] = _ts(start_time)
        if end_time is not None:
            q += " AND time <= :end"; p["end"] = _ts(end_time)
        df = pd.read_sql(text(q), con=engine, params=p)
        if df.empty:
            continue
        df["time"] = pd.to_datetime(df["time"])
        for c in df.columns.drop(["time", "staname"], errors="ignore"):
            s = df.loc[df[c].notna(), "time"]
            rows.append(dict(table=table, column=c, n=len(s),
                             first=s.min() if len(s) else None,
                             last=s.max() if len(s) else None))
    return pd.DataFrame(rows)


def _ts(t) -> str:
    return pd.Timestamp(t).strftime("%Y-%m-%d %H:%M:%S")


# --------- CLI ---------

def parse_args():
    p = argparse.ArgumentParser(description="Fetch station data from the CAIC database.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--site", help="site key in pipeline.SITES (uses its weather_stations)")
    src.add_argument("--station", action="append",
                     help="station sql_id (repeatable); fetches all its columns")
    p.add_argument("--start", default="2009-10-01", help="start, UTC (default 2009-10-01)")
    p.add_argument("--end", default=None, help="end, inclusive, UTC (default: now)")
    p.add_argument("--out", default=None,
                   help="output CSV (default: data/weather/<site or stations>.csv)")
    p.add_argument("--columns", action="store_true",
                   help="list each station's columns and coverage instead of fetching")
    return p.parse_args()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    args = parse_args()
    if args.site:
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from pipeline import SITES
        stations = stations_from_cfg(SITES[args.site])
        tag = args.site
    else:
        stations = {sid: ["*"] for sid in args.station}
        tag = "_".join(args.station)

    if args.columns:
        eng = get_engine()
        for sid in stations:
            print(f"\n== {sid}")
            print(station_columns(sid, eng, args.start, args.end).to_string(index=False))
        raise SystemExit(0)

    df = fetch_weather(stations, args.start, args.end)
    if df.empty:
        raise SystemExit("no data collected from any station; CSV not written")
    out = Path(args.out) if args.out else ROOT_PATH / "data" / "weather" / f"{tag}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out)
    logger.info(f"wrote {len(df)} rows × {df.shape[1]} columns to {out}")
