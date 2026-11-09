#!/usr/bin/env python3
"""Measure the whole search process, excluding dedispersion and follow-up."""
import argparse, datetime, hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("command", nargs=argparse.REMAINDER)
    a=p.parse_args(); cmd=a.command[1:] if a.command[:1]==["--"] else a.command
    if not cmd: p.error("Missing command")
    if a.out.exists(): raise FileExistsError(a.out)
    record=dict(command=cmd, hostname=platform.node(), python=sys.version,
        affinity=sorted(os.sched_getaffinity(0)), job_id=os.environ.get("SLURM_JOB_ID"),
        grid_index=os.environ.get("SEARCH_GRID_INDEX"),
        array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID"), task_id=os.environ.get("SLURM_ARRAY_TASK_ID"),
        utc_start=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        cpu_model=next((s for s in Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")),""),
        thread_environment={k:os.environ.get(k) for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"]},
        scope="whole search command: startup, input reading, spectrogram, Viterbi, lean products; no dedispersion or follow-up",
        file_sha256={s:hashlib.sha256(Path(s).read_bytes()).hexdigest() for s in cmd if Path(s).is_file()})
    started=time.perf_counter(); child=subprocess.Popen(cmd)
    _,status,usage=os.wait4(child.pid,0)
    child.returncode=os.waitstatus_to_exitcode(status)
    record.update(wall_s=time.perf_counter()-started,user_s=usage.ru_utime,system_s=usage.ru_stime,
        cpu_s=usage.ru_utime+usage.ru_stime,max_rss_kib=usage.ru_maxrss,returncode=child.returncode,
        utc_end=datetime.datetime.now(datetime.timezone.utc).isoformat())
    a.out.write_text(json.dumps(record,indent=2)+"\n")
    print("SEARCH_TIMING "+json.dumps(record),flush=True)
    sys.exit(child.returncode if child.returncode>=0 else 128-child.returncode)
if __name__=="__main__": main()
