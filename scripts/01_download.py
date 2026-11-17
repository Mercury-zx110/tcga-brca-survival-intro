# -*- coding: utf-8 -*-
"""01 下载 TCGA-BRCA 开放层数据（GDC API，无需认证）

- data/clinical.json    全部 BRCA 病例的临床字段（vital_status、年龄、分期、随访）
- data/expression/      100 例原发肿瘤 STAR - Counts 表达文件（gzip TSV）
- data/manifest.json    文件-病例对应表

注意：GDC 过滤语法中 "=" 操作的 content 必须是 {"field","value"} 对象，
"in" 数组形式在当前部署会返回 HTTP 500。
"""
import json, io, os, sys, tarfile, time, urllib.parse, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)
EXPR = os.path.join(DATA, "expression")
os.makedirs(EXPR, exist_ok=True)

UA = {"User-Agent": "tcga-brca-survival-intro/1.0 (educational; contact: student)",
      "Accept": "application/json"}
N_EXPR = 100  # 表达文件下载数量（练手规模）


def gdc_get(path, params, binary=False, retries=4):
    url = "https://api.gdc.cancer.gov" + path + "?" + urllib.parse.urlencode(params)
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            raw = urllib.request.urlopen(req, timeout=180).read()
            return raw if binary else json.loads(raw)
        except Exception as e:
            print(f"  retry {i+1}/{retries} ({e})", flush=True)
            time.sleep(8 * (i + 1))
    raise RuntimeError(f"GET failed: {url}")


def eq(field, value):
    return {"op": "=", "content": {"field": field, "value": value}}


# ---------- 1. 临床数据 ----------
print("[1/3] 下载临床数据 ...", flush=True)
flt = json.dumps(eq("project.project_id", "TCGA-BRCA"))
fields = ",".join([
    "case_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.age_at_diagnosis",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.days_to_death",
    "diagnoses.days_to_last_follow_up",
])
# follow_ups 是嵌套数组实体，必须用 expand 展开；
# 随访时长字段是 days_to_follow_up（不是 days_to_last_follow_up）
d = gdc_get("/cases", {"filters": flt, "fields": fields,
                       "expand": "follow_ups", "size": "1300"})
clinical = d["data"]["hits"]
json.dump(clinical, open(os.path.join(DATA, "clinical.json"), "w", encoding="utf-8"),
          ensure_ascii=False)
print(f"  临床病例数: {len(clinical)}", flush=True)

# ---------- 2. 表达文件清单（原发肿瘤，STAR - Counts） ----------
print("[2/3] 查询表达文件清单 ...", flush=True)
fflt = json.dumps({"op": "and", "content": [
    eq("cases.project.project_id", "TCGA-BRCA"),
    eq("files.data_category", "Transcriptome Profiling"),
    eq("files.data_type", "Gene Expression Quantification"),
    eq("files.analysis.workflow_type", "STAR - Counts"),
    eq("cases.samples.sample_type", "Primary Tumor"),
]})
d = gdc_get("/files", {"filters": fflt,
                       "fields": "file_id,file_size,cases.case_id",
                       "size": str(N_EXPR)})
hits = d["data"]["hits"]
manifest = []
for h in hits:
    case_id = h["cases"][0]["case_id"]
    manifest.append({"file_id": h["file_id"], "case_id": case_id,
                     "file_size": h.get("file_size")})
json.dump(manifest, open(os.path.join(DATA, "manifest.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
total_mb = sum(m["file_size"] for m in manifest) / 1e6
print(f"  表达文件: {len(manifest)} 个, 合计约 {total_mb:.0f} MB", flush=True)

# ---------- 3. 逐文件下载（curl.exe 并发 4 路；本网络下 GDC /data 会重置 Python TLS 连接） ----------
print("[3/3] 逐文件下载表达文件 ...", flush=True)
import subprocess
from concurrent.futures import ThreadPoolExecutor

LIMIT = int(os.environ.get("N_EXPR", "0"))  # 调试用：限制下载数量
todo = manifest[:LIMIT] if LIMIT > 0 else manifest

def fetch(m):
    out = os.path.join(EXPR, m["file_id"] + ".gz")
    if os.path.exists(out) and os.path.getsize(out) > 1e5:
        return m["file_id"], True
    url = "https://api.gdc.cancer.gov/data/" + m["file_id"]
    for i in range(6):
        try:
            r = subprocess.run(["curl.exe", "-s", "--retry", "3", "--retry-delay", "5",
                                "-o", out, "-w", "%{http_code}", url],
                               capture_output=True, text=True, timeout=180)
            code = r.stdout.strip()
        except subprocess.TimeoutExpired:
            code = "TIMEOUT"
        if code == "200" and os.path.exists(out) and os.path.getsize(out) > 1e5:
            return m["file_id"], True
        time.sleep(6 * (i + 1))
    return m["file_id"], False

ok, fail = 0, []
with ThreadPoolExecutor(max_workers=4) as ex:
    for k, (fid, success) in enumerate(ex.map(fetch, todo), 1):
        if success:
            ok += 1
        else:
            fail.append(fid)
        if k % 10 == 0:
            print(f"  进度 {k}/{len(todo)}", flush=True)
print(f"  完成: {ok}/{len(todo)} 个表达文件 -> {EXPR}", flush=True)
if fail:
    print("  失败 file_id:", fail, flush=True)
