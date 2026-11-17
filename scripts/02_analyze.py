# -*- coding: utf-8 -*-
"""02 临床—组学联合生存分析（TCGA-BRCA，100 例原发肿瘤）

流程：
  1. 临床表构建：生存时间（年）、事件（死亡=1）、年龄、分期
  2. STAR - Counts 解析：protein_coding 基因，CPM 标准化，log2(CPM+1)
  3. 目标基因（ESR1 / ERBB2 / MKI67）中位数分组 → KM 曲线 + 对数秩检验
  4. 多变量 Cox 比例风险模型（Efron 结近似，Newton-Raphson 求解）

输出：results/summary.json、figures/km_genes.png
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "Noto Sans SC", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

GENES = ["ESR1", "ERBB2", "MKI67"]
rng = np.random.default_rng(20261116)  # 固定随机种子（本项目随机性仅来自数据顺序）


# ---------- 1. 临床表 ----------
def first_not_null(vals):
    for v in vals:
        if v is not None:
            return v
    return None


def load_clinical():
    rows = []
    for h in json.load(open(os.path.join(DATA, "clinical.json"), encoding="utf-8")):
        demo = h.get("demographic") or {}
        diags = h.get("diagnoses") or []
        prim = next((d for d in diags if d.get("age_at_diagnosis") is not None), diags[0] if diags else {})
        event = 1 if demo.get("vital_status") == "Dead" else 0
        t_death = first_not_null([prim.get("days_to_death"), demo.get("days_to_death")])
        # 随访时长：优先 diagnoses.days_to_last_follow_up，其次 follow_ups 里最大的 days_to_follow_up
        fus = h.get("follow_ups") or []
        fu_days = [f["days_to_follow_up"] for f in fus if f.get("days_to_follow_up") is not None]
        t_fu = first_not_null([prim.get("days_to_last_follow_up"),
                               max(fu_days) if fu_days else None])
        t = t_death if (event == 1 and t_death is not None) else t_fu
        if t is None or t <= 0:
            continue
        stage = prim.get("ajcc_pathologic_stage") or ""
        # 注意 "Stage IIIA" 同时以 "Stage I" 和 "II" 开头，必须先判晚期再判早期
        s = stage.replace("Stage ", "").strip()
        if s.startswith(("III", "IV")):
            stage_grp = "late"
        elif s == "I" or s.startswith(("IA", "IB", "II")):
            stage_grp = "early"
        else:
            stage_grp = None
        rows.append({
            "case_id": h["case_id"],
            "time_years": t / 365.25,
            "event": event,
            "age_years": (prim.get("age_at_diagnosis") or np.nan) / 365.25,
            "stage_raw": stage or None,
            "stage_grp": stage_grp,
        })
    return pd.DataFrame(rows)


# ---------- 2. 表达矩阵 ----------
def load_expression(manifest):
    files = sorted(glob.glob(os.path.join(DATA, "expression", "*.gz")))
    if len(files) != len(manifest):
        print(f"警告: 文件数 {len(files)} != 清单数 {len(manifest)}，按可匹配者处理")
    name2case = {}
    for m in manifest:
        name2case[m["file_id"]] = m["case_id"]

    mat = {}
    lib_sizes = {}
    gene_type = {}
    used = 0
    for f in files:
        fid = os.path.basename(f).replace(".gz", "")
        case_id = name2case.get(fid)
        if case_id is None:  # tar 成员名可能带原始文件名，退化匹配
            case_id = next((v for k, v in name2case.items() if k in fid), None)
        if case_id is None:
            continue
        # GDC 的 STAR-Counts 经 /data 下发时为纯文本 TSV（首行 # 注释），
        # 虽以 .gz 结尾但并未压缩，须显式 compression=None 阻止 pandas 按后缀推断 gzip
        df = pd.read_csv(f, sep="\t", dtype=str, comment="#", compression=None)
        df = df[~df["gene_id"].str.startswith("N_")]
        counts = pd.to_numeric(df["unstranded"], errors="coerce").fillna(0.0)
        lib_sizes[case_id] = float(counts.sum())
        for g, c, t in zip(df["gene_name"], counts, df["gene_type"]):
            gene_type[g] = t
            mat.setdefault(g, {})[case_id] = mat.setdefault(g, {}).get(case_id, 0.0) + c
        used += 1
    print(f"表达文件匹配成功: {used}")

    expr = pd.DataFrame(mat).T  # mat 为 {基因: {病例: 计数}}，转置后基因=行、病例=列
    lib = pd.Series(lib_sizes)
    # CPM + log2 转换（仅保留蛋白编码基因）
    keep = [g for g in expr.index if gene_type.get(g) == "protein_coding"]
    expr = expr.loc[keep]
    cpm = expr.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


# ---------- 3. KM + 对数秩 ----------
def km_curve(time, event):
    order = np.argsort(time)
    t, e = time[order], event[order]
    ts = np.unique(t)
    surv, at_risk, n = [], len(t), []
    for x in ts:
        d = int(((t == x) & (e == 1)).sum())
        c = int(((t == x) & (e == 0)).sum())
        if at_risk > 0:
            surv.append((1 - d / at_risk))
        else:
            surv.append(1.0)
        n.append(at_risk)
        at_risk -= d + c
    return ts, np.cumprod(surv), n


def logrank(time, event, group):
    ts = np.unique(time)
    o1, e1, v = 0.0, 0.0, 0.0
    for x in ts:
        at = (time >= x)
        d = int(((time == x) & (event == 1)).sum())
        if d == 0 or at.sum() == 1:
            continue
        n1 = int((at & (group == 1)).sum())
        n = int(at.sum())
        o1 += int(((time == x) & (event == 1) & (group == 1)).sum())
        e1 += d * n1 / n
        v += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1) if n > 1 else 0.0
    chi2 = (o1 - e1) ** 2 / v if v > 0 else np.nan
    p = 1 - stats.chi2.cdf(chi2, 1) if not np.isnan(chi2) else np.nan
    return chi2, p


# ---------- 4. Cox PH（Efron 结近似，Newton-Raphson） ----------
def cox_efron(time, event, X, max_iter=50, tol=1e-8):
    """Efron 结近似 Cox 比例风险模型。

    每个死亡时刻 d（死亡数 m）：
      logL 贡献 = Σ_{i∈D} xb_i − Σ_{j=0}^{m-1} log S0_j
      S_k_j = Σ_{i∈R} exb_i·x_i^k − (j/m)·Σ_{i∈D} exb_i·x_i^k   (R 为风险集)
    得分与信息阵对 j 求和后做 Newton-Raphson 迭代。
    """
    t = np.asarray(time, float)
    ev = np.asarray(event, float)
    X = np.asarray(X, float)
    n, p = X.shape
    beta = np.zeros(p)
    uniq = np.unique(t[ev == 1])
    for _ in range(max_iter):
        xb = X @ beta
        exb = np.exp(xb - xb.max())
        g = np.zeros(p)
        H = np.zeros((p, p))
        for d in uniq:
            R = t >= d
            D = (t == d) & (ev == 1)
            m = int(D.sum())
            s0_R, s1_R, s2_R = exb[R].sum(), (exb[R, None] * X[R]).sum(0), (exb[R, None] * X[R]).T @ X[R]
            s0_D, s1_D, s2_D = exb[D].sum(), (exb[D, None] * X[D]).sum(0), (exb[D, None] * X[D]).T @ X[D]
            g += (X[D].sum(0) if m > 0 else 0.0)
            for j in range(m):
                w = j / m
                s0 = s0_R - w * s0_D
                s1 = s1_R - w * s1_D
                s2 = s2_R - w * s2_D
                g -= s1 / s0
                H += s2 / s0 - np.outer(s1, s1) / s0**2
        step = np.linalg.solve(H + 1e-10 * np.eye(p), g)
        beta += step
        if np.max(np.abs(step)) < tol:
            break
    se = np.sqrt(np.diag(np.linalg.inv(H + 1e-10 * np.eye(p))))
    z = beta / se
    pvals = 2 * (1 - stats.norm.cdf(np.abs(z)))
    return beta, se, pvals, H


def main():
    clin = load_clinical()
    manifest = json.load(open(os.path.join(DATA, "manifest.json"), encoding="utf-8"))
    expr = load_expression(manifest)

    merged = clin.set_index("case_id").join(expr.T, how="inner")
    n_all = len(merged)
    n_event = int(merged["event"].sum())
    print(f"分析样本: {n_all} 例, 死亡事件 {n_event} 个")
    print(f"蛋白编码基因数: {expr.shape[0]}")

    # ---------- KM + log-rank ----------
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    lr_results = {}
    for ax, g in zip(axes, GENES):
        if g not in merged.columns:
            ax.set_title(f"{g}（无数据）")
            continue
        med = merged[g].median()
        grp = (merged[g] >= med).astype(int)
        for label, gv, color in [("高表达", 1, "#1F6E7E"), ("低表达", 0, "#C0504D")]:
            sub = merged[grp == gv]
            ts, s, _ = km_curve(sub["time_years"].values, sub["event"].values)
            ax.step(ts, s, where="post", color=color, label=f"{label} (n={len(sub)})")
        chi2, p = logrank(merged["time_years"].values, merged["event"].values, grp.values)
        lr_results[g] = {"median_log2cpm": round(float(med), 3),
                         "logrank_chi2": round(float(chi2), 4),
                         "logrank_p": float(f"{p:.6g}")}
        ax.set_title(f"{g}  对数秩 p = {p:.4f}")
        ax.set_xlabel("生存时间（年）")
        ax.set_ylabel("总生存率")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("TCGA-BRCA 原发肿瘤（n={}，事件 {}）—— 基因表达中位数分组的 Kaplan-Meier 曲线".format(n_all, n_event), y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "km_genes.png"), bbox_inches="tight", dpi=180)
    print("KM 图 -> figures/km_genes.png")

    # ---------- Cox ----------
    # 仅 5 个死亡事件，按 EPV≥10 经验法则最多放 1 个协变量，
    # 多变量模型（3 协变量）严重过拟合且晚期组零事件导致完全分离（HR→∞），
    # 故改为 3 个基因各自的单变量 Cox，并在局限中说明。
    cox_df = merged.dropna(subset=["age_years"]).copy()
    cox_out = []
    for g in GENES:
        if g not in cox_df.columns:
            continue
        med = cox_df[g].median()
        X = (cox_df[g] >= med).astype(float).values.reshape(-1, 1)
        beta, se, pvals, H = cox_efron(cox_df["time_years"].values, cox_df["event"].values, X)
        b, s, pv = float(beta[0]), float(se[0]), float(pvals[0])
        ci_lo = np.exp(b - 1.96 * s) if b - 1.96 * s < 700 else np.inf
        ci_hi = np.exp(b + 1.96 * s) if b + 1.96 * s < 700 else np.inf
        cox_out.append({"变量": f"{g}高表达(≥中位)", "logHR": round(b, 4),
                        "HR": round(float(np.exp(np.clip(b, -700, 700))), 3),
                        "HR_95CI": [round(float(ci_lo), 3) if np.isfinite(ci_lo) else "inf",
                                    round(float(ci_hi), 3) if np.isfinite(ci_hi) else "inf"],
                        "p": float(f"{pv:.6g}")})
    print(pd.DataFrame(cox_out).to_string(index=False))

    # 晚期组事件数（供局限说明：是否完全分离）
    stage_ev = merged[merged["stage_grp"] == "late"]["event"].sum() if "stage_grp" in merged else None
    n_late = int((merged["stage_grp"] == "late").sum()) if "stage_grp" in merged else 0

    summary = {
        "数据": {"分析样本": int(n_all), "死亡事件": int(n_event),
                 "蛋白编码基因数": int(expr.shape[0]),
                 "Cox 子样本": int(len(cox_df)),
                 "晚期例数": n_late, "晚期事件数": int(stage_ev) if stage_ev is not None else None},
        "对数秩检验": lr_results,
        "单变量Cox": cox_out,
        "说明": "事件数=5，按EPV≥10仅够单变量模型；多变量Cox因晚期组零事件完全分离，不报告",
    }
    json.dump(summary, open(os.path.join(RES, "summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("结果 -> results/summary.json")


if __name__ == "__main__":
    main()
