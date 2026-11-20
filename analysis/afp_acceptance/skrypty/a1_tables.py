#!/usr/bin/env python3
"""a1: tabele z res_*.json (lokalnie). Zapis: tabele/afp_tables.md, tabele/summary.json."""
import json, numpy as np, os
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); R = W + '/runs/lxplus/'
mc = json.load(open(R + 'res_mc.json')); sd = mc['mc_mb_sd']
D = {r: json.load(open(R + f'res_{r}.json'))[f'data_{r}'] for r in ('427929', '428770')}
L = []; S = {}
f = lambda e: f"{e[0]:.4f} ± {e[1]:.4f}" if e and e[0] is not None else '–'
TB = ['0.03-0.035', '0.035-0.04', '0.04-0.05', '0.05-0.06', '0.06-0.07', '0.07-0.08', '0.08-0.09', '0.09-0.095', '0.095-0.1']
L.append('## T1. Tag-and-probe near/far (≥1 tor w stacji próbnej | dokładnie 1 tor w stacji tagu, ξ_est(x) w binie); dane vs fast-sim MC SD\n')
for tg in ('tagNEAR_probeFAR', 'tagFAR_probeNEAR'):
    L.append(f'### {tg} (wariant xical dla danych)\n\n| ξ_est | ' + ' | '.join(f'{r} {s}' for r in ('427929', '428770', 'MC') for s in 'AC') + ' |\n|---|' + '---|' * 6)
    for b in TB + ['0.035-0.09', '0.04-0.08', '0.03-0.10']:
        row = []
        for r in ('427929', '428770', 'MC'):
            for s in (0, 1):
                t = (D[r] if r != 'MC' else sd)['tp'][f's{s}|{tg}|xical']
                e = t['bins'][b] if b in t['bins'] else t[b]
                row.append(f"{e['probe>=1'][0]:.4f}" if e['probe>=1'][0] is not None else '–')
        L.append(f'| {b} | ' + ' | '.join(row) + ' |')
    L.append('')
L.append('## T2. P(dokładnie 1 proton po stronie | tag stacji), rdzeń 0,04–0,08\n\n| próbka | strona | near→: far≥1 | far==1 | prot==1 | far→: near≥1 | near==1 | prot==1 | ε_tag ≈ P(prot==1|near)·P(near==1|far) |\n|---|---|---|---|---|---|---|---|---|')
for r in ('427929', '428770', 'MC'):
    for s in (0, 1):
        t = (D[r] if r != 'MC' else sd)['tp']; a = t[f's{s}|tagNEAR_probeFAR|xical']['0.04-0.08']; b = t[f's{s}|tagFAR_probeNEAR|xical']['0.04-0.08']
        e = a['prot==1'][0] * b['probe==1'][0]; S[f'eps_tag_{r}_s{s}'] = e
        S[f'eps_parts_{r}_s{s}'] = dict(far_ge1=a['probe>=1'][0], far_eq1=a['probe==1'][0], prot1_near=a['prot==1'][0], near_ge1=b['probe>=1'][0], near_eq1=b['probe==1'][0], prot1_far=b['prot==1'][0])
        L.append(f"| {r} | {'AC'[s]} | {f(a['probe>=1'])} | {f(a['probe==1'])} | {f(a['prot==1'])} | {f(b['probe>=1'])} | {f(b['probe==1'])} | {f(b['prot==1'])} | {e:.4f} |")
L.append('\n## T3. Wydajność L1 AFP dla tagów offline T2 (odniesienie: 427929 MBTS, 428770 J12)\n\n| run | strona | N | L1_AFP_strona | FS | NS | L1_AFP_A_OR_C |\n|---|---|---|---|---|---|---|')
for r in ('427929', '428770'):
    for s in (0, 1):
        x = D[r]['l1'][f's{s}']['all']; L.append(f"| {r} | {'AC'[s]} | {x['N']} | {f(x['L1_AFP_side'])} | {f(x['L1_FS'])} | {f(x['L1_NS'])} | {f(x['L1_A_OR_C'])} |"); S[f'l1_{r}_s{s}'] = x['L1_A_OR_C'][0]
    L.append(f"| {r} | zależność od ξ (A) | | " + ', '.join(f"{b}: {v['L1_AFP_side'][0]:.3f}" for b, v in D[r]['l1']['s0']['bins'].items()) + ' | | | |')
    L.append(f"| {r} | zależność od ξ (C) | | " + ', '.join(f"{b}: {v['L1_AFP_side'][0]:.3f}" for b, v in D[r]['l1']['s1']['bins'].items()) + ' | | | |')
L.append(f"\nMBTS przy tagu T2 w zdarzeniach HLT_noalg_L1AFP (427929): L1_MBTS_1 {f(D['427929']['mbts_given_afp_tag']['L1_MBTS_1'])}, HLT_mb_mbts {f(D['427929']['mbts_given_afp_tag']['HLT_mb_mbts'])}, N = {D['427929']['mbts_given_afp_tag']['L1_MBTS_1'][2]}")
S['mbts_given_afp_tag_427929'] = D['427929']['mbts_given_afp_tag']
L.append('\n## T4. Krotność protonów po stronie przy ≥1 protonie w oknie (n_prot = 1, 2, 3, 4, 5, ≥6)\n')
for r in ('427929', '428770'):
    for s in (0, 1):
        v = np.array(D[r]['mult'][f's{s}']['nprot_given_nwin>=1'][1:], float); L.append(f"- {r} {'AC'[s]}: {v.astype(int).tolist()} → udział n_prot ≥ 2: {1 - v[0] / v.sum():.4f}"); S[f'multi_{r}_s{s}'] = 1 - v[0] / v.sum()
L.append('\n## T5. Akceptancja fast-sim (MC SD, 2022 twiss): A(ξ_prawda, pT_prawda), strona A\n')
XIB = [0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.095, 0.10, 0.105, 0.11, 0.12]; PTB = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 3.0]
L.append('| ξ \\ pT [GeV] | ' + ' | '.join(f'{PTB[j]}–{PTB[j + 1]}' for j in range(7)) + ' |\n|---|' + '---|' * 7)
for i, row in enumerate(sd['A_xi_pt_s0']): L.append(f'| {XIB[i]}–{XIB[i + 1]} | ' + ' | '.join(f"{e[0]:.3f}" if e[0] is not None else '–' for e in row) + ' |')
L.append('\n## T6. Poprawka fiducjalna C = N(tag rek.)/N(tag prawdziwy) z fast-sim; ε = P(rek. | prawda), f_nonfid = udział tagów rek. spoza obszaru prawdy\n\n| próbka | selekcja | okno ξ | cięcie pT | N_prawda | ε | f_nonfid | C |\n|---|---|---|---|---|---|---|---|')
for nm in ('mc_mb_sd', 'mc_mb_dd', 'mc_fd_nd'):
    for k, v in mc[nm]['fid'].items():
        if 'side' in k or (nm != 'mc_mb_sd' and not k.startswith(('all|0.03-0.1|', 'D_2_30|0.03-0.1|'))): continue
        a, b, c_ = k.split('|'); L.append(f"| {nm} | {a} | {b} | {c_} | {v['N_truth']:.0f} | {v['eff_reco_given_truth']:.4f} | {v['f_nonfid_in_reco']:.4f} | {v['C_reco_over_truth']:.4f} ± {v['C_reco_over_truth'] * v['rel_stat']:.4f} |")
L.append('\n## T7. Wrażliwość C na model (MC SD przeważone): nachylenie b → b + Δb [GeV⁻²], kształt ξ^δ\n\n| wariant | C (0,03–0,10) | C (0,035–0,09) |\n|---|---|---|')
for nm in ('nom', 'b-3', 'b-2', 'b+2', 'b+4', 'xi^-1', 'xi^-0.5', 'xi^+0.5', 'xi^+1'): L.append(f"| {nm} | {sd['reweight'][nm + '|0.03-0.1']:.4f} | {sd['reweight'][nm + '|0.035-0.09']:.4f} |")
L.append(f"\nMC SD: b_eff = 1/⟨pT²⟩ = {sd['b_slope_window_GeV-2']:.2f} GeV⁻² w oknie; udział pT > 0,5 / 0,8 / 1,0 / 1,5 GeV: " + ' / '.join(f"{sd['frac_window_pt_gt'][k]:.3f}" for k in ('0.5', '0.8', '1.0', '1.5')))
L.append(f"ND (fd_nd) proton wiodący w oknie: b_eff = {mc['mc_fd_nd']['b_slope_window_GeV-2']:.2f} GeV⁻², udział pT > 1 GeV {mc['mc_fd_nd']['frac_window_pt_gt']['1.0']:.2f}")
L.append('\n## T8. Krawędzie: rozkład ξ pojedynczego protonu (bin 0,0025), dane 427929 vs MC, unormowane do sumy w 0,04–0,08\n\n| ξ | dane A | dane C | MC A | MC C | dane/MC A | dane/MC C |\n|---|---|---|---|---|---|---|')
h = {}
for nm, src in (('dA', D['427929']['all1']['s0']), ('dC', D['427929']['all1']['s1']), ('mA', sd['all1']['s0']), ('mC', sd['all1']['s1'])):
    x = np.array(src['xi_1prot_0_0.2_80']); h[nm] = x / x[16:32].sum()
for i in range(6, 50):
    L.append(f"| {i * 0.0025:.4f} | {h['dA'][i]:.4f} | {h['dC'][i]:.4f} | {h['mA'][i]:.4f} | {h['mC'][i]:.4f} | {h['dA'][i] / max(h['mA'][i], 1e-9):.2f} | {h['dC'][i] / max(h['mC'][i], 1e-9):.2f} |")
L.append('\n## T9. Rozkład x pojedynczego toru (bin 0,25 mm): NEAR i FAR, dane 427929 (krawędzie detektora)\n\n| x [mm] | near A | near C | far A | far C | far A bez near | far C bez near | MC near A | MC far A |\n|---|---|---|---|---|---|---|---|---|')
a = D['427929']['all1']
for i in range(8, 76):
    L.append(f"| {-20 + i * 0.25:.2f} | " + ' | '.join(str(int(v)) for v in (a['s0']['x_near_1trk_-20_0_80'][i], a['s1']['x_near_1trk_-20_0_80'][i], a['s0']['x_far_1trk_-20_0_80'][i], a['s1']['x_far_1trk_-20_0_80'][i],
                                                                                   a['s0']['x_far_1trk_near0_-20_0_80'][i], a['s1']['x_far_1trk_near0_-20_0_80'][i], sd['all1']['s0']['x_near_1trk_-20_0_80'][i], sd['all1']['s0']['x_far_1trk_-20_0_80'][i])) + ' |')
for r in ('427929', '428770'):
    for s in (0, 1): S[f'xi_of_x_{r}_s{s}'] = D[r][f'xi_of_x_poly2_s{s}']
S['xi_of_x_MC_s0'] = sd['xi_of_x_poly2_s0']; S['fid'] = sd['fid']; S['reweight'] = sd['reweight']
os.makedirs(W + '/tabele', exist_ok=True)
open(W + '/tabele/afp_tables.md', 'w').write('# AFP: tabele akceptancji i wydajności (a1_tables.py; źródło runs/lxplus/res_*.json)\n\n' + '\n'.join(L) + '\n')
json.dump(S, open(W + '/tabele/summary.json', 'w'), indent=1)
print('\n'.join(L[:75])); print(json.dumps({k: v for k, v in S.items() if k.startswith(('eps_tag', 'l1_', 'multi', 'xi_of_x'))}, indent=0)[:3000])
