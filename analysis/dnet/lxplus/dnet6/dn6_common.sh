# dn6_common.sh — paczka DNET etap 6 (19.11.2026, projekt B): sieć D0 z torami otoczenia BEZ toru zgodnego z wolnym pionem D*→D0π
#   (maska bez prawdy MC: Δm = m(D0π) − m(D0) < DN6_DM_MAX, ładunek toru = −ładunek K; identyczna w danych i MC, sygnał i tło, trening i zastosowanie).
#   Warianty D0.tr.H+e04np.mc.d5.s1–s3 (bloki np1–np3; reszta jak tr.H+e04: ta sama pula, foldy, DisCo masy λ = 5, 20 epok).
#   Węzły główne: trening (GPU: 1 job na ziarno = 5 foldów + łączenie; CPU: fold × 5 + merge) → ana (P1, maska, analyze.py compute) → DONE.
#   Węzły miękkie PO ana (nie blokują DONE): app (100 % zdarzeń, drzewo H + ziarno 1, bramka P1) → anar (dn4_ana.py r --force-data: G_fit) → DONE_R.
# Na lxplus: ~/private/charm_q/Q_lxplus/dnet6. Wymaga ../dnet5, ../dnet4, ../dnet2, ../dnet1 (kod wspólny, BEZ ZMIAN). Kryteria: DNET/KRYTERIA_etap6.md.
DN6=${DN6:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DN6/.." && pwd)}
DN5=${DN5:-$Q/dnet5}
# shellcheck disable=SC1091
source "$DN5/dn5_common.sh"                              # DN4, DN2, DN1, DN_OUT, DN2_OUT, DN4_OUT, DN_FOLDS, dn_put, dn4_env, dn2_env, dn2_kv, dn4_files
: "${DN6_OUT:=$EOS_OUT/charm_q/dnet6}"                   # wyjście: models/<tag>/, folds/<tag>/, sc/, meta/, app/d150/, fail/, dnet6_*.json, DONE, DONE_R
: "${DN6_BLOCKS:=np1 np2 np3}"
: "${DN6_DM_MAX:=0.155}"; : "${DN6_SIGN:=1}"             # maska (dn6_mask.py); zmiana = inne modele → inny DN6_OUT
: "${DN6_EPOCHS:=20}"; : "${DN6_CAP_SIG:=300000}"; : "${DN6_CAP_BKG:=500000}"
: "${DN6_CPUS:=4}"; : "${DN6_FOLD_MEM:=12GB}"; : "${DN6_MERGE_MEM:=12GB}"; : "${DN6_ANA_MEM:=14GB}"      # ≤ 4 rdzenie, ≤ 16 GB (priorytet konta niski)
: "${DN6_MAXJOBS:=30}"
: "${DN6_TRAIN:=gpu}"                                    # gpu = 1 job GPU na ziarno (5 foldów + łączenie) | cpu = fold × K + merge na CPU (tylko tagi bez meta na EOS)
: "${DN6_GPU_WAIT:=2700}"                                # job GPU, który nie wystartuje w 45 min od wysłania, usuwa się sam → DAG pada → DN6_TRAIN=cpu (dn6_switch.sh / ORCH)
: "${DN6_APP:=1}"                                        # 1 = węzły miękkie app + anar po ana (G_fit z fitu piku; bramka P1 ziarna 1)
: "${DN6_APP150:=data_428770:20 data_427929:20 mb_nd:8 fd_nd:2}"
: "${DN6_APP_TAG:=D0.tr.H+e04np.mc.d5.s1}"; : "${DN6_APP_REF:=D0.xgb.H.sb.n.s1}"    # modele w węzłach app: sieć z maską (ziarno 1) + drzewo odniesienia etapu 4
: "${DN6_KEEP:=^D0\\.tr\\.H\\+e04(np)?\\.mc\\.d[0-9]+\\.s[123]\$}"                   # analyze.py: sieć z maską + pełna sieć etapu 2 (porównanie) + odniesienia
export DN6 DN5 DN4 DN2 DN1 DN Q DN6_OUT DN6_BLOCKS DN6_DM_MAX DN6_SIGN DN6_EPOCHS DN6_CAP_SIG DN6_CAP_BKG DN6_CPUS DN6_FOLD_MEM DN6_MERGE_MEM DN6_ANA_MEM \
       DN6_MAXJOBS DN6_TRAIN DN6_GPU_WAIT DN6_APP DN6_APP150 DN6_APP_TAG DN6_APP_REF DN6_KEEP
dn6_tag() { case $1 in np[123]) echo "D0.tr.H+e04np.mc.d${DN_LAM%%,*}.s${1#np}" ;; *) echo "" ;; esac; }
dn6_ready() { [ -s "$DN6_OUT/meta/$1.json" ] && [ -s "$DN6_OUT/sc/$1.npz" ] && [ -s "$DN6_OUT/models/$1/meta.json" ]; }     # tag gotowy (meta na końcu)
