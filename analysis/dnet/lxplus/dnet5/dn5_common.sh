# dn5_common.sh — paczka DNET etap 5a (19.11.2026, projekt B): trening zrównoleglony po foldach + warianty niezależne od wyników etapu 4.
#   (1) infrastruktura: fold = osobny job, węzeł łączący (score OOF + meta w formacie etapu 4), wagi foldów na EOS;
#   (2) krzywa uczenia D0.tr.H+e04 przy 25/50 % sygnału MC (100 % = istniejące score etapu 2, ziarna 1–3);
#   (3) cechy par w uwadze wg ParT (D0, 1 ziarno); (4) tagger otoczenia E (D0 i D+; trening na wszystkich kandydatach i po cięciu H WP70);
#   (5) sonda GPU (osobny job poza DAG, sam wycofuje się po DN5_GPU_WAIT s w kolejce).
# Na lxplus: ~/private/charm_q/Q_lxplus/dnet5. Wymaga ../dnet4, ../dnet2, ../dnet1 (kod wspólny, BEZ ZMIAN) oraz na EOS: $DN_OUT/{data,sc,meta}
# (etap 1), $DN2_OUT/{parts,sc} (etap 2), $DN4_OUT/parts_dp (tory otoczenia D+ z etapu 4). NIE czeka na wyniki etapu 4. Kryteria: DNET/KRYTERIA_etap5.md.
DN5=${DN5:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DN5/.." && pwd)}
DN4=${DN4:-$Q/dnet4}
# shellcheck disable=SC1091
source "$DN4/dn4_common.sh"                              # DN1, DN2, DN_OUT, DN2_OUT, DN4_OUT, DN4_EXTDP, DN_FOLDS, dn_put, dn4_env, dn2_env, dn2_kv
: "${DN5_OUT:=$EOS_OUT/charm_q/dnet5}"                   # wyjście: models/<tag>/, folds/<tag>/, sc/, meta/, test/, gpu/, fail/, dnet5_*.json, DONE
: "${DN5_BLOCKS:=lc25 lc50 pair eD0 eD0h eDp eDph}"
: "${DN5_EPOCHS:=20}"; : "${DN5_CAP_SIG:=300000}"; : "${DN5_CAP_BKG:=500000}"
: "${DN5_CPUS:=4}"                                       # rdzeni na job foldu (KRYTERIA_etap5.md §1). UWAGA: CERN przydziela rdzenie wg pamięci (~3 GB/rdzeń): 24 GB → 8 rdzeni → długa kolejka; trzymać ≤ 12 GB
: "${DN5_FOLD_MEM:=12GB}"; : "${DN5_MERGE_MEM:=12GB}"; : "${DN5_ANA_CPUS:=4}"
: "${DN5_MAXJOBS:=40}"                                   # condor_submit_dag -maxjobs
: "${DN5_GPU:=1}"; : "${DN5_GPU_WAIT:=21600}"            # sonda GPU: 1 = wyślij; czas w kolejce [s], po którym job sam się usuwa (periodic_remove)
: "${DN5_GPU_BLOCKS:=pair base}"                         # sonda liczy fold 0 tych bloków
: "${DN5_TRAIN:=gpu}"                                   # gpu = 1 job GPU na wariant (5 foldów + łączenie; od 19.11 07:20, po sondzie GPU) | cpu = fold × K + merge na CPU
: "${DN5_GPU_TRAIN_WAIT:=5400}"                          # job gtrain, który nie wystartuje w tym czasie [s], sam się usuwa → DAG pada → DN5_TRAIN=cpu ./submit_dn5.sh submit
: "${DN5_EQBENCH:=1}"                                    # węzeł miękki: test równoważności + skalowanie wątków na węźle roboczym
: "${DN5_KEEP:=^D0\\.(tr|trp)\\.H\\+e04\\.(mc|mcS25|mcS50)\\.d[0-9]+\\.s[123]\$}"      # konfiguracje do analyze.py (odniesienia dokłada sam)
export DN5 DN4 DN2 DN1 DN Q DN5_OUT DN5_BLOCKS DN5_EPOCHS DN5_CAP_SIG DN5_CAP_BKG DN5_CPUS DN5_FOLD_MEM DN5_MERGE_MEM DN5_ANA_CPUS DN5_MAXJOBS DN5_GPU \
       DN5_GPU_WAIT DN5_GPU_BLOCKS DN5_EQBENCH DN5_KEEP DN5_TRAIN DN5_GPU_TRAIN_WAIT
dn5_tag() { python3 "$DN5/dn5_train.py" list | awk -v b="$1" '$1 == b {print $2}'; }
