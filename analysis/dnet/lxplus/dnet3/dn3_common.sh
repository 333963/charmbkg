# dn3_common.sh — paczka DNET etap 3 (18.11.2026, projekt B): ablacje torów otoczenia i dekorelacja od aktywności dla transformera D0
# (tr.H+e04, DisCo masy λ = 5) na danych JUŻ wyekstrahowanych w etapie 2. source, nie uruchamiać.
# Na lxplus: ~/private/charm_q/Q_lxplus/dnet3. Wymaga ../dnet2 i ../dnet1 (kod wspólny, bez zmian) oraz wyników na EOS:
#   $DN_OUT/data (tabele etapu 1), $DN2_OUT/parts (tory otoczenia), $DN_OUT/{sc,meta} i $DN2_OUT/{sc,meta} (odniesienia do analizy).
DN3=${DN3:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DN3/.." && pwd)}
DN2=${DN2:-$Q/dnet2}
# shellcheck disable=SC1091
source "$DN2/dn2_common.sh"                              # DN1, DN_OUT, DN2_OUT, DN_SAMPLES, DN_FOLDS, DN_LAM, dn_put, dn2_env, dn2_kv, q_env
: "${DN3_OUT:=$EOS_OUT/charm_q/dnet3}"                   # wyjście: sc/, meta/, dnet3_results.json, DONE
: "${DN3_BLOCKS:=tr04g tr04pt1 actS2 actS10 actB10}"     # poza domyślną listą: tr04np tr04c ref (dn3_train.py list)
: "${DN3_GPU:=1}"                                        # 1 → trening z request_gpus = 1 (dn3_gpu.sub, Capability > 7.0); 0 → CPU (8 rdzeni, flavour tomorrow)
: "${DN3_SMOKE:=1}"                                      # 1 → węzeł smoke (dn3_train.py selftest, CPU) przed treningami
: "${DN3_EPOCHS:=$DN2_EPOCHS}"; : "${DN3_CAP_SIG:=$DN2_CAP_SIG}"; : "${DN3_CAP_BKG:=$DN2_CAP_BKG}"; : "${DN3_CPUS:=8}"   # jak etap 2 (ta sama pula)
export DN3 DN2 DN1 DN Q DN3_OUT DN3_BLOCKS DN3_GPU DN3_SMOKE DN3_EPOCHS DN3_CAP_SIG DN3_CAP_BKG DN3_CPUS
dn3_env() {   # środowisko jak w etapie 2 (GPU: widok LCG z CUDA gdy DN3_USE_GPU=1 i jest karta; inaczej q_env)
  DN2_USE_GPU=${DN3_USE_GPU:-0}; DN2_CPUS=$DN3_CPUS; export DN2_USE_GPU DN2_CPUS
  dn2_env
}
