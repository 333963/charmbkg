# dn2_common.sh — paczka DNET etap 2 (17.11.2026, projekt B): modele na zbiorze torów (córki + otoczenie z drzewa Tracks) dla D0. source, nie uruchamiać.
# Na lxplus: ~/private/charm_q/Q_lxplus/dnet2. Wymaga ../dnet1 (kod wspólny) i wyników etapu 1 na EOS ($DN_OUT/data, sc/, meta/). Tylko pojedyncze D.
DN2=${DN2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$DN2/.." && pwd)}
DN1=${DN1:-$Q/dnet1}; DN=$DN1
# shellcheck disable=SC1091
source "$DN1/dn_common.sh"                               # Q, EOS_OUT, DN_OUT, dn_env, dn_put, q_env
: "${RAWV20:=$EOS_OUT/charm_raw_v20}"                    # ROOT v2.0: <próbka>/*._NNNNNN.*.root (drzewo Tracks)
: "${DN2_OUT:=$EOS_OUT/charm_q/dnet2}"                   # wyjście: parts/dn2_<próbka>_<k>.npz, sc/, meta/, onnx/, dnet2_results.json, DONE
: "${DN2_NPARTS:=data_428770:10 data_427929:8 mb_nd:4 fd_nd:2}"     # węzłów ekstrakcji na próbkę → 24
: "${DN2_FRAC:=data_428770:0.5 data_427929:1.0 mb_nd:0.4 fd_nd:1.0}" # podpróbka zdarzeń tabel a_ (pary w oknie masy: 100 %)
: "${DN2_BLOCKS:=dau dauS ds04 tr04 tr04S2 tr04S3 tr10 tr04x tr04sb}"
: "${DN2_GPU:=1}"                                        # 1 → request_gpus = 1 dla treningu (widok LCG z CUDA); 0 → CPU (8 rdzeni, flavour tomorrow)
: "${DN2_EPOCHS:=20}"; : "${DN2_CAP_SIG:=300000}"; : "${DN2_CAP_BKG:=500000}"; : "${DN2_CPUS:=8}"
export DN2 DN1 DN Q RAWV20 DN2_OUT DN2_NPARTS DN2_FRAC DN2_BLOCKS DN2_GPU DN2_EPOCHS DN2_CAP_SIG DN2_CAP_BKG DN2_CPUS
dn2_kv() { local x; for x in $1; do [ "${x%%:*}" = "$2" ] && { echo "${x##*:}"; return; }; done; echo "$3"; }
dn2_env() {   # GPU: widok LCG z CUDA (jak s2_flow/s2_env.sh), inaczej q_env
  local cv; cv=$(ls -d /cvmfs/sft.cern.ch/lcg/views/LCG_1[01][0-9]*_cuda/x86_64-el9-gcc*-opt 2>/dev/null | sort -V | tail -1)
  if [ "${DN2_USE_GPU:-0}" = 1 ] && [ -n "$cv" ] && command -v nvidia-smi >/dev/null 2>&1; then
    set +u; source "$cv/setup.sh"; set -u; echo "== LCG CUDA: $cv"; nvidia-smi -L 2>/dev/null || true
  else
    q_env
  fi
  export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=$DN2_CPUS MKL_NUM_THREADS=$DN2_CPUS DN_CPUS=$DN2_CPUS
  python3 -c 'import numpy, torch, uproot, awkward, yaml, pandas; print("torch", torch.__version__, "cuda", torch.cuda.is_available())' || { echo "BŁĄD: środowisko"; return 1; }
}
