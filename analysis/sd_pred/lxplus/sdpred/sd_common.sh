# sd_common.sh — paczka sdpred (17.11.2026, projekt B; source). Wymaga ../common.sh (Q, EOS_OUT, q_env, put_output, eos_exists).
SD=${SD:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
Q=${Q:-$(cd "$SD/.." && pwd)}
# shellcheck disable=SC1091
source "$Q/common.sh"
: "${SD_OUT:=$EOS_OUT/charm_q/sdpred}"      # <model>/part_<seed>.{tsv.gz,tsv.json,npz}, merged/<model>.json, DONE
: "${SD_MODELS:=$SD/models.txt}"
: "${SD_ONLY:=}"                            # lista modeli po przecinku (puste = wszystkie)
: "${SD_NEV:=}"                             # nadpisanie liczby zdarzeń na job (test)
: "${SD_NJOBS:=}"                           # nadpisanie liczby jobów na model (test)
: "${MAXJOBS:=150}"
export SD Q SD_OUT SD_MODELS
sd_models() {   # → wiersze: indeks model njobs nev cmnd  (indeks = numer wiersza w models.txt → stałe ziarna)
  awk -v only="$SD_ONLY" -v nev="$SD_NEV" -v nj="$SD_NJOBS" '
    /^[ \t]*(#|$)/ {next} {i++; if (only != "" && index(","only",", ","$1",") == 0) next;
      print i, $1, (nj != "" ? nj : $2), (nev != "" ? nev : $3), $4}' "$SD_MODELS"
}
sd_build() {    # kompilacja bin/pythia_sd (Pythia 8 z LCG) — w podpowłoce, żeby LCG nie zostało w env submitu
  ( q_env >/dev/null 2>&1; command -v pythia8-config >/dev/null || { echo "BŁĄD: brak pythia8-config w LCG"; exit 1; }
    mkdir -p "$SD/bin"
    if [ ! -x "$SD/bin/pythia_sd" ] || [ "$SD/pythia_sd.cc" -nt "$SD/bin/pythia_sd" ]; then
      g++ -O2 -std=c++17 $(pythia8-config --cxxflags) "$SD/pythia_sd.cc" -o "$SD/bin/pythia_sd" $(pythia8-config --libs) -ldl || exit 1
      echo "  skompilowano bin/pythia_sd (Pythia $(pythia8-config --version))"
    fi )
}
