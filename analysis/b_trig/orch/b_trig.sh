# B_TRIG (17.11.2026, projekt B): odporność SF klasyfikatora H (tag D*), ε_B i izolacji na trigger (428770: J12only / MBTSbased …)
#   i niezależny run 427929. Paczka charm_work/B_trig/lx (+ *.py z charm_work/B_trig) → ~/private/charm_q/Q_lxplus/b_trig.
#   Ścieżka A (BT_PATH=A): tabele zdarzeń z ROOT v2.0 (tylko rucio_v20) → trg_<próbka>.npz. Ścieżka B (AB): + btrig_extract na h5_v20.
#   Wynik → charm_work/B_trig/runs/lxplus/ (+ diag.txt). Kopia w ORCH/stages/b_trig.sh (17.11); wpis do pipeline.yaml
#   (do zatwierdzenia):
#     - name: b_trig      # B_trig: SF/ε_B/izo per trigger + run 427929 (projekt B) — po h5_v20; ~1–2 h
#       needs: [h5_v20]
#       max_hours: 8
#   Wariant tylko-A przed h5_v20: needs: [rucio_v20], BT_PATH=A w stage_check/stage_submit (zmienna poniżej).
source "$HERE/stages/_common.sh"
: "${BT_PATH:=AB}"
SRC=$W/B_trig; R=$RQ/b_trig; EOSD=$EOSQ/b_trig; DST=$W/B_trig/runs/lxplus
stage_available() { [ -f "$SRC/lx/submit_bt.sh" ] && [ -f "$SRC/btrig_extract.py" ]; }
stage_sync()   { rssh "mkdir -p ~/$R ~/$RQ/logs/b_trig && touch ~/$R/.orch_start" \
                 && rpush --exclude lists/ --exclude 'bt.dag*' "$SRC/lx/" "$H:$R/" \
                 && rpush --include 'btrig_common.py' --include 'btrig_extract.py' --include 'trg_root_extract.py' --exclude '*' "$SRC/" "$H:$R/"; }
stage_check()  { rssh "cd ~/$R && chmod +x *.sh *.py && BT_PATH=$BT_PATH FORCE=${FORCE:-0} ./submit_bt.sh check" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_submit() { rssh "cd ~/$R && BT_PATH=$BT_PATH FORCE=${FORCE:-0} ./submit_bt.sh submit; ./submit_bt.sh status | tail -5" | tail -10 | tee -a "$LOG" >/dev/null; return "${PIPESTATUS[0]}"; }
stage_done()   { mkdir -p "$DST"; rssh "condor_q -batch 2>/dev/null | grep -E 'charm_btrig'; cd ~/$R; echo '## dagman.out'; tail -n 15 bt.dag.dagman.out 2>/dev/null; echo '## logs'; ls -lt ~/$RQ/logs/b_trig | head -6; for f in \$(ls -t ~/$RQ/logs/b_trig/*.err 2>/dev/null | head -2); do echo \"### \$f\"; tail -n 15 \$f; done; echo '## EOS'; ls -la $EOSD 2>/dev/null | head -20; ls $EOSD/parts 2>/dev/null | wc -l" > "$DST/diag.txt" 2>&1
                 if [ "$BT_PATH" = A ]; then rssh "[ $EOSD/trg_data_427929.npz -nt ~/$R/.orch_start ] && [ $EOSD/trg_data_428770.npz -nt ~/$R/.orch_start ]";
                 else rssh "[ -f ~/$R/.orch_start ] && [ -f $EOSD/DONE ] && [ $EOSD/DONE -nt ~/$R/.orch_start ]"; fi; }
stage_fetch()  { mkdir -p "$DST" && rpull --include 'bt_*.npz' --include 'trg_*.npz' --include '*.json' --include DONE --exclude 'parts/' --exclude '*' "$H:$EOSD/" "$DST/"; }
stage_result() { echo "B_trig/runs/lxplus/{bt_*,trg_*}.npz → analiza: bash charm_work/B_trig/run_analysis.sh → runs/lxplus/out/RAPORT_DANE.md"; }
