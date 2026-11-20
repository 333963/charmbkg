#!/bin/bash
# FONLL (formularz sieciowy; metoda jak XSEC_2026-11-19/xsec2/fonll/run.sh): D* i D+ przy 7 TeV (|eta| < 2,1) i 13,6 TeV (|eta| < 2,2), prompt i B -> D(*).
cd "$(dirname "$0")"
req(){ # collider eta quark final decay ptmin ptmax np out
 [ -s "$9" ] && grep -q '^ *[0-9]' "$9" && return
 curl -sS -m 900 -o "$9" -w "$9 %{http_code} %{time_total}\n" -e 'https://www.lpthe.jussieu.fr/~cacciari/fonll/fonllform.html' -A 'Mozilla/5.0' \
 --data-urlencode collider=$1 --data-urlencode pdf=260000 --data-urlencode order=1 --data-urlencode heavy_quark=$3 --data-urlencode "final_state=$4" --data-urlencode "decay=$5" --data-urlencode xsect_type=2 --data-urlencode band=3 --data-urlencode pdfunc=on --data-urlencode ptmin=$6 --data-urlencode ptmax=$7 --data-urlencode ymin=-$2 --data-urlencode ymax=$2 --data-urlencode use_eta=1 --data-urlencode npoints=$8 --data-urlencode npff=0 --data-urlencode nppar=default --data-urlencode BRcD=1 --data-urlencode BRbB=1 --data-urlencode BRBD=0.823 --data-urlencode BRBDs=0.173 --data-urlencode email=0 \
 https://www.lpthe.jussieu.fr/~cacciari/cgi-bin/call_band-18.pl
 sleep 2
}
for ce in "7 2.1" "19 2.2"; do
 set -- $ce; c=$1; e=$2
 for m in "D*" "D+"; do
  mm=${m/+/p}; mm=${mm/\*/st}
  req $c $e charm "$m" "-" 3.5 5 7  prompt_${mm}_c${c}_eta${e}_pt3.5-5.txt
  req $c $e charm "$m" "-" 5 8 13   prompt_${mm}_c${c}_eta${e}_pt5-8.txt
  req $c $e charm "$m" "-" 8 30 45  prompt_${mm}_c${c}_eta${e}_pt8-30.txt
  [ $c = 19 ] && req $c $e charm "$m" "-" 3 5 9 prompt_${mm}_c${c}_eta${e}_pt3-5.txt
 done
 for d in "B to D*" "B to D"; do
  dd=${d// /}; dd=${dd/\*/st}
  req $c $e bottom "B hadron" "$d" 3.5 5 7  ${dd}_c${c}_eta${e}_pt3.5-5.txt
  req $c $e bottom "B hadron" "$d" 5 8 13   ${dd}_c${c}_eta${e}_pt5-8.txt
  req $c $e bottom "B hadron" "$d" 8 30 45  ${dd}_c${c}_eta${e}_pt8-30.txt
  [ $c = 19 ] && req $c $e bottom "B hadron" "$d" 3 5 9 ${dd}_c${c}_eta${e}_pt3-5.txt
 done
done
echo DONE
