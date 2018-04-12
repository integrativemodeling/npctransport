#!/bin/bash
imp=$HOME/imp_git/fast/setup_environment.sh
imppy="$imp python"
npc_show_stats="$imppy $HOME/imp_git/repository/modules/npctransport/utility/show_statistics.py"

SKIP_ROWS=${1:-100}
for fg in  FSFG_generic GLFG_generic Nup1 Nup100 Nup49 Nsp1; do
    echo -n "$fg ";
    n_fgs=$(cat config_$fg.txt | awk 'BEGIN{is=0}{if(is==1){print $NF; is=2; } else { if($1 ~ /number_of_b/){is=1}}}')
    echo -n "n_fgs=$n_fgs n_res=$((n_fgs*20)) "
    # In next line, correct Rg from coarse grain to compare directly to SAXS data,
    # using Rg_real^2 = Rg_CG^2 + 0.8*Rbead^2, where 0.8Rbead^2 is the formula for a sphere's Rg
    Rbead=8.0
    $npc_show_stats output_$fg.pb | grep mean_radius_of | awk '{if(NR>'"$SKIP_ROWS"'){Rbead='"$Rbead"'; rg=sqrt($2*$2+0.8*Rbead*Rbead); s=s+rg; s2=s2+rg*rg; n=n+1}}END{e=s/n; e2=s2/n; var=e2-e*e; sem=sqrt(var/(n-1)); print e, "+-", 1.96*sem, "[95% conf]  std-dev", sqrt(var), "  n =", n, " For nres=120 at Flory=0.5: ", e*sqrt(120/20.0/'"$n_fgs"'), " at Flory=0.4: ", e*(120/20.0/'"$n_fgs"')^0.4}';
done;
